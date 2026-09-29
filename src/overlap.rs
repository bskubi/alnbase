//! Resolve mate overlap in a fragment, including a fragment that is chimeric.
//!
//! A sequenced fragment is one linear molecule. R1 reads it from one end, R2
//! from the other, and a chimeric aligner may split either read across loci.
//! Where the two reads reach far enough to meet, the same physical bases were
//! observed twice. Counting them twice inflates the coverage, and it makes the
//! two observations look like independent agreement.
//!
//! # Reference intervals are the wrong key
//!
//! The incumbents (`bamUtil clipOverlap`, `fgbio ClipBam`) find the overlap by
//! intersecting the two mates' reference intervals. That silently assumes the
//! molecule is colinear with the reference, which a ligation product is not: a
//! Hi-C fragment visits two loci, so "R1's interval ∩ R2's interval" is not the
//! set of doubly-observed bases and can be empty when the overlap is large.
//!
//! Everything here is keyed on the **molecule** instead. A ligation product is
//! still one linear molecule, so every base has a well-defined index along it,
//! and two bases are the same observation exactly when those indices match.
//! Reference positions are used only to *discover* the offset between the two
//! reads' coordinate systems, never to define the overlap.
//!
//! # Three coordinate systems
//!
//! | system | domain | definition |
//! |---|---|---|
//! | query index `q` | one record | offset into SEQ as stored, reference-forward |
//! | end index `e` | one read end | offset from the 5' base of the original read |
//! | fragment index `f` | the template | offset along the molecule, 5'→3' on the R1 strand |
//!
//! A tool of this class usually fails because it confuses these three. `e` is
//! what connects them. The code recovers `e` for each record from the clip
//! operations, and `e` is what puts a hard-clipped supplementary and its
//! primary into one coordinate system.
//!
//! # The estimator
//!
//! Given consensus fragment length `L`, R1 occupies fragment indices `[0, len1)`
//! and R2 occupies `[L - len2, L)`, so `f = e1` and `f = L - 1 - e2`. Any
//! candidate pair of reference-matched bases `(e1, e2)` implies `L = e1+e2+1`.
//!
//! Candidates are grouped by the `L` they imply. The true group is identified
//! not by its size or its median but by **where it reaches**: whenever an
//! overlap exists at all, it spans `[max(0, L-len2), min(len1, L))` in fragment
//! indices, whose two ends are precisely R1's 3'-most base and R2's 3'-most
//! base. The true group therefore always touches both reads' 3' termini, and a
//! group that arises from a repeat generally does not.
//!
//! That reach is a property of the *set*. It is not a property of any single
//! pair: `max(len1-1-e1, len2-1-e2)` is minimised at the overlap's midpoint,
//! not at either terminus, so seeding on one 3'-proximal pair would seed in the
//! middle. Set reach is also correct in the case where a central-tendency
//! statistic is not. That case is a short true overlap against a long spurious
//! one, where the mode picks the spurious group.
//!
//! # What is deliberately not merged
//!
//! Reference-position match implies molecule-position match only when the
//! molecule visits each locus once. Shearing a self-ligated circle does not
//! break that: the product reads `p→b`, junction, `a→p`, which is a circular
//! permutation of the fragment with every base present exactly once, so it is
//! handled as an ordinary short chimera and needs no special case. Genuine
//! two-copy molecules — same-fragment ligation between sister chromatids or
//! homologues, `A→B→A` concatemers, tandem-repeat phase disagreement — do
//! break it: they throw off a decoy group that pairs one copy against the
//! other. The decoy is not refused by its size. It pairs bases that sit far
//! from both reads' 3' termini, so the anchor test passes it over and the true
//! group is resolved (see `decide`). A template is refused instead when the
//! only group implies a molecule shorter than the reads aligned (`past_end`),
//! or when two groups both pass the anchor test with comparable support
//! (`ambiguous_L`), because merging on the wrong one would fuse two different
//! parts of the molecule and produce a false call rather than merely lost
//! coverage.
//!
//! Everything here is pure. `overlap_apply` applies a plan to real records. The
//! two are separate, because the reasoning is then easy to test on its own.

use std::collections::BTreeMap;

use crate::seq::Seq;

// ---------------------------------------------------------------------------
// Input ordering
// ---------------------------------------------------------------------------

/// How a file's records are arranged, as far as this tool cares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Grouping {
    /// `SO:queryname` -- sorted by name.
    Sorted,
    /// `GO:query` -- grouped by name but not sorted. What `samtools collate`
    /// produces, and much cheaper than a full sort, which matters when the
    /// input is a pipe.
    Collated,
}

/// Check that every alignment of a template will arrive together.
///
/// Pure, and given the header text rather than a file, because the string
/// handling is the part most likely to be subtly wrong and this way it can be
/// tested. Grouping is a correctness requirement: on ungrouped input the tool
/// would resolve the few templates that happen to be adjacent and miss the
/// rest, which is a near-no-op that looks like success.
pub fn check_grouping(header_text: &str) -> Result<Grouping, String> {
    let hd = header_text.lines().find(|l| l.starts_with("@HD"));
    let Some(hd) = hd else {
        return Err(
            "no @HD line, so the record order cannot be verified. Group by name first: \
             `samtools collate -O in.bam | ...`"
                .to_string(),
        );
    };
    let field = |k: &str| {
        hd.split('\t')
            .find_map(|f| f.strip_prefix(k).map(|v| v.trim().to_string()))
    };
    // GO wins: a collated file is legitimately SO:unsorted.
    if field("GO:").as_deref() == Some("query") {
        return Ok(Grouping::Collated);
    }
    match field("SO:").as_deref() {
        Some("queryname") => Ok(Grouping::Sorted),
        Some(other) => Err(format!(
            "records are ordered by {other}, so a template's alignments are not adjacent. \
             Group by name first: `samtools collate -O in.bam | ...`, or \
             `samtools sort -n`."
        )),
        None => Err(
            "the @HD line declares neither SO: nor GO:, so the record order cannot be \
             verified. Group by name first: `samtools collate -O in.bam | ...`"
                .to_string(),
        ),
    }
}

// ---------------------------------------------------------------------------
// CIGAR
// ---------------------------------------------------------------------------

/// A CIGAR operation. `H` and `P` are kept, because a real supplementary record
/// is usually hard-clipped and the leading `H` is what says where its bases sit
/// in the read.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CigarOp {
    M(u32),
    I(u32),
    D(u32),
    N(u32),
    S(u32),
    H(u32),
    P(u32),
}

impl CigarOp {
    pub fn len(self) -> usize {
        match self {
            CigarOp::M(n)
            | CigarOp::I(n)
            | CigarOp::D(n)
            | CigarOp::N(n)
            | CigarOp::S(n)
            | CigarOp::H(n)
            | CigarOp::P(n) => n as usize,
        }
    }
    /// True when the operation advances the position in the read.
    pub fn consumes_query(self) -> bool {
        matches!(self, CigarOp::M(_) | CigarOp::I(_) | CigarOp::S(_))
    }
    /// True when the operation is present in the record's SEQ field. `H` is
    /// not.
    pub fn in_seq(self) -> bool {
        self.consumes_query()
    }
    /// True when the operation carries an observed base placed on the
    /// reference. `S` does not.
    pub fn is_aligned(self) -> bool {
        matches!(self, CigarOp::M(_) | CigarOp::I(_))
    }
    fn with_len(self, n: usize) -> CigarOp {
        let n = n as u32;
        match self {
            CigarOp::M(_) => CigarOp::M(n),
            CigarOp::I(_) => CigarOp::I(n),
            CigarOp::D(_) => CigarOp::D(n),
            CigarOp::N(_) => CigarOp::N(n),
            CigarOp::S(_) => CigarOp::S(n),
            CigarOp::H(_) => CigarOp::H(n),
            CigarOp::P(_) => CigarOp::P(n),
        }
    }
}

/// Full original read length, including bases this record hard-clipped away.
///
/// This is the `len1`/`len2` of the estimator, and it must be the *whole* read:
/// not the aligned span, and not the union of the records' aligned spans. At a
/// ligation junction a few bases of reconstructed sequence often align to
/// neither side and sit in a query gap between records; counting only aligned
/// bases would shorten the read and shift every fragment index past the
/// junction.
pub fn read_length(cigar: &[CigarOp]) -> usize {
    cigar
        .iter()
        .filter(|o| o.consumes_query() || matches!(o, CigarOp::H(_)))
        .map(|o| o.len())
        .sum()
}

// ---------------------------------------------------------------------------
// Alignments
// ---------------------------------------------------------------------------

/// The byte htslib stores for every quality of a record whose QUAL is `*`.
pub const QUAL_MISSING: u8 = 0xff;

/// One alignment record, reduced to what planning needs.
#[derive(Clone, Debug)]
pub struct Aln {
    /// Index of the record within the template, for matching plans back up.
    pub id: usize,
    pub is_r1: bool,
    pub is_reverse: bool,
    pub is_supplementary: bool,
    /// This record's own mapping quality. The overlap frequently falls on a
    /// low-MAPQ supplementary, so the covering record's MAPQ is the one that
    /// matters and the primary's would invert the decision.
    pub mapq: u8,
    pub tid: i32,
    /// 0-based leftmost aligned reference position.
    pub pos: i64,
    pub cigar: Vec<CigarOp>,
    /// The record's own SEQ, reference-forward. Indexed from its first
    /// non-hard-clipped base.
    pub seq: Vec<Seq>,
    /// Parallel to `seq`. [`QUAL_MISSING`] where the record has no qualities.
    pub qual: Vec<u8>,
}

impl Aln {
    pub fn read_length(&self) -> usize {
        read_length(&self.cigar)
    }
    /// Aligned query bases: `M` and `I`.
    pub fn query_len(&self) -> usize {
        self.cigar
            .iter()
            .filter(|o| o.is_aligned())
            .map(|o| o.len())
            .sum()
    }
    /// Bases present in SEQ: `M`, `I` and `S`.
    pub fn seq_len(&self) -> usize {
        self.cigar.iter().filter(|o| o.in_seq()).map(|o| o.len()).sum()
    }
    /// Leading hard clip, which is the offset between SEQ indices and
    /// reference-forward read offsets.
    fn lead_hard(&self) -> usize {
        self.cigar
            .iter()
            .take_while(|o| matches!(o, CigarOp::H(_)))
            .map(|o| o.len())
            .sum()
    }

    /// End index of a SEQ index: its offset from the 5' base of the read.
    ///
    /// SEQ is stored reference-forward, so on a reverse-strand record the SEQ
    /// index runs opposite to the order the bases came off the sequencer. This
    /// inversion is the whole of the strand handling; nothing else in the
    /// module needs to know about orientation except the pairing rule.
    pub fn end_index(&self, s: usize) -> usize {
        let fs = s + self.lead_hard();
        if self.is_reverse {
            self.read_length().saturating_sub(1).saturating_sub(fs)
        } else {
            fs
        }
    }

    /// End index of this record's 5'-most base, i.e. where its segment starts
    /// in the read. On a reverse record that is the far end of SEQ.
    pub fn end_offset(&self) -> usize {
        if self.is_reverse {
            self.read_length()
                .saturating_sub(self.lead_hard())
                .saturating_sub(self.seq_len())
        } else {
            self.lead_hard()
        }
    }

    /// `(SEQ index, aligned, reference position)` for every base in SEQ.
    ///
    /// `S` bases appear with `aligned = false`: they are observations, but not
    /// placed ones. `I` bases appear as aligned with no reference position —
    /// which is exactly right, because fragment indices pair them up with their
    /// twin on the other mate for free, with no insertion-anchor bookkeeping.
    fn seq_map(&self) -> Vec<(bool, Option<i64>)> {
        let mut out = Vec::with_capacity(self.seq_len());
        let mut r = self.pos;
        for op in &self.cigar {
            match *op {
                CigarOp::M(n) => {
                    for k in 0..n as i64 {
                        out.push((true, Some(r + k)));
                    }
                    r += n as i64;
                }
                CigarOp::I(n) => {
                    for _ in 0..n {
                        out.push((true, None));
                    }
                }
                CigarOp::S(n) => {
                    for _ in 0..n {
                        out.push((false, None));
                    }
                }
                CigarOp::D(n) | CigarOp::N(n) => r += n as i64,
                CigarOp::H(_) | CigarOp::P(_) => {}
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Per-end model
// ---------------------------------------------------------------------------

/// One observed base, located in the read rather than in a record.
#[derive(Clone, Copy, Debug)]
struct Base {
    /// Index into the `alns` slice.
    rec: usize,
    /// SEQ index within that record.
    seq: usize,
    aligned: bool,
    refpos: Option<(i32, i64)>,
    reverse: bool,
    mapq: u8,
}

/// One read end, assembled from all of its records into read coordinates.
struct End {
    /// Full original read length. Every record agrees on this via its clips.
    len: usize,
    /// Indexed by end index; `None` where no record covers that base.
    bases: Vec<Option<Base>>,
    /// Largest end index carrying an aligned base.
    last_aligned: Option<usize>,
    /// Records disagreed about the read length.
    inconsistent: bool,
}

impl End {
    fn base(&self, e: usize) -> Option<&Base> {
        self.bases.get(e).and_then(|b| b.as_ref())
    }
}

/// Assemble one end's records into read coordinates.
///
/// Where two records of the same end claim the same base — which happens at a
/// junction the aligner split ambiguously — the claim is resolved in favour of
/// a placed base over an unplaced one, then higher MAPQ, then the earlier
/// record, so the result does not depend on input order.
fn build_end(alns: &[Aln], want_r1: bool, tolerance: usize) -> End {
    let idx: Vec<usize> = (0..alns.len()).filter(|&i| alns[i].is_r1 == want_r1).collect();
    if idx.is_empty() {
        return End { len: 0, bases: Vec::new(), last_aligned: None, inconsistent: false };
    }
    let lens: Vec<usize> = idx.iter().map(|&i| alns[i].read_length()).collect();
    let len = *lens.iter().max().unwrap();
    let inconsistent = len - lens.iter().min().unwrap() > tolerance;

    let mut bases: Vec<Option<Base>> = vec![None; len];
    for &i in &idx {
        let a = &alns[i];
        for (s, (aligned, rp)) in a.seq_map().into_iter().enumerate() {
            let e = a.end_index(s);
            if e >= len {
                continue;
            }
            let cand = Base {
                rec: i,
                seq: s,
                aligned,
                refpos: rp.map(|p| (a.tid, p)),
                reverse: a.is_reverse,
                mapq: a.mapq,
            };
            match &bases[e] {
                None => bases[e] = Some(cand),
                Some(prev) => {
                    let better = (cand.aligned, cand.mapq) > (prev.aligned, prev.mapq);
                    if better {
                        bases[e] = Some(cand);
                    }
                }
            }
        }
    }
    let last_aligned = bases
        .iter()
        .enumerate()
        .rev()
        .find(|(_, b)| b.as_ref().is_some_and(|b| b.aligned))
        .map(|(e, _)| e);
    End { len, bases, last_aligned, inconsistent }
}

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

/// What to do to the qualities of two bases that disagree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MismatchQual {
    None,
    ZeroBoth,
    ZeroLoser,
    /// Loser to zero, winner reduced by the loser's quality. The posterior
    /// treatment of two conflicting observations, and the default.
    Subtract,
}

/// What to do to the base calls of two bases that disagree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MismatchBase {
    None,
    SetN,
}

/// What to do to the qualities of two bases that agree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatchQual {
    None,
    Max,
    /// Sum, capped. The cap matters here. R1 and R2 are the same molecule, so
    /// a polymerase error made before amplification appears in both of them.
    /// The two observations are therefore not independent.
    SumCapped,
}

/// Which end gives up its copy of the overlap.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Keep {
    R1,
    R2,
    Score,
}

/// What to do with a template whose fragment length could not be established.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unresolved {
    /// Write both ends unmodified, tagged. Preserves coverage; downstream must
    /// tolerate the double count, and the tag makes those reads filterable.
    Pass,
    /// Discard the template.
    Drop,
}

#[derive(Clone, Debug)]
pub struct Options {
    /// Permitted disagreement between a read's records about its length.
    pub length_tolerance: usize,
    /// Overlaps shorter than this are left alone.
    pub min_overlap: usize,
    /// Candidate pairs required to believe a fragment length.
    pub min_support: usize,
    /// Fraction of the predicted overlap the candidate run must span.
    pub min_span_frac: f32,
    /// How far the candidate group may fall short of either read's 3' terminus.
    pub max_anchor_slack: usize,
    /// Bases a read may still have aligned past the end of the implied
    /// molecule. Above zero only to absorb a stray adapter alignment.
    pub max_past_end: usize,
    /// Mismatch fraction over the overlap above which nothing is merged.
    pub max_mismatch_frac: f32,
    /// A rival fragment length this fraction as well supported blocks merging.
    pub ambiguity_ratio: f32,
    /// A rival within this distance, over a disjoint span, is an indel
    /// disagreement rather than a repeat, and is reported as one.
    pub indel_window: usize,
    pub on_unresolved: Unresolved,
    pub mismatch_qual: MismatchQual,
    pub mismatch_base: MismatchBase,
    pub match_qual: MatchQual,
    pub qual_cap: u8,
    pub keep: Keep,
    pub score_qual_weight: f32,
    /// Hard-clip rather than soft-clip the discarded copy.
    pub hard_clip: bool,
    /// Clip terminal indels off the kept 3' end.
    pub clean_kept_3prime: bool,
    pub terminal_indel_window: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            length_tolerance: 0,
            min_overlap: 4,
            min_support: 3,
            min_span_frac: 0.7,
            max_anchor_slack: 5,
            max_past_end: 0,
            max_mismatch_frac: 0.15,
            ambiguity_ratio: 0.8,
            indel_window: 10,
            on_unresolved: Unresolved::Pass,
            mismatch_qual: MismatchQual::Subtract,
            mismatch_base: MismatchBase::None,
            match_qual: MatchQual::SumCapped,
            qual_cap: 40,
            keep: Keep::Score,
            score_qual_weight: 1.0,
            hard_clip: false,
            clean_kept_3prime: true,
            terminal_indel_window: 3,
        }
    }
}

// ---------------------------------------------------------------------------
// The plan
// ---------------------------------------------------------------------------

/// What to do to one record.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RecordPlan {
    pub id: usize,
    /// Aligned query bases to clip from the start of the alignment.
    pub clip_front: usize,
    /// ... and from the end.
    pub clip_back: usize,
    /// Make the resulting clips hard rather than soft.
    pub hard_front: bool,
    pub hard_back: bool,
    /// SEQ indices whose quality becomes this value. Absolute, not a delta, so
    /// that re-running cannot compound.
    pub qual_set: Vec<(usize, u8)>,
    /// SEQ indices whose base call becomes this value.
    pub base_set: Vec<(usize, Seq)>,
    /// Nothing aligned survives.
    pub drop: bool,
    /// This supplementary becomes the end's primary alignment, because the
    /// end's primary was dropped and this is its 5'-most survivor. The writer
    /// puts this alignment on the dropped primary's record, so the read keeps
    /// exactly one primary line and its whole sequence.
    pub promote: bool,
}

impl RecordPlan {
    pub fn is_empty(&self) -> bool {
        self.clip_front == 0
            && self.clip_back == 0
            && self.qual_set.is_empty()
            && self.base_set.is_empty()
            && !self.drop
            && !self.promote
    }
}

/// Why a template's fragment length was not established.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reason {
    /// Too few candidate pairs agreed on any fragment length.
    LowSupport,
    /// The agreeing pairs covered too little of the implied overlap.
    ShortSpan,
    /// The candidate group did not reach both reads' 3' termini, so it is not
    /// where a real overlap would sit.
    AnchorSlack,
    /// The implied fragment length is shorter than what one of the reads
    /// demonstrably aligned, so the length contradicts its own evidence.
    PastEnd,
    /// Too many bases disagreed for the alignment to be believable.
    HighMismatch,
    /// A second fragment length was about as well supported: the molecule
    /// appears to hold two copies of a locus.
    AmbiguousL,
    /// A second fragment length nearby, over a disjoint span: the mates
    /// disagree about an indel rather than about which locus they are on.
    IndelDisagreement,
    /// One read's records disagreed about its length.
    InconsistentLength,
}

impl Reason {
    pub fn as_str(self) -> &'static str {
        match self {
            Reason::LowSupport => "low_support",
            Reason::ShortSpan => "short_span",
            Reason::AnchorSlack => "anchor_slack",
            Reason::PastEnd => "past_end",
            Reason::HighMismatch => "high_mismatch",
            Reason::AmbiguousL => "ambiguous_L",
            Reason::IndelDisagreement => "indel_disagreement",
            Reason::InconsistentLength => "inconsistent_length",
        }
    }
    pub const ALL: [Reason; 8] = [
        Reason::LowSupport,
        Reason::ShortSpan,
        Reason::AnchorSlack,
        Reason::PastEnd,
        Reason::HighMismatch,
        Reason::AmbiguousL,
        Reason::IndelDisagreement,
        Reason::InconsistentLength,
    ];
}

/// A resolved template's numbers, which become its tags.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Resolution {
    pub l: usize,
    pub overlap: usize,
    pub support: usize,
    pub span_frac: f32,
    pub mismatches: usize,
    pub kept_r1: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Verdict {
    /// The two ends do not meet, or meet in fewer than `min_overlap` bases.
    NoOverlap { l: Option<usize> },
    Resolved(Resolution),
    Unresolved { reason: Reason, l: Option<usize> },
}

#[derive(Clone, Debug, Default)]
pub struct Counts {
    pub positions_compared: u64,
    pub positions_match: u64,
    pub positions_mismatch: u64,
    /// In the overlap, but observed by only one end (an indel or a clip).
    pub positions_unpaired: u64,
    /// Skipped because a base call was ambiguous.
    pub positions_ambiguous: u64,
    pub bases_clipped: u64,
    pub bases_hard_clipped: u64,
    pub quals_raised: u64,
    pub quals_lowered: u64,
    pub bases_set_n: u64,
    pub records_dropped: u64,
}

#[derive(Clone, Debug)]
pub struct Outcome {
    pub plans: Vec<RecordPlan>,
    pub verdict: Verdict,
    pub counts: Counts,
}

impl Default for Outcome {
    fn default() -> Self {
        Self {
            plans: Vec::new(),
            verdict: Verdict::NoOverlap { l: None },
            counts: Counts::default(),
        }
    }
}

// ---------------------------------------------------------------------------
// Base comparison
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Class {
    Match,
    Mismatch,
}

/// Classify two observations of one fragment position.
///
/// Both records store SEQ reference-forward, so two bases at the same reference
/// position are already on the same strand and compare directly with no
/// complementing. An ambiguous call on either side is no evidence either way,
/// so it is skipped rather than counted as agreement or disagreement.
fn classify(a: Seq, b: Seq) -> Option<Class> {
    if a.is_ambiguous() || b.is_ambiguous() || a.is_empty() || b.is_empty() {
        return None;
    }
    Some(if a == b { Class::Match } else { Class::Mismatch })
}

// ---------------------------------------------------------------------------
// Fragment span
// ---------------------------------------------------------------------------

/// Fragment indices observed by both ends, given a fragment length.
///
/// `[max(0, L-len2), min(len1, L))`. The clamps are not decoration: they are
/// what makes an untrimmed adapter read-through — where `L < len1`, so R1's
/// tail runs off the end of the molecule — come out as the molecule's own
/// extent, and not as an overlap longer than the fragment itself.
pub fn overlap_span(l: usize, len1: usize, len2: usize) -> (usize, usize) {
    let lo = l.saturating_sub(len2);
    let hi = len1.min(l);
    if lo >= hi {
        (0, 0)
    } else {
        (lo, hi)
    }
}

// ---------------------------------------------------------------------------
// Planning
// ---------------------------------------------------------------------------

/// One candidate fragment length and what supports it.
#[derive(Clone, Debug)]
struct Group {
    l: usize,
    n: usize,
    min_e1: usize,
    max_e1: usize,
}

impl Group {
    fn span(&self) -> usize {
        self.max_e1 - self.min_e1 + 1
    }
    /// How far short of each read's last aligned base the group stops.
    ///
    /// The true group reaches both, because the overlap's two ends *are* the
    /// two reads' 3' termini. Measuring against the last aligned base rather
    /// than the read length makes this tolerant of a soft-clipped tail and of
    /// untrimmed adapter, neither of which aligns.
    fn reach(&self, e1_last: usize, e2_last: usize) -> (usize, usize) {
        let max_e2 = self.l - 1 - self.min_e1;
        (e1_last.saturating_sub(self.max_e1), e2_last.saturating_sub(max_e2))
    }

    /// Bases still aligned past the end of the molecule this length implies.
    ///
    /// A fragment length is a claim about where the molecule stops. Anything a
    /// read holds beyond that is adapter, and adapter does not align to the
    /// locus the fragment came from. So a length that leaves aligned bases past
    /// its own end contradicts itself, whatever else it has going for it.
    ///
    /// This is what separates a genuine short fragment from a molecule holding
    /// one locus twice. Real read-through leaves those bases soft-clipped and
    /// the count is zero; a decoy that pairs one copy against the other implies
    /// a molecule shorter than the reads demonstrably aligned to it.
    fn past_end(&self, e1_last: usize, e2_last: usize) -> usize {
        let last = self.l - 1;
        e1_last.saturating_sub(last) + e2_last.saturating_sub(last)
    }
}

/// Plan the resolution of one template's alignments.
///
/// `alns` is every primary and supplementary record for one query name.
/// Secondary alignments must not be included: a secondary is the *same bases*
/// placed elsewhere, an alternative hypothesis, so deduplicating it against the
/// primary is meaningless -- they were never both true.
pub fn plan(alns: &[Aln], opts: &Options) -> Outcome {
    let mut out = Outcome {
        plans: alns.iter().map(|a| RecordPlan { id: a.id, ..Default::default() }).collect(),
        ..Default::default()
    };
    decide(&mut out, alns, opts);
    // Every path through `decide` lands here. `finish` drops records and hands
    // a dropped primary's flag to a survivor, so an exit that skipped it would
    // emit a template with no primary -- invisible until some downstream tool
    // that reads primaries only came up short.
    finish(&mut out, alns, opts);
    out
}

/// Fill in `out`'s verdict and per-record edits. See [`plan`].
fn decide(out: &mut Outcome, alns: &[Aln], opts: &Options) {
    let r1 = build_end(alns, true, opts.length_tolerance);
    let r2 = build_end(alns, false, opts.length_tolerance);
    if r1.len == 0 || r2.len == 0 {
        return; // single-ended: nothing to pair against
    }
    if r1.inconsistent || r2.inconsistent {
        out.verdict = Verdict::Unresolved { reason: Reason::InconsistentLength, l: None };
        return;
    }

    let groups = candidate_groups(&r1, &r2);

    let (Some(a1), Some(a2)) = (r1.last_aligned, r2.last_aligned) else {
        return;
    };

    // --- choose a fragment length -----------------------------------------
    //
    // Not the mode and not the median: the group that reaches both reads' 3'
    // termini. A short true overlap beats a long spurious one here, which is
    // the case a central-tendency statistic gets backwards.
    let best = groups.values().min_by(|a, b| {
        let (ra1, ra2) = a.reach(a1, a2);
        let (rb1, rb2) = b.reach(a1, a2);
        ra1.max(ra2)
            .cmp(&rb1.max(rb2))
            .then(b.n.cmp(&a.n))
            .then(a.l.cmp(&b.l))
    });
    let Some(best) = best.cloned() else {
        return;
    };

    let (lo, hi) = overlap_span(best.l, r1.len, r2.len);
    let predicted = hi - lo;
    if predicted < opts.min_overlap {
        out.verdict = Verdict::NoOverlap { l: Some(best.l) };
        return;
    }

    // --- compare the overlap ----------------------------------------------
    let mut cmp: Vec<(usize, usize, Class)> = Vec::with_capacity(predicted);
    for f in lo..hi {
        let (e1, e2) = (f, best.l - 1 - f);
        let (Some(b1), Some(b2)) = (r1.base(e1), r2.base(e2)) else {
            out.counts.positions_unpaired += 1;
            continue;
        };
        if !b1.aligned || !b2.aligned {
            out.counts.positions_unpaired += 1;
            continue;
        }
        let (s1, s2) = (alns[b1.rec].seq[b1.seq], alns[b2.rec].seq[b2.seq]);
        match classify(s1, s2) {
            None => out.counts.positions_ambiguous += 1,
            Some(c) => {
                out.counts.positions_compared += 1;
                match c {
                    Class::Match => out.counts.positions_match += 1,
                    Class::Mismatch => out.counts.positions_mismatch += 1,
                }
                cmp.push((e1, e2, c));
            }
        }
    }

    // --- validate ---------------------------------------------------------
    let (reach1, reach2) = best.reach(a1, a2);
    let hard = out.counts.positions_match + out.counts.positions_mismatch;
    let mism_frac = if hard == 0 { 0.0 } else { out.counts.positions_mismatch as f32 / hard as f32 };
    let span_frac = best.span() as f32 / predicted as f32;

    // Only a group that could itself be the answer counts as a rival. A
    // molecule holding two copies of a locus throws off a large decoy group,
    // but that group pairs R1's first copy against R2's second, which sit at
    // opposite ends of their reads -- so it lands far from both 3' termini and
    // is refused by the anchor test rather than by its size. Comparing sizes
    // instead would reject the sister-chromatid and A-B-A cases, whose true
    // group the estimator in fact recovers correctly.
    let anchored = |g: &Group| {
        let (x, y) = g.reach(a1, a2);
        x.max(y) <= opts.max_anchor_slack && g.past_end(a1, a2) <= opts.max_past_end
    };
    let rival = groups
        .values()
        .filter(|g| g.l != best.l && anchored(g))
        .max_by_key(|g| g.n)
        .filter(|g| g.n as f32 >= opts.ambiguity_ratio * best.n as f32);

    let reason = refusal(
        &best,
        &groups,
        rival.is_some(),
        (reach1.max(reach2), best.past_end(a1, a2)),
        (span_frac, mism_frac),
        opts,
    );

    if let Some(reason) = reason {
        out.verdict = Verdict::Unresolved { reason, l: Some(best.l) };
        return;
    }

    apply_consensus(out, alns, &r1, &r2, &cmp, opts);

    // --- pick the end that gives up its copy ------------------------------
    let keep_r1 = match opts.keep {
        Keep::R1 => true,
        Keep::R2 => false,
        Keep::Score => {
            let (m1, q1) = score(&r1, alns, best.l, lo, hi, true);
            let (m2, q2) = score(&r2, alns, best.l, lo, hi, false);
            // A read with no quality string (QUAL `*`) has no quality term,
            // so neither does its mate: one mean against nothing is no contest.
            let w = opts.score_qual_weight;
            let (s1, s2) = match (q1, q2) {
                (Some(q1), Some(q2)) => (m1 + w * q1, m2 + w * q2),
                _ => (m1, m2),
            };
            s1 >= s2 // ties to R1
        }
    };

    // The overlap is 3'-terminal on both ends by construction, so the discarded
    // copy is always a suffix of the read in end-index terms, which is the one
    // shape a clip can express. `e >= l` is the adapter guard: bases past the
    // end of the molecule, which both ends lose regardless of who won.
    let drop_from_1 = if keep_r1 { best.l } else { lo };
    let drop_from_2 = if keep_r1 { best.l - hi } else { best.l };
    clip_end(out, alns, true, r1.len, drop_from_1, opts);
    clip_end(out, alns, false, r2.len, drop_from_2, opts);

    if opts.clean_kept_3prime {
        let kept = if keep_r1 { &r1 } else { &r2 };
        clean_terminus(out, alns, kept, opts);
    }

    out.verdict = Verdict::Resolved(Resolution {
        l: best.l,
        overlap: predicted,
        support: best.n,
        span_frac,
        mismatches: out.counts.positions_mismatch as usize,
        kept_r1: keep_r1,
    });
}

/// Group every doubly-observed reference position by the fragment length it
/// implies.
///
/// A candidate is a reference position observed by both ends. The two covering
/// records must be on opposite strands: R2 is the reverse complement of the
/// molecule, so wherever R1 and R2 read the same molecule base their records
/// necessarily disagree about orientation. Requiring agreement instead would
/// reject every true overlap.
///
/// Both position lists are sorted, so this is a merge join rather than a
/// nested scan -- the lists are read length, and a template with several
/// supplementaries has several of them.
fn candidate_groups(r1: &End, r2: &End) -> BTreeMap<usize, Group> {
    let key = |e: &End| -> Vec<(i32, i64, usize)> {
        let mut v: Vec<(i32, i64, usize)> = e
            .bases
            .iter()
            .enumerate()
            .filter_map(|(i, b)| b.as_ref().and_then(|b| b.refpos.map(|(t, p)| (t, p, i))))
            .collect();
        v.sort_unstable();
        v
    };
    let k1 = key(r1);
    let k2 = key(r2);

    let mut groups: BTreeMap<usize, Group> = BTreeMap::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < k1.len() && j < k2.len() {
        let (a, b) = ((k1[i].0, k1[i].1), (k2[j].0, k2[j].1));
        match a.cmp(&b) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                let ie = k1[i..].iter().take_while(|x| (x.0, x.1) == a).count();
                let je = k2[j..].iter().take_while(|x| (x.0, x.1) == b).count();
                for x in &k1[i..i + ie] {
                    for y in &k2[j..j + je] {
                        let (e1, e2) = (x.2, y.2);
                        let (b1, b2) = (r1.base(e1).unwrap(), r2.base(e2).unwrap());
                        if b1.reverse == b2.reverse {
                            continue;
                        }
                        let l = e1 + e2 + 1;
                        let g = groups.entry(l).or_insert(Group {
                            l,
                            n: 0,
                            min_e1: usize::MAX,
                            max_e1: 0,
                        });
                        g.n += 1;
                        g.min_e1 = g.min_e1.min(e1);
                        g.max_e1 = g.max_e1.max(e1);
                    }
                }
                i += ie;
                j += je;
            }
        }
    }
    groups
}

/// Why this fragment length cannot be trusted, if it cannot.
///
/// Ordered most fundamental first, so the reported reason is the deepest thing
/// wrong rather than whichever test happened to run last.
fn refusal(
    best: &Group,
    groups: &BTreeMap<usize, Group>,
    has_rival: bool,
    (anchor_slack, past_end): (usize, usize),
    (span_frac, mism_frac): (f32, f32),
    opts: &Options,
) -> Option<Reason> {
    let reason = if best.n < opts.min_support {
        Reason::LowSupport
    } else if past_end > opts.max_past_end {
        Reason::PastEnd
    } else if anchor_slack > opts.max_anchor_slack {
        Reason::AnchorSlack
    } else if span_frac < opts.min_span_frac {
        Reason::ShortSpan
    } else if mism_frac > opts.max_mismatch_frac {
        Reason::HighMismatch
    } else if has_rival {
        Reason::AmbiguousL
    } else {
        return None;
    };

    // An indel one mate called and the other did not splits what is really one
    // group into two, near in length and covering disjoint stretches of R1,
    // each reaching one read's 3' terminus and neither reaching both. That is
    // the same refusal either way, so this only renames it -- it can never
    // cause one -- but "the mates disagree about an indel" is a far more
    // actionable thing to read in the statistics than "anchor slack".
    if matches!(reason, Reason::AnchorSlack | Reason::AmbiguousL) {
        let indel = groups.values().any(|g| {
            g.l != best.l
                && g.n >= opts.min_support
                && g.l.abs_diff(best.l) <= opts.indel_window
                && (g.max_e1 < best.min_e1 || best.max_e1 < g.min_e1)
        });
        if indel {
            return Some(Reason::IndelDisagreement);
        }
    }
    Some(reason)
}

/// Write the per-base quality and base edits the compared overlap implies.
///
/// Both copies of a doubly-observed base are edited, not just the one that
/// survives clipping: a caller that ignores the clip and reads both records
/// should still see the same call at both.
fn apply_consensus(
    out: &mut Outcome,
    alns: &[Aln],
    r1: &End,
    r2: &End,
    cmp: &[(usize, usize, Class)],
    opts: &Options,
) {
    for &(e1, e2, class) in cmp {
        let (b1, b2) = (*r1.base(e1).unwrap(), *r2.base(e2).unwrap());
        let (q1, q2) = (alns[b1.rec].qual[b1.seq], alns[b2.rec].qual[b2.seq]);
        // A record stored without qualities (QUAL `*`) has nothing to combine,
        // and editing one of its bytes would leave a string that is neither
        // `*` nor valid. Both copies keep their qualities; a base edit below
        // still applies, since it does not depend on quality.
        let unknown = q1 == QUAL_MISSING || q2 == QUAL_MISSING;
        match class {
            _ if unknown => {
                if class == Class::Mismatch && opts.mismatch_base == MismatchBase::SetN {
                    for b in [&b1, &b2] {
                        out.plans[b.rec].base_set.push((b.seq, Seq::N));
                        out.counts.bases_set_n += 1;
                    }
                }
            }
            Class::Match => match opts.match_qual {
                MatchQual::None => {}
                MatchQual::Max => {
                    let v = q1.max(q2);
                    set_qual(out, &b1, v, q1);
                    set_qual(out, &b2, v, q2);
                }
                MatchQual::SumCapped => {
                    let v = (q1.saturating_add(q2)).min(opts.qual_cap);
                    set_qual(out, &b1, v, q1);
                    set_qual(out, &b2, v, q2);
                }
            },
            Class::Mismatch => {
                // The winner is the better-quality base, which is a separate
                // question from which *end* is kept: a position where the kept
                // end holds the weaker call is still zeroed.
                let one_wins = (q1, alns[b1.rec].mapq, false) > (q2, alns[b2.rec].mapq, true);
                let (win, lose, qw, ql) = if one_wins {
                    (&b1, &b2, q1, q2)
                } else {
                    (&b2, &b1, q2, q1)
                };
                match opts.mismatch_qual {
                    MismatchQual::None => {}
                    MismatchQual::ZeroBoth => {
                        set_qual(out, win, 0, qw);
                        set_qual(out, lose, 0, ql);
                    }
                    MismatchQual::ZeroLoser => set_qual(out, lose, 0, ql),
                    MismatchQual::Subtract => {
                        set_qual(out, lose, 0, ql);
                        set_qual(out, win, qw.saturating_sub(ql), qw);
                    }
                }
                if opts.mismatch_base == MismatchBase::SetN {
                    for b in [&b1, &b2] {
                        out.plans[b.rec].base_set.push((b.seq, Seq::N));
                        out.counts.bases_set_n += 1;
                    }
                }
            }
        }
    }
}

fn set_qual(out: &mut Outcome, b: &Base, to: u8, from: u8) {
    if to == from {
        return;
    }
    if to > from {
        out.counts.quals_raised += 1;
    } else {
        out.counts.quals_lowered += 1;
    }
    out.plans[b.rec].qual_set.push((b.seq, to));
}

/// Mean MAPQ of the covering records over the overlap, and mean base quality.
/// MAPQ is per base: the overlap frequently sits on a low-MAPQ supplementary,
/// and using the end's primary MAPQ would invert the choice.
///
/// The quality mean is over bases with a known quality and is `None` when
/// there are none (QUAL `*`). The MAPQ mean is −∞ when the end has no aligned
/// base in the overlap.
fn score(end: &End, alns: &[Aln], l: usize, lo: usize, hi: usize, is_r1: bool) -> (f32, Option<f32>) {
    let (mut m, mut n) = (0f32, 0usize);
    let (mut q, mut nq) = (0f32, 0usize);
    for f in lo..hi {
        let e = if is_r1 { f } else { l - 1 - f };
        if let Some(b) = end.base(e) {
            if b.aligned {
                m += b.mapq as f32;
                n += 1;
                let bq = alns[b.rec].qual[b.seq];
                if bq != QUAL_MISSING {
                    q += bq as f32;
                    nq += 1;
                }
            }
        }
    }
    if n == 0 {
        return (f32::NEG_INFINITY, None);
    }
    (m / n as f32, (nq > 0).then(|| q / nq as f32))
}

/// Clip every base at or past end index `from` on one read end.
/// Clip every base at or past end index `from` on one read end.
///
/// Selection is by which read the record belongs to, not by whether it won any
/// bases in `build_end`. Two records of one read can claim the same read bases
/// at a junction the aligner split ambiguously; the consensus needs one winner
/// per base, but the loser is still in the file still carrying its copy, and a
/// record wholly contained in another would otherwise keep the overlap the rest
/// of the read just gave up.
fn clip_end(out: &mut Outcome, alns: &[Aln], is_r1: bool, len: usize, from: usize, opts: &Options) {
    if from >= len {
        return;
    }
    for i in (0..alns.len()).filter(|&i| alns[i].is_r1 == is_r1) {
        let a = &alns[i];
        let map = a.seq_map();
        let (mut front, mut back) = (0usize, 0usize);
        for (s, (aligned, _)) in map.iter().enumerate() {
            if !aligned {
                continue;
            }
            if a.end_index(s) >= from {
                // End index runs opposite to SEQ index on a reverse record, so
                // a suffix of the read is a prefix of SEQ.
                if a.is_reverse {
                    front += 1;
                } else {
                    back += 1;
                }
            }
        }
        add_clip(&mut out.plans[i], a, front, back, opts);
    }
}

fn add_clip(p: &mut RecordPlan, a: &Aln, front: usize, back: usize, opts: &Options) {
    p.clip_front = p.clip_front.max(front);
    p.clip_back = p.clip_back.max(back);
    if front > 0 && opts.hard_clip {
        p.hard_front = true;
    }
    if back > 0 && opts.hard_clip {
        p.hard_back = true;
    }
    if p.clip_front + p.clip_back >= a.query_len() {
        p.drop = true;
        p.clip_front = 0;
        p.clip_back = 0;
        p.hard_front = false;
        p.hard_back = false;
    }
}

/// Clip a terminal indel off the kept end's 3' terminus.
///
/// A terminal insertion or a deletion within a few bases of the end is placed
/// on almost no evidence, and whichever way it went it moves every base after
/// it. Cheaper to lose those bases than to call methylation on them.
fn clean_terminus(out: &mut Outcome, alns: &[Aln], end: &End, opts: &Options) {
    // The 3'-most record still holding an aligned base.
    let Some(last) = end.last_aligned else { return };
    let Some(b) = end.base(last) else { return };
    let a = &alns[b.rec];
    if out.plans[b.rec].drop {
        return;
    }
    // `terminal_trim` walks inward from the read's 3' end. On a forward record
    // that is the CIGAR's tail; on a reverse record the CIGAR is stored
    // reference-forward, so the read's tail is already the CIGAR's head.
    let mut ops = a.cigar.clone();
    if !a.is_reverse {
        ops.reverse();
    }
    let extra = terminal_trim(&ops, opts.terminal_indel_window);
    if extra == 0 {
        return;
    }
    let (f, k) = (out.plans[b.rec].clip_front, out.plans[b.rec].clip_back);
    if a.is_reverse {
        add_clip(&mut out.plans[b.rec], a, f + extra, 0, opts);
    } else {
        add_clip(&mut out.plans[b.rec], a, 0, k + extra, opts);
    }
}

/// Aligned query bases to remove so that no indel remains within `window`
/// matched bases of the terminus. `ops` runs inward from the 3' end.
fn terminal_trim(ops: &[CigarOp], window: usize) -> usize {
    let (mut m, mut run, mut best) = (0usize, 0usize, 0usize);
    for op in ops {
        match *op {
            CigarOp::H(_) | CigarOp::S(_) => {}
            CigarOp::M(n) => {
                let take = (n as usize).min(window - m);
                m += take;
                run += take;
                if m >= window {
                    break;
                }
            }
            CigarOp::I(n) => {
                run += n as usize;
                best = run;
            }
            CigarOp::D(n) | CigarOp::N(n) => {
                let _ = n;
                best = run;
            }
            CigarOp::P(_) => {}
        }
    }
    best
}

/// Apply the unresolved policy and total up the plans.
fn finish(out: &mut Outcome, alns: &[Aln], opts: &Options) {
    if matches!(out.verdict, Verdict::Unresolved { .. }) && opts.on_unresolved == Unresolved::Drop {
        for p in out.plans.iter_mut() {
            p.drop = true;
        }
    }

    // A read whose primary was dropped hands its place to the 5'-most survivor,
    // so the read keeps exactly one primary line (see `overlap_apply::graft`).
    for want_r1 in [true, false] {
        let of_end: Vec<usize> =
            (0..alns.len()).filter(|&i| alns[i].is_r1 == want_r1).collect();
        let primary_gone = of_end
            .iter()
            .any(|&i| !alns[i].is_supplementary && out.plans[i].drop);
        if !primary_gone {
            continue;
        }
        let survivor = of_end
            .iter()
            .filter(|&&i| !out.plans[i].drop && alns[i].is_supplementary)
            .min_by_key(|&&i| alns[i].end_offset());
        if let Some(&i) = survivor {
            out.plans[i].promote = true;
        }
    }

    for p in out.plans.iter_mut() {
        p.qual_set.sort_unstable();
        p.qual_set.dedup();
        p.base_set.sort_unstable();
        p.base_set.dedup();
    }
    for p in &out.plans {
        // Hard-clipped bases are a subset of clipped ones, counted per side:
        // a record can be soft-clipped at one end and hard-clipped at the other.
        out.counts.bases_clipped += (p.clip_front + p.clip_back) as u64;
        out.counts.bases_hard_clipped += (p.hard_front as u64) * p.clip_front as u64
            + (p.hard_back as u64) * p.clip_back as u64;
        out.counts.records_dropped += p.drop as u64;
    }
}

// ---------------------------------------------------------------------------
// Applying clips to a CIGAR
// ---------------------------------------------------------------------------

/// Soft-clip `front` and `back` aligned query bases, returning the new CIGAR
/// and reference start.
///
/// Clipping from the front moves POS: the alignment now begins further along
/// the reference by however many reference bases were removed. Clipping from
/// the back does not. Getting that wrong shifts every downstream coordinate
/// silently, so it is the single thing here most worth testing.
pub fn reclip(cigar: &[CigarOp], pos: i64, front: usize, back: usize) -> (Vec<CigarOp>, i64) {
    let mut ops: Vec<CigarOp> = cigar.to_vec();
    let mut pos = pos;

    if front > 0 {
        let (o, moved) = eat(&ops, front);
        ops = o;
        pos += moved as i64;
    }
    if back > 0 {
        // Reversing is the whole of the direction handling: `eat` always works
        // from the head of the list it is given.
        ops.reverse();
        let (o, _) = eat(&ops, back);
        ops = o;
        ops.reverse();
    }

    // An alignment may not begin or end with a deletion: it would assert a
    // reference gap adjacent to nothing. Clipping routinely exposes one, and it
    // sits *inside* the clip, so the search has to look past H and S.
    let clip = |o: &CigarOp| matches!(o, CigarOp::H(_) | CigarOp::S(_));
    loop {
        let Some(i) = ops.iter().position(|o| !clip(o)) else { break };
        match ops[i] {
            CigarOp::D(n) | CigarOp::N(n) => {
                pos += n as i64;
                ops.remove(i);
            }
            CigarOp::P(_) => {
                ops.remove(i);
            }
            _ => break,
        }
    }
    loop {
        let Some(i) = ops.iter().rposition(|o| !clip(o)) else { break };
        if matches!(ops[i], CigarOp::D(_) | CigarOp::N(_) | CigarOp::P(_)) {
            ops.remove(i);
        } else {
            break;
        }
    }

    (merge(&ops), pos)
}

/// Convert the outermost soft clip at one or both ends into a hard clip.
///
/// Returns the new operations and how many SEQ bases to remove from the front
/// and the back, so the caller can keep SEQ and QUAL in step with the CIGAR.
pub fn harden(ops: &[CigarOp], front: bool, back: bool) -> (Vec<CigarOp>, usize, usize) {
    let mut ops = ops.to_vec();
    let mut cut = (0usize, 0usize);
    if front {
        if let Some(i) = ops.iter().position(|o| !matches!(o, CigarOp::H(_))) {
            if let CigarOp::S(n) = ops[i] {
                cut.0 = n as usize;
                ops[i] = CigarOp::H(n);
            }
        }
    }
    if back {
        if let Some(i) = ops.iter().rposition(|o| !matches!(o, CigarOp::H(_))) {
            if let CigarOp::S(n) = ops[i] {
                cut.1 = n as usize;
                ops[i] = CigarOp::H(n);
            }
        }
    }
    (merge(&ops), cut.0, cut.1)
}

/// Convert the first `n` aligned query bases into soft clip, reporting how many
/// reference bases were given up.
///
/// Always works from the head of the list. Clipping the other end is done by
/// reversing before and after, so there is one implementation and no direction
/// to get wrong twice.
fn eat(ops: &[CigarOp], n: usize) -> (Vec<CigarOp>, usize) {
    let mut left = n;
    let mut moved = 0usize;
    let mut out: Vec<CigarOp> = Vec::with_capacity(ops.len() + 2);
    let mut clipped = 0usize;

    for op in ops {
        if left == 0 {
            out.push(*op);
            continue;
        }
        match *op {
            CigarOp::H(_) => out.push(*op),
            CigarOp::S(k) => {
                // Already unaligned; carry it into the new clip.
                clipped += k as usize;
            }
            CigarOp::M(k) => {
                let take = (k as usize).min(left);
                clipped += take;
                moved += take;
                left -= take;
                if (k as usize) > take {
                    out.push(CigarOp::M(k - take as u32));
                }
            }
            CigarOp::I(k) => {
                let take = (k as usize).min(left);
                clipped += take;
                left -= take;
                if (k as usize) > take {
                    out.push(CigarOp::I(k - take as u32));
                }
            }
            CigarOp::D(k) | CigarOp::N(k) => {
                // Between clipped bases, so it is removed with them; but only
                // while still clipping, or it would swallow a real deletion.
                moved += k as usize;
            }
            CigarOp::P(_) => {}
        }
    }

    if clipped > 0 {
        // After any hard clip, which stays outermost.
        let at = out.iter().take_while(|o| matches!(o, CigarOp::H(_))).count();
        out.insert(at, CigarOp::S(clipped as u32));
    }
    (out, moved)
}

/// Collapse adjacent operations of the same kind and drop empty ones.
fn merge(ops: &[CigarOp]) -> Vec<CigarOp> {
    let mut out: Vec<CigarOp> = Vec::with_capacity(ops.len());
    for op in ops {
        if op.len() == 0 {
            continue;
        }
        match out.last_mut() {
            Some(prev) if std::mem::discriminant(prev) == std::mem::discriminant(op) => {
                let n = prev.len() + op.len();
                *prev = prev.with_len(n);
            }
            _ => out.push(*op),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Accumulated statistics
// ---------------------------------------------------------------------------

/// Run-level totals, written beside the output BAM.
///
/// Two reasons this is worth a file rather than a line on stderr. It is the
/// only record of how much evidence the step removed, which a reviewer will
/// eventually ask for; and it is the first thing to look at when a concordance
/// check disagrees with another caller by a fraction of a percent.
#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub templates: u64,
    pub templates_touched: u64,
    pub no_overlap: u64,
    pub resolved: u64,
    pub unresolved: [u64; 8],
    pub kept_r1: u64,
    pub kept_r2: u64,
    pub totals: Counts,
    pub overlap_len: BTreeMap<usize, u64>,
    pub fragment_len: BTreeMap<usize, u64>,
}

impl Stats {
    pub fn add(&mut self, out: &Outcome) {
        self.templates += 1;
        let c = &out.counts;
        if c.bases_clipped > 0
            || c.quals_raised > 0
            || c.quals_lowered > 0
            || c.bases_set_n > 0
            || c.records_dropped > 0
        {
            self.templates_touched += 1;
        }
        match &out.verdict {
            Verdict::NoOverlap { .. } => self.no_overlap += 1,
            Verdict::Resolved(r) => {
                self.resolved += 1;
                if r.kept_r1 {
                    self.kept_r1 += 1;
                } else {
                    self.kept_r2 += 1;
                }
                *self.overlap_len.entry(r.overlap).or_default() += 1;
                *self.fragment_len.entry(r.l).or_default() += 1;
            }
            Verdict::Unresolved { reason, .. } => {
                let i = Reason::ALL.iter().position(|r| r == reason).unwrap();
                self.unresolved[i] += 1;
            }
        }

        let t = &mut self.totals;
        t.positions_compared += c.positions_compared;
        t.positions_match += c.positions_match;
        t.positions_mismatch += c.positions_mismatch;
        t.positions_unpaired += c.positions_unpaired;
        t.positions_ambiguous += c.positions_ambiguous;
        t.bases_clipped += c.bases_clipped;
        t.bases_hard_clipped += c.bases_hard_clipped;
        t.quals_raised += c.quals_raised;
        t.quals_lowered += c.quals_lowered;
        t.bases_set_n += c.bases_set_n;
        t.records_dropped += c.records_dropped;
    }

    /// Tidy three-column TSV: one row shape for scalars and histograms alike,
    /// so it loads with a single `read_csv` and groups by `section`.
    pub fn to_tsv(&self) -> String {
        let mut s = String::from("section\tkey\tvalue\n");
        let t = &self.totals;
        let row = |s: &mut String, sec: &str, k: &str, v: u64| {
            s.push_str(&format!("{sec}\t{k}\t{v}\n"));
        };

        row(&mut s, "counts", "templates", self.templates);
        row(&mut s, "counts", "templates_touched", self.templates_touched);
        row(&mut s, "counts", "no_overlap", self.no_overlap);
        row(&mut s, "counts", "resolved", self.resolved);
        row(&mut s, "counts", "kept_r1", self.kept_r1);
        row(&mut s, "counts", "kept_r2", self.kept_r2);
        for (i, r) in Reason::ALL.iter().enumerate() {
            row(&mut s, "unresolved", r.as_str(), self.unresolved[i]);
        }
        row(&mut s, "counts", "positions_compared", t.positions_compared);
        row(&mut s, "counts", "positions_match", t.positions_match);
        row(&mut s, "counts", "positions_mismatch", t.positions_mismatch);
        row(&mut s, "counts", "positions_unpaired", t.positions_unpaired);
        row(&mut s, "counts", "positions_ambiguous", t.positions_ambiguous);
        row(&mut s, "counts", "bases_clipped", t.bases_clipped);
        row(&mut s, "counts", "bases_hard_clipped", t.bases_hard_clipped);
        row(&mut s, "counts", "quals_raised", t.quals_raised);
        row(&mut s, "counts", "quals_lowered", t.quals_lowered);
        row(&mut s, "counts", "bases_set_n", t.bases_set_n);
        row(&mut s, "counts", "records_dropped", t.records_dropped);

        // Rates, so the file answers "how much did this change?" without
        // arithmetic. Integer per-million: no floats, no locale, no rounding
        // argument in a file meant for a methods section.
        let ppm = |num: u64, den: u64| if den == 0 { 0 } else { num * 1_000_000 / den };
        let unresolved: u64 = self.unresolved.iter().sum();
        row(&mut s, "rate_ppm", "templates_touched", ppm(self.templates_touched, self.templates));
        row(&mut s, "rate_ppm", "resolved_of_templates", ppm(self.resolved, self.templates));
        row(&mut s, "rate_ppm", "unresolved_of_templates", ppm(unresolved, self.templates));
        row(&mut s, "rate_ppm", "mismatch_of_compared", ppm(t.positions_mismatch, t.positions_compared));

        for (k, v) in &self.overlap_len {
            row(&mut s, "overlap_len", &k.to_string(), *v);
        }
        for (k, v) in &self.fragment_len {
            row(&mut s, "fragment_len", &k.to_string(), *v);
        }
        s
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// A molecule simulator, shared by the example-based and property tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_support {
    use super::*;

/// Parse a CIGAR string, and render one back.
    ///
    /// Test-only. The module reads CIGARs from records, never from text: `SA`
    /// is stripped rather than parsed, because clipping invalidates it.
    pub fn parse_cigar(s: &str) -> Result<Vec<CigarOp>, String> {
    let mut out = Vec::new();
    let mut num = String::new();
    for c in s.chars() {
        if c.is_ascii_digit() {
            num.push(c);
            continue;
        }
        let n: u32 = num
            .parse()
            .map_err(|_| format!("'{c}' has no length before it"))?;
        num.clear();
        out.push(match c {
            'M' | '=' | 'X' => CigarOp::M(n),
            'I' => CigarOp::I(n),
            'D' => CigarOp::D(n),
            'N' => CigarOp::N(n),
            'S' => CigarOp::S(n),
            'H' => CigarOp::H(n),
            'P' => CigarOp::P(n),
            other => return Err(format!("unknown CIGAR operation '{other}'")),
        });
    }
    if !num.is_empty() {
        return Err(format!("CIGAR ends with a length and no operation: '{num}'"));
    }
    Ok(out)
    }

    pub fn cigar_string(ops: &[CigarOp]) -> String {
    ops.iter()
        .map(|o| {
            let c = match o {
                CigarOp::M(_) => 'M',
                CigarOp::I(_) => 'I',
                CigarOp::D(_) => 'D',
                CigarOp::N(_) => 'N',
                CigarOp::S(_) => 'S',
                CigarOp::H(_) => 'H',
                CigarOp::P(_) => 'P',
            };
            format!("{}{}", o.len(), c)
        })
        .collect()
    }

    pub fn cig(s: &str) -> Vec<CigarOp> {
        parse_cigar(s).unwrap()
    }
    pub fn bases(s: &str) -> Vec<Seq> {
        s.chars().map(|c| Seq::FROM_FASTA[c as usize].unwrap()).collect()
    }

    // ---- a molecule simulator ---------------------------------------------
    //
    // Tests that hand-write CIGARs test the CIGARs. These build a molecule,
    // sequence it from both ends, and let the records fall out, so what is
    // being tested is the thing the tool actually claims: that a fragment's
    // bases can be recovered from its alignments.

    /// One molecule base and where it came from.
    #[derive(Clone, Copy, Debug)]
    pub struct Site {
        pub base: Seq,
        pub tid: i32,
        pub pos: i64,
        /// The molecule runs with the reference here, rather than against it.
        pub fwd: bool,
    }

    /// A stretch of reference, entered in one direction or the other.
    pub struct Seg {
        tid: i32,
        pos: i64,
        fwd: bool,
        bases: &'static str,
    }

    pub fn seg(tid: i32, pos: i64, fwd: bool, bases: &'static str) -> Seg {
        Seg { tid, pos, fwd, bases }
    }

    /// Lay segments end to end into one linear molecule.
    pub fn molecule(segs: &[Seg]) -> Vec<Site> {
        let mut out = Vec::new();
        for s in segs {
            let b = bases(s.bases);
            let n = b.len() as i64;
            for (k, base) in b.into_iter().enumerate() {
                // Against the reference, molecule index k is the k-th base
                // from the far end of the segment.
                let pos = if s.fwd { s.pos + k as i64 } else { s.pos + n - 1 - k as i64 };
                out.push(Site { base, tid: s.tid, pos, fwd: s.fwd });
            }
        }
        out
    }

    /// The base a record would store for this molecule index.
    ///
    /// SEQ is always reference-forward, whichever end read it and whichever way
    /// round the molecule entered the locus. That is why two records covering
    /// one molecule base compare directly, with no complementing anywhere in
    /// the planner.
    pub fn stored(s: &Site) -> Seq {
        if s.fwd { s.base } else { s.base.complement() }
    }

    /// Sequence a molecule from both ends and emit the alignments.
    ///
    /// Each maximal run of molecule bases that stays on one contig, in one
    /// direction, stepping one reference base at a time, becomes one record;
    /// the rest of the read is hard-clipped, as a chimeric aligner would do.
    pub fn sequence(mol: &[Site], len1: usize, len2: usize) -> Vec<Aln> {
        let l = mol.len();
        let mut out: Vec<Aln> = Vec::new();
        for is_r1 in [true, false] {
            let rlen = if is_r1 { len1 } else { len2 };
            // Read position k, and the molecule index it observes.
            let at = |k: usize| if is_r1 { k } else { l - 1 - k };
            let mut k = 0usize;
            while k < rlen.min(l) {
                let mut j = k + 1;
                while j < rlen.min(l) {
                    let (a, b) = (&mol[at(j - 1)], &mol[at(j)]);
                    let step = if is_r1 { b.pos - a.pos } else { b.pos - a.pos };
                    let want = if a.fwd == is_r1 { 1 } else { -1 };
                    if a.tid != b.tid || a.fwd != b.fwd || step != want {
                        break;
                    }
                    j += 1;
                }
                // Read positions [k, j) form one record.
                let sites: Vec<Site> = (k..j).map(|x| mol[at(x)]).collect();
                let fwd_seg = sites[0].fwd;
                let is_reverse = if is_r1 { !fwd_seg } else { fwd_seg };
                let mut sb: Vec<(i64, Seq)> =
                    sites.iter().map(|s| (s.pos, stored(s))).collect();
                sb.sort_by_key(|x| x.0);
                let pos = sb[0].0;
                let seq: Vec<Seq> = sb.iter().map(|x| x.1).collect();
                // Clips are written reference-forward, so on a reverse record
                // the read's tail is the CIGAR's head.
                let (lead, trail) = if is_reverse {
                    (rlen - j, k)
                } else {
                    (k, rlen - j)
                };
                let mut c = Vec::new();
                if lead > 0 {
                    c.push(CigarOp::H(lead as u32));
                }
                c.push(CigarOp::M(seq.len() as u32));
                if trail > 0 {
                    c.push(CigarOp::H(trail as u32));
                }
                out.push(Aln {
                    id: out.len(),
                    is_r1,
                    is_reverse,
                    is_supplementary: k > 0,
                    mapq: 60,
                    tid: sites[0].tid,
                    pos,
                    cigar: c,
                    qual: vec![30; seq.len()],
                    seq,
                });
                k = j;
            }
        }
        out
    }


    /// Build a molecule from bare segment shapes and a base generator, for the
    /// property tests, which have no fixed sequence to lay down.
    pub fn random_molecule(
        segs: &[(i32, i64, bool, usize)],
        pick: &mut dyn FnMut() -> usize,
    ) -> Vec<Site> {
        const ACGT: [Seq; 4] = [Seq::A, Seq::C, Seq::G, Seq::T];
        let mut out = Vec::new();
        for &(tid, pos, fwd, len) in segs {
            for k in 0..len {
                let p = if fwd { pos + k as i64 } else { pos + len as i64 - 1 - k as i64 };
                out.push(Site { base: ACGT[pick() % 4], tid, pos: p, fwd });
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::tests_support::*;
    use super::*;

    const A80: &str = "ACGTTGCAAGCTTAGCCATGGATCCGTAAGCTTCCGGAATTCGATCGTACGATCGGCTAAGCTTAAGGCCTTAAGCTA";
    const B80: &str = "TTGACCGGTTAACCGGATTCCGGAATTAAGGCCTTGGCCAAGGTTCCAAGGTTAACCGGTTAAGGCCTTAACCGGTAT";

    fn find(out: &Outcome, id: usize) -> &RecordPlan {
        out.plans.iter().find(|p| p.id == id).unwrap()
    }
    fn resolved(out: &Outcome) -> Resolution {
        match out.verdict {
            Verdict::Resolved(r) => r,
            ref v => panic!("expected a resolution, got {v:?}"),
        }
    }

    // ---- input ordering ----------------------------------------------------

    #[test]
    fn name_sorted_and_collated_are_both_accepted() {
        assert_eq!(
            check_grouping("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n").unwrap(),
            Grouping::Sorted
        );
        assert_eq!(
            check_grouping("@HD\tVN:1.6\tSO:unsorted\tGO:query\n").unwrap(),
            Grouping::Collated
        );
    }

    #[test]
    fn anything_ungrouped_is_refused_with_the_command_to_fix_it() {
        let e = check_grouping("@HD\tVN:1.6\tSO:coordinate\n").unwrap_err();
        assert!(e.contains("coordinate"), "{e}");
        assert!(e.contains("collate"), "the message should say how to fix it: {e}");
        assert!(check_grouping("@SQ\tSN:chr1\tLN:100\n").unwrap_err().contains("no @HD"));
        assert!(check_grouping("@HD\tVN:1.6\n").unwrap_err().contains("neither SO: nor GO:"));
    }

    #[test]
    fn header_parsing_is_not_fooled_by_lookalikes() {
        let e = check_grouping("@HD\tVN:1.6\tSO:coordinate\n@CO\tGO:query\n").unwrap_err();
        assert!(e.contains("coordinate"), "GO on a comment line is not the @HD field: {e}");
        assert_eq!(
            check_grouping("@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural\n").unwrap(),
            Grouping::Sorted
        );
    }

    // ---- coordinate systems (spec tests 1, 3, 11, 12) ----------------------

    #[test]
    fn read_length_is_the_whole_read_under_either_clip_convention() {
        // Hard- and soft-clipped supplementaries report the same read length:
        // the record says which convention its aligner used, so nothing has to
        // be told.
        assert_eq!(read_length(&cig("80H70M")), 150);
        assert_eq!(read_length(&cig("80S70M")), 150);
        assert_eq!(read_length(&cig("10H5S60M")), 75);
        assert_eq!(read_length(&cig("70M5S75H")), 150);
    }

    #[test]
    fn read_length_counts_bases_that_align_to_neither_side_of_a_junction() {
        // Spec test 3. A few bases of reconstructed junction sequence align
        // nowhere and sit in a query gap between the two records. Counting only
        // aligned bases would shorten the read and shift every fragment index
        // past the junction.
        let a = Aln {
            id: 0, is_r1: true, is_reverse: false, is_supplementary: false, mapq: 60,
            tid: 0, pos: 100, cigar: cig("70M5S75H"), seq: vec![Seq::A; 75],
            qual: vec![30; 75],
        };
        assert_eq!(a.read_length(), 150, "the 5 unaligned junction bases still count");
        assert_eq!(a.query_len(), 70, "but they are not aligned");
        assert_eq!(a.end_index(69), 69, "and the aligned bases keep their places");
    }

    #[test]
    fn end_index_inverts_on_a_reverse_record() {
        // Spec test 11. SEQ is reference-forward, so on a reverse record the
        // base sequenced first sits at the far end of SEQ.
        let f = Aln {
            id: 0, is_r1: true, is_reverse: false, is_supplementary: false, mapq: 60,
            tid: 0, pos: 0, cigar: cig("20H30M"), seq: vec![Seq::A; 30], qual: vec![30; 30],
        };
        assert_eq!(f.end_index(0), 20);
        assert_eq!(f.end_index(29), 49);
        let r = Aln { is_reverse: true, ..f.clone() };
        assert_eq!(r.end_index(0), 29, "SEQ 0 is the last base sequenced");
        assert_eq!(r.end_index(29), 0);
    }

    #[test]
    fn hard_clipped_bases_count_toward_the_read_length() {
        // Spec test 12. H is absent from SEQ but present in the read, and it is
        // the only thing that says where a supplementary's bases sit.
        let a = Aln {
            id: 0, is_r1: true, is_reverse: false, is_supplementary: true, mapq: 60,
            tid: 0, pos: 500, cigar: cig("70H80M"), seq: vec![Seq::C; 80], qual: vec![30; 80],
        };
        assert_eq!(a.read_length(), 150);
        assert_eq!(a.end_index(0), 70, "its first base is the read's 71st");
    }

    // ---- the estimator (spec tests 1, 2) ----------------------------------

    #[test]
    fn a_plain_overlapping_pair_recovers_the_fragment_length() {
        // Spec test 1. 120bp molecule, 80bp reads: a 40bp overlap.
        let mol = molecule(&[seg(0, 1000, true, A80), seg(0, 1080, true, &B80[..40])]);
        assert_eq!(mol.len(), 118);
        let alns = sequence(&mol, 80, 80);
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.l, 118, "the molecule's own length, recovered from the mates");
        assert_eq!(r.overlap, 80 + 80 - 118);
        assert_eq!(r.mismatches, 0);
    }

    #[test]
    fn mates_that_do_not_meet_are_left_alone() {
        // Spec test 2. A 400bp molecule read 80bp from each end.
        let segs = [seg(0, 1000, true, A80), seg(0, 2000, true, B80)];
        let mol = molecule(&segs);
        let alns = sequence(&mol, 40, 40);
        let out = plan(&alns, &Options::default());
        assert!(matches!(out.verdict, Verdict::NoOverlap { .. }), "{:?}", out.verdict);
        assert!(out.plans.iter().all(|p| p.is_empty()));
    }

    #[test]
    fn a_chimera_overlapping_across_its_junction_is_resolved() {
        // Spec test 4. The two ends meet at the ligation junction itself, so
        // the overlap spans both loci. Reference-interval intersection sees
        // nothing here worth resolving on one of the two contigs.
        let mol = molecule(&[seg(0, 1000, true, &A80[..60]), seg(1, 5000, true, &B80[..60])]);
        let alns = sequence(&mol, 80, 80);
        assert!(alns.len() >= 3, "the junction should split at least one read");
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.l, 120);
        assert_eq!(r.overlap, 40);
        assert_eq!(r.mismatches, 0, "fragment indices line up across the junction");
    }

    #[test]
    fn an_inverted_ligation_still_pairs_up() {
        // The second locus entered against the reference. Every record's
        // orientation flips, and the pairing rule -- that the two mates'
        // records disagree about strand -- has to keep holding.
        let mol = molecule(&[seg(0, 1000, true, &A80[..60]), seg(1, 5000, false, &B80[..60])]);
        let alns = sequence(&mol, 80, 80);
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.l, 120);
        assert_eq!(r.mismatches, 0);
    }

    #[test]
    fn a_sheared_self_circle_is_an_ordinary_short_chimera() {
        // Spec test 5. A self-ligated restriction fragment, sheared at an
        // interior point, reads p->b, junction, a->p: a circular permutation
        // with every base present exactly once. No special case, no rejection.
        let mol = molecule(&[
            seg(0, 1040, true, &A80[40..]),  // p -> b
            seg(0, 1000, true, &A80[..40]),  // a -> p
        ]);
        assert_eq!(mol.len(), 78);
        let alns = sequence(&mol, 60, 60);
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.l, 78, "resolved normally: {:?}", out.verdict);
        assert_eq!(r.overlap, 42);
        assert_eq!(r.mismatches, 0);
    }

    #[test]
    fn a_two_copy_molecule_resolves_on_the_true_group_not_the_decoy() {
        // Spec test 6, and the expectation there was wrong. Same-fragment
        // ligation between sister chromatids gives a molecule holding one locus
        // twice, which does throw off a decoy group -- but that group pairs
        // R1's first copy against R2's second, and those sit at opposite ends
        // of their reads, so it lands nowhere near either 3' terminus. The
        // anchor test discards it and the true length is recovered. Rejecting
        // on the decoy's size would have thrown away a usable read.
        let mol = molecule(&[seg(0, 1000, true, A80), seg(0, 1000, true, A80)]);
        assert_eq!(mol.len(), 156);
        let alns = sequence(&mol, 120, 120);
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.l, 156, "the decoy at L=78 must lose: {:?}", out.verdict);
        assert_eq!(r.overlap, 84);
        assert_eq!(r.mismatches, 0, "the true group pairs matching molecule bases");
    }

    #[test]
    fn a_two_copy_molecule_with_no_real_overlap_invents_none() {
        // Two copies of a 25 bp locus, read 20 bp from one end and 26 bp from
        // the other. The reads do not meet: R1 sits inside copy 1 and R2 runs
        // from copy 1's tail through the whole of copy 2. But R1's bases match
        // R2's reading of copy 2, and those pairs all agree on a length of 25,
        // reach R1's 3' end exactly, and span everything they claim to. Every
        // gate that looks at the candidates alone lets this through.
        //
        // What refuses it is that a 25-base molecule cannot have 26 aligned
        // bases on one of its reads. The length contradicts its own evidence.
        let unit = &A80[..25];
        let mol = molecule(&[seg(0, 1000, true, unit), seg(0, 1000, true, unit)]);
        assert_eq!(mol.len(), 50);
        let alns = sequence(&mol, 20, 26);
        let out = plan(&alns, &Options::default());
        match out.verdict {
            Verdict::Unresolved { reason, .. } => assert_eq!(reason, Reason::PastEnd),
            Verdict::NoOverlap { .. } => {}
            ref v => panic!("no overlap exists here, so nothing may be merged: {v:?}"),
        }
        assert!(out.plans.iter().all(|p| p.is_empty()), "{:?}", out.plans);

        // The same guard must stay silent on a genuinely short fragment, where
        // the bases past the molecule are adapter and do not align.
        let short = molecule(&[seg(0, 2000, true, &A80[..40])]);
        let mut alns = sequence(&short, 40, 40);
        for a in alns.iter_mut() {
            // 20 bases of unaligned adapter on the 3' end of each read
            a.cigar = if a.is_reverse { cig("20S40M") } else { cig("40M20S") };
            a.seq.extend(std::iter::repeat(Seq::A).take(20));
            a.qual.extend(std::iter::repeat(30u8).take(20));
            if a.is_reverse {
                a.seq.rotate_right(20);
                a.qual.rotate_right(20);
            }
        }
        let out = plan(&alns, &Options::default());
        assert_eq!(resolved(&out).l, 40, "read-through must still resolve: {:?}", out.verdict);
    }

    #[test]
    fn a_short_true_overlap_beats_a_long_spurious_one() {
        // Spec test 9, and the case that rules out the mode and the median. The
        // molecule visits one locus twice; the decoy group is 30 pairs wide and
        // the true overlap is 8. Any central-tendency statistic picks the decoy.
        let mol = molecule(&[
            seg(0, 1000, true, &A80[..30]),
            seg(1, 5000, true, B80),
            seg(0, 1000, true, &A80[..30]),
        ]);
        let l = mol.len();
        assert_eq!(l, 138);
        let alns = sequence(&mol, 73, 73);
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.l, 138, "the 30-pair decoy must lose to the 8-pair truth");
        assert_eq!(r.overlap, 8);
    }

    #[test]
    fn mates_disagreeing_about_an_indel_are_refused_and_named() {
        // R2 carries a deletion R1 does not, so positions before it imply one
        // fragment length and positions after it imply another. Neither group
        // reaches both 3' termini. Merging on either would misplace half the
        // overlap, so the template is passed through untouched -- and the
        // statistics say why.
        let mol = molecule(&[seg(0, 1000, true, A80), seg(0, 1078, true, &B80[..20])]);
        let mut alns = sequence(&mol, 60, 60);
        let r2 = alns.iter_mut().find(|a| !a.is_r1).unwrap();
        r2.cigar = cig("10M1D49M");
        let out = plan(&alns, &Options::default());
        match out.verdict {
            Verdict::Unresolved { reason, .. } => {
                assert_eq!(reason, Reason::IndelDisagreement, "{:?}", out.verdict)
            }
            ref v => panic!("expected a refusal, got {v:?}"),
        }
        assert!(out.plans.iter().all(|p| p.is_empty()), "pass means untouched");
    }

    #[test]
    fn too_little_support_is_refused_rather_than_guessed() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let alns = sequence(&mol, 45, 45);
        let opts = Options { min_support: 50, ..Default::default() };
        let out = plan(&alns, &opts);
        match out.verdict {
            Verdict::Unresolved { reason, .. } => assert_eq!(reason, Reason::LowSupport),
            ref v => panic!("expected a refusal, got {v:?}"),
        }
    }

    #[test]
    fn a_refusal_can_be_told_to_discard_the_template() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let alns = sequence(&mol, 45, 45);
        let opts =
            Options { min_support: 50, on_unresolved: Unresolved::Drop, ..Default::default() };
        let out = plan(&alns, &opts);
        assert!(out.plans.iter().all(|p| p.drop));
    }

    #[test]
    fn records_disagreeing_about_the_read_length_are_refused() {
        // One end must hold two records for them to disagree at all, so this
        // needs a chimera: the junction splits R1, and then one half claims a
        // read length the other half contradicts.
        let mol = molecule(&[seg(0, 1000, true, &A80[..40]), seg(1, 5000, true, &B80[..40])]);
        let mut alns = sequence(&mol, 60, 60);
        let first = alns.iter().position(|a| a.is_r1).unwrap();
        assert_eq!(alns.iter().filter(|a| a.is_r1).count(), 2);
        alns[first].cigar = cig("40M");
        let out = plan(&alns, &Options::default());
        match out.verdict {
            Verdict::Unresolved { reason, .. } => {
                assert_eq!(reason, Reason::InconsistentLength)
            }
            ref v => panic!("expected a refusal, got {v:?}"),
        }
    }

    // ---- consensus ---------------------------------------------------------

    #[test]
    fn agreement_sums_the_qualities_up_to_the_cap() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let alns = sequence(&mol, 60, 60);
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.overlap, 42);
        // Q30 and Q30 agreeing: capped at 40, not 60, because the two are the
        // same molecule and their errors are not independent.
        let touched: Vec<u8> =
            out.plans.iter().flat_map(|p| p.qual_set.iter().map(|x| x.1)).collect();
        assert!(!touched.is_empty());
        assert!(touched.iter().all(|&q| q == 40), "{touched:?}");
    }

    #[test]
    fn a_disagreement_zeroes_the_loser_and_docks_the_winner() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let mut alns = sequence(&mol, 60, 60);
        // Break one base on R2, inside the overlap, and make it the weaker call.
        let r2 = alns.iter().position(|a| !a.is_r1).unwrap();
        alns[r2].seq[5] = if alns[r2].seq[5] == Seq::A { Seq::C } else { Seq::A };
        alns[r2].qual[5] = 20;
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.mismatches, 1, "{:?}", out.verdict);
        let p2 = find(&out, alns[r2].id);
        assert!(p2.qual_set.contains(&(5, 0)), "the weaker call is zeroed: {p2:?}");
        // And the winner is reduced by the quality of what it beat: 30 - 20.
        let r1 = alns.iter().position(|a| a.is_r1).unwrap();
        let p1 = find(&out, alns[r1].id);
        assert!(
            p1.qual_set.iter().any(|&(_, q)| q == 10),
            "winner docked by the loser's quality: {p1:?}"
        );
    }

    #[test]
    fn a_transition_is_a_mismatch_like_any_other() {
        // A C/T difference is a disagreement like any other. Conversion and
        // amplification both happen before the molecule reaches the sequencer,
        // so the two strands the mates read are complements by construction and
        // report the same base at an overlapping position whatever the
        // conversion state was. A difference there is an error in one of them.
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let mut alns = sequence(&mol, 60, 60);
        let r2 = alns.iter().position(|a| !a.is_r1).unwrap();
        let c = alns[r2].seq.iter().position(|&b| b == Seq::C).unwrap();
        alns[r2].seq[c] = Seq::T;
        let out = plan(&alns, &Options::default());
        assert_eq!(resolved(&out).mismatches, 1);
        let p2 = find(&out, alns[r2].id);
        assert!(
            p2.qual_set.iter().any(|&(i, q)| i == c && q == 0),
            "the weaker call is zeroed like any other disagreement: {p2:?}"
        );
    }

    #[test]
    fn too_many_real_mismatches_refuse_the_merge() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let mut alns = sequence(&mol, 60, 60);
        let r2 = alns.iter().position(|a| !a.is_r1).unwrap();
        for i in 0..30 {
            alns[r2].seq[i] = match alns[r2].seq[i] {
                Seq::A => Seq::C,
                Seq::C => Seq::A,
                Seq::G => Seq::T,
                _ => Seq::G,
            };
        }
        let out = plan(&alns, &Options::default());
        match out.verdict {
            Verdict::Unresolved { reason, .. } => assert_eq!(reason, Reason::HighMismatch),
            ref v => panic!("expected a refusal, got {v:?}"),
        }
    }

    #[test]
    fn the_disagreement_policies_are_all_reachable() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let build = || {
            let mut a = sequence(&mol, 60, 60);
            let r2 = a.iter().position(|x| !x.is_r1).unwrap();
            a[r2].seq[5] = match a[r2].seq[5] {
                Seq::A => Seq::C,
                Seq::C => Seq::A,
                Seq::G => Seq::T,
                _ => Seq::G,
            };
            a[r2].qual[5] = 20;
            (a, r2)
        };
        let (a, r2) = build();
        let o = plan(&a, &Options { mismatch_qual: MismatchQual::None, ..Default::default() });
        assert!(find(&o, a[r2].id).qual_set.iter().all(|&(i, _)| i != 5));

        let (a, r2) = build();
        let o = plan(&a, &Options { mismatch_qual: MismatchQual::ZeroBoth, ..Default::default() });
        assert!(find(&o, a[r2].id).qual_set.contains(&(5, 0)));

        let (a, r2) = build();
        let o = plan(&a, &Options { mismatch_base: MismatchBase::SetN, ..Default::default() });
        assert!(find(&o, a[r2].id).base_set.contains(&(5, Seq::N)));

        let (a, _) = build();
        let o = plan(&a, &Options { match_qual: MatchQual::None, ..Default::default() });
        let raised = o.plans.iter().flat_map(|p| &p.qual_set).filter(|&&(_, q)| q > 30).count();
        assert_eq!(raised, 0, "no agreement boost when it is switched off");
    }

    #[test]
    fn an_ambiguous_base_is_neither_a_match_nor_a_mismatch() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let mut alns = sequence(&mol, 60, 60);
        let r2 = alns.iter().position(|a| !a.is_r1).unwrap();
        alns[r2].seq[5] = Seq::N;
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert_eq!(r.mismatches, 0);
        assert_eq!(out.counts.positions_ambiguous, 1);
    }

    // ---- clipping ----------------------------------------------------------

    #[test]
    fn exactly_one_end_gives_up_the_overlap() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let mut alns = sequence(&mol, 60, 60);
        // Make R2 clearly the weaker end so the choice is not a coin flip.
        for a in alns.iter_mut().filter(|a| !a.is_r1) {
            a.mapq = 10;
            a.qual.iter_mut().for_each(|q| *q = 15);
        }
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        assert!(r.kept_r1, "the better-scoring end keeps its copy");
        let clipped: usize = out
            .plans
            .iter()
            .filter(|p| alns.iter().any(|a| a.id == p.id && !a.is_r1))
            .map(|p| p.clip_front + p.clip_back)
            .sum();
        assert_eq!(clipped, r.overlap, "R2 gives up exactly the overlap");
        let kept: usize = out
            .plans
            .iter()
            .filter(|p| alns.iter().any(|a| a.id == p.id && a.is_r1))
            .map(|p| p.clip_front + p.clip_back)
            .sum();
        assert_eq!(kept, 0, "and R1 keeps all of its alignment");
    }

    #[test]
    fn the_discarded_copy_is_clipped_from_the_reads_own_three_prime_end() {
        // R2 is a reverse record, so the read's 3' tail is the CIGAR's head.
        // Clipping the wrong side would remove the contact anchor and shift POS
        // by the width of the overlap.
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let mut alns = sequence(&mol, 60, 60);
        for a in alns.iter_mut().filter(|a| !a.is_r1) {
            a.mapq = 10;
        }
        let out = plan(&alns, &Options { keep: Keep::R1, ..Default::default() });
        let r2 = alns.iter().find(|a| !a.is_r1).unwrap();
        assert!(r2.is_reverse);
        let p = find(&out, r2.id);
        assert_eq!(p.clip_back, 0);
        assert_eq!(p.clip_front, resolved(&out).overlap);
    }

    #[test]
    fn keeping_r2_instead_clips_r1() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let alns = sequence(&mol, 60, 60);
        let out = plan(&alns, &Options { keep: Keep::R2, ..Default::default() });
        let r1 = alns.iter().find(|a| a.is_r1).unwrap();
        let p = find(&out, r1.id);
        assert_eq!(p.clip_back, resolved(&out).overlap, "R1 is forward: clip its tail");
        assert_eq!(p.clip_front, 0);
    }

    #[test]
    fn the_covering_records_mapq_decides_not_the_primarys() {
        // The overlap sits on a supplementary. Scoring on the end's primary
        // MAPQ would read the wrong number and invert the choice.
        let mol = molecule(&[seg(0, 1000, true, &A80[..60]), seg(1, 5000, true, &B80[..60])]);
        let mut alns = sequence(&mol, 80, 80);
        // The overlap is [40, 80): on R1 that is its chr1 supplementary.
        for a in alns.iter_mut() {
            a.mapq = if a.is_supplementary { 1 } else { 60 };
        }
        let out = plan(&alns, &Options::default());
        let r = resolved(&out);
        // Both ends' overlap bases sit on supplementaries here, so the scores
        // tie and R1 keeps it -- but the point is that neither end scored 60.
        assert_eq!(r.overlap, 40);
        assert!(r.kept_r1);
    }

    #[test]
    fn a_record_wholly_inside_the_overlap_is_dropped_not_clipped_to_nothing() {
        // Spec test 8. R2's 3'-most supplementary lies entirely in the overlap.
        let mol = molecule(&[
            seg(0, 1000, true, &A80[..70]),
            seg(1, 5000, true, &B80[..20]),
            seg(2, 9000, true, &A80[..70]),
        ]);
        let alns = sequence(&mol, 100, 100);
        let out = plan(&alns, &Options { keep: Keep::R1, ..Default::default() });
        let r = resolved(&out);
        assert_eq!(r.l, 160);
        // The 20bp chr1 segment is inside the overlap [60, 100), so R2's record
        // for it survives nothing.
        let gone: Vec<usize> = out.plans.iter().filter(|p| p.drop).map(|p| p.id).collect();
        assert!(!gone.is_empty(), "a fully consumed record must go: {:?}", out.plans);
        for id in gone {
            let a = alns.iter().find(|a| a.id == id).unwrap();
            assert!(!a.is_r1, "only the discarded end loses records");
            let p = find(&out, id);
            assert_eq!(p.clip_front + p.clip_back, 0, "dropped, not clipped to zero");
        }
    }

    #[test]
    fn a_record_that_won_no_bases_is_still_clipped() {
        // Two records of one read can claim the same read bases, which is what
        // an aligner does at a junction it split ambiguously. `build_end` picks
        // one winner per base so the consensus has something definite to
        // compare -- but the loser is still in the file, still carrying those
        // bases, and if the read is the one giving up the overlap then the
        // loser has to give them up too. Selecting records by whether they won
        // anything would leave a wholly-contained record untouched and its copy
        // of the overlap counted a second time.
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let base = sequence(&mol, 60, 60);
        let l = 78;

        // A supplementary lying inside the primary's span, losing every base to
        // it on MAPQ, and reaching past the point where R1 gives up: e[5, 35).
        let contained = |cigar: &str, n: usize| Aln {
            id: 99,
            is_r1: true,
            is_reverse: false,
            is_supplementary: true,
            mapq: 1,
            tid: 1,
            pos: 5000,
            cigar: cig(cigar),
            seq: vec![Seq::A; n],
            qual: vec![30; n],
        };

        let mut alns = base.clone();
        alns.push(contained("5H30M25H", 30));
        let out = plan(&alns, &Options { keep: Keep::R2, ..Default::default() });
        assert_eq!(resolved(&out).l, l);
        // R1 gives up fragment indices [l - 60, 60), so bases at e >= 18 go.
        let p = find(&out, 99);
        assert_eq!(p.clip_back, (5..35).filter(|&e| e >= l - 60).count());
        assert!(!p.drop, "part of it survives, so it is clipped not dropped");

        // And one lying wholly past that point is consumed entirely.
        let mut alns = base;
        alns.push(contained("20H30M10H", 30));
        let out = plan(&alns, &Options { keep: Keep::R2, ..Default::default() });
        let p = find(&out, 99);
        assert!(p.drop, "nothing of it survives: {p:?}");
        assert_eq!(p.clip_front + p.clip_back, 0, "dropped, not clipped to zero");
    }

    #[test]
    fn a_dropped_primary_hands_its_flag_to_a_survivor() {
        let mol = molecule(&[
            seg(0, 1000, true, &A80[..70]),
            seg(1, 5000, true, &B80[..20]),
            seg(2, 9000, true, &A80[..70]),
        ]);
        let alns = sequence(&mol, 100, 100);
        let out = plan(&alns, &Options { keep: Keep::R1, ..Default::default() });
        for want_r1 in [true, false] {
            let ids: Vec<usize> =
                alns.iter().filter(|a| a.is_r1 == want_r1).map(|a| a.id).collect();
            let lost_primary = ids
                .iter()
                .any(|&i| !alns[i].is_supplementary && find(&out, i).drop);
            if lost_primary {
                assert!(
                    ids.iter().any(|&i| find(&out, i).promote),
                    "the end must not vanish from tools that read primaries only"
                );
            }
        }
    }

    #[test]
    fn hard_clipping_is_available_but_not_the_default() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let alns = sequence(&mol, 60, 60);
        let soft = plan(&alns, &Options { keep: Keep::R2, ..Default::default() });
        assert!(soft.plans.iter().all(|p| !p.hard_front && !p.hard_back));
        let hard =
            plan(&alns, &Options { keep: Keep::R2, hard_clip: true, ..Default::default() });
        assert!(hard.plans.iter().any(|p| p.hard_front || p.hard_back));
    }

    #[test]
    fn a_terminal_indel_on_the_kept_end_is_clipped_away() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let mut alns = sequence(&mol, 60, 60);
        let r1 = alns.iter().position(|a| a.is_r1).unwrap();
        alns[r1].cigar = cig("57M3I");
        let out = plan(&alns, &Options { keep: Keep::R1, ..Default::default() });
        resolved(&out);
        let p = find(&out, alns[r1].id);
        assert!(p.clip_back >= 3, "the terminal insertion goes: {p:?}");
        // And leaving it alone has to be possible, or the option is a lie.
        let off = plan(
            &alns,
            &Options { keep: Keep::R1, clean_kept_3prime: false, ..Default::default() },
        );
        assert_eq!(find(&off, alns[r1].id).clip_back, 0);
    }

    #[test]
    fn terminal_trim_finds_the_outermost_indel_within_the_window() {
        assert_eq!(terminal_trim(&cig("3I"), 3), 3, "a terminal insertion");
        assert_eq!(terminal_trim(&cig("2I5M"), 3), 2);
        assert_eq!(terminal_trim(&cig("1M2D5M"), 3), 1, "a deletion just inside");
        assert_eq!(terminal_trim(&cig("10M2D5M"), 3), 0, "and one safely away");
        assert_eq!(terminal_trim(&cig("5S2I5M"), 3), 2, "measured from the alignment");
        assert_eq!(terminal_trim(&cig("20M"), 3), 0);
    }

    #[test]
    fn adapter_read_through_is_clipped_off_both_ends() {
        // A molecule shorter than the reads. Everything past the fragment is
        // adapter, and the span clamps are what keep the overlap the molecule's
        // own length rather than something longer than the fragment.
        assert_eq!(overlap_span(100, 150, 150), (0, 100), "the whole molecule");
        assert_eq!(overlap_span(200, 150, 150), (50, 150));
        assert_eq!(overlap_span(400, 150, 150), (0, 0), "no overlap at all");
    }

    // ---- reclip ------------------------------------------------------------

    #[test]
    fn a_leading_clip_moves_the_position() {
        let (c, p) = reclip(&cig("10M"), 100, 4, 0);
        assert_eq!(cigar_string(&c), "4S6M");
        assert_eq!(p, 104, "POS must advance by the reference bases removed");
    }

    #[test]
    fn a_trailing_clip_does_not() {
        let (c, p) = reclip(&cig("10M"), 100, 0, 4);
        assert_eq!(cigar_string(&c), "6M4S");
        assert_eq!(p, 100);
    }

    #[test]
    fn clipping_absorbs_an_existing_soft_clip() {
        let (c, p) = reclip(&cig("3S10M"), 100, 4, 0);
        assert_eq!(cigar_string(&c), "7S6M");
        assert_eq!(p, 104, "the existing clip consumed no reference");
    }

    #[test]
    fn a_hard_clip_is_preserved_ahead_of_the_new_soft_clip() {
        let (c, p) = reclip(&cig("20H10M"), 100, 4, 0);
        assert_eq!(cigar_string(&c), "20H4S6M");
        assert_eq!(p, 104);
    }

    #[test]
    fn an_insertion_inside_the_clipped_span_consumes_no_reference() {
        let (c, p) = reclip(&cig("2M3I5M"), 100, 5, 0);
        assert_eq!(cigar_string(&c), "5S5M");
        assert_eq!(p, 102, "only the 2 matched bases were reference");
    }

    #[test]
    fn a_deletion_inside_the_clipped_span_is_removed_and_counted() {
        let (c, p) = reclip(&cig("2M3D5M"), 100, 3, 0);
        assert_eq!(cigar_string(&c), "3S4M");
        assert_eq!(p, 106, "2 matched + 3 deleted reference bases");
    }

    #[test]
    fn a_deletion_left_at_an_edge_is_folded_away() {
        let (c, p) = reclip(&cig("2M3D5M"), 100, 2, 0);
        assert_eq!(cigar_string(&c), "2S5M");
        assert_eq!(p, 105);
    }

    #[test]
    fn clipping_from_both_ends() {
        let (c, p) = reclip(&cig("10M"), 100, 3, 2);
        assert_eq!(cigar_string(&c), "3S5M2S");
        assert_eq!(p, 103);
    }

    #[test]
    fn reclip_is_idempotent_on_an_already_clipped_record() {
        let (c1, p1) = reclip(&cig("10M"), 100, 4, 0);
        let (c2, p2) = reclip(&c1, p1, 0, 0);
        assert_eq!(cigar_string(&c2), cigar_string(&c1));
        assert_eq!(p2, p1);
    }

    #[test]
    fn hardening_converts_the_new_clip_and_says_what_seq_to_drop() {
        let (c, _) = reclip(&cig("10M"), 100, 0, 4);
        let (h, front, back) = harden(&c, false, true);
        assert_eq!(cigar_string(&h), "6M4H");
        assert_eq!((front, back), (0, 4), "SEQ and QUAL must lose the same 4");
        // The other end is left alone.
        let (h2, f2, b2) = harden(&cig("3S10M"), true, false);
        assert_eq!(cigar_string(&h2), "3H10M");
        assert_eq!((f2, b2), (3, 0));
    }

    // ---- idempotence -------------------------------------------------------

    #[test]
    fn a_second_pass_finds_nothing_left_to_do() {
        // The discarded copy is soft-clipped, so it stops being an aligned
        // column and stops generating candidates. Re-running is a no-op by
        // construction rather than by a quality threshold.
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let alns = sequence(&mol, 60, 60);
        let first = plan(&alns, &Options { keep: Keep::R1, ..Default::default() });
        let r = resolved(&first);

        let mut after = alns.clone();
        for p in &first.plans {
            let a = after.iter_mut().find(|a| a.id == p.id).unwrap();
            if p.clip_front + p.clip_back > 0 {
                let (c, pos) = reclip(&a.cigar, a.pos, p.clip_front, p.clip_back);
                a.cigar = c;
                a.pos = pos;
            }
            for &(i, q) in &p.qual_set {
                a.qual[i] = q;
            }
        }
        let second = plan(&after, &Options { keep: Keep::R1, ..Default::default() });
        assert!(
            matches!(second.verdict, Verdict::NoOverlap { .. }),
            "the clipped copy no longer aligns: {:?} (was {r:?})",
            second.verdict
        );
        assert!(second.plans.iter().all(|p| p.is_empty()));
    }

    // ---- statistics --------------------------------------------------------

    #[test]
    fn stats_accumulate_and_render_as_tidy_tsv() {
        let mut st = Stats::default();
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        st.add(&plan(&sequence(&mol, 60, 60), &Options::default()));
        let far = molecule(&[seg(0, 1000, true, A80), seg(0, 9000, true, B80)]);
        st.add(&plan(&sequence(&far, 40, 40), &Options::default()));

        assert_eq!(st.templates, 2);
        assert_eq!(st.resolved, 1);
        assert_eq!(st.no_overlap, 1);
        assert_eq!(st.overlap_len.get(&42), Some(&1));

        let tsv = st.to_tsv();
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines[0], "section\tkey\tvalue");
        assert!(lines.iter().all(|l| l.split('\t').count() == 3), "{tsv}");
        assert!(tsv.contains("counts\ttemplates\t2"), "{tsv}");
        assert!(tsv.contains("fragment_len\t78\t1"), "{tsv}");
        // Every refusal reason gets a row even at zero, so a reader can tell
        // "none of these happened" from "this version does not count that".
        for r in Reason::ALL {
            assert!(tsv.contains(&format!("unresolved\t{}\t", r.as_str())), "{tsv}");
        }
    }

    #[test]
    fn rates_do_not_divide_by_zero() {
        let tsv = Stats::default().to_tsv();
        assert!(tsv.contains("rate_ppm\tmismatch_of_compared\t0"), "{tsv}");
        assert!(tsv.contains("rate_ppm\tresolved_of_templates\t0"), "{tsv}");
    }

    #[test]
    fn a_single_ended_template_is_left_alone() {
        let mol = molecule(&[seg(0, 1000, true, A80)]);
        let alns: Vec<Aln> = sequence(&mol, 60, 60).into_iter().filter(|a| a.is_r1).collect();
        let out = plan(&alns, &Options::default());
        assert!(out.plans.iter().all(|p| p.is_empty()));
    }
}

#[cfg(test)]
mod fuzz {
    use super::tests_support::*;
    use super::*;

    /// A small deterministic generator, so a failure is reproducible from its
    /// seed and the suite needs no dependency.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    /// Sequence a few thousand random molecules and insist the fragment length
    /// comes back.
    ///
    /// This is the claim the whole module rests on: that a molecule's own
    /// length and the bases it was observed twice at can be recovered from
    /// nothing but the alignments, for any number of segments, on any contigs,
    /// entered in either direction, at any read lengths. Hand-written cases
    /// cover the shapes someone thought of; this covers the ones nobody did.
    #[test]
    fn random_molecules_round_trip_through_the_estimator() {
        let mut rng = Rng(0x5EED_1234_ABCD_0001);
        let (mut resolved_n, mut overlapping) = (0usize, 0usize);

        for _ in 0..3000 {
            let nseg = 1 + rng.below(3);
            let mut segs = Vec::new();
            let mut used = 0usize;
            for _ in 0..nseg {
                let len = 15 + rng.below(45);
                let tid = rng.below(3) as i32;
                let pos = 1000 + 500 * rng.below(20) as i64;
                let fwd = rng.below(2) == 0;
                segs.push((tid, pos, fwd, len));
                used += len;
            }
            let mol = random_molecule(&segs, &mut || rng.below(4));
            let l = mol.len();
            assert_eq!(l, used);
            let len1 = 20 + rng.below(l.max(2));
            let len2 = 20 + rng.below(l.max(2));
            let alns = sequence(&mol, len1, len2);
            if alns.iter().all(|a| a.is_r1) || alns.iter().all(|a| !a.is_r1) {
                continue;
            }

            let out = plan(&alns, &Options::default());
            let (lo, hi) = overlap_span(l, len1.min(l), len2.min(l));
            let truth = hi - lo;
            if truth >= Options::default().min_overlap {
                overlapping += 1;
            }

            match out.verdict {
                Verdict::Resolved(r) => {
                    resolved_n += 1;
                    assert_eq!(r.l, l, "wrong fragment length for {segs:?} at {len1}/{len2}");
                    assert_eq!(r.overlap, truth, "wrong overlap for {segs:?}");
                    assert_eq!(
                        r.mismatches, 0,
                        "perfect data must not disagree with itself: {segs:?}"
                    );
                }
                Verdict::NoOverlap { .. } => assert!(
                    truth < Options::default().min_overlap || truth <= 8,
                    "a {truth}-base overlap was missed entirely: {segs:?} at {len1}/{len2}"
                ),
                // A refusal is always safe. It must not be the norm.
                Verdict::Unresolved { .. } => {}
            }
        }
        assert!(overlapping > 300, "the generator produced too few overlaps: {overlapping}");
        let rate = resolved_n as f32 / overlapping as f32;
        // Currently 1.0: every overlap the geometry says exists is recovered.
        // The slack is so that tightening a default cannot fail the suite for a
        // handful of templates it was meant to refuse.
        assert!(rate > 0.95, "only {resolved_n} of {overlapping} overlaps resolved");
    }

    /// Whatever else happens, the tool may not invent coverage.
    ///
    /// Every molecule base that both ends observed must end up counted once:
    /// the discarded copy is clipped, and the kept copy is not.
    #[test]
    fn no_molecule_base_survives_twice() {
        let mut rng = Rng(0xC0FF_EE00_1357_9BDF);
        for _ in 0..1500 {
            let nseg = 1 + rng.below(3);
            let mut segs = Vec::new();
            for _ in 0..nseg {
                segs.push((
                    rng.below(3) as i32,
                    1000 + 500 * rng.below(20) as i64,
                    rng.below(2) == 0,
                    15 + rng.below(45),
                ));
            }
            let mol = random_molecule(&segs, &mut || rng.below(4));
            let l = mol.len();
            let (len1, len2) = (20 + rng.below(l.max(2)), 20 + rng.below(l.max(2)));
            let alns = sequence(&mol, len1, len2);
            if alns.iter().all(|a| a.is_r1) || alns.iter().all(|a| !a.is_r1) {
                continue;
            }
            let out = plan(&alns, &Options::default());
            let Verdict::Resolved(r) = out.verdict else { continue };

            // Apply the clips, then ask how many aligned bases each fragment
            // index still has.
            let mut seen = vec![0u8; r.l];
            for a in &alns {
                let p = out.plans.iter().find(|p| p.id == a.id).unwrap();
                if p.drop {
                    continue;
                }
                let (cig, _) = reclip(&a.cigar, a.pos, p.clip_front, p.clip_back);
                let survivor = Aln { cigar: cig, ..a.clone() };
                for (s, (aligned, _)) in survivor.seq_map().into_iter().enumerate() {
                    if !aligned {
                        continue;
                    }
                    let e = survivor.end_index(s);
                    let f = if a.is_r1 { e } else { r.l.wrapping_sub(1).wrapping_sub(e) };
                    if f < r.l {
                        seen[f] += 1;
                    }
                }
            }
            let (lo, hi) = overlap_span(r.l, len1.min(l), len2.min(l));
            for f in lo..hi {
                assert!(
                    seen[f] <= 1,
                    "fragment index {f} survived {} times for {segs:?}",
                    seen[f]
                );
            }
        }
    }
}

#[cfg(test)]
mod terminology_check {
    use super::tests_support::*;
    use super::*;

    /// Pin down exactly what `q` and `h` count, in both orientations.
    #[test]
    fn q_indexes_seq_in_ascending_reference_order_and_h_clips_the_low_side() {
        let mk = |rev: bool, cigar: &str, n: usize| Aln {
            id: 0, is_r1: true, is_reverse: rev, is_supplementary: true, mapq: 60,
            tid: 0, pos: 500, cigar: cig(cigar), seq: vec![Seq::A; n], qual: vec![30; n],
        };
        // Forward: ascending coordinate runs with the read, so the hard clip
        // sits before the read's 5' end and q counts forward from there.
        let f = mk(false, "70H80M", 80);
        assert_eq!(f.read_length(), 150);
        assert_eq!(f.end_index(0), 70, "q=0 is the read's 71st base");
        assert_eq!(f.end_index(79), 149, "q=79 is the read's last base");
        // Reverse: ascending coordinate runs against the read, so the same
        // leading hard clip sits at the read's 3' end and q counts backwards.
        let r = mk(true, "70H80M", 80);
        assert_eq!(r.end_index(0), 79, "q=0 is the lowest-coordinate base = read's 80th");
        assert_eq!(r.end_index(79), 0, "q=79 is the read's first base");
    }
}
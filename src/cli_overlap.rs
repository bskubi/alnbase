//! Parse the arguments of the `overlap` subcommand. This module does nothing
//! else.
//!
//! Shares nothing with the scan CLI. The two are independent tools that
//! happen to ship in one binary: this one is a BAM-to-BAM filter meant to run
//! upstream of *any* caller, including callers that are not this one.
//!
//! The defaults are the conservative ones. A template whose fragment length
//! cannot be established is passed through untouched rather than guessed at,
//! the discarded copy is soft-clipped rather than destroyed, and qualities are
//! never raised past a cap that assumes the two mates are not independent
//! observations -- because they are not.

use std::path::PathBuf;

use clap::{Args, ValueEnum};

use crate::overlap::{Keep, MatchQual, MismatchBase, MismatchQual, Options, Unresolved};

/// What to do to the qualities of two bases that disagree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MismatchQualArg {
    /// Leave both alone.
    None,
    /// Zero both. Symmetric, and picks no winner.
    ZeroBoth,
    /// Zero the weaker call only.
    ZeroLoser,
    /// Zero the weaker call and reduce the stronger one by what it beat.
    Subtract,
}

/// What to do to the base calls of two bases that disagree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MismatchBaseArg {
    /// Keep both calls; the quality policy already marks the conflict.
    None,
    /// Replace both with N.
    SetN,
}

/// What to do to the qualities of two bases that agree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MatchQualArg {
    None,
    Max,
    /// Sum, capped by `--qual-cap`.
    SumCapped,
}

/// Which end keeps its copy of the overlap.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum KeepArg {
    R1,
    R2,
    /// Whichever scores higher on the covering records' MAPQ and base quality.
    Score,
}

/// What to do with a template whose fragment length could not be established.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum UnresolvedArg {
    /// Write it unmodified, tagged. Coverage is preserved and the tag makes the
    /// affected reads filterable downstream.
    Pass,
    /// Discard the template.
    Drop,
}

/// How the discarded copy is removed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ClipMode {
    /// Convert to soft clip. MethylDackel, biscuit and mpileup all ignore
    /// soft-clipped bases, so this is enough to keep them out of a call, and
    /// unlike a hard clip it can be undone.
    Soft,
    /// Convert to hard clip, discarding the bases.
    Hard,
}

/// What to do about stale alignment tags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum StaleTags {
    /// Recompute NM and MD from the reference; drop AS, which cannot be.
    Recompute,
    /// Drop NM, MD and AS. Absent beats wrong.
    Strip,
}

/// Resolve mate overlap within a sequenced fragment
#[derive(Args, Debug)]
pub struct OverlapArgs {
    /// Input BAM, or '-' for standard input.
    ///
    /// Must be grouped by name: either sorted (`samtools sort -n`) or collated
    /// (`samtools collate -O`, which is cheaper and enough). Anything else is
    /// rejected rather than warned about -- every alignment of a template has to
    /// be visible together, and ungrouped input would resolve the few that
    /// happen to be adjacent and miss the rest, which is a near-no-op that looks
    /// like success.
    ///
    ///     samtools collate -O in.bam | alnbase overlap - - | alnbase query --query-file q.toml - ...
    ///
    /// Adapter trimming should already have happened. Read-through is handled
    /// -- bases past the end of the molecule are clipped either way -- but an
    /// adapter that mismaps rather than soft-clipping will perturb the estimate.
    pub input: PathBuf,

    /// Output BAM, or '-' for standard output. Same order as the input.
    pub output: PathBuf,

    /// Reference index. Needed only by `--stale-tags recompute`, since NM and
    /// MD are defined against the reference and clipping invalidates them.
    #[arg(long, value_name = "FILE")]
    pub refr: Option<PathBuf>,

    // -- establishing the fragment length ----------------------------------
    /// Permitted disagreement between one read's records about its length.
    #[arg(long, value_name = "N", default_value_t = 0, help_heading = "Fragment length")]
    pub length_tolerance: usize,

    /// Overlaps shorter than this are left alone.
    #[arg(long, value_name = "N", default_value_t = 4, help_heading = "Fragment length")]
    pub min_overlap: usize,

    /// Reference-matched base pairs required to believe a fragment length.
    #[arg(long, value_name = "N", default_value_t = 3, help_heading = "Fragment length")]
    pub min_support: usize,

    /// Fraction of the implied overlap the agreeing pairs must span.
    #[arg(long, value_name = "F", default_value_t = 0.7, help_heading = "Fragment length")]
    pub min_span_frac: f32,

    /// How far the supporting pairs may stop short of either read's last
    /// aligned base.
    ///
    /// This is the test that does most of the work. An overlap's two ends are
    /// the two reads' 3' termini, so the true set of matched pairs reaches both
    /// by construction; a set thrown off by a repeat, or by the same locus
    /// appearing twice in the molecule, generally does not.
    #[arg(long, value_name = "N", default_value_t = 5, help_heading = "Fragment length")]
    pub max_anchor_slack: usize,

    /// Bases a read may still have aligned past the end of the implied molecule.
    ///
    /// A fragment length is a claim about where the molecule stops, so anything
    /// beyond it is adapter and should not align. Genuine read-through leaves
    /// those bases soft-clipped and this count is zero. Raise it only to absorb
    /// a stray adapter alignment.
    #[arg(long, value_name = "N", default_value_t = 0, help_heading = "Fragment length")]
    pub max_past_end: usize,

    /// Mismatch fraction over the overlap above which nothing is merged.
    ///
    #[arg(long, value_name = "F", default_value_t = 0.15, help_heading = "Fragment length")]
    pub max_mismatch_frac: f32,

    /// A rival fragment length this well supported blocks the merge.
    #[arg(long, value_name = "F", default_value_t = 0.8, help_heading = "Fragment length")]
    pub ambiguity_ratio: f32,

    /// A rival this close in length, over a disjoint stretch of the read, is
    /// reported as an indel disagreement rather than as an ambiguous length.
    #[arg(long, value_name = "N", default_value_t = 10, help_heading = "Fragment length")]
    pub indel_window: usize,

    /// What to do with a template whose fragment length could not be
    /// established.
    #[arg(long, value_enum, default_value = "pass", help_heading = "Fragment length")]
    pub on_unresolved: UnresolvedArg,

    // -- consensus ---------------------------------------------------------
    /// What to do to the qualities of two bases that disagree.
    #[arg(long, value_enum, default_value = "subtract", help_heading = "Consensus")]
    pub mismatch_qual: MismatchQualArg,

    /// What to do to the base calls of two bases that disagree.
    #[arg(long, value_enum, default_value = "none", help_heading = "Consensus")]
    pub mismatch_base: MismatchBaseArg,

    /// What to do to the qualities of two bases that agree.
    #[arg(long, value_enum, default_value = "sum-capped", help_heading = "Consensus")]
    pub match_qual: MatchQualArg,

    /// Ceiling on a combined quality.
    ///
    /// The cap is not a formality. Two mates are the same molecule, so a
    /// pre-amplification polymerase error appears in both and the observations
    /// are not independent; summing to the SAM maximum would assert a
    /// confidence the evidence does not support, and several callers mishandle
    /// anything above 41 in any case.
    #[arg(long, value_name = "N", default_value_t = 40, help_heading = "Consensus")]
    pub qual_cap: u8,

    // -- resolution --------------------------------------------------------
    /// Which end keeps its copy of the overlap.
    #[arg(long, value_enum, default_value = "score", help_heading = "Resolution")]
    pub keep: KeepArg,

    /// Weight on mean base quality relative to mean MAPQ when scoring the ends.
    #[arg(long, value_name = "F", default_value_t = 1.0, help_heading = "Resolution")]
    pub score_qual_weight: f32,

    /// How the discarded copy is removed.
    #[arg(long, value_enum, default_value = "soft", help_heading = "Resolution")]
    pub clip_mode: ClipMode,

    /// Do not clip terminal indels off the kept 3' end.
    #[arg(long, help_heading = "Resolution")]
    pub keep_terminal_indels: bool,

    /// Matched bases within which an indel counts as terminal.
    #[arg(long, value_name = "N", default_value_t = 3, help_heading = "Resolution")]
    pub terminal_indel_window: usize,

    // -- output ------------------------------------------------------------
    /// Handling of NM, MD and AS, which describe the pre-clip alignment.
    #[arg(long, value_enum, default_value = "strip", help_heading = "Output")]
    pub stale_tags: StaleTags,

    /// Tag to carry the inferred fragment length.
    ///
    /// This is the one number another tool is likely to want, so it gets a tag
    /// of its own rather than sharing the diagnostic prefix below. Two
    /// characters, written only when a fragment resolves; a fragment whose
    /// length could not be established gets no length tag rather than a guess.
    #[arg(long, value_name = "TAG", default_value = "XL", help_heading = "Output")]
    pub length_tag: String,

    /// Prefix for the per-template tags this step writes.
    ///
    /// Two characters, of which this is the first. The default is lower case
    /// deliberately: `XM` and `XG` are Bismark's methylation-call string and
    /// genome-conversion tag, and this is a methyl pipeline, so writing there
    /// would corrupt the very thing downstream is about to read.
    #[arg(long, value_name = "C", default_value = "o", help_heading = "Output")]
    pub tag_prefix: String,

    /// Do not write the per-template tags (the length tag and the prefix
    /// tags).
    ///
    /// Tags that describe an alignment this step changed are still kept
    /// correct, because leaving them stale would make the record wrong: NM, MD,
    /// AS and XA of a clipped record (see --stale-tags), SA of a read whose
    /// alignments changed, and MC/MQ of a template whose mate fields changed.
    #[arg(long, help_heading = "Output")]
    pub no_tag: bool,

    /// Write per-run statistics to this file, or '-' for standard error.
    ///
    /// A tidy three-column TSV: `section`, `key`, `value`, so scalars and
    /// histograms share one row shape and the whole file loads with a single
    /// `read_csv`. It is the only record of how much evidence this step
    /// removed, which is both what a reviewer eventually asks for and the first
    /// thing to check when a concordance test disagrees by a fraction of a
    /// percent.
    #[arg(long, value_name = "FILE", help_heading = "Output")]
    pub stats: Option<PathBuf>,

    /// Report each refused template to standard error.
    #[arg(long, help_heading = "Output")]
    pub report_unresolved: bool,

    /// Compression threads for reading and writing the BAM.
    #[arg(short = '@', long, default_value_t = 1, help_heading = "Output")]
    pub threads: usize,

    /// BGZF compression level of the output, 0 (uncompressed) to 9.
    #[arg(long, value_name = "N", default_value_t = 3, help_heading = "Output")]
    pub compression_level: i32,
}

impl OverlapArgs {
    pub fn options(&self) -> Options {
        Options {
            length_tolerance: self.length_tolerance,
            min_overlap: self.min_overlap,
            min_support: self.min_support,
            min_span_frac: self.min_span_frac,
            max_anchor_slack: self.max_anchor_slack,
            max_past_end: self.max_past_end,
            max_mismatch_frac: self.max_mismatch_frac,
            ambiguity_ratio: self.ambiguity_ratio,
            indel_window: self.indel_window,
            on_unresolved: match self.on_unresolved {
                UnresolvedArg::Pass => Unresolved::Pass,
                UnresolvedArg::Drop => Unresolved::Drop,
            },
            mismatch_qual: match self.mismatch_qual {
                MismatchQualArg::None => MismatchQual::None,
                MismatchQualArg::ZeroBoth => MismatchQual::ZeroBoth,
                MismatchQualArg::ZeroLoser => MismatchQual::ZeroLoser,
                MismatchQualArg::Subtract => MismatchQual::Subtract,
            },
            mismatch_base: match self.mismatch_base {
                MismatchBaseArg::None => MismatchBase::None,
                MismatchBaseArg::SetN => MismatchBase::SetN,
            },
            match_qual: match self.match_qual {
                MatchQualArg::None => MatchQual::None,
                MatchQualArg::Max => MatchQual::Max,
                MatchQualArg::SumCapped => MatchQual::SumCapped,
            },
            qual_cap: self.qual_cap,
            keep: match self.keep {
                KeepArg::R1 => Keep::R1,
                KeepArg::R2 => Keep::R2,
                KeepArg::Score => Keep::Score,
            },
            score_qual_weight: self.score_qual_weight,
            hard_clip: self.clip_mode == ClipMode::Hard,
            clean_kept_3prime: !self.keep_terminal_indels,
            terminal_indel_window: self.terminal_indel_window,
            // Only worth collecting when there is a reference to score them
            // against.
        }
    }

    /// The two-character tag for a one-character suffix.
    pub fn tag(&self, suffix: char) -> [u8; 2] {
        [self.tag_prefix.as_bytes()[0], suffix as u8]
    }

    /// The tag carrying the inferred fragment length.
    pub fn length_tag(&self) -> [u8; 2] {
        let b = self.length_tag.as_bytes();
        [b[0], b[1]]
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.stale_tags == StaleTags::Recompute && self.refr.is_none() {
            return Err(
                "--stale-tags recompute needs --refr: NM and MD are defined against the \
                 reference. AS is stripped either way, because it is the aligner's own \
                 score and cannot be honestly regenerated."
                    .to_string(),
            );
        }
        if self.tag_prefix.len() != 1 || !self.tag_prefix.as_bytes()[0].is_ascii_alphabetic() {
            return Err(format!(
                "--tag-prefix must be a single letter, not {:?}",
                self.tag_prefix
            ));
        }
        // X, Y and Z are the SAM local-use space, but a methyl pipeline has
        // Bismark's XM/XR/XG in it and colliding with the methylation call
        // string would be a silent, total corruption of the output.
        if self.tag_prefix.eq_ignore_ascii_case("x") && !self.no_tag {
            return Err(
                "--tag-prefix x collides with Bismark's XM, XR and XG tags, which carry \
                 the methylation call string. Use a lower-case prefix such as the \
                 default 'o', or pass --no-tag."
                    .to_string(),
            );
        }
        // Reserved by the SAM spec for tools to define locally, but a few
        // two-letter names are spoken for by convention and one of them would
        // be a silent corruption rather than a mere collision.
        let t = self.length_tag.as_bytes();
        if t.len() != 2 || !t[0].is_ascii_alphabetic() || !t[1].is_ascii_alphanumeric() {
            return Err(format!(
                "--length-tag must be two characters, a letter then a letter or digit, \
                 not {:?}",
                self.length_tag
            ));
        }
        for (taken, whose) in [
            ("XM", "Bismark's methylation call string"),
            ("XR", "Bismark's read-conversion tag"),
            ("XG", "Bismark's genome-conversion tag"),
            ("NM", "the edit distance"),
            ("MD", "the mismatch string"),
            ("AS", "the aligner's score"),
            ("SA", "the supplementary alignment list"),
        ] {
            if self.length_tag.eq_ignore_ascii_case(taken) {
                return Err(format!(
                    "--length-tag {taken} is {whose}. Choose another two-letter tag; \
                     the default XL is free."
                ));
            }
        }
        if !self.no_tag && self.length_tag.as_bytes()[0] == self.tag_prefix.as_bytes()[0] {
            return Err(format!(
                "--length-tag {} shares its first letter with --tag-prefix {}, so it \
                 could collide with the diagnostic tags. Choose one or the other.",
                self.length_tag, self.tag_prefix
            ));
        }
        if self.qual_cap > 93 {
            return Err("--qual-cap above 93 is not representable in SAM".to_string());
        }
        if !(0.0..=1.0).contains(&self.min_span_frac) {
            return Err("--min-span-frac is a fraction, so it must be in [0, 1]".to_string());
        }
        if !(0.0..=1.0).contains(&self.max_mismatch_frac) {
            return Err("--max-mismatch-frac is a fraction, so it must be in [0, 1]".to_string());
        }
        if self.min_support == 0 {
            return Err(
                "--min-support 0 would let a fragment length be believed on no evidence \
                 at all"
                    .to_string(),
            );
        }
        if !(0..=9).contains(&self.compression_level) {
            return Err(format!(
                "--compression-level must be 0 (uncompressed) to 9, not {}",
                self.compression_level
            ));
        }
        if self.min_overlap == 0 {
            return Err(
                "--min-overlap 0 would try to resolve a zero-length overlap".to_string()
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::{CommandFactory, Parser};

    #[derive(Parser, Debug)]
    struct Harness {
        #[command(flatten)]
        args: OverlapArgs,
    }

    fn parse(extra: &[&str]) -> OverlapArgs {
        let mut argv = vec!["overlap"];
        argv.extend_from_slice(extra);
        argv.extend(["in.bam", "out.bam"]);
        Harness::try_parse_from(argv).expect("parse").args
    }

    #[test]
    fn args_are_well_formed() {
        Harness::command().debug_assert();
    }

    #[test]
    fn the_defaults_are_the_conservative_ones() {
        let a = parse(&[]);
        let o = a.options();
        assert_eq!(o.on_unresolved, Unresolved::Pass, "never guess a fragment length");
        assert!(!o.hard_clip, "reversible beats irreversible");
        assert_eq!(o.mismatch_base, MismatchBase::None, "never destroy a base call");
        assert_eq!(o.qual_cap, 40, "the mates are not independent observations");
        assert_eq!(a.stale_tags, StaleTags::Strip);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn recomputing_tags_needs_a_reference() {
        let e = parse(&["--stale-tags", "recompute"]).validate().unwrap_err();
        assert!(e.contains("--refr"), "{e}");
        assert!(parse(&["--stale-tags", "recompute", "--refr", "r.mm"]).validate().is_ok());
    }

    #[test]
    fn the_bismark_tag_space_is_refused() {
        // XM is the methylation call string. Writing a mismatch count there
        // would not be a collision so much as a silent corruption of the one
        // field the rest of the pipeline exists to read.
        let e = parse(&["--tag-prefix", "X"]).validate().unwrap_err();
        assert!(e.contains("Bismark"), "{e}");
        assert!(e.contains("XM"), "the message should name the tag: {e}");
        // Unless no tags are written at all.
        assert!(parse(&["--tag-prefix", "X", "--no-tag"]).validate().is_ok());
        assert!(parse(&[]).validate().is_ok());
    }

    #[test]
    fn the_fragment_length_gets_a_tag_of_its_own() {
        // It is the one number another tool is likely to read, so it does not
        // share the diagnostic prefix and can be pointed anywhere.
        assert_eq!(&parse(&[]).length_tag(), b"XL");
        assert_eq!(&parse(&["--length-tag", "ZF"]).length_tag(), b"ZF");
        assert!(parse(&["--length-tag", "ZF"]).validate().is_ok());
    }

    #[test]
    fn a_length_tag_that_is_spoken_for_is_refused() {
        // XM is Bismark's methylation call string. Writing a fragment length
        // there would replace the field a methyl pipeline exists to read.
        let e = parse(&["--length-tag", "XM"]).validate().unwrap_err();
        assert!(e.contains("Bismark"), "{e}");
        for taken in ["NM", "MD", "AS", "SA", "xm"] {
            assert!(
                parse(&["--length-tag", taken]).validate().is_err(),
                "{taken} should be refused"
            );
        }
        // Malformed names are refused too.
        for bad in ["X", "XLL", "1L"] {
            assert!(parse(&["--length-tag", bad]).validate().is_err(), "{bad}");
        }
    }

    #[test]
    fn the_length_tag_may_not_collide_with_the_diagnostic_prefix() {
        let e = parse(&["--length-tag", "oL"]).validate().unwrap_err();
        assert!(e.contains("--tag-prefix"), "{e}");
        // Unless the diagnostic tags are not being written at all.
        assert!(parse(&["--length-tag", "oL", "--no-tag"]).validate().is_ok());
    }

    #[test]
    fn tags_are_built_from_the_prefix() {
        assert_eq!(&parse(&[]).tag('L'), b"oL");
        assert_eq!(&parse(&["--tag-prefix", "q"]).tag('V'), b"qV");
    }

    #[test]
    fn nonsense_thresholds_are_refused() {
        assert!(parse(&["--min-support", "0"]).validate().is_err());
        assert!(parse(&["--min-overlap", "0"]).validate().is_err());
        assert!(parse(&["--qual-cap", "94"]).validate().is_err());
        assert!(parse(&["--min-span-frac", "1.5"]).validate().is_err());
        assert!(parse(&["--tag-prefix", "ov"]).validate().is_err());
    }

    #[test]
    fn a_dash_is_accepted_for_either_end() {
        let a = Harness::try_parse_from(["overlap", "-", "-"]).unwrap().args;
        assert_eq!(a.input.as_os_str(), "-");
        assert_eq!(a.output.as_os_str(), "-");
        assert!(a.validate().is_ok());
    }

    #[test]
    fn stats_are_off_unless_asked_for() {
        assert!(parse(&[]).stats.is_none());
        // Stderr, so a stats file does not collide with a BAM on stdout.
        assert_eq!(parse(&["--stats", "-"]).stats.unwrap().as_os_str(), "-");
    }

    #[test]
    fn policies_map_through() {
        assert_eq!(
            parse(&["--mismatch-qual", "zero-both"]).options().mismatch_qual,
            MismatchQual::ZeroBoth
        );
        assert_eq!(
            parse(&["--match-qual", "none"]).options().match_qual,
            MatchQual::None
        );
        assert_eq!(parse(&["--keep", "r2"]).options().keep, Keep::R2);
        assert!(parse(&["--clip-mode", "hard"]).options().hard_clip);
        assert!(!parse(&["--keep-terminal-indels"]).options().clean_kept_3prime);
        assert_eq!(
            parse(&["--on-unresolved", "drop"]).options().on_unresolved,
            Unresolved::Drop
        );
    }

}
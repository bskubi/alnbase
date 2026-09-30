//! Match a BAM tid to a reference contig by name.
//!
//! A BAM's `tid` is an index into *that BAM's* `@SQ` list. A [`Aref`] tid is
//! an index into the order the contigs appeared in the FASTA. Nothing keeps the
//! two lists in the same order, and for the common case of a karyotype-sorted
//! BAM against a lexicographically sorted UCSC FASTA they are not: `chr2` is
//! tid 1 in the BAM and tid 16 in the reference, and reference tid 1 is
//! `chr10`. Using one where the other is meant reads real sequence from the
//! wrong chromosome, which produces plausible-looking output — the bases are
//! real bases, the coordinates are real coordinates — and is close to
//! undetectable downstream.
//!
//! So the mapping is built once, by name, and validated. A name in the BAM that
//! the reference does not have is not an error by itself (a BAM aligned against
//! a larger assembly is ordinary), but records on it cannot be scanned and are
//! skipped and counted. A name both have with *different lengths* is an error:
//! same name, different sequence means the BAM was aligned against a different
//! assembly, and every coordinate is then suspect.
//!
//! # Why a `Vec` rather than a `HashMap` at the call site
//!
//! The lookup runs once per record. `tid2name` returns bytes that must be
//! validated as UTF-8 and hashed; indexing a `Vec<Option<usize>>` by tid does
//! neither. The names are resolved once, at construction, where the cost is one
//! hash per contig instead of one per record.


use anyhow::{bail, Result};
use rust_htslib::bam::HeaderView;

use crate::aref::Aref;

/// A contig in the BAM header that the reference does not have.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unresolved {
    pub name: String,
    /// The BAM's length for it, for the diagnostic.
    pub len: u64,
}

/// Maps BAM tids onto reference tids by contig name.
#[derive(Clone, Debug)]
pub struct ContigMap {
    /// Indexed by BAM tid; `None` for a contig the reference does not have.
    to_refr: Vec<Option<usize>>,
    unresolved: Vec<Unresolved>,
    #[cfg_attr(not(test), allow(dead_code))] // read by the contig-mismatch tests
    n_resolved: usize,
    /// Indexed by BAM tid: true for a contig the reference has whose `@SQ` line
    /// carries no M5. Empty unless [`require_m5`](Self::require_m5) was called.
    missing_m5: Vec<bool>,
}

impl ContigMap {
    /// Resolve every contig in `header` against `refr`.
    ///
    /// Fails if the two share no contig at all, or if a shared name has two
    /// different lengths.
    pub fn build(header: &HeaderView, refr: &Aref) -> Result<Self> {
        let n = header.target_count() as usize;
        let mut to_refr = Vec::with_capacity(n);
        let mut unresolved = Vec::new();
        let mut n_resolved = 0usize;
        let mut mismatched: Vec<String> = Vec::new();
        let m5 = header_m5(header);
        let mut md5_mismatched: Vec<String> = Vec::new();

        for tid in 0..n {
            let raw = header.tid2name(tid as u32);
            let name = match std::str::from_utf8(raw) {
                Ok(s) => s,
                Err(_) => bail!(
                    "BAM contig at tid {tid} is not valid UTF-8: {:?}",
                    String::from_utf8_lossy(raw)
                ),
            };
            let bam_len = header.target_len(tid as u32).unwrap_or(0);

            match refr.tid(name) {
                Some(rt) => {
                    // Same name, different length: not a naming difference, a
                    // different assembly. Coordinates would be meaningless.
                    let refr_len = refr.len_of(rt).unwrap_or(0) as u64;
                    if refr_len != bam_len {
                        if mismatched.len() < 10 {
                            mismatched.push(format!(
                                "{name}: {bam_len} in the BAM, {refr_len} in the reference"
                            ));
                        }
                    }
                    // Same name and length can still be different sequence (a
                    // masked analysis set, a patched assembly). The header's M5,
                    // when the aligner or `samtools dict` wrote one, settles it.
                    if let (Some(bam_md5), Some(refr_md5)) = (m5.get(name), refr.md5_hex(rt)) {
                        if !bam_md5.eq_ignore_ascii_case(&refr_md5) && md5_mismatched.len() < 10 {
                            md5_mismatched.push(format!(
                                "{name}: M5 {bam_md5} in the BAM header, {refr_md5} in the reference"
                            ));
                        }
                    }
                    to_refr.push(Some(rt));
                    n_resolved += 1;
                }
                None => {
                    to_refr.push(None);
                    if unresolved.len() < 1000 {
                        unresolved.push(Unresolved { name: name.to_string(), len: bam_len });
                    }
                }
            }
        }

        if !mismatched.is_empty() {
            bail!(
                "the BAM and the reference disagree about contig lengths, so they are \
                 different assemblies:\n  {}\n\
                 Re-index the FASTA the BAM was aligned against.",
                mismatched.join("\n  ")
            );
        }

        if !md5_mismatched.is_empty() {
            bail!(
                "the BAM header's @SQ M5 checksums disagree with the reference, so a \
                 contig with the same name and length holds different sequence:\n  {}\n\
                 Re-index the FASTA the BAM was aligned against.",
                md5_mismatched.join("\n  ")
            );
        }

        if n > 0 && n_resolved == 0 {
            let sample: Vec<&str> = (0..n.min(3)).filter_map(|i| refr.name(i)).collect();
            bail!(
                "no contig in the BAM header is present in the reference by name.\n\
                 The BAM's first contigs are {:?}; the reference's are {:?}.\n\
                 These are different references, or the same contigs under different \
                 naming conventions (`chr1` vs `1`, `chr1_GL456210_random` vs \
                 `GL456210.1`).",
                (0..n.min(3))
                    .map(|i| String::from_utf8_lossy(header.tid2name(i as u32)).into_owned())
                    .collect::<Vec<_>>(),
                sample
            );
        }

        Ok(Self { to_refr, unresolved, n_resolved, missing_m5: Vec::new() })
    }

    /// Refuse records on contigs whose `@SQ` line has no M5 (`--require-m5`).
    ///
    /// Without an M5 the only reference check is name and length, which a
    /// masked analysis set or a patched assembly passes. Decided per record
    /// rather than from the header, like an off-reference contig: a header
    /// routinely lists contigs no record is on, and those do no harm.
    pub fn require_m5(mut self, header: &HeaderView) -> Self {
        let m5 = header_m5(header);
        self.missing_m5 = (0..self.to_refr.len())
            .map(|tid| {
                self.to_refr[tid].is_some()
                    && !m5.contains_key(String::from_utf8_lossy(header.tid2name(tid as u32)).as_ref())
            })
            .collect();
        self
    }

    /// Whether a record on this BAM tid is refused for lacking an M5.
    #[inline]
    pub fn lacks_required_m5(&self, bam_tid: i32) -> bool {
        bam_tid >= 0 && self.missing_m5.get(bam_tid as usize).copied().unwrap_or(false)
    }


    /// The reference tid for a BAM tid, or `None` if the reference lacks it.
    ///
    /// `None` for an unmapped record's `-1`, and for any tid past the header.
    #[inline]
    pub fn refr_tid(&self, bam_tid: i32) -> Option<usize> {
        if bam_tid < 0 {
            return None;
        }
        self.to_refr.get(bam_tid as usize).copied().flatten()
    }

    /// True when every contig in the header was found in the reference.
    #[inline]
    #[cfg_attr(not(test), allow(dead_code))] // read by the tests that pin the contig-mismatch reporting
    pub fn is_complete(&self) -> bool {
        self.unresolved.is_empty()
    }

    #[inline]
    #[cfg_attr(not(test), allow(dead_code))] // read by the tests that pin the contig-mismatch reporting
    pub fn n_resolved(&self) -> usize {
        self.n_resolved
    }

    /// Contigs in the header the reference does not have, in header order.
    #[cfg_attr(not(test), allow(dead_code))] // read by the tests that pin the contig-mismatch reporting
    pub fn unresolved(&self) -> &[Unresolved] {
        &self.unresolved
    }

    /// One line per unresolved contig, capped, for a startup warning.
    pub fn warning(&self) -> Option<String> {
        if self.unresolved.is_empty() {
            return None;
        }
        let shown = self.unresolved.len().min(10);
        let mut s = format!(
            "{} of {} BAM contigs are not in the reference; records on them will be skipped:",
            self.unresolved.len(),
            self.to_refr.len()
        );
        for u in &self.unresolved[..shown] {
            s.push_str(&format!("\n  {} ({} bp)", u.name, u.len));
        }
        if self.unresolved.len() > shown {
            s.push_str(&format!("\n  ... and {} more", self.unresolved.len() - shown));
        }
        Some(s)
    }

}

/// `@SQ SN` -> `M5` for every contig whose header line carries an M5 tag.
pub fn header_m5(header: &HeaderView) -> std::collections::HashMap<String, String> {
    let text = String::from_utf8_lossy(header.as_bytes());
    let mut out = std::collections::HashMap::new();
    for line in text.lines().filter(|l| l.starts_with("@SQ\t")) {
        let field = |tag: &str| line.split('\t').find_map(|f| f.strip_prefix(tag));
        if let (Some(name), Some(m5)) = (field("SN:"), field("M5:")) {
            out.insert(name.to_string(), m5.to_string());
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aref::write_from_fasta;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::Header;

    /// A reference whose contig order deliberately differs from the header's.
    fn refr(tag: &str, contigs: &[(&str, &str)]) -> Aref {
        let dir = crate::test_support::temp_dir();
        let fa = dir.join(format!("contigmap_{tag}.fa"));
        let mm = dir.join(format!("contigmap_{tag}.refrmm"));
        let mut body = String::new();
        for (n, s) in contigs {
            body.push_str(&format!(">{n}\n{s}\n"));
        }
        std::fs::write(&fa, body).unwrap();
        write_from_fasta(&fa, &mm).unwrap();
        Aref::open(&mm).unwrap()
    }

    fn head(contigs: &[(&str, usize)]) -> HeaderView {
        let mut h = Header::new();
        for (n, l) in contigs {
            let mut sq = HeaderRecord::new(b"SQ");
            sq.push_tag(b"SN", n);
            sq.push_tag(b"LN", l);
            h.push_record(&sq);
        }
        HeaderView::from_header(&h)
    }

    /// The bug this module exists for: karyotype-ordered BAM, lexicographically
    /// ordered FASTA. Ordinal lookup sends chr2 to chr10.
    #[test]
    fn resolves_across_a_reordering() {
        let r = refr("reorder", &[("chr1", "AAAA"), ("chr10", "CCCC"), ("chr2", "GGGG")]);
        let h = head(&[("chr1", 4), ("chr2", 4), ("chr10", 4)]);
        let m = ContigMap::build(&h, &r).unwrap();

        assert_eq!(m.refr_tid(0), Some(0), "chr1 -> chr1");
        assert_eq!(m.refr_tid(1), Some(2), "chr2 -> chr2, not reference tid 1");
        assert_eq!(m.refr_tid(2), Some(1), "chr10 -> chr10");
        assert!(m.is_complete());
        assert_eq!(m.n_resolved(), 3);
        // The whole point: the identity map would have been wrong here.
        assert_ne!(m.refr_tid(1), Some(1));
    }

    #[test]
    fn unmapped_and_out_of_range_tids_resolve_to_none() {
        let r = refr("range", &[("chr1", "AAAA")]);
        let m = ContigMap::build(&head(&[("chr1", 4)]), &r).unwrap();
        assert_eq!(m.refr_tid(-1), None, "unmapped record");
        assert_eq!(m.refr_tid(7), None, "past the end of the header");
    }

    /// A BAM aligned against a superset assembly is ordinary; those records are
    /// skipped rather than failing the run.
    #[test]
    fn a_contig_missing_from_the_reference_is_reported_not_fatal() {
        let r = refr("missing", &[("chr1", "AAAA")]);
        let h = head(&[("chr1", 4), ("chrUn_scaffold", 9)]);
        let m = ContigMap::build(&h, &r).unwrap();

        assert_eq!(m.refr_tid(0), Some(0));
        assert_eq!(m.refr_tid(1), None);
        assert!(!m.is_complete());
        assert_eq!(m.unresolved().len(), 1);
        assert_eq!(m.unresolved()[0].name, "chrUn_scaffold");
        assert_eq!(m.unresolved()[0].len, 9);
        assert!(m.warning().unwrap().contains("chrUn_scaffold"));
    }

    /// Same name, different length. This is the check that would have caught
    /// an ordinal mapping at startup rather than after a full run.
    #[test]
    fn a_length_disagreement_is_fatal() {
        let r = refr("lens", &[("chr1", "AAAA")]);
        let h = head(&[("chr1", 999)]);
        let err = ContigMap::build(&h, &r).unwrap_err().to_string();
        assert!(err.contains("different assemblies"), "{err}");
        assert!(err.contains("chr1"), "{err}");
    }

    /// Same name and length, different sequence: only the header's M5 can tell,
    /// and when it is there it must agree. MD5("AAAA") is
    /// 098890dde069e9abad63f19a0d9e1f32, as `samtools dict` writes it.
    #[test]
    fn an_m5_disagreement_is_fatal_and_agreement_or_absence_is_not() {
        let r = refr("m5", &[("chr1", "AAAA")]);
        let with_m5 = |m5: &str| {
            let text = format!("@SQ\tSN:chr1\tLN:4\tM5:{m5}\n");
            HeaderView::from_header(&Header::from_template(&HeaderView::from_bytes(text.as_bytes())))
        };
        let err = ContigMap::build(&with_m5("0123456789abcdef0123456789abcdef"), &r).unwrap_err();
        assert!(err.to_string().contains("M5"), "{err}");
        assert!(ContigMap::build(&with_m5("098890DDE069E9ABAD63F19A0D9E1F32"), &r).is_ok(), "case-insensitive");
        assert!(ContigMap::build(&head(&[("chr1", 4)]), &r).is_ok(), "no M5, nothing to compare");
    }

    /// `--require-m5` refuses a contig the reference has and the header gives no
    /// M5 for; a contig with one, and a contig the reference lacks, are not
    /// refused by it.
    #[test]
    fn required_m5_marks_only_resolved_contigs_without_one() {
        let r = refr("require_m5", &[("chr1", "AAAA"), ("chr2", "AAAA")]);
        let text = "@SQ\tSN:chr1\tLN:4\tM5:098890dde069e9abad63f19a0d9e1f32\n\
                    @SQ\tSN:chr2\tLN:4\n@SQ\tSN:chrUn\tLN:9\n";
        let h = HeaderView::from_bytes(text.as_bytes());
        let plain = ContigMap::build(&h, &r).unwrap();
        assert!(!(0..3).any(|t| plain.lacks_required_m5(t)), "not required, never refused");
        let m = ContigMap::build(&h, &r).unwrap().require_m5(&h);
        assert_eq!((0..3).map(|t| m.lacks_required_m5(t)).collect::<Vec<_>>(), [false, true, false]);
        assert!(!m.lacks_required_m5(-1), "unmapped");
    }

    #[test]
    fn no_overlap_at_all_is_fatal() {
        let r = refr("disjoint", &[("1", "AAAA")]);
        let h = head(&[("chr1", 4)]);
        let err = ContigMap::build(&h, &r).unwrap_err().to_string();
        assert!(err.contains("no contig in the BAM header"), "{err}");
    }

    #[test]
    fn an_empty_header_is_not_an_error() {
        let r = refr("empty", &[("chr1", "AAAA")]);
        let m = ContigMap::build(&head(&[]), &r).unwrap();
        assert_eq!(m.n_resolved(), 0);
        assert_eq!(m.refr_tid(0), None);
    }
}
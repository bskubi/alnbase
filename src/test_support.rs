//! Hold the fixtures that the test modules share. A test needs a record, or a
//! genome, or both.
//!
//! Compiled only under `cfg(test)`. It exists because the alternative is each
//! test module growing its own record builder, and two builders that disagree
//! about, say, which flag bits mean the bottom strand would let a test pass
//! against a record the scanner would never see.
//!
//! [`built`] is trusted because `alignment::tests::built_records_walk_like_parsed_ones`
//! checks it against the same record round-tripped through a real BAM file.

use std::path::PathBuf;
use std::sync::Arc;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};

use crate::contig_map::ContigMap;
use crate::aref::{write_from_fasta, Aref};

/// SAM flag 0x10, read reverse strand.
pub const REVERSE: u16 = 0x10;
/// SAM flag 0x80, last segment in the template.
pub const LAST_IN_TEMPLATE: u16 = 0x80;

/// A distinct temp path per test, so tests that write files can run in
/// parallel. Callers pass a tag unique within the crate.
pub fn temp_path(tag: &str, ext: &str) -> PathBuf {
    std::env::temp_dir().join(format!("alnbase_test_{tag}.{ext}"))
}

/// A one-contig reference index on disk, opened the way the scanner opens one.
pub fn reference(tag: &str, contig: &str, seq: &str) -> Aref {
    let fa = temp_path(tag, "fa");
    let mm = temp_path(tag, "aref");
    std::fs::write(&fa, format!(">{contig}\n{seq}\n")).unwrap();
    write_from_fasta(&fa, &mm).unwrap();
    Aref::open(&mm).unwrap()
}

/// A reference holding several contigs, in the order given.
///
/// The order is the point: it is what a real FASTA fixes and what a BAM header
/// is free to disagree with.
pub fn reference_many(tag: &str, contigs: &[(&str, &str)]) -> Aref {
    let fa = temp_path(tag, "fa");
    let mm = temp_path(tag, "aref");
    let mut body = String::new();
    for (n, s) in contigs {
        body.push_str(&format!(">{n}\n{s}\n"));
    }
    std::fs::write(&fa, body).unwrap();
    write_from_fasta(&fa, &mm).unwrap();
    Aref::open(&mm).unwrap()
}

/// A header naming several contigs, in the order given.
pub fn header_many(contigs: &[(&str, usize)]) -> HeaderView {
    let mut h = Header::new();
    for (n, l) in contigs {
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", n);
        sq.push_tag(b"LN", l);
        h.push_record(&sq);
    }
    HeaderView::from_header(&h)
}

/// A header naming one contig, for the record columns that resolve tids.
pub fn header(contig: &str, len: usize) -> HeaderView {
    let mut h = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", contig);
    sq.push_tag(b"LN", len);
    h.push_record(&sq);
    HeaderView::from_header(&h)
}

/// A mapped record assembled in memory.
///
/// `flags` carries the strand. `walk_alignment` reads the *fragment* strand as
/// `is_last_in_template() == is_reverse()`, so a lone [`REVERSE`] is the bottom
/// strand and `REVERSE | LAST_IN_TEMPLATE` is the top one.
pub fn built(name: &[u8], seq: &[u8], cigar: &[Cigar], pos: i64, flags: u16) -> Record {
    let mut rec = Record::new();
    let qual = vec![37u8; seq.len()];
    rec.set(name, Some(&CigarString(cigar.to_vec())), seq, &qual);
    rec.set_tid(0);
    rec.set_pos(pos);
    rec.set_mapq(60);
    rec.set_mtid(-1);
    rec.set_mpos(-1);
    rec.set_insert_size(0);
    rec.set_flags(flags);
    rec
}

/// Write `records` to a real BAM file and return its path.
pub fn write_bam(tag: &str, contig: &str, len: usize, records: &[Record]) -> PathBuf {
    let path = temp_path(tag, "bam");
    let mut h = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", contig);
    sq.push_tag(b"LN", len);
    h.push_record(&sq);
    let mut w = Writer::from_path(&path, &h, Format::Bam).unwrap();
    for r in records {
        w.write(r).unwrap();
    }
    path
}

/// A BAM whose header lists several contigs, for the cases where the header
/// is broader than the reference or a record sits on a contig the reference
/// does not have.
pub fn write_bam_many(tag: &str, contigs: &[(&str, usize)], records: &[Record]) -> PathBuf {
    let path = temp_path(tag, "bam");
    let mut h = Header::new();
    for (n, l) in contigs {
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", n);
        sq.push_tag(b"LN", l);
        h.push_record(&sq);
    }
    let mut w = Writer::from_path(&path, &h, Format::Bam).unwrap();
    for r in records {
        w.write(r).unwrap();
    }
    path
}

/// The same record after a round trip through a real BAM file.
///
/// Used to check that [`built`] is a faithful stand-in for a parsed record.
pub fn round_tripped(rec: &Record, tag: &str, contig: &str, len: usize) -> Record {
    use rust_htslib::bam::Read as _;

    let path = write_bam(tag, contig, len, std::slice::from_ref(rec));
    let mut r = rust_htslib::bam::Reader::from_path(&path).unwrap();
    r.records().next().unwrap().unwrap()
}

/// A reference behind an `Arc`, which is how the scanner takes one.
pub fn shared_reference(tag: &str, contig: &str, seq: &str) -> Arc<Aref> {
    Arc::new(reference(tag, contig, seq))
}

/// The tid map for a one-contig fixture, resolved the way the real path does.
///
/// Built through [`ContigMap::build`] rather than by hand so the fixtures
/// exercise the same resolution the scanner uses. With one contig named the
/// same on both sides the result is the identity, which is exactly why a test
/// that only ever uses one contig cannot catch an ordinal lookup -- see
/// `contig_map::tests::resolves_across_a_reordering` for the one that can.
pub fn contig_map(refr: &Aref, contig: &str, len: usize) -> Arc<ContigMap> {
    Arc::new(ContigMap::build(&header(contig, len), refr).unwrap())
}
/// Queries for a test, written as a query file.
///
/// The one way tests build queries. They used to use the terse one-liner
/// parser, which is no longer an input format: a query written here is written
/// the way a user writes one, so a test cannot pass against a syntax nobody
/// can type.
///
/// ```ignore
/// let specs = queries("[query.m]\nread = \"Y~\"\nrefr = \"CG\"\n");
/// ```
pub fn queries(toml: &str) -> Vec<crate::dsl::QuerySpec> {
    crate::query_toml::parse_file(toml)
        .unwrap_or_else(|e| panic!("test query file:\n{toml}\n{e}"))
        .queries
}

/// A single query, for the tests that render or explain exactly one.
pub fn query(toml: &str) -> crate::dsl::QuerySpec {
    let mut qs = queries(toml);
    assert_eq!(qs.len(), 1, "expected exactly one query:\n{toml}");
    qs.pop().expect("checked above")
}

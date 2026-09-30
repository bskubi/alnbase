//! Pin the output of the whole pipeline, byte for byte.
//!
//! The unit tests elsewhere check one decision each. This checks the thing a
//! user actually gets: take a reference and a BAM, run the real query
//! pipeline, and compare every tag and every extracted row against text
//! recorded here.
//!
//! It exists for refactoring. A deletion that was supposed to be dead code, a
//! parser swapped for another, a pipeline rearranged -- each is safe to do
//! quickly exactly when something would notice if the output moved. Unit tests
//! do not notice: they test the parts that were kept.
//!
//! # When it fails
//!
//! Either the change was not supposed to alter output, and this is the bug
//! report, or it was, and the fixture is regenerated **in the same commit** so
//! the diff is reviewable. Run with `ALNBASE_GOLDEN=show` to print the current
//! output in the form the constants below expect.
//!
//! # What the fixture covers
//!
//! One record per thing that has ever been got wrong: each of the four
//! conversion strands, an insertion, a deletion, a soft clip, a read ending
//! mid-context so the flank supplies the rest, a reference `N` (which must
//! reach `u`/`U` rather than a context code), and an unmapped record.

#![cfg(test)]

use std::path::Path;
use std::sync::Arc;

use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{Read, Reader, Record};

use crate::bam_out::{BamOptions, BamStream};
use crate::hit_writer::{OutputFormat, ParquetOpts};
use crate::ordered::{scan_ordered, OrderedConfig};
use crate::parallel::{scan_sharded, ShardConfig};
use crate::query::QuerySet;
use crate::query_toml;
use crate::record_field::RecordField;
use crate::scanner::{OutputConfig, WalkConfig};
use crate::tag_extract::{self, TagExtract};
use crate::test_support::{built, shared_reference, temp_path, write_bam, LAST_IN_TEMPLATE, REVERSE};

const CONTIG: &str = "chr1";

/// Positions           0         1         2         3
///                     0123456789012345678901234567890123456789
/// The `N` at 23 leaves the cytosine at 22 with an unknowable context, which
/// must land on `u`/`U` rather than on a context code.
const GENOME: &str = "TTTCGTTTCGTTTTACGTACCACNTAACGGTTAACCGTAA";

const UNMAPPED: u16 = 0x4;

/// Bismark's eight codes, including `u`/`U` as the complement of the three
/// determinable contexts. The same file shipped as `docs/bismark-xm.toml`.
const QUERIES: &str = r#"
[pattern.cg]
read = "~~~"
refr = "CG~"
[pattern.chg]
read = "~~~"
refr = "CHG"
[pattern.chh]
read = "~~~"
refr = "CHH"

[query.Z]
read = "C~~"
refr = "CG~"
[query.z]
read = "T~~"
refr = "CG~"
[query.X]
read = "C~~"
refr = "CHG"
[query.x]
read = "T~~"
refr = "CHG"
[query.H]
read = "C~~"
refr = "CHH"
[query.h]
read = "T~~"
refr = "CHH"

[query.U]
read  = "C~~"
refr  = "C~~"
where = "U and not cg and not chg and not chh"
[query.u]
read  = "T~~"
refr  = "C~~"
where = "u and not cg and not chg and not chh"

[tag.XM.bases]
fill = "."
z = "z"
Z = "Z"
x = "x"
X = "X"
h = "h"
H = "H"
u = "u"
U = "U"

[tag.XR.strand]
CT = ["OT", "OB"]
GA = ["CTOT", "CTOB"]

[tag.XG.strand]
CT = ["OT", "CTOT"]
GA = ["OB", "CTOB"]
"#;

/// `(name, seq, cigar, pos, flags)`, one per case named in the module comment.
fn records() -> Vec<Record> {
    let g = GENOME.as_bytes();
    let sub = |a: usize, b: usize| g[a..b].to_vec();
    // A top-strand read with the cytosines converted: every C becomes T.
    let converted = |a: usize, b: usize| -> Vec<u8> {
        sub(a, b).iter().map(|&c| if c == b'C' { b'T' } else { c }).collect()
    };

    vec![
        // Unconverted top strand: every cytosine protected.
        built(b"ot_protected", &sub(0, 14), &[Cigar::Match(14)], 0, 0),
        // Fully converted top strand.
        built(b"ot_converted", &converted(0, 14), &[Cigar::Match(14)], 0, 0),
        // The other three strands, same stretch, so the tags differ only by strand.
        built(b"ob", &sub(0, 14), &[Cigar::Match(14)], 0, REVERSE),
        built(b"ctot", &sub(0, 14), &[Cigar::Match(14)], 0, LAST_IN_TEMPLATE | REVERSE),
        built(b"ctob", &sub(0, 14), &[Cigar::Match(14)], 0, LAST_IN_TEMPLATE),
        // An insertion of two bases after 6 aligned: by default no column is
        // emitted for them, so a context can span the insertion.
        built(
            b"insertion",
            &[sub(0, 6), b"GG".to_vec(), sub(6, 12)].concat(),
            &[Cigar::Match(6), Cigar::Ins(2), Cigar::Match(6)],
            0,
            0,
        ),
        // A deletion of two reference bases.
        built(
            b"deletion",
            &[sub(0, 6), sub(8, 14)].concat(),
            &[Cigar::Match(6), Cigar::Del(2), Cigar::Match(6)],
            0,
            0,
        ),
        // Soft clip at each end: skipped by default, but read_off still counts them.
        built(
            b"softclip",
            &[b"AA".to_vec(), sub(2, 12), b"GG".to_vec()].concat(),
            &[Cigar::SoftClip(2), Cigar::Match(10), Cigar::SoftClip(2)],
            2,
            0,
        ),
        // Ends *on* the cytosine at 3, so both context bases come from the
        // flank past the end of the read. Without the flank it is uncallable.
        built(b"ends_on_c", &sub(0, 4), &[Cigar::Match(4)], 0, 0),
        // Spans the N at 24: the cytosines at 22 and 23 are unknown-context.
        built(b"near_n", &sub(18, 30), &[Cigar::Match(12)], 18, 0),
        // Written through untagged.
        built(b"unmapped", b"ACGTACGTACGT", &[Cigar::Match(12)], 0, UNMAPPED),
    ]
}

fn render_bam(path: &Path) -> String {
    let mut reader = Reader::from_path(path).unwrap();
    let mut out = String::new();
    for r in reader.records() {
        let r = r.unwrap();
        let tag = |t: &[u8]| match r.aux(t) {
            Ok(Aux::String(v)) => v.to_string(),
            _ => "-".to_string(),
        };
        out.push_str(&format!(
            "{:<13} {:<14} XM={:<14} XR={} XG={}\n",
            String::from_utf8_lossy(r.qname()),
            String::from_utf8_lossy(&r.seq().as_bytes()),
            tag(b"XM"),
            tag(b"XR"),
            tag(b"XG"),
        ));
    }
    out
}

fn render_parquet(path: &Path) -> String {
    use arrow::util::display::{ArrayFormatter, FormatOptions};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path).unwrap();
    let reader = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap();
    let mut rows: Vec<String> = Vec::new();
    for batch in reader {
        let batch = batch.unwrap();
        let opts = FormatOptions::default().with_null("-");
        let cols: Vec<(String, ArrayFormatter)> = batch
            .schema()
            .fields()
            .iter()
            .zip(batch.columns())
            .map(|(f, c)| (f.name().clone(), ArrayFormatter::try_new(c.as_ref(), &opts).unwrap()))
            .collect();
        for r in 0..batch.num_rows() {
            rows.push(
                cols.iter()
                    .map(|(n, f)| format!("{n}={}", f.value(r)))
                    .collect::<Vec<_>>()
                    .join(" "),
            );
        }
    }
    rows.join("\n") + "\n"
}

/// Run the real pipelines over the fixture: tag a BAM, then extract from it.
fn run() -> (String, String) {
    let tag = "golden";
    let file = query_toml::parse_file(QUERIES).unwrap_or_else(|e| panic!("{e}"));
    let queries = Arc::new(QuerySet::compile(&file.queries).unwrap());
    let recs = records();
    let bam = write_bam(tag, CONTIG, GENOME.len(), &recs);
    let refr = shared_reference(tag, CONTIG, GENOME);

    let tagged = temp_path("golden_tagged", "bam");
    let stream = BamStream::new(
        &file.tags,
        &queries,
        BamOptions {
            walk: WalkConfig::default(),
            strand: crate::strand_rule::directional(),
            overwrite_tags: false,
            bad_data: crate::batch::OnBadData::Stop,
            // Fixed, so the @PG line cannot drift with the test runner's argv.
            command_line: "alnbase query --query-file bismark-xm.toml in.bam ref.aref out.bam"
                .to_string(),
            query_sources: Vec::new(),
        },
    )
    .unwrap();
    scan_ordered(
        &bam,
        refr,
        Arc::clone(&queries),
        &OrderedConfig {
            n_workers: 2,
            reader_threads: 1,
            writer_threads: 1,
            out: tagged.clone(), bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        },
        &stream,
    )
    .unwrap();

    // One worker and one shard, so the row order is the record order rather
    // than a hash of the partition key.
    let header = String::from_utf8_lossy(Reader::from_path(&tagged).unwrap().header().as_bytes())
        .into_owned();
    let choice = tag_extract::choose(&header, Some((file.queries.clone(), file.tags.clone())), None)
        .unwrap_or_else(|e| panic!("{e}"));
    let chosen = Arc::new(QuerySet::compile(&choice.queries).unwrap());
    let extract = TagExtract::new(
        &choice,
        &chosen,
        OutputConfig {
            bad_data: crate::batch::OnBadData::Stop,
            fields: vec![RecordField::RefName, RecordField::RefrStrand],
            batch_rows: 64,
            parquet: ParquetOpts::default(),
            keep_hitless: false,
            format: OutputFormat::Parquet,
            walk: WalkConfig::default(),
            strand: crate::strand_rule::directional(),
        },
        crate::strand_rule::directional(),
    )
    .unwrap();
    let shard = ShardConfig {
        strand: crate::strand_rule::directional(),
        n_workers: 1,
        shards_per_worker: 1,
        partition_by: vec![RecordField::Qname],
        reader_threads: 1,
        out_prefix: temp_path("golden_rows", "parquet"),
        bad_data: crate::batch::OnBadData::Stop,
        require_m5: false,
    };
    scan_sharded(&tagged, None, chosen, &shard, &extract).unwrap();

    (render_bam(&tagged), render_parquet(&shard.out_path(0, 0)))
}

#[test]
fn golden_output_is_unchanged() {
    let (bam, rows) = run();
    if std::env::var("ALNBASE_GOLDEN").as_deref() == Ok("show") {
        panic!("\n--- BAM ---\n{bam}\n--- ROWS ---\n{rows}");
    }
    assert_eq!(bam, GOLDEN_BAM, "\nthe tagged BAM changed; see the module comment");
    assert_eq!(rows, GOLDEN_ROWS, "\nthe extracted rows changed; see the module comment");
}

/// Every `src/*.rs` is declared in `main.rs`.
///
/// An undeclared file is compiled by nothing, linted by nothing and tested by
/// nothing, so it is invisible to every other check in the build --
/// `record_view.rs` sat there for months. This is the only thing that notices.
#[test]
fn every_source_file_is_part_of_the_crate() {
    let main = std::fs::read_to_string("src/main.rs").expect("run from the crate root");
    let mut orphans: Vec<String> = Vec::new();
    for entry in std::fs::read_dir("src").unwrap() {
        let path = entry.unwrap().path();
        let Some(stem) = path.file_stem().and_then(|s| s.to_str()) else { continue };
        if path.extension().and_then(|e| e.to_str()) != Some("rs") || stem == "main" {
            continue;
        }
        // Line-based: a `mod` on the first line of the file has no newline
        // before it, which a substring search would miss.
        let declared = main
            .lines()
            .any(|l| l.trim() == format!("mod {stem};") || l.trim() == format!("pub mod {stem};"));
        if !declared {
            orphans.push(format!("src/{stem}.rs"));
        }
    }
    assert!(orphans.is_empty(), "not declared in main.rs, so never compiled: {orphans:?}");
}

const GOLDEN_BAM: &str = "\
ot_protected  TTTCGTTTCGTTTT XM=...Z....Z..... XR=CT XG=CT
\
ot_converted  TTTTGTTTTGTTTT XM=...z....z..... XR=CT XG=CT
\
ob            TTTCGTTTCGTTTT XM=....Z....Z.... XR=CT XG=GA
\
ctot          TTTCGTTTCGTTTT XM=...Z....Z..... XR=GA XG=CT
\
ctob          TTTCGTTTCGTTTT XM=....Z....Z.... XR=GA XG=GA
\
insertion     TTTCGTGGTTCGTT XM=...Z......Z... XR=CT XG=CT
\
deletion      TTTCGTCGTTTT   XM=...Z..Z.....   XR=CT XG=CT
\
softclip      AATCGTTTCGTTGG XM=...Z....Z..... XR=CT XG=CT
\
ends_on_c     TTTC           XM=...Z           XR=CT XG=CT
\
near_n        ACCACNTAACGG   XM=.HH.U....Z..   XR=CT XG=CT
\
unmapped      ACGTACGTACGT   XM=-              XR=- XG=-
\
";
const GOLDEN_ROWS: &str = "\
shard=0 record_id=0 ref_name=chr1 strand=+ name=Z read_5p=3 read_3p=10 refr_pos=3 qual=37 read_base=C refr_base=C
\
shard=0 record_id=0 ref_name=chr1 strand=+ name=Z read_5p=8 read_3p=5 refr_pos=8 qual=37 read_base=C refr_base=C
\
shard=0 record_id=1 ref_name=chr1 strand=+ name=z read_5p=3 read_3p=10 refr_pos=3 qual=37 read_base=T refr_base=C
\
shard=0 record_id=1 ref_name=chr1 strand=+ name=z read_5p=8 read_3p=5 refr_pos=8 qual=37 read_base=T refr_base=C
\
shard=0 record_id=2 ref_name=chr1 strand=- name=Z read_5p=4 read_3p=9 refr_pos=9 qual=37 read_base=C refr_base=C
\
shard=0 record_id=2 ref_name=chr1 strand=- name=Z read_5p=9 read_3p=4 refr_pos=4 qual=37 read_base=C refr_base=C
\
shard=0 record_id=3 ref_name=chr1 strand=+ name=Z read_5p=10 read_3p=3 refr_pos=3 qual=37 read_base=C refr_base=C
\
shard=0 record_id=3 ref_name=chr1 strand=+ name=Z read_5p=5 read_3p=8 refr_pos=8 qual=37 read_base=C refr_base=C
\
shard=0 record_id=4 ref_name=chr1 strand=- name=Z read_5p=9 read_3p=4 refr_pos=9 qual=37 read_base=C refr_base=C
\
shard=0 record_id=4 ref_name=chr1 strand=- name=Z read_5p=4 read_3p=9 refr_pos=4 qual=37 read_base=C refr_base=C
\
shard=0 record_id=5 ref_name=chr1 strand=+ name=Z read_5p=3 read_3p=10 refr_pos=3 qual=37 read_base=C refr_base=C
\
shard=0 record_id=5 ref_name=chr1 strand=+ name=Z read_5p=10 read_3p=3 refr_pos=8 qual=37 read_base=C refr_base=C
\
shard=0 record_id=6 ref_name=chr1 strand=+ name=Z read_5p=3 read_3p=8 refr_pos=3 qual=37 read_base=C refr_base=C
\
shard=0 record_id=6 ref_name=chr1 strand=+ name=Z read_5p=6 read_3p=5 refr_pos=8 qual=37 read_base=C refr_base=C
\
shard=0 record_id=7 ref_name=chr1 strand=+ name=Z read_5p=3 read_3p=10 refr_pos=3 qual=37 read_base=C refr_base=C
\
shard=0 record_id=7 ref_name=chr1 strand=+ name=Z read_5p=8 read_3p=5 refr_pos=8 qual=37 read_base=C refr_base=C
\
shard=0 record_id=8 ref_name=chr1 strand=+ name=Z read_5p=3 read_3p=0 refr_pos=3 qual=37 read_base=C refr_base=C
\
shard=0 record_id=9 ref_name=chr1 strand=+ name=H read_5p=1 read_3p=10 refr_pos=19 qual=37 read_base=C refr_base=C
\
shard=0 record_id=9 ref_name=chr1 strand=+ name=H read_5p=2 read_3p=9 refr_pos=20 qual=37 read_base=C refr_base=C
\
shard=0 record_id=9 ref_name=chr1 strand=+ name=U read_5p=4 read_3p=7 refr_pos=22 qual=37 read_base=C refr_base=C
\
shard=0 record_id=9 ref_name=chr1 strand=+ name=Z read_5p=9 read_3p=2 refr_pos=27 qual=37 read_base=C refr_base=C
\
";

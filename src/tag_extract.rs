//! Read parquet hits back out of a `bases` tag.
//!
//! A tagged BAM already holds its calls: one character per base of SEQ. Given
//! the tag's definition -- which query each character stands for -- those
//! calls can be written as the same hit rows a direct parquet run writes,
//! without a reference and without walking a single alignment.
//!
//! # Which definition
//!
//! By default, the query file alnbase stored in the BAM's header when it wrote
//! the tag ([`crate::header_query`]). When the BAM has been through alnbase more
//! than once, the most recent run that defines the tag is used. When no tag is
//! named, the one `bases` tag the stored runs define is used, and naming one
//! is required only when there are several.
//!
//! Query files given on the command line replace the stored ones. That is also
//! how a tag written by another tool -- Bismark's `XM` -- is read: a file that
//! maps its characters to queries is all the definition needs.
//!
//! # What a row holds
//!
//! One row per character that is not the fill, in the order a walk would have
//! found it:
//!
//! - `name` of the query the character stands for.
//! - `off_5p` / `off_3p` from the character's position in SEQ, counted from
//!   the ends of the read as sequenced, exactly as a direct `query --parquet`
//!   run reports them.
//! - `refr_pos` of that base, by the walk's conventions: an inserted or
//!   soft-clipped base carries the nearest 5' reference coordinate.
//! - `qual`, and the read base as the flat capture, complemented on the bottom
//!   strands as a walk reports it.
//! - The reference base when the query pins it down -- `T~~@CHH` can only
//!   fire over a C, and `=` means the read's own base -- and null when it does
//!   not, since no reference is read; see [`QuerySpec::anchor_refr`].
//!
//! Captures are always flat: a tag records the anchor and nothing else a query
//! may have captured.
//!
//! Rows come only from the queries the tag's characters stand for; a direct
//! run reports every query. Each query has its own character, and two
//! characters cannot mark one base, so there is never more than one row per
//! base.

use std::path::Path;
use std::sync::Arc;

use anyhow::{bail, Result};
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{HeaderView, Record};
use rust_htslib::htslib;

use crate::contig_map::ContigMap;
use crate::dsl::{AnchorRefr, QuerySpec};
use crate::header_query::{self, QuerySource};
use crate::hit_writer::HitWriter;
use crate::hits::HitBuilder;
use crate::parallel::{ShardSink, SinkCounts, SinkFactory};
use crate::query::QuerySet;
use crate::aref::Aref;
use crate::scanner::OutputConfig;
use crate::seq::Seq;
use crate::tags::{BasesTag, TagConfig, TagKind, TagName};
use crate::strand_rule::StrandRule;

/// A tag's definition, chosen from the stored runs or the given files.
#[derive(Debug, Clone)]
pub struct Choice {
    pub tag: TagName,
    pub bases: BasesTag,
    /// The queries of the files that define the tag.
    pub queries: Vec<QuerySpec>,
    /// Where the definition came from, for telling the user.
    pub source: String,
    /// The text of the query files the definition came from, when it came from
    /// the BAM header. Empty for query files given on the command line, whose
    /// text the caller already holds.
    pub stored_files: Vec<QuerySource>,
}

/// Choose which tag to extract and with which definition.
///
/// `given` is the queries and tags of query files named on the command line;
/// when present, the header is not consulted. `want` is `--tag`.
/// Every tag a BAM's stored runs define, each resolved to its most recent
/// definition.
///
/// Tags are reconciled one at a time, not run at a time. A BAM tagged twice
/// may have `XM` redefined by the second run while `XE` is only defined by the
/// first; walking the runs newest-first and taking each tag from the first run
/// that defines it gives the right answer for both, where picking one run
/// would be wrong for one of them.
type Defined = (TagName, TagKind, Vec<QuerySpec>, String, Vec<QuerySource>);

fn reconcile(header_text: &str) -> Result<Vec<Defined>, String> {
    let mut runs = header_query::stored_runs(header_text)?;
    runs.reverse();
    let mut out: Vec<Defined> = Vec::new();
    for run in runs {
        let (files, tags) = run.resolve()?;
        let paths: Vec<&str> = run.sources.iter().map(|s| s.path.as_str()).collect();
        let desc = format!("the query file stored under @PG {} ({})", run.pg_id, paths.join(", "));
        let queries: Vec<QuerySpec> = files.into_iter().flat_map(|f| f.queries).collect();
        for t in &tags.tags {
            if out.iter().any(|(n, ..)| *n == t.name) {
                continue; // a newer run already defined it
            }
            out.push((t.name, t.kind.clone(), queries.clone(), desc.clone(), run.sources.clone()));
        }
    }
    Ok(out)
}

/// Choose which tag to extract and with which definition.
///
/// `given` is the queries and tags of query files named on the command line;
/// when present, the header is not consulted. `want` is `--tag`.
pub fn choose(
    header_text: &str,
    given: Option<(Vec<QuerySpec>, TagConfig)>,
    want: Option<TagName>,
) -> Result<Choice, String> {
    let defined: Vec<Defined> = match given {
        Some((queries, tags)) => {
            // Files that declare only a strand rule, say, are a plausible
            // command line and define nothing to extract; the header is not
            // consulted once --query-file is given, so there is nothing else
            // to fall back on.
            if tags.tags.is_empty() {
                return Err("the given query files declare no tag to extract; a file whose \
                            [tag.XX.bases] table maps the tag's characters to queries is its \
                            definition"
                    .into());
            }
            tags.tags
                .iter()
                .map(|t| {
                    (t.name, t.kind.clone(), queries.clone(), "the given query files".to_string(), Vec::new())
                })
                .collect()
        }
        None => {
            let d = reconcile(header_text)?;
            if d.is_empty() {
                return Err("the BAM has no stored alnbase query file; pass --query-file to say \
                            how to read its tags"
                    .into());
            }
            d
        }
    };

    let name = match want {
        Some(name) => name,
        None => {
            let bases: Vec<TagName> = defined
                .iter()
                .filter(|(_, k, ..)| matches!(k, TagKind::Bases(_)))
                .map(|(n, ..)| *n)
                .collect();
            match bases.as_slice() {
                [one] => *one,
                [] => {
                    let src = &defined.first().expect("non-empty: both arms above are").3;
                    return Err(format!("{src} defines no bases tag"));
                }
                many => {
                    let list: Vec<String> = many.iter().map(|n| n.to_string()).collect();
                    return Err(format!(
                        "several bases tags are defined ({}); choose one with --tag",
                        list.join(", ")
                    ));
                }
            }
        }
    };

    match defined.iter().find(|(n, ..)| *n == name) {
        Some((_, TagKind::Bases(b), queries, desc, stored)) => Ok(Choice {
            tag: name,
            bases: b.clone(),
            queries: queries.clone(),
            source: desc.clone(),
            stored_files: stored.clone(),
        }),
        Some((_, TagKind::Strand(_), _, desc, _)) => Err(format!(
            "{name} is a strand tag in {desc}: it holds one value per read, not calls to extract"
        )),
        None => {
            let known: Vec<String> = defined.iter().map(|(n, ..)| n.to_string()).collect();
            Err(format!(
                "no definition of tag {name}; the BAM defines {}",
                if known.is_empty() { "none".to_string() } else { known.join(", ") }
            ))
        }
    }
}


/// What one character of the tag writes.
#[derive(Debug, Clone)]
struct Label {
    name: String,
    /// The reference symbol under the anchor, when the query fixes it.
    refr: AnchorRefr,
}

/// Makes a [`TagExtractShard`] per output file.
pub struct TagExtract {
    tag: [u8; 2],
    tag_name: TagName,
    fill: u8,
    labels: Arc<Vec<Option<Label>>>,
    cfg: OutputConfig,
    strand: Arc<StrandRule>,
}

impl TagExtract {
    /// `queries` is `choice.queries` compiled, which supplies each query's
    /// base-level expression.
    pub fn new(choice: &Choice, queries: &QuerySet, cfg: OutputConfig, strand: Arc<StrandRule>) -> Result<Self> {
        let mut labels: Vec<Option<Label>> = vec![None; 256];
        for code in &choice.bases.codes {
            let q = &code.query;
            let Some(compiled) = queries.queries.iter().find(|c| &c.name == q) else {
                bail!("tag {}: code '{}' names no query '{q}'", choice.tag, code.code as char);
            };
            let refr = choice
                .queries
                .iter()
                .find(|s| &s.name == q)
                .map(|s| s.anchor_refr())
                .unwrap_or(AnchorRefr::Unknown);
            labels[code.code as usize] =
                Some(Label { name: compiled.name.clone(), refr });
        }
        Ok(Self {
            tag: choice.tag.bytes(),
            tag_name: choice.tag,
            fill: choice.bases.fill,
            labels: Arc::new(labels),
            cfg,
            strand,
        })
    }
}

impl SinkFactory for TagExtract {
    type Sink<'q> = TagExtractShard;

    fn keeps_unscanned(&self) -> bool {
        true
    }

    /// The position of every base, its quality and its read base, the tag
    /// itself, and the columns this output writes.
    fn cram_decoding(&self) -> Option<(u32, bool)> {
        let fields = &self.cfg.fields;
        let own = htslib::sam_fields_SAM_FLAG
            | htslib::sam_fields_SAM_RNAME
            | htslib::sam_fields_SAM_POS
            | htslib::sam_fields_SAM_CIGAR
            | htslib::sam_fields_SAM_SEQ
            | htslib::sam_fields_SAM_QUAL
            | htslib::sam_fields_SAM_AUX;
        let required = fields.iter().fold(own, |m, f| m | f.sam_fields());
        Some((required, fields.iter().any(|f| f.needs_md_nm())))
    }

    fn open<'q>(
        &self,
        _queries: &'q QuerySet,
        _refr: Option<Arc<Aref>>,
        _contigs: Option<Arc<ContigMap>>,
        _header: &HeaderView,
        path: &Path,
        slot: u32,
    ) -> Result<TagExtractShard> {
        let cap = self.cfg.batch_rows.min(4096);
        let builder = HitBuilder::new(&self.cfg.fields, cap, true, slot);
        let writer =
            HitWriter::new(path, builder, self.cfg.batch_rows, self.cfg.format, &self.cfg.parquet)?;
        Ok(TagExtractShard {
            writer,
            tag: self.tag,
            tag_name: self.tag_name,
            fill: self.fill,
            labels: Arc::clone(&self.labels),
            keep_hitless: self.cfg.keep_hitless,
            strand: self.strand.clone(),
            positions: Vec::new(),
            scanned: 0,
            unknown_strand: 0,
        })
    }
}

/// One output file of decoded hits.
pub struct TagExtractShard {
    writer: HitWriter<HitBuilder>,
    tag: [u8; 2],
    tag_name: TagName,
    fill: u8,
    labels: Arc<Vec<Option<Label>>>,
    keep_hitless: bool,
    strand: Arc<StrandRule>,
    /// Reference position of each SEQ base, for the current record. Reused.
    positions: Vec<Option<i64>>,
    scanned: u64,
    /// Records the run's strand rule declined to call; see [`ShardSink::scan`].
    unknown_strand: u64,
}

impl ShardSink for TagExtractShard {
    fn scan(&mut self, rid: u64, record: &mut Record, header: &HeaderView) -> Result<()> {
        // Decoding a tag needs the same call the tag was written with: the
        // read's own 5' end is at the far end of SEQ on a minus-strand walk.
        // A record the rule declines is skipped and counted, as in the scan.
        let Some(call) = self.strand.call(record)? else {
            self.unknown_strand += 1;
            return Ok(());
        };
        let builder = self.writer.builder_mut();
        builder.begin_record(rid, record, header, Some(call));
        let len = record.seq_len();
        let qname = || String::from_utf8_lossy(record.qname()).into_owned();

        let value = match record.aux(&self.tag) {
            Ok(Aux::String(v)) => Some(v.as_bytes()),
            Ok(_) => bail!("record {}: tag {} is not a string", qname(), self.tag_name),
            Err(_) => None,
        };

        let mut any = false;
        if let Some(v) = value {
            if v.len() != len {
                bail!(
                    "record {}: tag {} has {} characters, but the read has {len} bases",
                    qname(),
                    self.tag_name,
                    v.len()
                );
            }
            let bottom = call.walk_reversed();
            ref_positions(record, bottom, &mut self.positions);
            let seq = record.seq();
            let quals = record.qual();
            // Emitted order: along the strand, as a walk would find them.
            for k in 0..len {
                let i = if bottom { len - 1 - k } else { k };
                let c = v[i];
                if c == self.fill {
                    continue;
                }
                let Some(label) = &self.labels[c as usize] else {
                    bail!(
                        "record {}: tag {} has '{}' at SEQ position {i} (0-based), which its definition \
                         does not include",
                        qname(),
                        self.tag_name,
                        c as char
                    );
                };
                let base = Seq(seq.encoded_base(i));
                let base = if bottom { base.complement() } else { base };
                let read = base.name();
                let refr = match label.refr {
                    AnchorRefr::Fixed(s) => Some(s.name()),
                    AnchorRefr::SameAsRead => Some(read.clone()),
                    AnchorRefr::Unknown => None,
                };
                builder.push_decoded(
                    &label.name,
                    k as i64,
                    self.positions[i],
                    quals.get(i).copied().filter(|&q| q != 0xff), // 0xff: QUAL is '*'
                    &read,
                    refr.as_deref(),
                );
                any = true;
            }
        }
        if !any && self.keep_hitless {
            builder.push_hitless();
        }
        self.writer.maybe_flush()?;
        self.scanned += 1;
        Ok(())
    }

    /// A record with no tag to decode: the same all-null row the parquet
    /// output writes, so both accounts of a BAM have the same shape.
    fn pass(&mut self, rid: u64, record: &mut Record, header: &HeaderView) -> Result<()> {
        if !self.keep_hitless {
            return Ok(());
        }
        // Not held to the rule, for the reason `Scanner::pass_through` gives.
        let builder = self.writer.builder_mut();
        builder.begin_record(rid, record, header, self.strand.call(record).unwrap_or(None));
        builder.push_hitless();
        self.writer.maybe_flush()
    }

    fn finish(self) -> Result<SinkCounts> {
        let rows = self.writer.close()?;
        Ok(SinkCounts {
            scanned: self.scanned,
            rows,
            unknown_strand: self.unknown_strand,
            ..SinkCounts::default()
        })
    }
}

/// The reference coordinate of every base of SEQ, by the walk's conventions:
/// a base aligned to the reference has its own coordinate; an inserted or
/// soft-clipped base has none and carries the nearest 5' one along the strand
/// -- the coordinate before it on the top strands, after it on the bottom.
///
/// Held to [`crate::alignment::walk_alignment`] by
/// `tests::positions_agree_with_the_walk`.
fn ref_positions(record: &Record, bottom: bool, out: &mut Vec<Option<i64>>) {
    let len = record.seq_len();
    out.clear();
    out.resize(len, None);
    let mut r = record.pos();
    let mut q = 0usize;
    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for _ in 0..n {
                    if q < len {
                        out[q] = Some(r);
                    }
                    q += 1;
                    r += 1;
                }
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => {
                let nearest = if bottom { r } else { r - 1 };
                for _ in 0..n {
                    if q < len {
                        out[q] = Some(nearest);
                    }
                    q += 1;
                }
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => r += n as i64,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::tags::Library;
    use super::*;
    use crate::alignment::{walk_alignment, Insertions};
    use crate::header_query::{header_lines, QuerySource};
    use crate::hit_writer::{OutputFormat, ParquetOpts};
    use crate::query_toml;
    use crate::record_field::RecordField;
    use crate::scanner::WalkConfig;
    use crate::test_support::{built, shared_reference, write_bam, LAST_IN_TEMPLATE, REVERSE};
    use rust_htslib::bam::{Read, Reader};

    const CONTIG: &str = "chr1";
    const GENOME: &str = "TTTCGTTTCGTTTTACGTACGTACCGGTTAAC";

    #[test]
    fn positions_agree_with_the_walk() {
        let refr = shared_reference("extract_pos", CONTIG, GENOME);
        let bam = write_bam("extract_pos", CONTIG, GENOME.len(), &[built(b"x", b"T", &[Cigar::Match(1)], 0, 0)]);
        let reader = Reader::from_path(&bam).unwrap();
        let contigs = ContigMap::build(reader.header(), &refr).unwrap();
        let walk = WalkConfig { insertions: Insertions::Emit, ..WalkConfig::default() };

        let cigars: Vec<Vec<Cigar>> = vec![
            vec![Cigar::Match(10)],
            vec![Cigar::SoftClip(2), Cigar::Match(6), Cigar::SoftClip(3)],
            vec![Cigar::Match(3), Cigar::Ins(2), Cigar::Match(4)],
            vec![Cigar::Match(3), Cigar::Del(2), Cigar::Match(4)],
            vec![Cigar::SoftClip(1), Cigar::Match(2), Cigar::Ins(1), Cigar::Del(1), Cigar::Match(3), Cigar::SoftClip(2)],
            vec![Cigar::HardClip(4), Cigar::Match(4), Cigar::RefSkip(5), Cigar::Match(3), Cigar::HardClip(2)],
            vec![Cigar::Ins(2), Cigar::Match(5), Cigar::Ins(1)],
        ];
        for cigar in &cigars {
            let len: u32 = cigar
                .iter()
                .map(|c| match c {
                    Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) => *n,
                    _ => 0,
                })
                .sum();
            let seq: Vec<u8> = (0..len).map(|i| b"ACGT"[i as usize % 4]).collect();
            for flags in [0, REVERSE, LAST_IN_TEMPLATE, LAST_IN_TEMPLATE | REVERSE] {
                let rec = built(b"x", &seq, cigar, 4, flags);
                let bottom = Library::Directional.strand(&rec).is_bottom();
                let mut mine = Vec::new();
                ref_positions(&rec, bottom, &mut mine);

                let mut cols = Vec::new();
                let mut read_seq = Vec::new();
                walk_alignment(&rec, &refr, &contigs, &mut read_seq, walk.to_opts(1), Library::Directional.call(&rec), &mut |c| cols.push(c))
                    .unwrap();
                let mut checked = 0;
                for c in cols.iter().filter(|c| c.qual >= 0) {
                    let i = if bottom { len as i64 - 1 - c.read_off } else { c.read_off } as usize;
                    assert_eq!(mine[i], Some(c.refr_pos), "{cigar:?} flags {flags:#x} SEQ {i}");
                    checked += 1;
                }
                // Soft-clipped bases produce no columns, so every other base is checked.
                let clipped: u32 = cigar.iter().map(|c| if let Cigar::SoftClip(n) = c { *n } else { 0 }).sum();
                assert_eq!(checked, len - clipped, "{cigar:?}: every aligned base checked");
            }
        }
    }

    const TWO_TAGS: &str = r#"
[query.TG]
mark = "+."
where = "tg"
[query.CG]
mark = "+."
where = "cg"
[pattern.tg]
read = "T~"
refr = "CG"
[pattern.cg]
read = "C~"
refr = "CG"
[tag.XM.bases]
fill = "."
z = "TG"
Z = "CG"
[tag.XR.strand]
CT = ["OT", "OB"]
GA = ["CTOT", "CTOB"]
"#;

    fn header_of(runs: &[(&str, &str)]) -> String {
        let mut h = String::from("@HD\tVN:1.6\n");
        for (id, _) in runs {
            h.push_str(&format!("@PG\tID:{id}\tPN:alnbase\n"));
        }
        for (id, text) in runs {
            h.push_str(&header_lines(id, &[QuerySource { path: format!("{id}.toml"), text: text.to_string() }]));
        }
        h
    }

    #[test]
    fn the_one_chars_tag_is_taken_without_asking() {
        let c = choose(&header_of(&[("alnbase", TWO_TAGS)]), None, None).unwrap();
        assert_eq!(c.tag.to_string(), "XM");
        assert!(c.source.contains("@PG alnbase"), "{}", c.source);
    }

    #[test]
    fn the_most_recent_run_defining_the_tag_wins() {
        let newer = TWO_TAGS.replace("fill = \".\"", "fill = \"-\"");
        let other = "[query.q]\n[pattern.p]\nread = \"C\"\nrefr = \"C\"\n[tag.XE.bases]\nfill = \".\"\ne = \"q\"\n";
        let h = header_of(&[("alnbase", TWO_TAGS), ("alnbase.1", &newer), ("alnbase.2", other)]);

        let e = choose(&h, None, None).unwrap_err();
        assert!(e.contains("several bases tags") && e.contains("XM") && e.contains("XE"), "{e}");

        let c = choose(&h, None, Some(TagName::parse("XM").unwrap())).unwrap();
        assert_eq!(c.bases.fill, b'-', "the newer of the two XM runs");
        assert!(c.source.contains("alnbase.1"), "{}", c.source);

        let c = choose(&h, None, Some(TagName::parse("XE").unwrap())).unwrap();
        assert!(c.source.contains("alnbase.2"));

        let e = choose(&h, None, Some(TagName::parse("XR").unwrap())).unwrap_err();
        assert!(e.contains("strand tag"), "{e}");
        let e = choose(&h, None, Some(TagName::parse("YY").unwrap())).unwrap_err();
        assert!(e.contains("no definition of tag YY"), "{e}");
    }

    #[test]
    fn given_files_replace_the_stored_ones() {
        let h = header_of(&[("alnbase", TWO_TAGS)]);
        let given = query_toml::parse_file(&TWO_TAGS.replace("fill = \".\"", "fill = \"_\"")).unwrap();
        let c = choose(&h, Some((given.queries, given.tags)), None).unwrap();
        assert_eq!(c.bases.fill, b'_');
        assert_eq!(c.source, "the given query files");

        let e = choose("@HD\tVN:1.6\n", None, None).unwrap_err();
        assert!(e.contains("no stored alnbase query file"), "{e}");
    }

    fn parquet_cfg() -> OutputConfig {
        OutputConfig {
            bad_data: crate::batch::OnBadData::Stop,
            fields: vec![RecordField::Qname, RecordField::ConvStrand],
            batch_rows: 64,
            parquet: ParquetOpts::default(),
            keep_hitless: false,
            format: OutputFormat::Parquet,
            walk: WalkConfig::default(),
            strand: crate::strand_rule::directional(),
        }
    }

    /// Every row of every file, as text, sorted: order across shards is not
    /// part of either output's contract.
    fn rows(prefix: &Path, workers: usize, shards: usize) -> Vec<Vec<String>> {
        use arrow::util::display::{ArrayFormatter, FormatOptions};
        use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
        let mut out = Vec::new();
        for w in 0..workers {
            for s in 0..shards {
                let path = crate::parallel::ShardConfig {
                    strand: crate::strand_rule::directional(),
                    n_workers: workers,
                    shards_per_worker: shards,
                    partition_by: vec![RecordField::Qname],
                    reader_threads: 1,
                    out_prefix: prefix.to_path_buf(), bad_data: crate::batch::OnBadData::Stop,
                    require_m5: false,
                }
                .out_path(w, s);
                let file = std::fs::File::open(&path).unwrap();
                for batch in ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap() {
                    let batch = batch.unwrap();
                    let opts = FormatOptions::default().with_null("<null>");
                    let cols: Vec<(String, ArrayFormatter)> = batch
                        .schema()
                        .fields()
                        .iter()
                        .zip(batch.columns())
                        .filter(|(f, _)| f.name() != "record_id" && f.name() != "shard")
                        .map(|(f, c)| (f.name().clone(), ArrayFormatter::try_new(c.as_ref(), &opts).unwrap()))
                        .collect();
                    for r in 0..batch.num_rows() {
                        out.push(cols.iter().map(|(n, f)| format!("{n}={}", f.value(r))).collect());
                    }
                }
            }
        }
        out.sort();
        out
    }

    /// The same reads, two ways: a direct parquet run, and a BAM run whose tag
    /// is then extracted using the query file stored in its header. The rows
    /// agree in every column -- offsets, positions, the conversion strand, and
    /// the reference base, which these queries fix -- for all four strands and
    /// across an insertion, a deletion, a soft clip and a hard clip.
    #[test]
    fn extraction_matches_a_direct_parquet_run() {
        use crate::bam_out::{BamOptions, BamStream};
        use crate::ordered::{scan_ordered, OrderedConfig};
        use crate::parallel::{scan_sharded, ShardConfig};
        use crate::scanner::ParquetOutput;
        use crate::test_support::temp_path;

        let file = query_toml::parse_file(TWO_TAGS).unwrap();
        let queries = Arc::new(QuerySet::compile(&file.queries).unwrap());
        let m = |n| [Cigar::Match(n)];
        let recs = vec![
            built(b"ot", b"TTTTGTTTCGTTTTACG", &m(17), 0, 0),
            built(b"ob", b"TTTCATTTCGTTTTACA", &m(17), 0, REVERSE),
            built(b"ctot", b"GTTTTACGTACGTA", &m(14), 7, LAST_IN_TEMPLATE | REVERSE),
            built(b"ctob", b"TTCATTTTACGTAC", &m(14), 2, LAST_IN_TEMPLATE),
            built(b"ins", b"TTTTGTAATTTCG", &[Cigar::Match(6), Cigar::Ins(2), Cigar::Match(5)], 0, 0),
            built(b"clip", b"AATTTCGTTTTG", &[Cigar::SoftClip(2), Cigar::Match(10)], 0, 0),
            built(b"del", b"TTTTGTCGTTTTA", &[Cigar::Match(5), Cigar::Del(2), Cigar::Match(8)], 0, 0),
            // A supplementary alignment: offsets count the hard-clipped bases.
            built(b"hard", b"TTTCGTTTTA", &[Cigar::HardClip(5), Cigar::Match(10)], 5, REVERSE | 0x800),
            built(b"none", b"TTTTTTTT", &m(8), 20, 0),
        ];
        let refr = shared_reference("extract_eq", CONTIG, GENOME);
        let bam = write_bam("extract_eq", CONTIG, GENOME.len(), &recs);
        let shard = |tag: &str| ShardConfig {
            strand: crate::strand_rule::directional(),
            n_workers: 2,
            shards_per_worker: 2,
            partition_by: vec![RecordField::Qname],
            reader_threads: 1,
            out_prefix: temp_path(tag, "parquet"),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };

        let direct = shard("extract_eq_direct");
        scan_sharded(&bam, Some(Arc::clone(&refr)), Arc::clone(&queries), &direct, &ParquetOutput::new(parquet_cfg()))
            .unwrap();

        let tagged = temp_path("extract_eq_tagged", "bam");
        let stream = BamStream::new(
            &file.tags,
            &queries,
            BamOptions {
                walk: WalkConfig::default(),
                strand: crate::strand_rule::directional(),
                overwrite_tags: false,
                bad_data: crate::batch::OnBadData::Stop,
                command_line: "alnbase query".into(),
                query_sources: vec![QuerySource { path: "two.toml".into(), text: TWO_TAGS.into() }],
            },
        )
        .unwrap();
        let ocfg = OrderedConfig {
            n_workers: 1,
            reader_threads: 1,
            writer_threads: 1,
            out: tagged.clone(),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        scan_ordered(&bam, refr, Arc::clone(&queries), &ocfg, &stream).unwrap();

        let header = String::from_utf8_lossy(Reader::from_path(&tagged).unwrap().header().as_bytes()).into_owned();
        let choice = choose(&header, None, None).unwrap();
        let chosen = Arc::new(QuerySet::compile(&choice.queries).unwrap());
        let extract = TagExtract::new(&choice, &chosen, parquet_cfg(), crate::strand_rule::directional()).unwrap();
        let extracted = shard("extract_eq_extracted");
        let stats = scan_sharded(&tagged, None, chosen, &extracted, &extract).unwrap();
        assert_eq!(stats.records_scanned, recs.len() as u64);

        let a = rows(&direct.out_prefix, 2, 2);
        let b = rows(&extracted.out_prefix, 2, 2);
        assert!(a.len() >= 8, "enough calls to mean something: {a:?}");
        assert_eq!(a, b);
        assert!(b.iter().flatten().any(|c| c == "refr_base=C"), "the reference base is filled in");
        for (read, strand) in [("ot", "OT"), ("ob", "OB"), ("ctot", "CTOT"), ("del", "OT"), ("hard", "OB")] {
            let row = a.iter().find(|r| r.contains(&format!("qname={read}")));
            let row = row.unwrap_or_else(|| panic!("{read} has a call: {a:?}"));
            assert!(row.contains(&format!("conv_strand={strand}")), "{row:?}");
        }
    }

    #[test]
    fn a_character_the_definition_lacks_is_an_error() {
        use crate::parallel::{scan_sharded, ShardConfig};
        use crate::test_support::temp_path;

        let mut rec = built(b"odd", b"TTTTGTTTCG", &[Cigar::Match(10)], 0, 0);
        rec.push_aux(b"XM", Aux::String("...u......")).unwrap();
        let bam = write_bam("extract_odd", CONTIG, GENOME.len(), &[rec]);
        let file = query_toml::parse_file(TWO_TAGS).unwrap();
        let choice = choose("", Some((file.queries.clone(), file.tags.clone())), None).unwrap();
        let queries = Arc::new(QuerySet::compile(&choice.queries).unwrap());
        let extract = TagExtract::new(&choice, &queries, parquet_cfg(), crate::strand_rule::directional()).unwrap();
        let cfg = ShardConfig {
            strand: crate::strand_rule::directional(),
            n_workers: 1,
            shards_per_worker: 1,
            partition_by: vec![RecordField::Qname],
            reader_threads: 1,
            out_prefix: temp_path("extract_odd", "parquet"),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        let e = scan_sharded(&bam, None, queries, &cfg, &extract).unwrap_err();
        let msg = format!("{e:#}");
        assert!(msg.contains("record odd") && msg.contains("'u' at SEQ position 3"), "{msg}");
    }
    /// Tags are reconciled one at a time. Run 1 defines XM and XE; run 2
    /// redefines XM only. Asking for XE must give run 1's XE without dragging
    /// in run 1's stale XM, which picking a whole run would do.
    #[test]
    fn each_tag_resolves_to_its_own_most_recent_definition() {
        let run1 = "[query.a]\nread = \"C~\"\nrefr = \"CG\"\n\
                    [tag.XM.bases]\nfill = \".\"\nZ = \"a\"\n\
                    [tag.XE.bases]\nfill = \"-\"\ne = \"a\"\n";
        let run2 = "[query.b]\nread = \"T~\"\nrefr = \"CG\"\n\
                    [tag.XM.bases]\nfill = \"_\"\nz = \"b\"\n";
        let h = header_of(&[("alnbase", run1), ("alnbase.1", run2)]);

        let xm = choose(&h, None, Some(TagName::parse("XM").unwrap())).unwrap();
        assert_eq!(xm.bases.fill, b'_', "the newer XM");
        assert_eq!(xm.bases.codes[0].query, "b");
        assert!(xm.source.contains("alnbase.1"), "{}", xm.source);

        let xe = choose(&h, None, Some(TagName::parse("XE").unwrap())).unwrap();
        assert_eq!(xe.bases.fill, b'-', "XE only ever came from run 1");
        assert_eq!(xe.bases.codes[0].query, "a");
        assert!(xe.source.contains("@PG alnbase "), "{}", xe.source);

        // Both are defined, so the tag has to be named.
        let e = choose(&h, None, None).unwrap_err();
        assert!(e.contains("several bases tags") && e.contains("XM") && e.contains("XE"), "{e}");

        // And an undefined tag says what the BAM does define.
        let e = choose(&h, None, Some(TagName::parse("ZZ").unwrap())).unwrap_err();
        assert!(e.contains("no definition of tag ZZ") && e.contains("XM") && e.contains("XE"), "{e}");
    }

    /// The parquet path reads aux tags too, so a damaged block is the same
    /// kind of problem there: it stops the run unless the run is permissive,
    /// because a tag after the damage reads as absent rather than as an error.
    #[test]
    fn a_damaged_aux_block_stops_a_parquet_run() {
        use crate::parallel::{scan_sharded, ShardConfig};
        use crate::scanner::ParquetOutput;
        use crate::test_support::temp_path;
        use rust_htslib::htslib;

        let mut rec = built(b"broken", b"TTTTGTTTCG", &[Cigar::Match(10)], 0, 0);
        rec.push_aux(b"RG", Aux::String("g1")).unwrap();
        // A field nothing can size, so RG is readable and anything after is not.
        let r = unsafe {
            htslib::bam_aux_append(
                rec.inner_mut() as *mut htslib::bam1_t,
                b"XX".as_ptr() as *const std::os::raw::c_char,
                b'?' as std::os::raw::c_char,
                3,
                [1u8, 2, 3].as_ptr(),
            )
        };
        assert_eq!(r, 0);
        let bam = write_bam("extract_damaged", CONTIG, GENOME.len(), &[rec]);
        let refr = shared_reference("extract_damaged", CONTIG, GENOME);

        let file = query_toml::parse_file(TWO_TAGS).unwrap();
        let queries = Arc::new(QuerySet::compile(&file.queries).unwrap());
        let cfg = |tag: &str| ShardConfig {
            strand: crate::strand_rule::directional(),
            n_workers: 1,
            shards_per_worker: 1,
            partition_by: vec![RecordField::Qname],
            reader_threads: 1,
            out_prefix: temp_path(tag, "parquet"),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        let out = |bad| OutputConfig {
            bad_data: bad,
            fields: vec![RecordField::Aux(*b"RG", crate::record_field::AuxKind::Str)],
            ..parquet_cfg()
        };

        let e = scan_sharded(
            &bam,
            Some(Arc::clone(&refr)),
            Arc::clone(&queries),
            &cfg("extract_damaged_stop"),
            &ParquetOutput::new(out(crate::batch::OnBadData::Stop)),
        )
        .unwrap_err();
        let msg = format!("{e:#}");
        assert!(msg.contains("record broken") && msg.contains("damaged aux field"), "{msg}");
        assert!(msg.contains("--permissive"), "{msg}");

        let stats = scan_sharded(
            &bam,
            Some(refr),
            queries,
            &cfg("extract_damaged_warn"),
            &ParquetOutput::new(out(crate::batch::OnBadData::WarnAndCount)),
        )
        .unwrap();
        assert_eq!((stats.damaged_aux, stats.records_scanned), (1, 1));
    }

}

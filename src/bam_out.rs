//! Write a tagged BAM that holds every input record, in the order it went in.
//! The tags on each record are the ones that a query file declares.
//!
//! This is one of the two outputs of a scan. The threads and the ordering
//! belong to [`crate::ordered`]. This module holds what a worker does to each
//! record ([`Tagger`]) and how the result gets to disk ([`BamWriter`]).
//!
//! # What reaches the output
//!
//! Every record that the scan reads reaches the output, in the order in which
//! the scan read it. A record that the scan can walk gets every declared tag. A
//! record that the scan cannot walk goes to the output unchanged. An unmapped
//! record is one example, and a record on a contig that the reference does not
//! have is another. The output is therefore the whole input, and the scan drops
//! or reorders nothing without telling you.
//!
//! # A `bases` tag
//!
//! A `bases` tag has one character per base of SEQ. Where no code applies, it
//! has the `fill` character. When a query that a code names fires, the tag gets
//! that code at the SEQ position of the anchor column of the query.
//!
//! The walk reports a record along its conversion strand. On the bottom strands
//! the offsets of the walk therefore run backwards through SEQ. This module
//! flips the position back, which puts the characters in SEQ order, as the `XM`
//! tag of Bismark has them.
//!
//! A query that a bases tag uses must not be able to fire with its anchor on a
//! deletion column or an intron column, because such a column has no base to
//! mark. alnbase refuses the query file when it parses it. An anchor past either
//! end of the read also has no position. alnbase counts that hit and does not
//! write it.
//!
//! Two queries under different codes of one tag must not mark the same position.
//! This is an error and not a choice between the two queries. It means that the
//! queries overlap, and a bases tag cannot express that.
//!
//! # A `strand` tag
//!
//! A `strand` tag has one value per read. alnbase looks up the value by the
//! strand that the run decided the read came from. The run makes that decision
//! with the rule that the user declares, and not with a built-in reading of the
//! FLAG. See [`crate::strand_rule`] for the reason, and for what happens to a
//! record that the rule does not classify.
//!
//! # The header
//!
//! The output header is the header of the input, plus an `@PG` record for the
//! run, plus a copy of every query file that the run read. You can therefore
//! interpret the tags from the BAM alone. See [`crate::header_query`].
//!
//! # Why this module writes the file through htslib directly
//!
//! `rust_htslib::bam::Writer` closes its file in `Drop` and discards the result.
//! If the flush of the last block fails, for example because the disk is full,
//! you get a truncated BAM and no error. [`BamWriter`] therefore closes the file
//! itself and reports the result, because that file is the result of the run.

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{HeaderView, Record};
use rust_htslib::htslib;

use crate::contig_map::ContigMap;
use crate::header_query::{self, QuerySource};
use crate::query::QuerySet;
use crate::aref::Aref;
use crate::scanner::{Matcher, WalkConfig};
use crate::seq::Seq;
use crate::batch::OnBadData;
use crate::tags::{StrandTag, TagConfig, TagKind, TagName};
use crate::strand_rule::StrandRule;

/// Settings for [`BamStream`] that come from the command line.
#[derive(Debug, Clone)]
pub struct BamOptions {
    pub walk: WalkConfig,
    /// The rule this run calls each record's strand with; see [`crate::strand_rule`].
    pub strand: Arc<StrandRule>,
    /// Replace a tag the input record already has, rather than stop.
    pub overwrite_tags: bool,
    /// What to do about a record whose aux block is damaged.
    pub bad_data: OnBadData,
    /// Recorded in the `@PG` line of every output file.
    pub command_line: String,
    /// The query files the run read, copied into the output header.
    pub query_sources: Vec<QuerySource>,
}

/// The tags, resolved against the compiled queries.
#[derive(Debug)]
struct Plan {
    tags: Vec<Planned>,
    /// Query names by index, for error messages.
    query_names: Vec<String>,
}

#[derive(Debug)]
struct Planned {
    name: TagName,
    kind: PlannedKind,
}

#[derive(Debug)]
enum PlannedKind {
    Bases {
        fill: u8,
        /// The code each query writes, by query index; `None` for a query
        /// this tag does not use.
        code_by_query: Vec<Option<u8>>,
    },
    Strand(StrandTag),
}

/// The tags resolved against the queries, and how to write them: everything a
/// run needs to make [`Tagger`]s and its [`BamWriter`].
pub struct BamStream {
    plan: Arc<Plan>,
    opts: BamOptions,
    version: CString,
    command_line: CString,
    warnings: Vec<String>,
    /// Shared by every tagger; see [`crate::scanner::SharedConcordance`].
    concordance: Arc<crate::scanner::SharedConcordance>,
}

impl BamStream {
    /// Resolve `tags` against `queries`. Fails when there are no tags -- the
    /// output would be a copy of the input -- or a tag names a query the set
    /// does not have.
    pub fn new(tags: &TagConfig, queries: &QuerySet, opts: BamOptions) -> Result<Self> {
        if tags.is_empty() {
            bail!(
                "no tags to write: BAM output writes the tags a TOML query file declares, \
                 e.g. [tag.XM.bases]"
            );
        }
        let query_names: Vec<String> = queries.queries.iter().map(|q| q.name.clone()).collect();
        let index = |name: &str| query_names.iter().position(|n| n == name);

        let mut used = vec![false; query_names.len()];
        let mut planned = Vec::with_capacity(tags.tags.len());
        for t in &tags.tags {
            let kind = match &t.kind {
                TagKind::Bases(c) => {
                    let mut code_by_query = vec![None; query_names.len()];
                    for code in &c.codes {
                        let q = &code.query;
                        let i = index(q).ok_or_else(|| {
                            anyhow!("tag {}: code '{}' names no query '{q}'", t.name, code.code as char)
                        })?;
                        code_by_query[i] = Some(code.code);
                        used[i] = true;
                    }
                    PlannedKind::Bases { fill: c.fill, code_by_query }
                }
                TagKind::Strand(s) => PlannedKind::Strand(s.clone()),
            };
            planned.push(Planned { name: t.name, kind });
        }

        // A query no tag uses still runs, and contributes nothing to a BAM.
        let warnings = query_names
            .iter()
            .zip(&used)
            .filter(|(_, u)| !**u)
            .map(|(n, _)| {
                let shown = if n.is_empty() { "an unnamed query".to_string() } else { format!("query '{n}'") };
                format!("warning: {shown} is used by no tag, so it writes nothing to the BAM output")
            })
            .collect();

        let pg_text = |s: &str| CString::new(s.replace(['\t', '\n', '\r', '\0'], " "));
        Ok(Self {
            plan: Arc::new(Plan { tags: planned, query_names }),
            version: pg_text(env!("CARGO_PKG_VERSION"))?,
            command_line: pg_text(&opts.command_line)?,
            concordance: crate::scanner::SharedConcordance::new(opts.walk.max_discordance),
            opts,
            warnings,
        })
    }

    /// Notes about the configuration worth printing once, before the run.
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// The run's read/reference concordance, shared by every tagger.
    pub fn concordance(&self) -> &crate::scanner::SharedConcordance {
        &self.concordance
    }

    /// A tagger for one worker thread.
    pub fn tagger<'q>(
        &self,
        queries: &'q QuerySet,
        refr: Arc<Aref>,
        contigs: Arc<ContigMap>,
    ) -> Tagger<'q> {
        let n = self.plan.tags.len();
        Tagger {
            matcher: Matcher::new(queries, refr, contigs, self.opts.walk, self.opts.strand.clone(), Arc::clone(&self.concordance)),
            plan: Arc::clone(&self.plan),
            strand: self.opts.strand.clone(),
            overwrite: self.opts.overwrite_tags,
            bad_data: self.opts.bad_data,
            bufs: vec![Vec::new(); n],
            owners: vec![Vec::new(); n],
            counts: TagCounts::default(),
        }
    }

    /// The output file, with `template`'s header and this run's `@PG` line.
    /// `threads` above 1 compresses on that many htslib threads.
    pub fn writer(&self, path: &Path, template: &HeaderView, threads: usize) -> Result<BamWriter> {
        BamWriter::create(
            path,
            template,
            &self.version,
            &self.command_line,
            &self.opts.query_sources,
            threads,
        )
    }
}

/// Where the BAM goes: `-` for standard output, or a path ending in `.bam` or
/// with no extension at all (a named pipe, `/dev/stdout`). Anything else is
/// refused rather than written as a BAM under a name that says otherwise.
pub fn output_path(path: &Path) -> Result<PathBuf> {
    if path.as_os_str() == "-" {
        return Ok(path.to_path_buf());
    }
    match path.extension().and_then(|e| e.to_str()) {
        None => Ok(path.to_path_buf()),
        Some(ext) if ext.eq_ignore_ascii_case("bam") => Ok(path.to_path_buf()),
        Some(ext) => bail!(
            "the output is BAM, so its name should end in .bam, not .{ext}: {}",
            path.display()
        ),
    }
}

/// What a [`Tagger`] has done so far.
#[derive(Debug, Default, Clone, Copy)]
pub struct TagCounts {
    pub scanned: u64,
    /// Records whose aux block was damaged, whose tags were written at the
    /// front of the block; see [`prepend_tags`].
    pub damaged_aux: u64,
    /// Hits anchored past either end of the read, so with no SEQ position.
    pub unplaced_hits: u64,
    /// Records the run's strand rule declined to call, which are written
    /// through untagged and counted.
    pub unknown_strand: u64,
}

/// No query has written this position yet.
const NOBODY: u32 = u32::MAX;

/// Adds the declared tags to one record at a time. One per worker thread.
pub struct Tagger<'q> {
    matcher: Matcher<'q>,
    plan: Arc<Plan>,
    strand: Arc<StrandRule>,
    overwrite: bool,
    bad_data: OnBadData,
    /// Per tag, the value being built for the current record. Reused.
    bufs: Vec<Vec<u8>>,
    /// Per bases tag, which query wrote each position, for reporting an
    /// overlap by name.
    owners: Vec<Vec<u32>>,
    counts: TagCounts,
}

/// Two queries marking one position of one tag with different codes.
struct Clash {
    tag: usize,
    pos: usize,
    first: usize,
    second: usize,
}

impl Tagger<'_> {
    pub fn counts(&self) -> TagCounts {
        self.counts
    }

    /// Walk `record` and add every declared tag to it. Only for a record the
    /// scan can walk: mapped, on a contig the reference has.
    pub fn tag(&mut self, record: &mut Record) -> Result<()> {
        let Self { matcher, plan, strand, overwrite, bad_data, bufs, owners, counts } = self;
        let len = record.seq_len();
        // A record the rule declines is written through untagged and counted,
        // like one the scan cannot walk: a tag written on a guessed strand
        // would be indistinguishable from one the input supported.
        let Some(call) = strand.call(record)? else {
            counts.unknown_strand += 1;
            return Ok(());
        };
        let bottom = call.walk_reversed();

        for (t, (buf, own)) in plan.tags.iter().zip(bufs.iter_mut().zip(owners.iter_mut())) {
            if let PlannedKind::Bases { fill, .. } = t.kind {
                buf.clear();
                buf.resize(len, fill);
                own.clear();
                own.resize(len, NOBODY);
            }
        }

        let mut clash: Option<Clash> = None;
        let mut no_base: Option<usize> = None;
        matcher.run(record, call, &mut |qi, q, ring| {
            let anchor = ring.get(q.anchor_back);
            // Refused when the tags were parsed; see
            // `QuerySpec::anchor_can_lack_read_base`. Checked again here so a
            // gap in that reasoning is an error, not a mark on the wrong base.
            if anchor.read.0 & (Seq::GAP.0 | Seq::SKIP.0) != 0 {
                no_base.get_or_insert(qi);
                return;
            }
            let off = anchor.read_off;
            for (ti, t) in plan.tags.iter().enumerate() {
                let PlannedKind::Bases { code_by_query, .. } = &t.kind else { continue };
                let Some(code) = code_by_query[qi] else { continue };
                // A pad has no base in SEQ; a clip has one, but it is not
                // aligned, so neither is marked.
                if anchor.read == Seq::PAD || anchor.read == Seq::CLIP || off < 0 || off as usize >= len {
                    counts.unplaced_hits += 1;
                    continue;
                }
                // The walk's offset counts along the strand; SEQ is stored
                // along the reference.
                let pos = if bottom { len - 1 - off as usize } else { off as usize };
                let prev = owners[ti][pos];
                if prev == NOBODY {
                    bufs[ti][pos] = code;
                    owners[ti][pos] = qi as u32;
                } else if bufs[ti][pos] != code && clash.is_none() {
                    clash = Some(Clash { tag: ti, pos, first: prev as usize, second: qi });
                }
            }
        })?;

        if let Some(qi) = no_base {
            bail!(
                "record {}: query '{}' fired with its anchor on a deletion or intron, which a \
                 bases tag cannot mark",
                qname_of(record),
                plan.query_names[qi]
            );
        }
        if let Some(c) = clash {
            let t = &plan.tags[c.tag];
            let PlannedKind::Bases { code_by_query, .. } = &t.kind else { unreachable!() };
            let code = |q: usize| code_by_query[q].map(|b| b as char).unwrap_or('?');
            bail!(
                "record {}: queries '{}' and '{}' both mark SEQ position {} (0-based) of tag {}, as '{}' \
                 and '{}'. The queries one bases tag uses must never match at the same \
                 position; tighten one of them.",
                qname_of(record),
                plan.query_names[c.first],
                plan.query_names[c.second],
                c.pos,
                t.name,
                code(c.first),
                code(c.second),
            );
        }

        // One walk answers both questions: is the block intact, and which of
        // our tags does the record already carry. `bam_aux_get` would walk it
        // once per tag and could answer neither on a damaged block.
        let names: Vec<[u8; 2]> = plan.tags.iter().map(|t| t.name.bytes()).collect();
        let scan = crate::aux::scan(record, &names);
        if scan.damaged {
            if *bad_data == OnBadData::Stop {
                bail!(
                    "record {}: its aux block is damaged -- a field partway through cannot be \
                     read, so every tag after it is unreachable to any reader. Pass \
                     --permissive to tag such records anyway (the new tags go at the front of \
                     the block, where readers will find them) and count them.",
                    qname_of(record)
                );
            }
            counts.damaged_aux += 1;
        }

        for (ti, t) in plan.tags.iter().enumerate() {
            let tag = t.name.bytes();
            if scan.present[ti] {
                if !*overwrite {
                    bail!(
                        "record {} already carries tag {}; pass --overwrite-tags to replace it",
                        qname_of(record),
                        t.name
                    );
                }
                if let Err(e) = record.remove_aux(&tag) {
                    bail!("record {}: removing its {} tag: {e}", qname_of(record), t.name);
                }
            }
        }

        // Values first, then one write, because a damaged record needs all of
        // its tags placed at the front of the block in a single rebuild.
        let values: Vec<([u8; 2], String)> = plan
            .tags
            .iter()
            .enumerate()
            .map(|(ti, t)| {
                let v = match &t.kind {
                    PlannedKind::Bases { .. } => {
                        std::str::from_utf8(&bufs[ti]).expect("codes and fill are ASCII").to_string()
                    }
                    PlannedKind::Strand(s) => {
                        // A two-way rule cannot say which of OT/CTOT a read came
                        // from, and a strand tag names all four.
                        // `tags::strand_tags_need_origin` refuses that pair
                        // while the run is still reading its query files, so
                        // this is a backstop that should never speak; it stays a
                        // refusal rather than a guessed value, because a wrong
                        // strand written into a BAM outlives the run.
                        let Some(origin) = call.origin else {
                            bail!(
                                "record {}: tag {} names a strand of origin, but this run's \
                                 strand rule distinguishes only the conversion strand, so the \
                                 record's origin is unknown. Use a rule with a [strand.origin] \
                                 table, or drop the tag.",
                                qname_of(record),
                                t.name
                            );
                        };
                        s.value(origin).to_string()
                    }
                };
                Ok((t.name.bytes(), v))
            })
            .collect::<Result<Vec<_>>>()?;

        if scan.damaged {
            prepend_tags(record, &values)
                .with_context(|| format!("record {}: writing its tags", qname_of(record)))?;
        } else {
            for (tag, value) in &values {
                if let Err(e) = record.push_aux(tag, Aux::String(value)) {
                    bail!("record {}: adding a tag: {e}", qname_of(record));
                }
            }
        }

        counts.scanned += 1;
        Ok(())
    }
}

fn qname_of(record: &Record) -> String {
    String::from_utf8_lossy(record.qname()).into_owned()
}

/// Write `tags` at the *front* of a record's aux block.
///
/// For a damaged block only. Appending would put the new tag after the
/// damage, where no reader reaches it -- the record would look tagged and read
/// as untagged. Written first, it is the first thing any walk sees, and the
/// unreadable tail stays exactly as it was.
///
/// A tag that already exists past the damage cannot be removed, because
/// nothing can find it, so the block may end up holding two of that tag. Ours
/// comes first and is what every reader returns; the caller counts these
/// records so the run can say how many.
fn prepend_tags(record: &mut Record, tags: &[([u8; 2], String)]) -> Result<()> {
    let b = record.inner();
    if b.data.is_null() || b.l_data <= 0 {
        bail!("the record has no data block");
    }
    // SAFETY: htslib keeps `data` valid for `l_data` bytes while the record is
    // borrowed; this is the same slice `Record::data` builds.
    let data = unsafe { std::slice::from_raw_parts(b.data, b.l_data as usize) }.to_vec();
    let l_qseq = b.core.l_qseq.max(0) as usize;
    let aux_at =
        b.core.l_qname as usize + b.core.n_cigar as usize * 4 + l_qseq.div_ceil(2) + l_qseq;
    if aux_at > data.len() {
        bail!("the record's fixed-length fields overrun its data block");
    }

    let mut out = Vec::with_capacity(data.len() + tags.len() * 32);
    out.extend_from_slice(&data[..aux_at]);
    for (tag, value) in tags {
        out.extend_from_slice(tag);
        out.push(b'Z');
        out.extend_from_slice(value.as_bytes());
        out.push(0);
    }
    out.extend_from_slice(&data[aux_at..]);
    record.set_data(&out);
    Ok(())
}


/// A BAM file written through htslib, closed explicitly so that a failure to
/// finish it is an error. See the module comment.
pub struct BamWriter {
    fp: *mut htslib::htsFile,
    hdr: *mut htslib::sam_hdr_t,
    path: PathBuf,
}

// SAFETY: the handles belong to this value alone, and htslib file and header
// handles are not tied to the thread that opened them. The writer is created
// on one thread and used and closed on another, never by two at once.
unsafe impl Send for BamWriter {}

impl BamWriter {
    /// Create `path` (`-` for standard output) with `template`'s header plus an
    /// `@PG` line for this run, and `sources` stored under that line's ID.
    /// htslib gives the line a unique ID and chains it to the previous program.
    fn create(
        path: &Path,
        template: &HeaderView,
        version: &CString,
        cl: &CString,
        sources: &[QuerySource],
        threads: usize,
    ) -> Result<Self> {
        let before = header_query::pg_ids(&String::from_utf8_lossy(template.as_bytes()));
        // SAFETY: the template is a live header; the copy is owned by `w` from
        // here, so every early return below frees it.
        let hdr = unsafe { htslib::sam_hdr_dup(template.inner_ptr()) };
        if hdr.is_null() {
            bail!("copying the header for {}", path.display());
        }
        let mut w = BamWriter { fp: std::ptr::null_mut(), hdr, path: path.to_path_buf() };

        let added = unsafe {
            htslib::sam_hdr_add_pg(
                w.hdr,
                c"alnbase".as_ptr(),
                c"VN".as_ptr(),
                version.as_ptr(),
                c"CL".as_ptr(),
                cl.as_ptr(),
                std::ptr::null::<c_char>(),
            )
        };
        if added < 0 {
            bail!("adding the @PG line to the header of {}", path.display());
        }

        // The ID htslib chose. It adds one record per program chain in the
        // input, so the first new ID is the one the stored files refer to.
        if !sources.is_empty() {
            // SAFETY: `sam_hdr_str` returns the header's own NUL-terminated
            // text, valid until the header is next changed; it is copied out
            // before that.
            let text = unsafe {
                let p = htslib::sam_hdr_str(w.hdr);
                if p.is_null() {
                    bail!("reading back the header of {}", path.display());
                }
                CStr::from_ptr(p).to_string_lossy().into_owned()
            };
            let pg_id = header_query::pg_ids(&text)
                .into_iter()
                .find(|id| !before.contains(id))
                .ok_or_else(|| anyhow!("the @PG line added to {} has no ID", path.display()))?;
            let lines = header_query::header_lines(&pg_id, sources);
            let stored = unsafe {
                htslib::sam_hdr_add_lines(w.hdr, lines.as_ptr() as *const c_char, lines.len())
            };
            if stored < 0 {
                bail!("storing the query files in the header of {}", path.display());
            }
        }

        let cpath = CString::new(path.as_os_str().as_encoded_bytes())
            .with_context(|| format!("output path {} contains a NUL", path.display()))?;
        w.fp = unsafe { htslib::hts_open(cpath.as_ptr(), c"wb".as_ptr()) };
        if w.fp.is_null() {
            let os = std::io::Error::last_os_error();
            bail!("creating {}: {os}", path.display());
        }
        if threads > 1 && unsafe { htslib::hts_set_threads(w.fp, threads as c_int) } < 0 {
            bail!("starting {threads} compression threads for {}", path.display());
        }
        if unsafe { htslib::sam_hdr_write(w.fp, w.hdr) } < 0 {
            bail!("writing the header of {}", path.display());
        }
        Ok(w)
    }

    pub fn write(&mut self, record: &Record) -> Result<()> {
        // SAFETY: both handles are open for the writer's lifetime; the record
        // is a valid bam1_t for the length of the call.
        let r = unsafe {
            htslib::sam_write1(self.fp, self.hdr, record.inner() as *const htslib::bam1_t)
        };
        if r < 0 {
            bail!("writing a record to {}", self.path.display());
        }
        Ok(())
    }

    /// Flush the last block, write the EOF marker and close.
    pub fn close(mut self) -> Result<()> {
        let fp = std::mem::replace(&mut self.fp, std::ptr::null_mut());
        // SAFETY: `fp` was open and is closed exactly once; `Drop` sees null.
        if unsafe { htslib::hts_close(fp) } != 0 {
            bail!("finishing {}: the file is probably incomplete", self.path.display());
        }
        Ok(())
    }
}

impl Drop for BamWriter {
    fn drop(&mut self) {
        // SAFETY: each handle is freed at most once; `close` nulls `fp`.
        unsafe {
            if !self.fp.is_null() {
                htslib::hts_close(self.fp);
            }
            if !self.hdr.is_null() {
                htslib::sam_hdr_destroy(self.hdr);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::tags::Library;
    use super::*;
    use crate::ordered::{scan_ordered, OrderedConfig};
    use crate::query_toml;
    use crate::tags::Strand;
    use crate::test_support::{built, shared_reference, temp_path, write_bam, LAST_IN_TEMPLATE, REVERSE};
    use rust_htslib::bam::record::Cigar;
    use rust_htslib::bam::{Read, Reader};

    const CONTIG: &str = "chr1";
    //                    0         1
    //                    01234567890123
    const GENOME: &str = "TTTCGTTTCGTTTT";

    /// Bismark's CpG calls and conversion tags, for both strands at once: the
    /// walk orients a bottom-strand read so the same patterns apply.
    const BISMARK_CG: &str = r#"
[query.TG]
mark = "+."
where = "tg"

[query.CG]
mark = "+."
where = "cg"

[pat.tg]
read = "T~"
refr = "CG"

[pat.cg]
read = "C~"
refr = "CG"

[tag.XM.bases]
fill = "."
z = "TG"
Z = "CG"

[tag.XR.strand]
CT = ["OT", "OB"]
GA = ["CTOT", "CTOB"]

[tag.XG.strand]
CT = ["OT", "CTOT"]
GA = ["OB", "CTOB"]
"#;

    struct Run {
        stats: crate::batch::Stats,
        /// The tagged BAM, for tests that read it back themselves.
        path: std::path::PathBuf,
        /// The output's records, as `(qname, tags)`.
        records: Vec<(String, Vec<(String, String)>)>,
        header: String,
    }

    fn opts(overwrite: bool) -> BamOptions {
        BamOptions {
            walk: WalkConfig::default(),
            strand: crate::strand_rule::directional(),
            overwrite_tags: overwrite,
            bad_data: crate::batch::OnBadData::Stop,
            command_line: "alnbase query\ttest".to_string(),
            query_sources: Vec::new(),
        }
    }

    fn run_with(toml: &str, recs: &[Record], tag: &str, overwrite: bool) -> Result<Run> {
        run_with_opts(toml, recs, tag, |o| BamOptions { overwrite_tags: overwrite, ..o })
    }

    fn run_with_opts(
        toml: &str,
        recs: &[Record],
        tag: &str,
        tweak: impl FnOnce(BamOptions) -> BamOptions,
    ) -> Result<Run> {
        let overwrite = false;
        let _ = overwrite;
        let file = query_toml::parse_file(toml).unwrap_or_else(|e| panic!("{e}"));
        let queries = QuerySet::compile(&file.queries).unwrap();
        let output = BamStream::new(&file.tags, &queries, tweak(opts(false)))?;
        let bam = write_bam(tag, CONTIG, GENOME.len(), recs);
        let refr = shared_reference(tag, CONTIG, GENOME);
        let out = temp_path(&format!("{tag}_out"), "bam");
        let cfg = OrderedConfig {
            n_workers: 2,
            reader_threads: 1,
            writer_threads: 1,
            out: out.clone(),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        let stats = scan_ordered(&bam, refr, Arc::new(queries), &cfg, &output)?;

        let mut reader = Reader::from_path(&out).unwrap();
        let header = String::from_utf8_lossy(reader.header().as_bytes()).into_owned();
        let mut records = Vec::new();
        for r in reader.records() {
            let r = r.unwrap();
            let tags = ["XM", "XR", "XG", "XE"]
                .iter()
                .filter_map(|t| match r.aux(t.as_bytes()) {
                    Ok(Aux::String(v)) => Some((t.to_string(), v.to_string())),
                    _ => None,
                })
                .collect();
            records.push((String::from_utf8_lossy(r.qname()).into_owned(), tags));
        }
        Ok(Run { stats, path: out, records, header })
    }

    fn get<'a>(run: &'a Run, qname: &str, t: &str) -> Option<&'a str> {
        let (_, tags) = run
            .records
            .iter()
            .find(|(n, _)| n == qname)
            .unwrap_or_else(|| panic!("{qname} was not written"));
        tags.iter().find(|(k, _)| k == t).map(|(_, v)| v.as_str())
    }

    /// A read from each strand. On the top strands a CpG is called at its C;
    /// on the bottom strands at its G, which the read shows as A (converted)
    /// or G (protected). Either way the call sits at that base's SEQ index.
    #[test]
    fn bismark_xm_xr_xg_on_all_four_strands() {
        let m = [Cigar::Match(14)];
        let recs = vec![
            // C at 3 converted, C at 8 protected.
            built(b"ot", b"TTTTGTTTCGTTTT", &m, 0, 0),
            // G at 4 converted (A), G at 9 protected.
            built(b"ob", b"TTTCATTTCGTTTT", &m, 0, REVERSE),
            built(b"ctot", b"TTTTGTTTCGTTTT", &m, 0, LAST_IN_TEMPLATE | REVERSE),
            built(b"ctob", b"TTTCATTTCGTTTT", &m, 0, LAST_IN_TEMPLATE),
        ];
        let run = run_with(BISMARK_CG, &recs, "bamout_strands", false).unwrap();

        let want = [
            ("ot", "...z....Z.....", "CT", "CT"),
            ("ob", "....z....Z....", "CT", "GA"),
            ("ctot", "...z....Z.....", "GA", "CT"),
            ("ctob", "....z....Z....", "GA", "GA"),
        ];
        for (name, xm, xr, xg) in want {
            assert_eq!(get(&run, name, "XM"), Some(xm), "{name} XM");
            assert_eq!(get(&run, name, "XR"), Some(xr), "{name} XR");
            assert_eq!(get(&run, name, "XG"), Some(xg), "{name} XG");
        }
        assert_eq!(run.stats.records_scanned, 4);
        assert!(run.header.contains("@PG\tID:alnbase"), "{}", run.header);
        assert!(run.header.contains("CL:alnbase query test"), "tabs in CL are replaced");
    }

    /// The tag writer flips offsets for exactly the strands the walk runs
    /// backwards on; if the two ever disagreed, bottom-strand calls would land
    /// mirrored.
    #[test]
    fn the_library_and_the_walk_agree_on_orientation() {
        use crate::alignment::walk_alignment;
        let refr = shared_reference("bamout_orient", CONTIG, GENOME);
        let bam = write_bam(
            "bamout_orient",
            CONTIG,
            GENOME.len(),
            &[built(b"x", b"TTTT", &[Cigar::Match(4)], 0, 0)],
        );
        let reader = Reader::from_path(&bam).unwrap();
        let contigs = ContigMap::build(reader.header(), &refr).unwrap();
        for flags in [0, REVERSE, LAST_IN_TEMPLATE, LAST_IN_TEMPLATE | REVERSE] {
            let rec = built(b"x", b"TTTT", &[Cigar::Match(4)], 2, flags);
            let mut cols = Vec::new();
            let mut seq = Vec::new();
            walk_alignment(&rec, &refr, &contigs, &mut seq, WalkConfig { end_context: Some(0), ..WalkConfig::default() }.to_opts(1), Library::Directional.call(&rec), &mut |c| {
                cols.push(c)
            })
            .unwrap();
            let backwards = cols[0].refr_pos > cols[3].refr_pos;
            assert_eq!(Library::Directional.strand(&rec).is_bottom(), backwards, "flags {flags:#x}");
        }
    }

    #[test]
    fn a_tag_already_on_the_record_stops_the_run_unless_overwriting() {
        let mut rec = built(b"old", b"TTTTGTTTCGTTTT", &[Cigar::Match(14)], 0, 0);
        rec.push_aux(b"XM", Aux::String("stale")).unwrap();

        let err = run_with(BISMARK_CG, std::slice::from_ref(&rec), "bamout_exists", false)
            .err()
            .expect("an existing XM is an error");
        let msg = format!("{err:#}");
        assert!(msg.contains("already carries tag XM") && msg.contains("--overwrite-tags"), "{msg}");

        let run = run_with(BISMARK_CG, &[rec], "bamout_overwrite", true).unwrap();
        assert_eq!(get(&run, "old", "XM"), Some("...z....Z....."));
    }

    #[test]
    fn overlapping_queries_under_different_codes_are_an_error() {
        let toml = r#"
[query.TG]
mark = "+."
where = "tg"
[query.anyCG]
mark = "+."
where = "any"
[pat.tg]
read = "T~"
refr = "CG"
[pat.any]
read = "N~"
refr = "CG"
[tag.XM.bases]
fill = "."
z = "TG"
x = "anyCG"
"#;
        let rec = built(b"both", b"TTTTGTTTCGTTTT", &[Cigar::Match(14)], 0, 0);
        let err = run_with(toml, &[rec], "bamout_clash", false).err().expect("clash");
        let msg = format!("{err:#}");
        assert!(msg.contains("record both") && msg.contains("SEQ position 3"), "{msg}");
        assert!(msg.contains("'TG'") && msg.contains("'anyCG'"), "{msg}");
    }

    /// An anchor past the end of the read has no SEQ position: nothing is
    /// written for it, and it is counted.
    #[test]
    fn a_hit_anchored_off_the_read_is_counted_not_written() {
        let toml = r#"
[query.edge]
mark = ".+"
where = "p"
[pat.p]
read = "~_"
refr = "~~"
[tag.XE.bases]
fill = "-"
e = "edge"
"#;
        let rec = built(b"edge", b"TTTTGTTTCGTTT", &[Cigar::Match(13)], 0, 0);
        let run = run_with(toml, &[rec], "bamout_edge", false).unwrap();
        // Two pads per end (the span): the first and second pad past the 3' end,
        // and the second pad before the 5' end, whose window is all pads.
        assert_eq!(run.stats.unplaced_hits, 3);
        assert_eq!(get(&run, "edge", "XE"), Some("-".repeat(13).as_str()));
    }

    /// A clip column stands for a soft-clipped base, which is in SEQ but not
    /// aligned, so a hit on one is counted, not marked. Up to 0.1.8 the flank
    /// there was a pad carrying the clipped base's offset, and the clipped base
    /// was marked.
    #[test]
    fn a_hit_on_the_flank_beside_a_soft_clip_marks_nothing() {
        let toml = r#"
[alias]
f = "{_:}"
[query.flank]
read = "f"
refr = "~"
[tag.XE.bases]
fill = "."
P = "flank"
"#;
        let clip = [Cigar::SoftClip(2), Cigar::Match(6), Cigar::SoftClip(2)];
        for flags in [0, crate::test_support::REVERSE] {
            let rec = built(b"clip", b"TTTCGTTTTT", &clip, 2, flags);
            let run = run_with(toml, &[rec], "bamout_pad_clip", false).unwrap();
            assert_eq!(run.stats.unplaced_hits, 2, "flags {flags:#x}: one flank column at each end");
            assert_eq!(get(&run, "clip", "XE"), Some(".........."), "flags {flags:#x}");
        }
    }

    #[test]
    fn construction_checks() {
        let file = query_toml::parse_file(BISMARK_CG).unwrap();
        let queries = QuerySet::compile(&file.queries).unwrap();

        let err = BamStream::new(&TagConfig::default(), &queries, opts(false)).err().unwrap();
        assert!(err.to_string().contains("no tags to write"), "{err}");

        let with_spare = format!("{BISMARK_CG}\n[query.spare]\nwhere = \"cg\"\n");
        let file = query_toml::parse_file(&with_spare).unwrap();
        let queries = QuerySet::compile(&file.queries).unwrap();
        let out = BamStream::new(&file.tags, &queries, opts(false)).unwrap();
        assert_eq!(out.warnings().len(), 1);
        assert!(out.warnings()[0].contains("query 'spare' is used by no tag"));
        let _ = Strand::ALL;
    }

    #[test]
    fn output_names() {
        assert_eq!(output_path(Path::new("out.bam")).unwrap(), Path::new("out.bam"));
        assert_eq!(output_path(Path::new("-")).unwrap(), Path::new("-"));
        assert_eq!(output_path(Path::new("/dev/stdout")).unwrap(), Path::new("/dev/stdout"));
        assert!(output_path(Path::new("out.parquet")).is_err());
    }

    /// The query file is stored in the output header and comes back from the
    /// BAM alone, byte for byte, parsing to the same tags. Running over that
    /// output again stores a second run beside the first, and the latest is
    /// the second.
    #[test]
    fn the_query_file_travels_in_the_header() {
        use crate::header_query::{stored_runs, QuerySource};
        let latest_run = |h: &str| stored_runs(h).map(|mut r| r.pop());

        let first_text = format!("# with a\ttab and a backslash \\\n{BISMARK_CG}");
        let file = query_toml::parse_file(&first_text).unwrap();
        let queries = QuerySet::compile(&file.queries).unwrap();
        let sources = vec![QuerySource { path: "calls/bismark.toml".into(), text: first_text.clone() }];
        let output = BamStream::new(
            &file.tags,
            &queries,
            BamOptions { query_sources: sources.clone(), ..opts(false) },
        )
        .unwrap();
        let rec = built(b"r", b"TTTTGTTTCGTTTT", &[Cigar::Match(14)], 0, 0);
        let bam = write_bam("bamout_hdr", CONTIG, GENOME.len(), &[rec]);
        let refr = shared_reference("bamout_hdr", CONTIG, GENOME);
        let out1 = temp_path("bamout_hdr_out1", "bam");
        let cfg = OrderedConfig {
            n_workers: 1,
            reader_threads: 1,
            writer_threads: 1,
            out: out1.clone(),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        scan_ordered(&bam, Arc::clone(&refr), Arc::new(queries), &cfg, &output).unwrap();

        let header1 = String::from_utf8_lossy(Reader::from_path(&out1).unwrap().header().as_bytes()).into_owned();
        let run = latest_run(&header1).unwrap().expect("a stored run");
        assert_eq!(run.pg_id, "alnbase");
        assert_eq!(run.sources, sources, "the file comes back byte for byte");
        let (_, tags) = run.resolve().unwrap();
        assert_eq!(tags, file.tags, "and parses to the same tags");
        assert!(header1.contains("@CO\talnbase:v1:line\tPG:alnbase\tfile:0\tline:2\t[query.TG]"), "{header1}");

        // Again, over the first run's output, with a different file.
        let second_text = BISMARK_CG.replace("fill = \".\"", "fill = \"-\"");
        let file2 = query_toml::parse_file(&second_text).unwrap();
        let queries2 = QuerySet::compile(&file2.queries).unwrap();
        let sources2 = vec![QuerySource { path: "dash.toml".into(), text: second_text }];
        let output2 = BamStream::new(
            &file2.tags,
            &queries2,
            BamOptions { query_sources: sources2.clone(), ..opts(true) },
        )
        .unwrap();
        let out2 = temp_path("bamout_hdr_out2", "bam");
        let cfg2 = OrderedConfig { out: out2.clone(), ..cfg };
        scan_ordered(&out1, refr, Arc::new(queries2), &cfg2, &output2).unwrap();

        let header2 = String::from_utf8_lossy(Reader::from_path(&out2).unwrap().header().as_bytes()).into_owned();
        let runs = stored_runs(&header2).unwrap();
        let ids: Vec<&str> = runs.iter().map(|r| r.pg_id.as_str()).collect();
        assert_eq!(ids, ["alnbase", "alnbase.1"], "{header2}");
        assert_eq!(runs[0].sources, sources, "the first run is still there");
        assert_eq!(runs[1].sources, sources2);
        let (_, latest_tags) = latest_run(&header2).unwrap().unwrap().resolve().unwrap();
        assert_eq!(latest_tags, file2.tags);
    }
    /// Append a raw aux field, bypassing rust-htslib's checks. The only way to
    /// build the damaged block this is all about.
    fn raw_append(rec: &mut Record, tag: &[u8; 2], ty: u8, bytes: &[u8]) {
        let r = unsafe {
            htslib::bam_aux_append(
                rec.inner_mut() as *mut htslib::bam1_t,
                tag.as_ptr() as *const c_char,
                ty as c_char,
                bytes.len() as i32,
                bytes.as_ptr(),
            )
        };
        assert_eq!(r, 0, "bam_aux_append failed");
    }

    /// A record whose aux block breaks partway through: `NM` is readable, the
    /// `XX` field has a type nothing can size, and `ZZ` sits past it where no
    /// reader will ever reach.
    fn damaged_record() -> Record {
        let mut r = built(b"broken", b"TTTTGTTTCGTTTT", &[Cigar::Match(14)], 0, 0);
        r.push_aux(b"NM", Aux::I32(3)).unwrap();
        raw_append(&mut r, b"XX", b'?', &[1, 2, 3]);
        raw_append(&mut r, b"ZZ", b'Z', b"unreachable\0");
        r
    }

    /// The claim the whole design rests on: a tag written at the front of a
    /// damaged block is found by every reader, where an appended one would be
    /// past the damage and invisible.
    #[test]
    fn a_tag_prepended_to_a_damaged_block_is_readable() {
        let rec = damaged_record();
        // htslib itself cannot reach a tag past the damage.
        assert!(rec.aux(b"ZZ").is_err(), "the fixture is actually damaged");
        assert!(rec.aux(b"NM").is_ok(), "and readable up to the damage");

        let run = run_with(BISMARK_CG, &[rec], "bamout_damaged", false);
        let err = run.as_ref().err().map(|e| format!("{e:#}")).unwrap_or_default();
        assert!(err.contains("aux block is damaged") && err.contains("--permissive"), "{err}");

        let run = run_with_opts(BISMARK_CG, &[damaged_record()], "bamout_damaged_ok", |o| {
            BamOptions { bad_data: crate::batch::OnBadData::WarnAndCount, ..o }
        })
        .unwrap();
        assert_eq!(get(&run, "broken", "XM"), Some("...z....Z....."));
        assert_eq!(get(&run, "broken", "XR"), Some("CT"));
        assert_eq!(run.stats.damaged_aux, 1);

        // And the record is otherwise untouched: what was readable before is
        // readable still.
        let mut reader = Reader::from_path(&run.path).unwrap();
        let r = reader.records().next().unwrap().unwrap();
        assert_eq!(r.qname(), b"broken");
        assert_eq!(&r.seq().as_bytes(), b"TTTTGTTTCGTTTT");
        assert!(matches!(r.aux(b"NM"), Ok(Aux::I32(3))), "the readable prefix survives");
    }

}

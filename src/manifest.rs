//! Write the run manifest, which says which files make up a hit table. It also
//! records what produced them.
//!
//! A sharded output is a set of files, and no one file says whether the set is
//! complete, which run wrote it, or what its coordinates mean. The manifest
//! does, in two places:
//!
//! - **In every file**, before any row is written: parquet footer key-value
//!   metadata (Arrow schema metadata for an IPC stream) under
//!   [`MANIFEST_KEY`], beside `format_version` and `coordinate_base`. It
//!   describes the run's settings, not its outcome.
//! - **Beside the output**, as `{stem}.manifest.json`, written last and only
//!   when the run succeeded: the same content plus what the run found (record
//!   counts, the concordance totals, the rows in each file). Its presence is the
//!   completeness marker. A run removes any manifest already at that path before
//!   it creates a file, so a manifest left by an earlier run cannot vouch for a
//!   partial output.
//!
//! Both carry `run_id`, so a reader can tell this run's files from ones an
//! earlier run with more threads or shards left under the same naming pattern.
//! Such files are not removed (they were not written by this run), but they are
//! reported, and the manifest's file list is what a reader should open.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use rust_htslib::bam::HeaderView;
use serde::Serialize;

use crate::aref::Aref;
use crate::batch::Stats;

/// Coordinates in hit tables (`refr_pos`, `pos`, `end_pos`, ...) count from 0.
pub const COORDINATE_BASE: u8 = 0;

/// Footer and schema metadata key holding the manifest as JSON.
pub const MANIFEST_KEY: &str = "alnbase_manifest";

#[derive(Serialize, Clone, Debug)]
pub struct Manifest {
    pub alnbase_version: &'static str,
    pub format_version: &'static str,
    pub coordinate_base: u8,
    /// Distinct per run; the same in every file of one run.
    pub run_id: String,
    /// `query` or `extract`.
    pub command: &'static str,
    pub command_line: Vec<String>,
    /// The input BAM as given (`-` for standard input).
    pub input: String,
    /// The reference index as given; `null` for `extract`, which reads none.
    pub reference: Option<String>,
    /// Where this run's strand rule came from: the query file that declared
    /// it, or `--library directional` for the rule still compiled into the
    /// program. The rule itself is not repeated here -- a declared one is the
    /// `[strand.*]` tables of a file whose whole text is in `query_files`.
    ///
    /// Recorded because the rule is a property of the BAM rather than of the
    /// queries, so a file carrying both can be pointed at the wrong input and
    /// still run; the rows then say how their strands were called.
    pub strand: String,
    pub output: Output,
    /// The query files the rows' `name`s are defined in, with their text.
    pub query_files: Vec<QueryFile>,
    /// For `extract`: which tag, and where its definition came from.
    pub extract: Option<Extract>,
    /// For `query`: the walk the rows came from, with defaults resolved.
    pub walk: Option<Walk>,
    /// Every `@SQ` contig, in header order.
    pub contigs: Vec<Contig>,
}

#[derive(Serialize, Clone, Debug)]
pub struct Output {
    /// `parquet` or `ipc`.
    pub format: &'static str,
    /// OUT as given; the files are named after it.
    pub prefix: String,
    pub n_workers: usize,
    pub shards_per_worker: usize,
    pub partition_by: Vec<String>,
    /// Record columns, in file order.
    pub fields: Vec<String>,
    pub only_hits: bool,
    /// `flat` (`read_base`/`refr_base`) or `lists` (`capture_*`).
    pub capture_layout: &'static str,
    /// File names (beside OUT), indexed by slot: the `shard` column's value.
    pub files: Vec<String>,
}

#[derive(Serialize, Clone, Debug)]
pub struct QueryFile {
    pub path: String,
    pub text: String,
}

#[derive(Serialize, Clone, Debug)]
pub struct Extract {
    pub tag: String,
    pub definition: String,
}

#[derive(Serialize, Clone, Debug)]
pub struct Walk {
    /// `skip` or `emit`.
    pub insertions: &'static str,
    /// Pad columns past each read end.
    pub end_context: usize,
    /// Whether `--end-context` was given, rather than derived from `widest_span`.
    pub end_context_given: bool,
    /// Reference bases shown at each end of an intron.
    pub splice_context: usize,
    pub splice_context_given: bool,
    /// The widest query's span, which the contexts default to.
    pub widest_span: usize,
    pub max_discordance: f64,
    pub require_m5: bool,
}

#[derive(Serialize, Clone, Debug)]
pub struct Contig {
    pub name: String,
    pub length: u64,
    /// The reference's MD5 for a contig it has; otherwise the header's M5, if any.
    pub md5: Option<String>,
    /// Whether the reference has the contig; `null` when the run read no reference.
    pub in_reference: Option<bool>,
}

/// What a finished run adds to the manifest file.
#[derive(Serialize, Debug)]
struct RunResult {
    records_read: u64,
    records_scanned: u64,
    records_skipped: u64,
    records_off_reference: u64,
    records_without_seq: u64,
    damaged_aux: u64,
    non_utf8_aux: u64,
    /// Bases compared with the reference and how many disagreed; `null` for
    /// `extract`, which compares none.
    concordance: Option<Concordance>,
    /// Rows written to each file, by slot.
    rows: Vec<u64>,
}

#[derive(Serialize, Debug)]
struct Concordance {
    compared: u64,
    discordant: u64,
}

#[derive(Serialize)]
struct Complete<'a> {
    #[serde(flatten)]
    manifest: &'a Manifest,
    result: RunResult,
}

/// A run id: the time and the process id, in hex. Unique enough to tell two
/// runs' files apart, which is all it is for.
pub fn run_id() -> String {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    format!("{nanos:x}-{:x}", std::process::id())
}

/// Every `@SQ` contig in header order, with the reference's view of it when
/// there is a reference.
pub fn contigs(header: &HeaderView, refr: Option<&Aref>) -> Vec<Contig> {
    let m5 = crate::contig_map::header_m5(header);
    (0..header.target_count())
        .map(|tid| {
            let name = String::from_utf8_lossy(header.tid2name(tid)).into_owned();
            let length = header.target_len(tid).unwrap_or(0);
            let in_refr = refr.and_then(|r| r.tid(&name));
            let md5 = match (refr, in_refr) {
                (Some(r), Some(rt)) => r.md5_hex(rt),
                _ => m5.get(&name).map(|s| s.to_ascii_lowercase()),
            };
            Contig {
                name,
                length,
                md5,
                in_reference: refr.map(|_| in_refr.is_some()),
            }
        })
        .collect()
}

/// `calls.parquet` -> `calls.manifest.json`; `calls` -> `calls.manifest.json`.
pub fn path(prefix: &Path) -> PathBuf {
    let stem = match prefix.extension() {
        Some(_) => prefix.file_stem(),
        None => prefix.file_name(),
    };
    let stem = stem.unwrap_or_default().to_string_lossy();
    prefix.with_file_name(format!("{stem}.manifest.json"))
}

impl Manifest {
    /// The key-value metadata every file of the run carries.
    pub fn footer(&self) -> Vec<(String, String)> {
        vec![
            ("format_version".to_string(), self.format_version.to_string()),
            ("coordinate_base".to_string(), self.coordinate_base.to_string()),
            (
                MANIFEST_KEY.to_string(),
                serde_json::to_string(self).expect("a manifest serializes"),
            ),
        ]
    }

    /// Remove a manifest an earlier run left at this output's path, so that
    /// only a run that finishes can leave one there.
    pub fn remove_stale(&self) -> Result<()> {
        let p = path(Path::new(&self.output.prefix));
        match std::fs::remove_file(&p) {
            Ok(()) => Ok(()),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(e) => Err(e).with_context(|| format!("removing the earlier run's manifest {}", p.display())),
        }
    }

    /// Write the manifest file for a run that succeeded. Written to a scratch
    /// name and renamed, so a reader never sees half of one.
    pub fn write_complete(&self, stats: &Stats) -> Result<PathBuf> {
        let complete = Complete {
            manifest: self,
            result: RunResult {
                records_read: stats.records_read,
                records_scanned: stats.records_scanned,
                records_skipped: stats.records_skipped,
                records_off_reference: stats.records_off_reference,
                records_without_seq: stats.records_without_seq,
                damaged_aux: stats.damaged_aux,
                non_utf8_aux: stats.non_utf8_aux,
                concordance: self.walk.as_ref().map(|_| Concordance {
                    compared: stats.concordance.compared,
                    discordant: stats.concordance.discordant,
                }),
                rows: stats.file_rows.clone(),
            },
        };
        let p = path(Path::new(&self.output.prefix));
        let partial = p.with_extension("json.partial");
        let text = serde_json::to_string_pretty(&complete).expect("a manifest serializes");
        std::fs::write(&partial, text + "\n").with_context(|| format!("writing {}", partial.display()))?;
        std::fs::rename(&partial, &p).with_context(|| format!("renaming {} into place", partial.display()))?;
        Ok(p)
    }

    /// Files beside the output that follow its naming (`{stem}_W_S.{ext}`) but
    /// are not this run's: most likely left by an earlier run that wrote more
    /// files. Sorted.
    pub fn foreign_files(&self) -> Vec<String> {
        let prefix = Path::new(&self.output.prefix);
        let dir = match prefix.parent() {
            Some(d) if !d.as_os_str().is_empty() => d.to_path_buf(),
            _ => PathBuf::from("."),
        };
        let (stem, ext) = match prefix.extension() {
            Some(e) => (prefix.file_stem(), Some(e.to_string_lossy().into_owned())),
            None => (prefix.file_name(), None),
        };
        let stem = stem.unwrap_or_default().to_string_lossy().into_owned();
        let Ok(entries) = std::fs::read_dir(&dir) else {
            return Vec::new();
        };
        let mut out: Vec<String> = entries
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .filter(|name| follows_naming(name, &stem, ext.as_deref()))
            .filter(|name| !self.output.files.contains(name))
            .collect();
        out.sort();
        out
    }
}

/// Whether `name` is `{stem}_{digits}_{digits}` plus `.{ext}` when there is one.
fn follows_naming(name: &str, stem: &str, ext: Option<&str>) -> bool {
    let body = match ext {
        Some(ext) => name.strip_suffix(&format!(".{ext}")),
        None => Some(name),
    };
    let Some(ids) = body
        .and_then(|b| b.strip_prefix(stem))
        .and_then(|b| b.strip_prefix('_'))
    else {
        return false;
    };
    let parts: Vec<&str> = ids.split('_').collect();
    parts.len() == 2
        && parts
            .iter()
            .all(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_manifest_sits_beside_the_output_under_its_stem() {
        assert_eq!(
            path(Path::new("out/calls.parquet")),
            Path::new("out/calls.manifest.json")
        );
        assert_eq!(path(Path::new("calls")), Path::new("calls.manifest.json"));
    }

    #[test]
    fn only_the_output_naming_counts_as_a_run_file() {
        assert!(follows_naming("calls_3_12.parquet", "calls", Some("parquet")));
        assert!(follows_naming("calls_0_0", "calls", None));
        for other in [
            "calls.manifest.json",
            "calls_0.parquet",
            "calls_a_0.parquet",
            "calls_0_0_1.parquet",
            "calls_x_0_0.parquet",
        ] {
            assert!(!follows_naming(other, "calls", Some("parquet")), "{other}");
        }
    }
}

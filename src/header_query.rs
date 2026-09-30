//! Store the query files that a run used in the header of the BAM that the run
//! wrote.
//!
//! A tagged BAM is only interpretable alongside the definitions of its tags:
//! what `z` in `XM` means is whatever query the file mapped it to. So every
//! `--query-file` a run reads is copied into the output header, and can be
//! recovered from the BAM alone -- as text, or parsed back into the same
//! queries and tags with [`StoredRun::resolve`].
//!
//! # Format, version 1
//!
//! Plain `@CO` lines, one per line of each file, so the definitions read
//! naturally in `samtools view -H`:
//!
//! ```text
//! @PG  ID:alnbase    PN:alnbase    VN:0.1.0      CL:alnbase query --query-file bismark.toml ...
//! @CO  alnbase:v1:file  PG:alnbase  file:0  lines:23  path:bismark.toml
//! @CO  alnbase:v1:line  PG:alnbase  file:0  line:0    [query.TG]
//! @CO  alnbase:v1:line  PG:alnbase  file:0  line:1    mark = "+."
//! ```
//!
//! (Fields are separated by single tabs; spaces above are for alignment.)
//!
//! - `PG` is the ID of the `@PG` record the run added, which is what ties a
//!   stored file to one run when a BAM has been through alnbase more than
//!   once. htslib picks the ID -- `alnbase`, `alnbase.1`, ... -- and adds one
//!   `@PG` record per program chain in the input, so the first one added is
//!   the one used.
//! - `file` numbers the files in command-line order; `lines` is how many
//!   `line` records file has, so a lost comment is detected rather than
//!   silently read around; `path` is the path as given.
//! - `line` numbers a file's lines from 0, so the file is rebuilt correctly
//!   even if a tool reorders comments. The last field is the line's text. A
//!   file ending in a newline has a final, empty line, so the text comes back
//!   byte for byte.
//! - Text is escaped so that it stays on one line and inside one field:
//!   `\\` for a backslash, `\t` `\r` `\n` for those characters, `\xHH` for any
//!   other ASCII control character. Everything else, including UTF-8, is
//!   written as is.
//!
//! The version is in the key, so a later format can sit beside this one and a
//! reader of this one skips what it does not understand.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use crate::query_toml;
use crate::tags::TagConfig;
use crate::lower::QueryFile;

const FILE_KEY: &str = "alnbase:v1:file";
const LINE_KEY: &str = "alnbase:v1:line";

/// One query file as a run read it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuerySource {
    /// The path as given on the command line.
    pub path: String,
    pub text: String,
}

/// The `@CO` lines that store `sources` under the `@PG` ID `pg_id`, each
/// ending in a newline, ready to append to a header.
pub fn header_lines(pg_id: &str, sources: &[QuerySource]) -> String {
    let mut out = String::new();
    for (k, src) in sources.iter().enumerate() {
        let lines: Vec<&str> = src.text.split('\n').collect();
        let _ = writeln!(
            out,
            "@CO\t{FILE_KEY}\tPG:{pg_id}\tfile:{k}\tlines:{}\tpath:{}",
            lines.len(),
            escape(&src.path)
        );
        for (i, line) in lines.iter().enumerate() {
            let _ = writeln!(out, "@CO\t{LINE_KEY}\tPG:{pg_id}\tfile:{k}\tline:{i}\t{}", escape(line));
        }
    }
    out
}

/// The IDs of a header's `@PG` records, in header order.
pub fn pg_ids(header_text: &str) -> Vec<String> {
    header_text
        .lines()
        .filter_map(|l| l.strip_prefix("@PG\t"))
        .filter_map(|rest| rest.split('\t').find_map(|f| f.strip_prefix("ID:")))
        .map(str::to_string)
        .collect()
}

/// The query files one run stored.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StoredRun {
    /// The `@PG` ID the run's files are stored under.
    pub pg_id: String,
    pub sources: Vec<QuerySource>,
}

impl StoredRun {
    /// Parse the stored files, as TOML like every query file, and merge their
    /// tags.
    pub fn resolve(&self) -> Result<(Vec<QueryFile>, TagConfig), String> {
        let mut files = Vec::with_capacity(self.sources.len());
        let mut tags = TagConfig::default();
        for src in &self.sources {
            let file = query_toml::parse_file(&src.text)
            .map_err(|e| format!("{} (stored under @PG {}): {e}", src.path, self.pg_id))?;
            tags.merge(file.tags.clone()).map_err(|e| format!("{}: {e}", src.path))?;
            files.push(file);
        }
        Ok((files, tags))
    }
}

/// Every run's stored query files, oldest first: in the order of the runs'
/// `@PG` records, which htslib appends. A run whose `@PG` record has gone
/// missing -- some tool rewrote the header -- comes after the rest.
pub fn stored_runs(header_text: &str) -> Result<Vec<StoredRun>, String> {
    struct FileRec {
        lines: usize,
        path: String,
        text: BTreeMap<usize, String>,
    }
    // Keyed by PG ID, then file index. Lines may arrive before their file
    // record if a tool reordered the comments.
    let mut runs: Vec<(String, BTreeMap<usize, FileRec>)> = Vec::new();
    let mut lines_seen: Vec<(String, usize, usize, String)> = Vec::new();

    for (n, raw) in header_text.lines().enumerate() {
        let Some(rest) = raw.strip_prefix("@CO\t") else { continue };
        let mut f = rest.splitn(5, '\t');
        let key = f.next().unwrap_or("");
        if key != FILE_KEY && key != LINE_KEY {
            continue;
        }
        let bad = |what: &str| format!("header line {}: malformed {key} comment: {what}", n + 1);
        let field = |v: Option<&str>, name: &str| -> Result<String, String> {
            v.and_then(|v| v.strip_prefix(&format!("{name}:")))
                .map(str::to_string)
                .ok_or_else(|| bad(&format!("no {name}")))
        };
        let pg = field(f.next(), "PG")?;
        let file: usize = field(f.next(), "file")?.parse().map_err(|_| bad("file is not a number"))?;
        let run = match runs.iter_mut().position(|(id, _)| *id == pg) {
            Some(i) => &mut runs[i].1,
            None => {
                runs.push((pg.clone(), BTreeMap::new()));
                &mut runs.last_mut().expect("just pushed").1
            }
        };

        if key == FILE_KEY {
            let lines: usize =
                field(f.next(), "lines")?.parse().map_err(|_| bad("lines is not a number"))?;
            let path = unescape(&field(f.next(), "path")?).map_err(|e| bad(&e))?;
            if run.contains_key(&file) {
                return Err(bad(&format!("file {file} of @PG {pg} is described twice")));
            }
            run.insert(file, FileRec { lines, path, text: BTreeMap::new() });
        } else {
            let line: usize =
                field(f.next(), "line")?.parse().map_err(|_| bad("line is not a number"))?;
            // A tool that trims trailing whitespace can drop the tab before
            // an empty line's text; an absent field is an empty line.
            let text = unescape(f.next().unwrap_or("")).map_err(|e| bad(&e))?;
            lines_seen.push((pg, file, line, text));
        }
    }

    for (pg, file, line, text) in lines_seen {
        let run = &mut runs.iter_mut().find(|(id, _)| *id == pg).expect("created above").1;
        let rec = run
            .get_mut(&file)
            .ok_or_else(|| format!("@PG {pg} stores lines of file {file} but no record of the file"))?;
        match rec.text.get(&line) {
            Some(prev) if *prev != text => {
                return Err(format!(
                    "@PG {pg}, {}: line {line} is stored twice, differently",
                    rec.path
                ))
            }
            _ => {
                rec.text.insert(line, text);
            }
        }
    }

    let order = pg_ids(header_text);
    let rank = |id: &str| order.iter().position(|o| o == id).unwrap_or(usize::MAX);
    runs.sort_by_key(|(id, _)| rank(id));

    let mut out = Vec::with_capacity(runs.len());
    for (pg_id, files) in runs {
        let mut sources = Vec::with_capacity(files.len());
        for (expect, (k, rec)) in files.into_iter().enumerate() {
            if k != expect {
                return Err(format!("@PG {pg_id}: file {expect} is missing from the header"));
            }
            if rec.text.len() != rec.lines || rec.text.keys().next_back().is_some_and(|&l| l >= rec.lines) {
                return Err(format!(
                    "@PG {pg_id}, {}: {} of its {} lines are in the header",
                    rec.path,
                    rec.text.keys().filter(|&&l| l < rec.lines).count(),
                    rec.lines
                ));
            }
            let text = rec.text.into_values().collect::<Vec<_>>().join("\n");
            sources.push(QuerySource { path: rec.path, text });
        }
        out.push(StoredRun { pg_id, sources });
    }
    Ok(out)
}

/// A stored query file as text that can be used again: its run's `@PG` record
/// as a comment on top, then the file exactly as the run read it.
///
/// `pg` picks the run by its `@PG` ID, defaulting to the most recent. `file`
/// picks among the files that run stored, and is needed only when it stored
/// more than one.
pub fn dump(header_text: &str, pg: Option<&str>, file: Option<usize>) -> Result<String, String> {
    let runs = stored_runs(header_text)?;
    let run = match pg {
        Some(id) => runs.iter().find(|r| r.pg_id == id).ok_or_else(|| {
            let ids: Vec<&str> = runs.iter().map(|r| r.pg_id.as_str()).collect();
            if ids.is_empty() {
                "the BAM has no stored alnbase query file".to_string()
            } else {
                format!("no query file is stored under @PG {id}; stored under: {}", ids.join(", "))
            }
        })?,
        None => runs.last().ok_or("the BAM has no stored alnbase query file")?,
    };
    let k = match (file, run.sources.len()) {
        (Some(k), n) if k < n => k,
        (Some(k), n) => return Err(format!("@PG {} stored {n} file(s); there is no file {k}", run.pg_id)),
        (None, 1) => 0,
        (None, _) => {
            let list: Vec<String> =
                run.sources.iter().enumerate().map(|(i, s)| format!("{i}: {}", s.path)).collect();
            return Err(format!(
                "@PG {} stored {} files; choose one with --file ({})",
                run.pg_id,
                run.sources.len(),
                list.join(", ")
            ));
        }
    };

    let pg_line = header_text
        .lines()
        .find(|l| {
            l.strip_prefix("@PG\t")
                .is_some_and(|rest| rest.split('\t').any(|f| f.strip_prefix("ID:") == Some(run.pg_id.as_str())))
        })
        .map(str::to_string)
        .unwrap_or_else(|| format!("@PG\tID:{} (no @PG record with this ID remains in the header)", run.pg_id));
    let src = &run.sources[k];
    Ok(format!("# {pg_line}\n# file {k} of {}: {}\n{}", run.sources.len(), src.path, src.text))
}

fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            '\n' => out.push_str("\\n"),
            c if c.is_ascii_control() => {
                let _ = write!(out, "\\x{:02X}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

fn unescape(s: &str) -> Result<String, String> {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('\\') => out.push('\\'),
            Some('t') => out.push('\t'),
            Some('r') => out.push('\r'),
            Some('n') => out.push('\n'),
            Some('x') => {
                let hex: String = chars.by_ref().take(2).collect();
                let b = (hex.len() == 2)
                    .then(|| u8::from_str_radix(&hex, 16).ok())
                    .flatten()
                    .filter(|b| b.is_ascii_control())
                    .ok_or_else(|| format!("bad escape \\x{hex}"))?;
                out.push(b as char);
            }
            other => {
                return Err(format!("bad escape \\{}", other.map(String::from).unwrap_or_default()))
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn latest_run(header: &str) -> Result<Option<StoredRun>, String> {
        Ok(stored_runs(header)?.pop())
    }

    fn src(path: &str, text: &str) -> QuerySource {
        QuerySource { path: path.to_string(), text: text.to_string() }
    }

    fn header_with(pg: &[&str], co: &str) -> String {
        let mut h = String::from("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:14\n");
        for id in pg {
            h.push_str(&format!("@PG\tID:{id}\tPN:alnbase\n"));
        }
        h.push_str(co);
        h
    }

    #[test]
    fn escaping_round_trips_awkward_text() {
        for s in ["", "plain", "a\tb", "back\\slash", "crlf\r", "\\t literal", "ctrl\x01\x7f", "µ café ✓", "\\x41"] {
            assert_eq!(unescape(&escape(s)).unwrap(), s, "{s:?}");
            assert!(!escape(s).contains(['\t', '\n', '\r']), "{s:?} stays on one field");
        }
        assert!(unescape("\\q").is_err());
        assert!(unescape("\\x4").is_err());
        assert!(unescape("\\x41").is_err(), "only control characters are \\x-escaped");
        assert!(unescape("trailing\\").is_err());
    }

    #[test]
    fn files_come_back_byte_for_byte() {
        let sources = vec![
            src("a.toml", "[query.x]\n\n# a\tcomment with \\ and µ\r\nwhere = \"p\"\n"),
            src("dir/b with space.txt", "no trailing newline"),
            src("empty.toml", ""),
        ];
        let header = header_with(&["alnbase"], &header_lines("alnbase", &sources));
        let runs = stored_runs(&header).unwrap();
        assert_eq!(runs, vec![StoredRun { pg_id: "alnbase".into(), sources }]);
    }

    #[test]
    fn reordered_comments_still_rebuild_the_file() {
        let sources = vec![src("q.toml", "one\ntwo\nthree\n")];
        let lines = header_lines("alnbase", &sources);
        let mut co: Vec<&str> = lines.lines().collect();
        co.reverse();
        let header = header_with(&["alnbase"], &(co.join("\n") + "\n"));
        assert_eq!(latest_run(&header).unwrap().unwrap().sources, sources);
    }

    #[test]
    fn a_lost_line_is_an_error_not_a_shorter_file() {
        let lines = header_lines("alnbase", &[src("q.toml", "one\ntwo\nthree")]);
        let kept: Vec<&str> = lines.lines().filter(|l| !l.contains("line:1\t")).collect();
        let header = header_with(&["alnbase"], &(kept.join("\n") + "\n"));
        let e = stored_runs(&header).unwrap_err();
        assert!(e.contains("2 of its 3 lines"), "{e}");
    }

    #[test]
    fn runs_are_ordered_by_their_pg_records() {
        let co = header_lines("alnbase.1", &[src("second.toml", "2")])
            + &header_lines("alnbase", &[src("first.toml", "1")]);
        let header = header_with(&["alnbase", "alnbase.1"], &co);
        let runs = stored_runs(&header).unwrap();
        let ids: Vec<&str> = runs.iter().map(|r| r.pg_id.as_str()).collect();
        assert_eq!(ids, ["alnbase", "alnbase.1"]);
        assert_eq!(latest_run(&header).unwrap().unwrap().sources[0].path, "second.toml");
    }

    #[test]
    fn other_comments_are_ignored_and_no_runs_is_none() {
        let header = header_with(&["bwa"], "@CO\tuser comment\n@CO\talnbase:v2:file\tsomething new\n");
        assert!(stored_runs(&header).unwrap().is_empty());
        assert!(latest_run(&header).unwrap().is_none());
    }

    #[test]
    fn pg_ids_are_read_in_order() {
        let h = "@HD\tVN:1.6\n@PG\tID:bwa\tPN:bwa\n@PG\tPN:alnbase\tID:alnbase\tPP:bwa\n";
        assert_eq!(pg_ids(h), ["bwa", "alnbase"]);
    }

    #[test]
    fn a_dump_carries_its_pg_line_and_the_file_verbatim() {
        let text = "[query.x]\n[pattern.p]\nread = \"C\"\nrefr = \"C\"\n";
        let mut h = String::from("@HD\tVN:1.6\n@PG\tID:alnbase\tPN:alnbase\tCL:alnbase query a b c\n");
        h.push_str(&header_lines("alnbase", &[src("q.toml", text)]));
        let d = dump(&h, None, None).unwrap();
        assert_eq!(d, format!("# @PG\tID:alnbase\tPN:alnbase\tCL:alnbase query a b c\n# file 0 of 1: q.toml\n{text}"));
        // The comments leave it a query file that parses as the stored one did.
        assert!(crate::query_toml::parse_file(&d).is_ok());
    }

    #[test]
    fn a_dump_picks_the_run_and_file_or_says_how() {
        let co = header_lines("alnbase", &[src("a.toml", "A")])
            + &header_lines("alnbase.1", &[src("b.toml", "B"), src("c.toml", "C")]);
        let h = header_with(&["alnbase", "alnbase.1"], &co);
        let e = dump(&h, None, None).unwrap_err();
        assert!(e.contains("stored 2 files") && e.contains("0: b.toml, 1: c.toml"), "{e}");
        assert!(dump(&h, None, Some(1)).unwrap().ends_with("\nC"));
        assert!(dump(&h, Some("alnbase"), None).unwrap().ends_with("\nA"));
        let e = dump(&h, Some("nope"), None).unwrap_err();
        assert!(e.contains("stored under: alnbase, alnbase.1"), "{e}");
        assert!(dump(&h, None, Some(5)).unwrap_err().contains("no file 5"));
        assert!(dump("@HD\tVN:1.6\n", None, None).unwrap_err().contains("no stored"));
    }
}

//! Read a query file, which is written in TOML.
//!
//! ```toml
//! [query.mCpG]
//! mark  = ".......+...."
//! where = "cpg and not junc_end"
//!
//! [pat.cpg]
//! read = "~~~~~~~Yo~~~"
//! refr = "~~~~~~~CG~~~"
//!
//! [alias]
//! o = "{N.}"
//!
//! [tag.XM.bases]
//! fill = "."
//! Z = "mCpG"
//! ```
//!
//! Aliases, patterns, queries and `tag` tables. This module only reads the
//! file: once read, aliases, patterns and queries go through
//! [`crate::lower::resolve`], which decides what they mean, and tags go through
//! [`crate::tags`].
//!
//! # Why the reading is typed, and strict
//!
//! The file is deserialized straight into structs, never into a generic TOML
//! value, so every row and name arrives as exactly the text written. TOML has
//! no implicit typing to begin with -- an unquoted row is a syntax error, not a
//! number -- and it rejects a repeated table or key by itself.
//!
//! What the format does not guard against, the structs do: every table denies
//! unknown keys, so `marks = "+."` is an error rather than a query that
//! silently anchors on column 0.
//!
//! # Errors
//!
//! A syntax or schema error comes from the TOML crate with a position and is
//! reported by line. A problem found afterwards -- a row of the wrong width, a
//! code naming a query that does not exist -- is reported by table and key,
//! e.g. `[pat.cpg] read`, since typed structs keep no positions and every name
//! in a file is unique.

use indexmap::IndexMap;
use serde::de::{self, Deserializer, SeqAccess, Visitor};
use serde::Deserialize;

use crate::strand_rule::StrandRule;
use crate::tags::{self, Names, TagKind, TagName, TagSpec};
use crate::lower::{self, err, FileError, Loc, PendingPat, PendingQuery, QueryFile};

/// A worked example, shown by `--help` and parsed by a test so it cannot drift
/// from the parser.
pub const EXAMPLE: &str = r#"# Methyl-Hi-C. Tables may come in any order: the queries come first, the
# patterns they share come next, and the alias they use comes last.

[query.mCpG]                # a real CpG, not a ligation artifact
mark  = ".......+...."
where = "cpg and not junc_end and not junc_mid"

[query.mCpG_artifact]       # the artifacts themselves, for QC
mark  = ".......+...."
where = "cpg and (junc_end or junc_mid)"

[pat.cpg]                   # read has C or T where the reference has CG
read = "~~~~~~~Yo~~~"
refr = "~~~~~~~CG~~~"

[pat.junc_end]              # DpnII junction on columns 0-7: its final C is
read = "~~~~~~~~~~~~"       # column 7, so a genomic G at 8 fakes a CpG
refr = "GATCGATC~~~~"

[pat.junc_mid]              # junction on columns 4-11: its own internal CG
read = "~~~~~~~~~~~~"       # lands on columns 7-8
refr = "~~~~GATCGATC"

[alias]
o = "{N.}"                  # any base or a gap: on the read, not past its end
"#;

/// `--help` for `--query-file`: the format, the codes a row may hold, and
/// [`EXAMPLE`].
pub fn help() -> String {
    format!(
        "\
TOML file of queries and tags. Repeatable.

Every value is text in double quotes, or a list of them. Tables may come in any
order, and a query may name a pattern declared below it.

  [alias]           one-letter names for base sets: `o = \"{{N.}}\"`. A row holds one
                    character per column and most base sets have no one-character
                    code, so an alias names one. Lowercase letters only.
  [pat.NAME]        a pattern: `read` and `refr` rows of equal width, column i of
                    one against column i of the other.
  [query.NAME]      a query: `where`, which patterns it uses combined with
                    and / or / not and parentheses (optional when the file has
                    one pattern), and `mark`, a row as wide as its patterns:
                    '.' nothing, '+' the anchor, '^' also record this column. The
                    anchor defaults to column 0 and is always recorded.
                    `read` and `refr` in a query declare a pattern of the query's
                    own name, which it uses when it has no `where`:
                      [query.CG]
                      read = \"C~\"
                      refr = \"CG\"
  [tag.XX.bases]    a BAM tag with one character per base: `fill`, the character
                    written where nothing matched, and one entry per code naming
                    the query that writes it: `z = \"TG\"`. One query per code, and
                    a query's anchor must always be on a read base.
  [tag.XX.strand]   a BAM tag with one value per read, by strand of origin:
                    `CT = [\"OT\", \"CTOT\"]`. Every strand needs a value.
  [strand.original] the reference strand a read's sequence came from: one
                    condition per key, `forward = 'XG == \"CT\"'`, over named
                    fields, `flags == 99` and tags.
  [strand.aligned]  whether the read as sequenced aligns forward or reverse:
                    `forward = \"not is_reverse\"`. Both tables take the keys
                    forward, reverse and unknown; exactly one key must match
                    each record, and `unknown` skips it. The strand of origin
                    follows from the pair. Rules for each surveyed aligner ship
                    in `query/strand/`.

Codes in a row are uppercase IUPAC, plus:
  .  -  Z     gap                 =   both sides the same unambiguous base
  _  X        pad (off the read)  /   both sides unambiguous and different
  ,  J        junction (intron)   ~   anything at all
  a-z         an alias
Put '=' or '/' on one side; the other side supplies the base set. Run with
--list-codes for the full table.

Example:

{EXAMPLE}"
    )
}

// ---------------------------------------------------------------- schema --

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawFile {
    #[serde(default)]
    alias: IndexMap<String, String>,
    #[serde(default)]
    pat: IndexMap<String, RawPat>,
    #[serde(default)]
    query: IndexMap<String, RawQuery>,
    #[serde(default)]
    tag: IndexMap<String, RawTag>,
    strand: Option<RawStrand>,
}

/// How to recover each record's strand: the reference strand its sequence
/// came from, and the direction it aligns as sequenced. Both are required;
/// [`crate::strand_rule::StrandRule::from_tables`] says so.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawStrand {
    original: Option<IndexMap<String, String>>,
    aligned: Option<IndexMap<String, String>>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPat {
    read: String,
    refr: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawQuery {
    mark: Option<String>,
    r#where: Option<String>,
    /// The query's own pattern, named after it; see [`parse_file`].
    read: Option<String>,
    refr: Option<String>,
}

/// One table per kind. Exactly one must be present; that is checked after
/// reading, where the message can say so plainly.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawTag {
    bases: Option<IndexMap<String, Names>>,
    strand: Option<IndexMap<String, Names>>,
}

impl<'de> Deserialize<'de> for Names {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = Names;
            fn expecting(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str("a string or a list of strings")
            }
            fn visit_str<E: de::Error>(self, s: &str) -> Result<Names, E> {
                Ok(Names { items: vec![s.to_owned()], list: false })
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Names, A::Error> {
                let mut items = Vec::new();
                while let Some(s) = seq.next_element::<String>()? {
                    items.push(s);
                }
                Ok(Names { items, list: true })
            }
        }
        d.deserialize_any(V)
    }
}

// ----------------------------------------------------------------- parse --

/// Parse a TOML query file.
pub fn parse_file(src: &str) -> Result<QueryFile, FileError> {
    let raw: RawFile = toml::from_str(src).map_err(|e| toml_error(src, &e))?;

    if raw.query.is_empty() && raw.tag.is_empty() && raw.strand.is_none() {
        return Err(err(
            Loc::Toml("[query]".into()),
            "no queries, tags or strand rule in file",
        ));
    }

    let mut alias_lines = Vec::with_capacity(raw.alias.len());
    for (name, set) in &raw.alias {
        let at = Loc::Toml(format!("[alias] {}", key(name)));
        let (c, seq) = lower::alias_decl(name, set).map_err(|e| err(at.clone(), e))?;
        alias_lines.push((c, seq, at));
    }

    let mut pats = Vec::with_capacity(raw.pat.len());
    for (name, p) in &raw.pat {
        let table = format!("[pat.{}]", key(name));
        check_name(name, "pattern", &table)?;
        for (row, which) in [(&p.read, "read"), (&p.refr, "refr")] {
            if row.is_empty() {
                return Err(err(Loc::Toml(format!("{table} {which}")), format!("'{which}' row is empty")));
            }
        }
        pats.push(PendingPat {
            name: name.clone(),
            line: Loc::Toml(table.clone()),
            read: Some((p.read.clone(), Loc::Toml(format!("{table} read")))),
            refr: Some((p.refr.clone(), Loc::Toml(format!("{table} refr")))),
        });
    }

    let mut queries = Vec::with_capacity(raw.query.len());
    for (name, q) in &raw.query {
        let table = format!("[query.{}]", key(name));
        check_name(name, "query", &table)?;
        if q.r#where.as_deref().is_some_and(|w| w.trim().is_empty()) {
            return Err(err(Loc::Toml(format!("{table} where")), "'where' is empty"));
        }

        // The query's own pattern: `read` and `refr` written in the query
        // declare a pattern named after it, which it uses when it has no
        // `where` and which other queries may use like any other.
        let own = match (&q.read, &q.refr) {
            (Some(read), Some(refr)) => Some((read, refr)),
            (None, None) => None,
            (Some(_), None) => {
                return Err(err(
                    Loc::Toml(format!("{table} read")),
                    "has 'read' but no 'refr'; a query's own pattern needs both",
                ))
            }
            (None, Some(_)) => {
                return Err(err(
                    Loc::Toml(format!("{table} refr")),
                    "has 'refr' but no 'read'; a query's own pattern needs both",
                ))
            }
        };
        if let Some((read, refr)) = own {
            if raw.pat.contains_key(name) {
                return Err(err(
                    Loc::Toml(table.clone()),
                    format!(
                        "its 'read' and 'refr' declare a pattern named '{name}', and \
                         [pat.{}] declares another; a pattern name must mean one thing",
                        key(name)
                    ),
                ));
            }
            for (row, which) in [(read, "read"), (refr, "refr")] {
                if row.is_empty() {
                    return Err(err(Loc::Toml(format!("{table} {which}")), format!("'{which}' row is empty")));
                }
            }
            pats.push(PendingPat {
                name: name.clone(),
                line: Loc::Toml(table.clone()),
                read: Some((read.clone(), Loc::Toml(format!("{table} read")))),
                refr: Some((refr.clone(), Loc::Toml(format!("{table} refr")))),
            });
        }

        let where_ = match (&q.r#where, own) {
            (Some(w), _) => Some((w.clone(), Loc::Toml(format!("{table} where")))),
            (None, Some(_)) => Some((name.clone(), Loc::Toml(format!("{table} read")))),
            (None, None) => None,
        };
        queries.push(PendingQuery {
            name: name.clone(),
            line: Loc::Toml(table.clone()),
            mark: q.mark.clone().map(|m| (m, Loc::Toml(format!("{table} mark")))),
            where_,
        });
    }

    let mut file = lower::resolve(alias_lines, &pats, &queries)?;

    let query_names: Vec<&str> = raw.query.keys().map(String::as_str).collect();
    for (name, t) in &raw.tag {
        let table = format!("[tag.{}]", key(name));
        let tag_name = TagName::parse(name).map_err(|e| err(Loc::Toml(table.clone()), e))?;
        let kind = match (&t.bases, &t.strand) {
            (Some(entries), None) => {
                let at = Loc::Toml(format!("[tag.{}.bases]", key(name)));
                let tag = tags::bases_tag(&entries_of(entries), &query_names)
                    .map_err(|e| err(at.clone(), e))?;
                // A bases tag marks a base of the read, so a query that could
                // fire on a deletion or intron column has nothing to mark.
                tags::check_anchors(&tag, &file.queries).map_err(|e| err(at, e))?;
                TagKind::Bases(tag)
            }
            (None, Some(entries)) => {
                let at = Loc::Toml(format!("[tag.{}.strand]", key(name)));
                TagKind::Strand(tags::strand_tag(&entries_of(entries)).map_err(|e| err(at, e))?)
            }
            (None, None) => {
                return Err(err(
                    Loc::Toml(table),
                    "says nothing about what to write; give it a 'bases' or 'strand' table",
                ))
            }
            (Some(_), Some(_)) => {
                return Err(err(
                    Loc::Toml(table),
                    "has both a 'bases' and a 'strand' table; a tag holds one kind of value",
                ))
            }
        };
        file.tags.tags.push(TagSpec { name: tag_name, kind });
    }

    if let Some(st) = &raw.strand {
        let pairs = |t: &Option<IndexMap<String, String>>| -> Vec<(String, String)> {
            t.iter().flat_map(|m| m.iter().map(|(k, v)| (k.clone(), v.clone()))).collect()
        };
        let rule = StrandRule::from_tables(&pairs(&st.original), &pairs(&st.aligned))
        .map_err(|e| err(Loc::Toml("[strand]".into()), e.to_string()))?;
        file.strand = Some(rule);
    }

    Ok(file)
}

fn entries_of(m: &IndexMap<String, Names>) -> Vec<(String, Names)> {
    m.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
}

/// Pattern and query names follow the line syntax's rule, so a name valid in
/// one syntax is valid in the other and `where` can always tokenize it.
fn check_name(name: &str, what: &str, table: &str) -> Result<(), FileError> {
    let at = || Loc::Toml(table.to_string());
    if name.is_empty() {
        return Err(err(at(), format!("{what} name is empty")));
    }
    if name.contains(char::is_whitespace) {
        return Err(err(at(), format!("{what} name must not contain whitespace")));
    }
    Ok(())
}

/// A key as it would be written in a table header: bare when TOML allows it,
/// quoted otherwise, so an error names the table the way the file does.
fn key(k: &str) -> String {
    let bare = !k.is_empty() && k.chars().all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-');
    if bare {
        k.to_string()
    } else {
        format!("\"{}\"", k.replace('\\', "\\\\").replace('"', "\\\""))
    }
}

/// The TOML crate's error, reported by line like every other file error.
fn toml_error(src: &str, e: &toml::de::Error) -> FileError {
    let msg = e.message().trim().to_string();
    match e.span() {
        Some(span) => {
            let line = src[..span.start.min(src.len())].matches('\n').count() + 1;
            let text = src.lines().nth(line - 1).unwrap_or("");
            match unquoted_hint(text) {
                Some(hint) => err(line, format!("{hint} (TOML's own complaint: {msg})")),
                None => err(line, msg),
            }
        }
        None => err(Loc::Toml("file".into()), msg),
    }
}

/// When `line` is `key = value` with the value not in quotes, what to write
/// instead.
///
/// TOML reads an unquoted value as a number, a boolean or a date, and reports
/// failing to -- `z = TG` is "invalid boolean, expected `true`", because the T
/// might have started `true`. Every value in a query file is text or a list of
/// text, so on the line TOML stopped at, an unquoted value that is not a list
/// is the mistake, whatever TOML guessed it was.
fn unquoted_hint(line: &str) -> Option<String> {
    let (k, v) = line.split_once('=')?;
    let k = k.trim();
    if k.is_empty() || k.starts_with('[') {
        return None;
    }
    // An unquoted value has no string for a '#' to hide in.
    let v = v.split('#').next().unwrap_or("").trim();
    let first = v.chars().next()?;
    if matches!(first, '"' | '\'' | '[' | '{') {
        return None;
    }
    let quoted = v.replace('\\', "\\\\").replace('"', "\\\"");
    Some(format!("text must be in quotes in TOML: {k} = \"{quoted}\""))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tags::Strand;

    fn bad(src: &str) -> FileError {
        parse_file(src).unwrap_err()
    }

    const MIN: &str = "[pat.p]\nread = \"C~\"\nrefr = \"CG\"\n[query.x]\n";

    /// The example in `--help` parses, and means what it says it means: four
    /// patterns, two queries anchored on column 7, and the alias they use.
    #[test]
    fn the_help_example_parses() {
        let f = parse_file(EXAMPLE).unwrap_or_else(|e| panic!("{e}"));
        let names: Vec<&str> = f.queries.iter().map(|q| q.name.as_str()).collect();
        assert_eq!(names, ["mCpG", "mCpG_artifact"]);
        for q in &f.queries {
            assert_eq!(q.anchor, 7, "{}", q.name);
            assert_eq!(q.span, 12, "{}", q.name);
        }
        assert_eq!(f.queries[0].operands.len(), 3);
        assert!(f.aliases.get('o').is_some(), "the alias is declared");
        assert!(f.warnings.is_empty(), "{:?}", f.warnings);
        assert!(f.tags.tags.is_empty());
    }

    #[test]
    fn tables_interleave_and_queries_keep_file_order() {
        let src = "[query.b]\nwhere = \"p\"\n[pat.p]\nread = \"C\"\nrefr = \"C\"\n\
                   [query.a]\nwhere = \"q\"\n[pat.q]\nread = \"T\"\nrefr = \"C\"\n";
        let f = parse_file(src).unwrap_or_else(|e| panic!("{e}"));
        let names: Vec<_> = f.queries.iter().map(|q| q.name.as_str()).collect();
        assert_eq!(names, ["b", "a"]);
    }

    #[test]
    fn a_single_pattern_needs_no_where() {
        let f = parse_file(MIN).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(f.queries[0].span, 2);
    }

    #[test]
    fn schema_and_syntax_errors_come_with_a_line() {
        let e = bad(&format!("{MIN}marks = \"+.\"\n"));
        assert_eq!(e.at, Loc::Line(5));
        assert!(e.msg.contains("unknown field `marks`"), "{e}");

        let e = bad("[pat.p]\nread = C~\nrefr = \"CG\"\n[query.x]\n");
        assert_eq!(e.at, Loc::Line(2), "{e}");

        let e = bad(&format!("{MIN}[pat.p]\nread = \"C\"\nrefr = \"C\"\n"));
        assert!(e.msg.contains("duplicate"), "{e}");

        let e = bad("[pat.p]\nread = \"C\"\n[query.x]\n");
        assert!(e.msg.contains("missing field `refr`"), "{e}");

        let e = bad("[query.x]\nwhere = \"p\"\n[pat.p]\nread = \"C\"\nrefr = \"C\"\n[paterns.q]\n");
        assert!(e.msg.contains("unknown field `paterns`"), "{e}");
    }

    #[test]
    fn later_errors_name_the_table_and_key() {
        let e = bad("[pat.p]\nread = \"C~\"\nrefr = \"CGA\"\n[query.x]\n");
        assert_eq!(e.at, Loc::Toml("[pat.p] read".into()));
        assert!(e.msg.contains("read row is 2 columns, refr row is 3"), "{e}");

        let e = bad(&format!("{MIN}mark = \"+..\"\n"));
        assert_eq!(e.at, Loc::Toml("[query.x] mark".into()));

        let e = bad("[pat.p]\nread = \"C1\"\nrefr = \"CG\"\n[query.x]\n");
        assert_eq!(e.at, Loc::Toml("[pat.p] read".into()));
        assert!(e.msg.contains("digit"), "{e}");

        let e = bad("[pat.\"a b\"]\nread = \"C\"\nrefr = \"C\"\n[query.x]\n");
        assert_eq!(e.at, Loc::Toml("[pat.\"a b\"]".into()));
        assert!(e.msg.contains("whitespace"), "{e}");

        let e = bad("[pat.p]\nread = \"\"\nrefr = \"C\"\n[query.x]\n");
        assert_eq!(e.at, Loc::Toml("[pat.p] read".into()));

        let e = bad(&format!("{MIN}[alias]\nab = \"C\"\n"));
        assert_eq!(e.at, Loc::Toml("[alias] ab".into()));
        assert!(e.msg.contains("one character"), "{e}");

        let e = bad(&format!("{MIN}[pat.q]\nread = \"C\"\nrefr = \"C\"\n"));
        assert_eq!(e.at, Loc::Toml("[query.x]".into()));
        assert!(e.msg.contains("needs a 'where'"), "{e}");

        let e = bad("[pat.p]\nread = \"C\"\nrefr = \"C\"\n");
        assert!(e.msg.contains("no queries, tags or strand rule"), "{e}");
    }

    #[test]
    fn an_unused_pattern_is_a_warning_named_by_table() {
        let f = parse_file(&format!("{MIN}where = \"p\"\n[pat.q]\nread = \"C\"\nrefr = \"C\"\n"))
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(f.warnings, ["[pat.q]: pattern 'q' is defined but no query uses it"]);
    }

    /// Bismark's tags, for the top strands: XM from four queries, and the read
    /// and genome conversion tags from the strand a read came from.
    const BISMARK: &str = r#"
[query.TG]
mark = "+."
where = "tg"

[query.CG]
mark = "+."
where = "cg"

[query.THH]
mark = "+.."
where = "thh"

[query.CHH]
mark = "+.."
where = "chh"

[pat.tg]
read = "T~"
refr = "CG"

[pat.cg]
read = "C~"
refr = "CG"

[pat.thh]
read = "T~~"
refr = "CHH"

[pat.chh]
read = "C~~"
refr = "CHH"

[tag.XM.bases]
fill = "."
z = "TG"
Z = "CG"
h = "THH"
H = "CHH"

[tag.XR.strand]
CT = ["OT", "OB"]
GA = ["CTOT", "CTOB"]

[tag.XG.strand]
CT = ["OT", "CTOT"]
GA = ["OB", "CTOB"]
"#;

    #[test]
    fn bismark_tags_parse_into_the_config() {
        let f = parse_file(BISMARK).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(f.queries.len(), 4);
        let names: Vec<String> = f.tags.tags.iter().map(|t| t.name.to_string()).collect();
        assert_eq!(names, ["XM", "XR", "XG"]);

        let TagKind::Bases(xm) = &f.tags.tags[0].kind else { panic!("XM is chars") };
        assert_eq!(xm.fill, b'.');
        let codes: Vec<(char, &str)> = xm.codes.iter().map(|c| (c.code as char, c.query.as_str())).collect();
        assert_eq!(codes, [('z', "TG"), ('Z', "CG"), ('h', "THH"), ('H', "CHH")]);

        let TagKind::Strand(xr) = &f.tags.tags[1].kind else { panic!("XR is strand") };
        let TagKind::Strand(xg) = &f.tags.tags[2].kind else { panic!("XG is strand") };
        assert_eq!(xr.value(Strand::Ob), "CT");
        assert_eq!(xr.value(Strand::Ctot), "GA");
        assert_eq!(xg.value(Strand::Ctot), "CT");
        assert_eq!(xg.value(Strand::Ctob), "GA");
    }

    #[test]
    fn a_file_of_strand_tags_alone_is_fine() {
        let f = parse_file("[tag.XG.strand]\nCT = [\"OT\", \"CTOT\"]\nGA = [\"OB\", \"CTOB\"]\n")
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(f.queries.is_empty());
        assert_eq!(f.tags.tags.len(), 1);
    }

    #[test]
    fn tag_errors() {
        let e = bad(&format!("{MIN}[tag.XMM.bases]\nfill = \".\"\nz = \"x\"\n"));
        assert_eq!(e.at, Loc::Toml("[tag.XMM]".into()));

        let e = bad(&format!("{MIN}[tag.XM.char]\nfill = \".\"\n"));
        assert!(e.msg.contains("unknown field `char`"), "{e}");

        let e = bad(&format!("{MIN}[tag.XM]\n"));
        assert_eq!(e.at, Loc::Toml("[tag.XM]".into()));
        assert!(e.msg.contains("'bases' or 'strand'"), "{e}");

        let e = bad(&format!(
            "{MIN}[tag.XM.bases]\nfill = \".\"\nz = \"x\"\n[tag.XM.strand]\nA = [\"OT\"]\n"
        ));
        assert!(e.msg.contains("both"), "{e}");

        let e = bad(&format!("{MIN}[tag.XM.bases]\nfill = \".\"\nz = \"nope\"\n"));
        assert_eq!(e.at, Loc::Toml("[tag.XM.bases]".into()));
        assert!(e.msg.contains("no query named 'nope' (defined: x)"), "{e}");

        let e = bad(&format!("{MIN}[tag.XM.bases]\nfill = \".\"\nz = 5\n"));
        assert_eq!(e.at, Loc::Line(7), "{e}");
        assert!(e.msg.contains("a string or a list of strings"), "{e}");

        let e = bad(&format!("{MIN}[tag.XR.strand]\nCT = [\"OT\"]\n"));
        assert_eq!(e.at, Loc::Toml("[tag.XR.strand]".into()));
        assert!(e.msg.contains("no value for CTOT, OB, CTOB"), "{e}");

        // A code that TOML only accepts quoted.
        let f = parse_file(&format!("{MIN}[tag.XM.bases]\nfill = \"-\"\n\".\" = \"x\"\n"))
            .unwrap_or_else(|e| panic!("{e}"));
        let TagKind::Bases(t) = &f.tags.tags[0].kind else { panic!() };
        assert_eq!(t.codes[0].code, b'.');
    }

    #[test]
    fn keys_are_shown_the_way_toml_spells_them() {
        assert_eq!(key("cpg_1-x"), "cpg_1-x");
        assert_eq!(key("a.b"), "\"a.b\"");
        assert_eq!(key("say \"hi\""), "\"say \\\"hi\\\"\"");
    }

    /// A bases tag marks a base of the read, so its queries must not be able
    /// to fire with the anchor on a deletion or intron column.
    #[test]
    fn a_bases_tag_refuses_anchors_that_can_lack_a_read_base() {
        let file = |read: &str, mark: &str, extra: &str| {
            format!(
                "[query.q]\nmark = \"{mark}\"\nwhere = \"p{extra}\"\n[pat.p]\nread = \"{read}\"\nrefr = \"CG\"\n\
                 [pat.c]\nread = \"C~\"\nrefr = \"~~\"\n[tag.XM.bases]\nfill = \".\"\nz = \"q\"\n\
                 [alias]\nj = \"{{C-}}\"\nk = \"{{CT}}\"\n"
            )
        };
        // Anchor columns that exclude a gap or skip: accepted.
        for read in ["C~", "N~", "Y~", "=~", "k~"] {
            let src = file(read, "+.", "");
            assert!(parse_file(&src).is_ok(), "{read}: {:?}", parse_file(&src).err());
        }
        // A gap-admitting anchor column is refused ...
        // (`j` is an alias for C or a gap.)
        for read in ["~~", "-~", "j~"] {
            let e = bad(&file(read, "+.", ""));
            assert_eq!(e.at, Loc::Toml("[tag.XM.bases]".into()));
            assert!(e.msg.contains("query 'q' can fire with its anchor on a deletion"), "{e}");
        }
        // ... unless the anchor is somewhere else, or the query needs a pattern
        // that excludes a gap there anyway.
        assert!(parse_file(&file("~N", ".+", "")).is_ok());
        assert!(parse_file(&file("~~", ".+", "")).is_err(), "column 1 is ~ too");
        assert!(parse_file(&file("~~", "+.", " and c")).is_ok());
        // Or-ing one in does not help: the query can still fire without it.
        assert!(parse_file(&file("~~", "+.", " or c")).is_err());
        assert!(parse_file(&file("~~", "+.", " and not c")).is_err());
    }

    #[test]
    fn an_unquoted_value_is_reported_as_one() {
        let e = bad("[query.q]\n[pat.p]\nread = \"C\"\nrefr = \"C\"\n[tag.XM.bases]\nfill = \".\"\nz = TG\n");
        assert_eq!(e.at, Loc::Line(7));
        assert!(e.msg.starts_with("text must be in quotes in TOML: z = \"TG\""), "{e}");
        for (line, want) in [
            ("read = C~", "read = \"C~\""),
            ("mark = +.", "mark = \"+.\""),
            ("fill = .", "fill = \".\""),
            ("where = cg and not x   # comment", "where = \"cg and not x\""),
        ] {
            assert!(unquoted_hint(line).is_some_and(|h| h.ends_with(want)), "{line}: {:?}", unquoted_hint(line));
        }
        for line in ["read = \"C~\"", "z = [\"a\"]", "[tag.XM.bases]", "no equals sign", "k = 'x'"] {
            assert_eq!(unquoted_hint(line), None, "{line}");
        }
    }

    /// A query's own read and refr: the same query as a named pattern and a
    /// `where` naming it, down to the operand's name.
    #[test]
    fn a_query_can_hold_its_own_pattern() {
        let long = parse_file("[pat.CG]\nread = \"C~\"\nrefr = \"CG\"\n[query.CG]\nwhere = \"CG\"\n").unwrap();
        let short = parse_file("[query.CG]\nread = \"C~\"\nrefr = \"CG\"\n").unwrap();
        assert_eq!(format!("{:?}", short.queries), format!("{:?}", long.queries));
        assert_eq!(short.queries[0].anchor, 0, "no mark: the first column is the anchor");
        assert!(short.warnings.is_empty());
    }

    /// Bismark's XM contexts written the short way, against the long way.
    #[test]
    fn the_short_form_makes_the_same_tagged_queries() {
        let contexts = [("TG", "T~", "CG"), ("CG", "C~", "CG"), ("THH", "T~~", "CHH"), ("CHH", "C~~", "CHH")];
        let tag = "[tag.XM.bases]\nfill = \".\"\nz = \"TG\"\nZ = \"CG\"\nh = \"THH\"\nH = \"CHH\"\n";
        let mut short = String::new();
        let mut long = String::new();
        for (name, read, refr) in contexts {
            short.push_str(&format!("[query.{name}]\nread = \"{read}\"\nrefr = \"{refr}\"\n"));
            long.push_str(&format!("[pat.{name}]\nread = \"{read}\"\nrefr = \"{refr}\"\n[query.{name}]\nwhere = \"{name}\"\n"));
        }
        let a = parse_file(&(short + tag)).unwrap();
        let b = parse_file(&(long + tag)).unwrap();
        assert_eq!(format!("{:?}", a.queries), format!("{:?}", b.queries));
        assert_eq!(a.tags, b.tags);
    }

    #[test]
    fn a_query_with_its_own_pattern_can_still_say_where() {
        let src = "[query.cg]\nread = \"C~\"\nrefr = \"CG\"\nwhere = \"cg and not junc\"\n\
                   [query.other]\nwhere = \"cg\"\n\
                   [pat.junc]\nread = \"~~\"\nrefr = \"GA\"\n";
        let f = parse_file(src).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(f.queries[0].operands.len(), 2);
        assert_eq!(f.queries[1].operands.len(), 1, "another query can use it by name");
        assert!(f.warnings.is_empty(), "{:?}", f.warnings);

        // Its own pattern unused by its own where is still a pattern nothing
        // is obliged to use -- reported like any other.
        let f = parse_file("[query.a]\nread = \"C\"\nrefr = \"C\"\nwhere = \"b\"\n[pat.b]\nread = \"G\"\nrefr = \"G\"\n").unwrap();
        assert_eq!(f.warnings, ["[query.a]: pattern 'a' is defined but no query uses it"]);
    }

    #[test]
    fn own_pattern_errors() {
        let e = bad("[query.q]\nread = \"C~\"\n");
        assert_eq!(e.at, Loc::Toml("[query.q] read".into()));
        assert!(e.msg.contains("no 'refr'"), "{e}");
        let e = bad("[query.q]\nrefr = \"CG\"\n");
        assert!(e.msg.contains("no 'read'"), "{e}");
        let e = bad("[query.q]\nread = \"C~\"\nrefr = \"CG\"\n[pat.q]\nread = \"T~\"\nrefr = \"CG\"\n");
        assert_eq!(e.at, Loc::Toml("[query.q]".into()));
        assert!(e.msg.contains("[pat.q] declares another"), "{e}");
        let e = bad("[query.q]\nread = \"C~\"\nrefr = \"CGA\"\n");
        assert_eq!(e.at, Loc::Toml("[query.q] read".into()));
        assert!(e.msg.contains("read row is 2 columns, refr row is 3"), "{e}");

        // A `{...}` group is one column, so these rows are the same width.
        for (read, refr) in [("{CT}~", "CG"), ("C~", "{CT}G"), ("{C.}{G_}", "CG")] {
            let f = parse_file(&format!("[query.q]\nread = \"{read}\"\nrefr = \"{refr}\"\n"));
            assert!(f.is_ok(), "{read} / {refr}: {:?}", f.err());
        }
        let e = bad("[query.q]\nread = \"{CT}~\"\nrefr = \"CGA\"\n");
        assert!(e.msg.contains("read row is 2 columns, refr row is 3"), "{e}");
        let e = bad("[query.q]\nread = \"C~\"\nrefr = \"CG\"\nmark = \"+..\"\n");
        assert_eq!(e.at, Loc::Toml("[query.q] mark".into()));
    }

    /// A strand rule is a query file like any other, so a file that holds
    /// nothing else is still a file: `query/strand/*.toml` are all like this.
    #[test]
    fn a_file_of_only_a_strand_rule_parses() {
        let f = parse_file(
            "[strand.original]\n\
             forward = \"true\"\n\
             [strand.aligned]\n\
             forward = \"not is_reverse\"\n\
             reverse = \"is_reverse\"\n",
        )
        .unwrap_or_else(|e| panic!("{e}"));
        assert!(f.queries.is_empty());
        assert!(f.strand.is_some());
    }

    /// Every rule alnbase ships loads through the real parser, not only through
    /// the engine's own tests. The files are the interface.
    #[test]
    fn every_shipped_strand_rule_parses() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("query/strand");
        let mut seen = 0;
        for entry in std::fs::read_dir(&dir).unwrap() {
            let path = entry.unwrap().path();
            if path.extension().and_then(|e| e.to_str()) != Some("toml") {
                continue;
            }
            let name = path.file_name().unwrap().to_string_lossy().into_owned();
            let src = std::fs::read_to_string(&path).unwrap();
            let f = parse_file(&src).unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(f.strand.is_some(), "{name} declared no rule");
            seen += 1;
        }
        assert_eq!(seen, 12);
    }

    #[test]
    fn a_broken_strand_rule_is_reported_at_the_table() {
        let e = bad("[strand.original]\nforward = \"not is_reverse\"\n");
        assert_eq!(e.at, Loc::Toml("[strand]".into()));
        assert!(e.msg.contains("[strand.aligned]"), "{e}");

        let e = bad(
            "[strand.original]\nforward = \"not is_revarse\"\n[strand.aligned]\nforward = \"true\"\n",
        );
        assert!(e.msg.contains("unknown field"), "{e}");
        assert!(e.msg.contains("`forward`"), "{e}");

        // An unknown table under [strand] is a typo, not a new feature. That
        // includes the tables earlier versions took.
        for table in ["orignal", "origin", "conversion", "sequenced"] {
            let e = bad(&format!("[strand.{table}]\nforward = \"true\"\n"));
            assert!(e.msg.contains("unknown field"), "{table}: {e}");
        }
    }

    #[test]
    fn an_empty_file_says_what_a_file_may_hold() {
        let e = bad("");
        assert!(e.msg.contains("no queries, tags or strand rule"), "{e}");
    }
}

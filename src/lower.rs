//! Lower a query file into queries.
//!
//! The step between a parsed file and a compiled query set: aliases are
//! resolved, patterns are lowered through [`crate::dsl::parse_pattern_text`],
//! `where` expressions are resolved against pattern names, and mark rows
//! become anchors and captures.
//!
//! [`crate::query_toml`] reads the file and calls [`resolve`]. The two were
//! separate because there were once two surface syntaxes -- a line-oriented
//! one lived here, which is where the name comes from -- and they are still
//! separate because reading a file and deciding what it means are different
//! jobs with different failure modes: a TOML error has a line and a column,
//! while a width mismatch has a table and a key.


use crate::dsl::{self, Aliases, ParseError, PatternSpec, QuerySpec};
use crate::predicate::Expr;
use crate::seq::Seq;
use crate::tags::TagConfig;

/// Filler in a `mark` row. Not a space: editors strip trailing whitespace, and
/// a mark on the last column would vanish.
pub const MARK_NONE: char = '.';
pub const MARK_ANCHOR: char = '+';
pub const MARK_CAPTURE: char = '^';

/// Where in a query file something was declared.
///
/// A TOML error has a line and a column, and says so. A rule this module
/// checks -- a row of the wrong width, a `where` naming a pattern that does
/// not exist -- is found after the file has been read into typed structs,
/// which keep no positions, so it names the table and key instead:
/// `[pat.cpg] read`. Every name in a file is unique, so that is unambiguous.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Loc {
    Line(usize),
    Toml(String),
}

impl std::fmt::Display for Loc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Loc::Line(n) => write!(f, "line {n}"),
            Loc::Toml(key) => f.write_str(key),
        }
    }
}

impl From<usize> for Loc {
    fn from(line: usize) -> Self {
        Loc::Line(line)
    }
}

/// A parse failure, with where it happened.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileError {
    pub at: Loc,
    pub msg: String,
}

impl std::fmt::Display for FileError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.at, self.msg)
    }
}

impl std::error::Error for FileError {}

pub(crate) fn err(at: impl Into<Loc>, msg: impl Into<String>) -> FileError {
    FileError { at: at.into(), msg: msg.into() }
}

/// Everything one file declares.
#[derive(Debug, Clone, Default)]
pub struct QueryFile {
    pub aliases: Aliases,
    pub queries: Vec<QuerySpec>,
    /// BAM tags to write, by query. Only the TOML syntax can declare them;
    /// a file in this syntax always leaves this empty.
    pub tags: TagConfig,
    /// How to recover each record's original strand, when the file declares it.
    /// Only the TOML syntax can; a file in the line syntax always leaves this
    /// `None`, and a run with no rule at all keeps the historical behaviour.
    pub strand: Option<crate::strand_rule::StrandRule>,
    /// Non-fatal notes, e.g. a pattern nothing references.
    pub warnings: Vec<String>,
}

/// A `pat` block, before its rows are lowered.
#[derive(Debug, Clone)]
pub(crate) struct PendingPat {
    pub(crate) name: String,
    pub(crate) line: Loc,
    pub(crate) read: Option<(String, Loc)>,
    pub(crate) refr: Option<(String, Loc)>,
}

/// A `query` block, before its patterns are resolved.
#[derive(Debug, Clone)]
pub(crate) struct PendingQuery {
    pub(crate) name: String,
    pub(crate) line: Loc,
    pub(crate) mark: Option<(String, Loc)>,
    pub(crate) where_: Option<(String, Loc)>,
}

/// Lower a file's declarations into queries.
///
/// Shared by whatever reads a query file, so how a declaration is spelled and
/// what it means stay separate concerns.
pub(crate) fn resolve(
    alias_lines: Vec<(char, Seq, Loc)>,
    pats: &[PendingPat],
    queries: &[PendingQuery],
) -> Result<QueryFile, FileError> {
    let mut aliases = Aliases::new();
    for (c, s, at) in alias_lines {
        aliases.insert(c, s).map_err(|e| err(at, e))?;
    }

    // Every pattern is lowered once, against the complete alias table -- which
    // is why an alias may be declared after the pattern that uses it.
    let lowered: Vec<(String, PatternSpec, Loc)> = pats
        .iter()
        .map(|p| lower_pattern(p, &aliases))
        .collect::<Result<_, _>>()?;

    let mut used: Vec<bool> = vec![false; lowered.len()];
    let mut out = Vec::with_capacity(queries.len());
    for q in queries {
        out.push(build(q, &lowered, &mut used, pats.len())?);
    }

    let warnings = lowered
        .iter()
        .zip(used.iter())
        .filter(|(_, u)| !**u)
        .map(|((name, _, at), _)| {
            format!("{at}: pattern '{name}' is defined but no query uses it")
        })
        .collect();

    Ok(QueryFile { aliases, queries: out, tags: TagConfig::default(), strand: None, warnings })
}

/// Check one alias declaration: a single character naming a non-empty base
/// set. Whether the character is free is `Aliases::insert`'s job.
pub(crate) fn alias_decl(name: &str, set: &str) -> Result<(char, Seq), String> {
    let name = name.trim();
    let mut cs = name.chars();
    let (Some(c), None) = (cs.next(), cs.next()) else {
        return Err(format!("alias name must be one character, got '{name}'"));
    };
    let set = set.trim();
    let seq = Seq::from_name(set).ok_or_else(|| format!("'{set}' is not a base set"))?;
    if seq.is_empty() {
        return Err(format!("'{set}' names the empty set"));
    }
    Ok((c, seq))
}

/// How many columns a row describes. A `{...}` group is one column however
/// many characters it takes to write, so counting characters would reject
/// every row that uses one. An unclosed group counts to the end of the row;
/// the pattern parser reports it.
fn row_columns(row: &str) -> usize {
    let mut n = 0;
    let mut in_group = false;
    for c in row.chars() {
        match c {
            '{' if !in_group => in_group = true,
            '}' if in_group => {
                in_group = false;
                n += 1;
            }
            _ if in_group => {}
            _ => n += 1,
        }
    }
    n + usize::from(in_group)
}

/// What one pass over the file's lines yields, before anything is resolved.
#[cfg_attr(not(test), allow(dead_code))] // the line syntax is test scaffolding now; see the module comment

fn lower_pattern(
    p: &PendingPat,
    aliases: &Aliases,
) -> Result<(String, PatternSpec, Loc), FileError> {
    let (read, rl) = p
        .read
        .clone()
        .ok_or_else(|| err(p.line.clone(), format!("pattern '{}' has no 'read' row", p.name)))?;
    let (refr, _) = p
        .refr
        .clone()
        .ok_or_else(|| err(p.line.clone(), format!("pattern '{}' has no 'refr' row", p.name)))?;

    for (row, which) in [(&read, "read"), (&refr, "refr")] {
        if let Some(pos) = row.find(|c: char| c.is_ascii_digit()) {
            return Err(err(
                rl.clone(),
                format!(
                    "pattern '{}': digit in the {which} row at position {pos}. Repeat \
                     counts belong to the terse form; a grid holds one character per \
                     column, so write the columns out.",
                    p.name
                ),
            ));
        }
    }
    let (read_cols, refr_cols) = (row_columns(&read), row_columns(&refr));
    if read_cols != refr_cols {
        return Err(err(
            rl,
            format!(
                "pattern '{}': read row is {read_cols} columns, refr row is {refr_cols}. \
                 A misaligned grid shows up here rather than as shifted output.",
                p.name,
            ),
        ));
    }

    let spec = dsl::parse_pattern_text(&format!("{read}@{refr}"), aliases)
        .map_err(|e: ParseError| err(rl, format!("pattern '{}': {}", p.name, e.msg)))?;
    Ok((p.name.clone(), spec, p.line.clone()))
}

/// Resolve one query against the file's patterns.
fn build(
    q: &PendingQuery,
    pats: &[(String, PatternSpec, Loc)],
    used: &mut [bool],
    n_pats: usize,
) -> Result<QuerySpec, FileError> {
    // A query names the patterns it wants. With exactly one pattern in the file
    // there is nothing to choose, so `where` may be left off.
    let (where_text, wl) = match &q.where_ {
        Some((t, l)) => (t.clone(), l.clone()),
        None if n_pats == 1 => (pats[0].0.clone(), q.line.clone()),
        None => {
            return Err(err(
                q.line.clone(),
                format!(
                    "query '{}' needs a 'where' saying which patterns it uses \
                     (the file defines {n_pats})",
                    q.name
                ),
            ))
        }
    };

    // Only the patterns this query mentions become operands, numbered in order
    // of first appearance so the grid and the tree agree.
    let mut operands: Vec<PatternSpec> = Vec::new();
    let mut chosen: Vec<usize> = Vec::new();
    let expr = dsl::parse_bool_expr(&where_text, |tok, _| {
        let Some(gi) = pats.iter().position(|(n, _, _)| n == tok) else {
            return Err(ParseError {
                at: 0,
                msg: format!(
                    "no pattern named '{tok}' (defined: {})",
                    pats.iter().map(|(n, _, _)| n.as_str()).collect::<Vec<_>>().join(", ")
                ),
            });
        };
        used[gi] = true;
        Ok(Expr::Pattern(match chosen.iter().position(|&c| c == gi) {
            Some(local) => local,
            None => {
                chosen.push(gi);
                operands.push(pats[gi].1.clone());
                operands.len() - 1
            }
        }))
    })
    .map_err(|e| err(wl.clone(), e.msg))?;

    let span = operands[0].len();
    if let Some(bad) = operands.iter().position(|o| o.len() != span) {
        return Err(err(
            wl,
            format!(
                "query '{}' uses '{}' ({} columns) alongside '{}' ({span} columns); \
                 every pattern in a query describes the same window",
                q.name,
                pats[chosen[bad]].0,
                operands[bad].len(),
                pats[chosen[0]].0
            ),
        ));
    }

    // ---- markers ----------------------------------------------------------
    let mut anchor: Option<usize> = None;
    let mut captures: Vec<usize> = Vec::new();
    if let Some((row, ml)) = &q.mark {
        if row.chars().count() != span {
            return Err(err(
                ml.clone(),
                format!(
                    "mark row is {} columns, the patterns are {span}",
                    row.chars().count()
                ),
            ));
        }
        for (i, c) in row.chars().enumerate() {
            match c {
                MARK_NONE => {}
                MARK_CAPTURE => captures.push(i),
                MARK_ANCHOR => {
                    if anchor.is_some() {
                        return Err(err(ml.clone(), "more than one '+' in the mark row"));
                    }
                    anchor = Some(i);
                }
                other => {
                    return Err(err(
                        ml.clone(),
                        format!(
                            "'{other}' in a mark row; use '{MARK_NONE}' for nothing, \
                             '{MARK_ANCHOR}' for the anchor, '{MARK_CAPTURE}' to capture"
                        ),
                    ))
                }
            }
        }
    }

    // The names the user wrote, in the order the operands were numbered, so the
    // explainer can say `cpg` where the terse form can only say `[1]`.
    let operand_names = chosen.iter().map(|&gi| pats[gi].0.clone()).collect();

    dsl::finish(
        q.name.clone(),
        &where_text,
        0,
        expr,
        operands,
        operand_names,
        captures,
        anchor,
    )
    .map_err(|e| err(q.line.clone(), e.msg))
}

// The rules lowered here -- row widths, mark rows, alias declarations, `where`
// resolution, the unused-pattern warning -- are tested through the front end
// that reaches them, in `query_toml`. Testing them again against a syntax
// nobody can type is what the deleted tests here were doing.

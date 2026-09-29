//! Hold the query language that sits below the file format. The language has
//! three levels: columns, patterns, and boolean expressions over patterns.
//!
//! [`crate::query_toml`] reads the query files, and [`crate::lower`] lowers
//! them. `lower` gives each pattern to this module as `read@refr` text, and each
//! `where` expression to [`parse_bool_expr`].
//!
//! There are three layers. The innermost layer comes first.
//!
//! **Column** — one pair of a read value and a reference value. Each side is a
//! base code, a `{...}` group, or a relational code. A relational code names the
//! other side instead of a set of bases. At most one of the two sides can be
//! relational, because `=` on both sides would say nothing. A count after a
//! column repeats that column, so `~8` is eight columns in two characters.
//!
//! **Pattern** — two sides that `@` joins. The two sides describe the same
//! columns from each end, so they must agree on the number of columns. This is a
//! count of columns after the code expands the counts, and not a count of
//! characters. `~8@GATCGATC` is therefore well formed.
//!
//! **Query** — a boolean expression over patterns. All of the operands span the
//! same window, so column *i* means the same column in every operand. This is
//! what makes the anchor and the capture set properties of the query, and not of
//! one operand.
//!
//! The anchor (`+`) and the captures (`^`) therefore come from the `mark` row of
//! a query, and not from a pattern. The grammar here still recognises them
//! inside a pattern, but only so that it can refuse them by name.
//! [`parse_pattern_text`] is the entry point that the TOML front end uses. It
//! parses a pattern, then rejects any marker that it found and points the user
//! at the mark row. If it accepted the markers quietly, there would be two
//! places to write the same thing, and a pattern that two queries share could
//! carry the anchor of one of them.
//!
//! This module turns strings into structures and does nothing else. It never
//! touches a BAM, and every check that it makes is a compile-time property of
//! the text. The checks that need the compiled automaton live in `query.rs`.
//! Those are the check for a dead column and the check for satisfaction on an
//! empty window.

use std::fmt;

use crate::predicate::Expr;
use crate::seq::Seq;

/// Order in which auto-generated aliases are handed out.
///
/// Any lowercase letter is legal; this only decides which ones a machine picks
/// when converting a query to the matrix form. `f i j l o p q` come first
/// because they are the only letters that are not a base code in either case --
/// a minted `a` would read as a soft-masked A in a grid, which is exactly the
/// confusion aliases are supposed to avoid.
pub const ALIAS_MINT_ORDER: &str = "fijlopqbdehkmnrsuvwxyzacgt";

/// True when `c` can be used as an alias.
///
/// The whole rule: **uppercase letters are built-in codes, lowercase letters
/// are yours.** Twenty-six of them, no list to consult, and nothing to collide
/// with -- a pattern token always contains `@`, so it can never be mistaken for
/// the `and` / `or` / `not` keywords.
pub fn is_alias_char(c: char) -> bool {
    c.is_ascii_lowercase()
}

/// User-declared single-character names for base sets.
///
/// The matrix DSL needs exactly one character per column, but only 19 of the
/// 127 non-empty `Seq` values have a single-character form — anything mixing
/// bases with a flag, like `{C.}`, has none. An alias gives a one-character
/// form to the few values that an analysis needs. The alias is also declared at
/// the top of the file, and not written inline in a pattern.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Aliases {
    map: Vec<(char, Seq)>,
}

impl Aliases {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (char, Seq)> + '_ {
        self.map.iter().copied()
    }

    /// Lowercase letters not yet taken, in mint order.
    pub fn available(&self) -> Vec<char> {
        ALIAS_MINT_ORDER.chars().filter(|&c| self.get(c).is_none()).collect()
    }

    /// Fails on anything that is not a free lowercase letter.
    pub fn insert(&mut self, c: char, s: Seq) -> Result<(), String> {
        if !is_alias_char(c) {
            let why = if c.is_ascii_uppercase() {
                "uppercase letters are base codes"
            } else {
                "only letters can be aliases"
            };
            return Err(format!(
                "'{c}' cannot be an alias: {why}. Alias names are lowercase letters, \
                 a through z."
            ));
        }
        if let Some((_, prev)) = self.map.iter().find(|(k, _)| *k == c) {
            return Err(format!("'{c}' is already an alias for {}", prev.name()));
        }
        self.map.push((c, s));
        Ok(())
    }

    /// Fold another table in. Re-declaring a character with the same meaning is
    /// fine -- two query files may reasonably share a legend -- but a conflict
    /// is an error, because whichever won would silently change one file.
    pub fn merge(&mut self, other: &Aliases) -> Result<(), String> {
        for (c, s) in other.iter() {
            match self.get(c) {
                Some(prev) if prev == s => {}
                Some(prev) => {
                    return Err(format!(
                        "'{c}' is declared as {} in one place and {} in another",
                        prev.name(),
                        s.name()
                    ))
                }
                None => self.insert(c, s)?,
            }
        }
        Ok(())
    }

    pub fn get(&self, c: char) -> Option<Seq> {
        self.map.iter().find(|(k, _)| *k == c).map(|(_, v)| *v)
    }

    /// The alias standing for `s`, if one was declared. Used when rendering back
    /// to text so a round trip produces what the user wrote.
    pub fn name_of(&self, s: Seq) -> Option<char> {
        self.map.iter().find(|(_, v)| *v == s).map(|(k, _)| *k)
    }
}

/// Separator between the read and reference sides of a pattern.
pub const SIDE_SEP: char = '@';
/// Marks a column for capture: record what was observed there.
pub const CAPTURE: char = '^';
/// Marks the column whose coordinates are reported for the query.
pub const ANCHOR: char = '+';

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    /// Byte offset into whatever string was being parsed, for the caret.
    pub at: usize,
    pub msg: String,
}

impl ParseError {
    fn new(at: usize, msg: impl Into<String>) -> Self {
        Self { at: at, msg: msg.into() }
    }

    /// The message with a caret under the offending byte.
    #[cfg_attr(not(test), allow(dead_code))] // the terse syntax is no longer an input; kept for tests
    pub fn render(&self, src: &str) -> String {
        let mut s = String::new();
        s.push_str(src);
        s.push('\n');
        for _ in 0..self.at.min(src.len()) {
            s.push(' ');
        }
        s.push_str("^ ");
        s.push_str(&self.msg);
        s
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "at byte {}: {}", self.at, self.msg)
    }
}

impl std::error::Error for ParseError {}

// ---------------------------------------------------------------------------
// Columns
// ---------------------------------------------------------------------------

/// Constrains the two sides of a column against each other rather than against
/// a fixed base set.
///
/// Both require unambiguous A/C/G/T on *both* sides, so neither fires on a gap,
/// a pad, or an ambiguity code. They are therefore not complements: their union
/// is "both sides real and unambiguous", not "everything". If `Ne` were the
/// strict negation of `Eq` it would be true at every padded position past the
/// end of every read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Rel {
    /// `=` — both sides unambiguous ACGT and the same base.
    Eq,
    /// `/` — both sides unambiguous ACGT and different bases.
    Ne,
}

impl Rel {
    pub fn to_char(self) -> char {
        match self {
            Rel::Eq => '=',
            Rel::Ne => '/',
        }
    }

    /// True when an observed symbol satisfies this relation.
    ///
    /// `read` and `refr` are concrete observed values, not pattern sets.
    #[inline]
    pub fn holds(self, read: Seq, refr: Seq) -> bool {
        let a = read.0 & Seq::BASES;
        let b = refr.0 & Seq::BASES;
        // No flag bits, and exactly one base bit on each side.
        if read.0 != a || refr.0 != b {
            return false;
        }
        if a.count_ones() != 1 || b.count_ones() != 1 {
            return false;
        }
        match self {
            Rel::Eq => a == b,
            Rel::Ne => a != b,
        }
    }
}

/// One resolved column of a pattern.
///
/// An observed symbol matches iff the read observation is a subset of `read`,
/// the reference observation is a subset of `refr`, and `rel` holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ColumnSpec {
    pub read: Seq,
    pub refr: Seq,
    pub rel: Option<Rel>,
}

impl ColumnSpec {
    /// Render back to DSL text, as `(read_side, refr_side)`.
    pub fn to_text(self) -> (String, String) {
        self.to_text_with(&Aliases::new())
    }

    /// As [`to_text`](Self::to_text), preferring a declared alias for either
    /// side's base set.
    ///
    /// Every rendering of a column goes through here, so which side carries a
    /// relational code is decided once. A caller that rendered the sides itself
    /// and then substituted aliases afterwards would have to recognise `=` and
    /// `/` by their text to avoid aliasing them, and would be wrong the moment
    /// a base set rendered as one of those characters.
    pub fn to_text_with(self, aliases: &Aliases) -> (String, String) {
        let code = |s: Seq| match aliases.name_of(s) {
            Some(c) => c.to_string(),
            None => s.name(),
        };
        match self.rel {
            None => (code(self.read), code(self.refr)),
            // Whichever side carried the relational code got no base set of its
            // own, so it is the side that renders as '=' or '/'. When both are
            // unconstrained (`~@=`) the reference side takes it, matching how
            // it is almost always written.
            Some(r) if self.refr == Seq::N_GAP_PAD || self.refr == Seq::ANY => {
                (code(self.read), r.to_char().to_string())
            }
            Some(r) => (r.to_char().to_string(), code(self.refr)),
        }
    }
}

/// One column of one side, before the two sides are zipped together.
#[derive(Debug, Clone, Copy)]
struct SideCol {
    seq: Seq,
    rel: Option<Rel>,
    count: usize,
    capture: bool,
    anchor: bool,
    at: usize,
    /// Set once a count is seen, so markers after a count can be rejected and a
    /// second count can be caught.
    counted: bool,
}

impl SideCol {
    /// A bare column: one wide, unmarked, uncounted.
    ///
    /// Every column starts here; counts and markers are applied afterwards by
    /// the branches that scan them, which is why they are not parameters.
    fn new(seq: Seq, rel: Option<Rel>, at: usize) -> Self {
        SideCol { seq, rel, count: 1, capture: false, anchor: false, at, counted: false }
    }
}

/// Scan one side of a pattern into columns, counts expanded.
///
/// Returns one entry per emitted column; a `C3` produces three, all sharing the
/// original byte offset for error reporting.
fn parse_side_with(s: &str, aliases: &Aliases) -> Result<Vec<SideCol>, ParseError> {
    if !s.is_ascii() {
        return Err(ParseError::new(0, "patterns must be ASCII"));
    }
    let b = s.as_bytes();
    let mut out: Vec<SideCol> = Vec::new();
    let mut i = 0usize;

    while i < b.len() {
        let c = b[i] as char;
        match c {
            c if c.is_ascii_whitespace() => {
                return Err(ParseError::new(i, "whitespace is not allowed inside a pattern"))
            }

            '{' => {
                let close = b[i..]
                    .iter()
                    .position(|&x| x == b'}')
                    .map(|p| i + p)
                    .ok_or_else(|| ParseError::new(i, "unclosed '{'"))?;
                let body = &s[i + 1..close];
                if body.is_empty() {
                    return Err(ParseError::new(i, "empty group"));
                }
                if body.contains('{') {
                    return Err(ParseError::new(i + 1, "nested '{' is not allowed"));
                }
                // `,` is SKIP outside a group, but a comma inside one is
                // almost always a guess at separator syntax. Allowing it would
                // make `{A,C,G}` parse as A, C, G *and a junction* -- a set
                // that behaves identically on DNA and diverges only on spliced
                // data, which is the worst way to be wrong. `J` is SKIP here.
                if let Some(at) = body.find(',') {
                    return Err(ParseError::new(
                        i + 1 + at,
                        "commas are not allowed in a group; write the codes together, as \
                         {ACG}, and spell a junction 'J', as {NJ}",
                    ));
                }
                let seq = Seq::from_name(body)
                    .ok_or_else(|| ParseError::new(i + 1, format!("unknown code '{body}'")))?;
                if seq.is_empty() {
                    return Err(ParseError::new(i + 1, format!("'{body}' names the empty set")));
                }
                out.push(SideCol::new(seq, None, i));
                i = close + 1;
            }

            '=' | '/' => {
                let rel = if c == '=' { Rel::Eq } else { Rel::Ne };
                // The relational code supplies no base set of its own; the
                // relation does all the constraining.
                out.push(SideCol::new(Seq::N_GAP_PAD, Some(rel), i));
                i += 1;
            }

            '0'..='9' => {
                let start = i;
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
                let n: usize = s[start..i]
                    .parse()
                    .map_err(|_| ParseError::new(start, "repeat count is too large"))?;
                let Some(last) = out.last_mut() else {
                    return Err(ParseError::new(start, "repeat count with no column before it"));
                };
                if last.counted {
                    return Err(ParseError::new(start, "two repeat counts on one column"));
                }
                if last.capture || last.anchor {
                    return Err(ParseError::new(
                        start,
                        "repeat count must come before the markers, e.g. 'C3' not 'C^3'",
                    ));
                }
                if n == 0 {
                    return Err(ParseError::new(start, "repeat count must be at least 1"));
                }
                last.count = n;
                last.counted = true;
            }

            CAPTURE | ANCHOR => {
                let Some(last) = out.last_mut() else {
                    return Err(ParseError::new(i, format!("'{c}' with no column before it")));
                };
                if last.count > 1 {
                    return Err(ParseError::new(
                        i,
                        format!(
                            "'{c}' cannot mark a repeated column: it is ambiguous how many \
                             columns it applies to. Write the columns out separately."
                        ),
                    ));
                }
                let flag = if c == CAPTURE { &mut last.capture } else { &mut last.anchor };
                if *flag {
                    return Err(ParseError::new(i, format!("repeated '{c}'")));
                }
                *flag = true;
                i += 1;
            }

            '}' => return Err(ParseError::new(i, "unmatched '}'")),

            _ if aliases.get(c).is_some() => {
                out.push(SideCol::new(aliases.get(c).expect("just checked"), None, i));
                i += 1;
            }

            _ => {
                let body = &s[i..i + 1];
                if body.chars().next().is_some_and(|ch| ch.is_ascii_lowercase()) {
                    return Err(ParseError::new(
                        i,
                        format!(
                            "'{body}' is lowercase; base codes are uppercase so that 'and', \
                             'or' and 'not' stay unambiguous"
                        ),
                    ));
                }
                let seq = Seq::from_name(body)
                    .ok_or_else(|| ParseError::new(i, format!("unknown code '{body}'")))?;
                if seq.is_empty() {
                    return Err(ParseError::new(i, format!("'{body}' names the empty set")));
                }
                out.push(SideCol::new(seq, None, i));
                i += 1;
            }
        }
    }

    if out.is_empty() {
        return Err(ParseError::new(0, "empty pattern side"));
    }

    // Expand counts.
    let mut expanded = Vec::with_capacity(out.iter().map(|c| c.count).sum());
    for c in out {
        for _ in 0..c.count {
            expanded.push(SideCol { count: 1, ..c });
        }
    }
    Ok(expanded)
}

// ---------------------------------------------------------------------------
// Patterns
// ---------------------------------------------------------------------------

/// One operand of a query: a `read@refr` pattern, columns resolved.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PatternSpec {
    pub columns: Box<[ColumnSpec]>,
    /// As the user wrote it, for error messages and the output `expr` column.
    pub text: String,
}

impl PatternSpec {
    pub fn len(&self) -> usize {
        self.columns.len()
    }
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }
}

/// Markers found on a pattern, as column indices into the expanded pattern.
#[derive(Debug, Clone, Default)]
struct Markers {
    captures: Vec<usize>,
    anchor: Option<usize>,
}

/// Parse one `read@refr` pattern, discarding markers.
///
/// The matrix front end takes its markers from a dedicated row, so it wants the
/// columns without them -- but it must lower through exactly this code, or the
/// two surfaces could disagree about what a column means.
pub fn parse_pattern_text(s: &str, aliases: &Aliases) -> Result<PatternSpec, ParseError> {
    let (spec, markers) = parse_pattern(s, aliases)?;
    if markers.anchor.is_some() || !markers.captures.is_empty() {
        return Err(ParseError::new(
            0,
            format!(
                "'{CAPTURE}' and '{ANCHOR}' go in the mark row here, not in the pattern"
            ),
        ));
    }
    Ok(spec)
}

fn parse_pattern(s: &str, aliases: &Aliases) -> Result<(PatternSpec, Markers), ParseError> {
    let sep = s.find(SIDE_SEP);
    let Some(sep) = sep else {
        return Err(ParseError::new(
            s.len(),
            format!(
                "pattern has no '{SIDE_SEP}' separating the read and reference sides. \
                 Note that whitespace splits patterns, so a stray space inside one shows up \
                 here; quote the whole argument to keep the shell out of it."
            ),
        ));
    };
    if s[sep + 1..].contains(SIDE_SEP) {
        return Err(ParseError::new(
            sep + 1 + s[sep + 1..].find(SIDE_SEP).unwrap(),
            format!("more than one '{SIDE_SEP}' in a pattern"),
        ));
    }

    let read = parse_side_with(&s[..sep], aliases)?;
    // Reference-side offsets are relative to that substring; shift them so the
    // caret lands in the right place in the whole pattern.
    let refr = parse_side_with(&s[sep + 1..], aliases)
        .map_err(|e| ParseError { at: e.at + sep + 1, ..e })?;

    if read.len() != refr.len() {
        return Err(ParseError::new(
            sep,
            format!(
                "read side is {} columns, reference side is {}. Both sides describe the same \
                 columns, so they must be the same length.",
                read.len(),
                refr.len()
            ),
        ));
    }

    let mut columns = Vec::with_capacity(read.len());
    let mut markers = Markers::default();

    for (i, (r, f)) in read.iter().zip(refr.iter()).enumerate() {
        if f.capture || f.anchor {
            return Err(ParseError::new(
                f.at + sep + 1,
                format!(
                    "'{CAPTURE}' and '{ANCHOR}' go on the read side. They mark a column, and a \
                     column is the same on both sides."
                ),
            ));
        }
        if r.rel.is_some() && f.rel.is_some() {
            return Err(ParseError::new(
                f.at + sep + 1,
                "a relational code goes on one side only; the other side supplies the base set",
            ));
        }
        if r.capture {
            markers.captures.push(i);
        }
        if r.anchor {
            if markers.anchor.is_some() {
                return Err(ParseError::new(r.at, format!("more than one '{ANCHOR}'")));
            }
            markers.anchor = Some(i);
        }
        columns.push(ColumnSpec {
            read: r.seq,
            refr: f.seq,
            rel: r.rel.or(f.rel),
        });
    }

    Ok((
        PatternSpec { columns: columns.into_boxed_slice(), text: s.to_string() },
        markers,
    ))
}

// ---------------------------------------------------------------------------
// Query lexer
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
enum Tok {
    And,
    Or,
    Not,
    Open,
    Close,
    Pattern(String),
}

fn lex(s: &str) -> Result<Vec<(Tok, usize)>, ParseError> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_whitespace() {
            i += 1;
            continue;
        }
        match c {
            b'(' => {
                out.push((Tok::Open, i));
                i += 1;
            }
            b')' => {
                out.push((Tok::Close, i));
                i += 1;
            }
            _ => {
                let start = i;
                while i < b.len()
                    && !b[i].is_ascii_whitespace()
                    && b[i] != b'('
                    && b[i] != b')'
                {
                    i += 1;
                }
                let word = &s[start..i];
                out.push((
                    match word {
                        "and" => Tok::And,
                        "or" => Tok::Or,
                        "not" => Tok::Not,
                        _ => Tok::Pattern(word.to_string()),
                    },
                    start,
                ));
            }
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Query
// ---------------------------------------------------------------------------

impl ColumnSpec {
    /// Whether this column can match a column with no read base -- a deletion
    /// or a skipped intron. Such a column's read symbol is the gap or skip
    /// flag, which matches only a read set containing it, and never a
    /// relational column, which needs a plain base on each side.
    pub fn admits_no_read_base(self) -> bool {
        self.rel.is_none() && (self.read.0 & (Seq::GAP.0 | Seq::SKIP.0)) != 0
    }
}

/// One parsed query, before the automaton is built.
///
/// `expr` indexes `operands` — pattern ids are local to the query here, and are
/// remapped to global ids when the whole query set is compiled.
/// Operands beyond which [`QuerySpec::anchor_can_lack_read_base`] stops
/// enumerating and answers yes: 2^20 evaluations is the most it will spend.
const ANCHOR_CHECK_FREE_OPERANDS: usize = 20;

impl QuerySpec {
    /// Whether this query could fire with its anchor on a column that has no
    /// read base: a deletion or a skipped intron.
    ///
    /// An operand whose anchor column excludes those cannot be true on such a
    /// column. So if the expression cannot be true with all of those operands
    /// false -- every other operand free -- no firing has its anchor there.
    /// That is checked by trying every assignment of the free operands. It can
    /// answer yes for a query that could never actually fire that way, when
    /// operands overlap in ways the check does not model; it never answers no
    /// for one that could.
    pub fn anchor_can_lack_read_base(&self) -> bool {
        let free: Vec<usize> = (0..self.operands.len())
            .filter(|&i| self.operands[i].columns[self.anchor].admits_no_read_base())
            .collect();
        if free.len() > ANCHOR_CHECK_FREE_OPERANDS {
            return true;
        }
        let mut val = vec![false; self.operands.len()];
        (0u64..1u64 << free.len()).any(|mask| {
            for (bit, &i) in free.iter().enumerate() {
                val[i] = (mask >> bit) & 1 == 1;
            }
            eval_expr(&self.expr, &val)
        })
    }
}

/// What a query firing says about the reference symbol under its anchor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnchorRefr {
    /// Always this symbol, in the orientation the walk reports it.
    Fixed(Seq),
    /// Always the same base as the read's.
    SameAsRead,
    /// The query does not determine it.
    Unknown,
}

impl QuerySpec {
    /// The reference symbol under the anchor whenever this query fires, as far
    /// as the query alone determines it.
    ///
    /// A firing makes some set of operands true, and the observed reference
    /// symbol must then lie inside every true operand's reference set at the
    /// anchor column. So every assignment of the operands that makes the
    /// expression true is tried: if each pins the symbol to one and the same
    /// single-bit set, that symbol is [`Fixed`](AnchorRefr::Fixed); if each has
    /// an `=` column there, it is [`SameAsRead`](AnchorRefr::SameAsRead).
    /// Anything else -- an assignment with no true operand, a set of several
    /// symbols, assignments that disagree -- is [`Unknown`](AnchorRefr::Unknown).
    ///
    /// Sound: the actual truth values of any real firing are one of the
    /// assignments tried, so a determined answer is the one a walk would have
    /// observed. It can answer `Unknown` where overlapping operands would in
    /// fact pin the symbol; never the reverse.
    pub fn anchor_refr(&self) -> AnchorRefr {
        let n = self.operands.len();
        if n > ANCHOR_CHECK_FREE_OPERANDS {
            return AnchorRefr::Unknown;
        }
        let mut val = vec![false; n];
        let mut answer: Option<AnchorRefr> = None;
        for mask in 0u64..1u64 << n {
            for (i, v) in val.iter_mut().enumerate() {
                *v = (mask >> i) & 1 == 1;
            }
            if !eval_expr(&self.expr, &val) {
                continue;
            }
            let mut set = 0xFFu8;
            let (mut any, mut eq, mut ne) = (false, false, false);
            for (i, _) in val.iter().enumerate().filter(|(_, v)| **v) {
                let c = self.operands[i].columns[self.anchor];
                any = true;
                set &= c.refr.0;
                match c.rel {
                    Some(Rel::Eq) => eq = true,
                    Some(Rel::Ne) => ne = true,
                    None => {}
                }
            }
            // A relational column only matches a plain base on each side.
            if ne {
                set &= Seq::BASES;
            }
            let this = if !any {
                AnchorRefr::Unknown
            } else if eq {
                AnchorRefr::SameAsRead
            } else if set.count_ones() == 1 {
                AnchorRefr::Fixed(Seq(set))
            } else {
                AnchorRefr::Unknown
            };
            match answer {
                _ if this == AnchorRefr::Unknown => return AnchorRefr::Unknown,
                Some(prev) if prev != this => return AnchorRefr::Unknown,
                _ => answer = Some(this),
            }
        }
        answer.unwrap_or(AnchorRefr::Unknown)
    }
}

fn eval_expr(e: &Expr, val: &[bool]) -> bool {
    match e {
        Expr::Const(b) => *b,
        Expr::Pattern(p) => val[*p],
        Expr::Not(x) => !eval_expr(x, val),
        Expr::And(xs) => xs.iter().all(|x| eval_expr(x, val)),
        Expr::Or(xs) => xs.iter().any(|x| eval_expr(x, val)),
        Expr::Count { lo, hi, patterns } => {
            let n = patterns.iter().filter(|&&p| val[p]).count();
            *lo <= n && n <= *hi
        }
    }
}

#[derive(Debug, Clone)]
pub struct QuerySpec {
    pub name: String,
    /// The expression text, name prefix stripped.
    pub text: String,
    pub expr: Expr,
    pub operands: Vec<PatternSpec>,
    /// Column whose coordinates are reported. Defaults to 0.
    pub anchor: usize,
    /// Columns to record, ascending. Union across all operands.
    pub captures: Box<[usize]>,
    /// Columns spanned. Every operand has this length.
    pub span: usize,
    /// Names for the operands, parallel to `operands`.
    ///
    /// Empty for the terse DSL, which has nowhere to write one; the matrix form
    /// names every pattern. When a name exists it is what the explainer shows,
    /// because `cpg and not junc_end` says what the query does and
    /// `[1] and not [2]` makes the reader hold a lookup table in their head.
    pub operand_names: Box<[String]>,
}

impl QuerySpec {
    /// How to refer to operand `i`: its name, or `[n]` when it has none.
    pub fn operand_label(&self, i: usize) -> String {
        match self.operand_names.get(i) {
            Some(n) if !n.is_empty() => n.clone(),
            _ => format!("[{}]", i + 1),
        }
    }

    /// True when every operand carries a name.
    pub fn is_named(&self) -> bool {
        self.operand_names.len() == self.operands.len()
            && self.operand_names.iter().all(|n| !n.is_empty())
    }
}

/// Precedence-climbing parser over the operator tokens, with leaf handling
/// supplied by the caller.
///
/// Both front ends share this: the terse DSL's leaves are `read@refr` patterns,
/// the matrix DSL's are pattern names. Two copies of precedence climbing would
/// be two chances for `and` and `or` to bind differently depending on which
/// file the query came from.
struct Parser<'a, F> {
    toks: &'a [(Tok, usize)],
    pos: usize,
    leaf: F,
    src_len: usize,
}

impl<'a, F> Parser<'a, F>
where
    F: FnMut(&str, usize) -> Result<Expr, ParseError>,
{
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos).map(|(t, _)| t)
    }
    fn at(&self) -> usize {
        self.toks.get(self.pos).map(|(_, a)| *a).unwrap_or(self.src_len)
    }
    fn bump(&mut self) {
        self.pos += 1;
    }

    fn parse_or(&mut self) -> Result<Expr, ParseError> {
        let mut lhs = self.parse_and()?;
        while self.peek() == Some(&Tok::Or) {
            self.bump();
            let rhs = self.parse_and()?;
            lhs = match lhs {
                Expr::Or(mut v) => {
                    v.push(rhs);
                    Expr::Or(v)
                }
                other => Expr::Or(vec![other, rhs]),
            };
        }
        Ok(lhs)
    }

    fn parse_and(&mut self) -> Result<Expr, ParseError> {
        let mut lhs = self.parse_unary()?;
        while self.peek() == Some(&Tok::And) {
            self.bump();
            let rhs = self.parse_unary()?;
            lhs = match lhs {
                Expr::And(mut v) => {
                    v.push(rhs);
                    Expr::And(v)
                }
                other => Expr::And(vec![other, rhs]),
            };
        }
        Ok(lhs)
    }

    fn parse_unary(&mut self) -> Result<Expr, ParseError> {
        if self.peek() == Some(&Tok::Not) {
            self.bump();
            return Ok(Expr::not(self.parse_unary()?));
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<Expr, ParseError> {
        let at = self.at();
        match self.peek().cloned() {
            Some(Tok::Open) => {
                self.bump();
                let e = self.parse_or()?;
                if self.peek() != Some(&Tok::Close) {
                    return Err(ParseError::new(self.at(), "expected ')'"));
                }
                self.bump();
                Ok(e)
            }
            Some(Tok::Pattern(p)) => {
                self.bump();
                (self.leaf)(&p, at)
            }
            Some(Tok::Close) => Err(ParseError::new(at, "unexpected ')'")),
            Some(Tok::And) | Some(Tok::Or) => {
                Err(ParseError::new(at, "expected a pattern, found an operator"))
            }
            Some(Tok::Not) => unreachable!("handled in parse_unary"),
            None => Err(ParseError::new(at, "expected a pattern, found end of input")),
        }
    }
}

/// Parse a boolean expression whose leaves the caller resolves.
///
/// `leaf` receives the token text and its byte offset. It returns the `Expr`
/// that the token denotes. A front end does its own work inside `leaf`. There
/// it keeps one copy of each operand, records a marker, or looks up a name.
pub fn parse_bool_expr<F>(src: &str, mut leaf: F) -> Result<Expr, ParseError>
where
    F: FnMut(&str, usize) -> Result<Expr, ParseError>,
{
    let toks = lex(src)?;
    let mut p = Parser { toks: &toks, pos: 0, leaf: &mut leaf, src_len: src.len() };
    let e = p.parse_or()?;
    if p.pos != toks.len() {
        return Err(ParseError::new(
            p.at(),
            "trailing input; missing 'and' / 'or' between patterns?",
        ));
    }
    Ok(e)
}


/// Shared tail for both front ends: the query-level checks and the anchor rule.
pub(crate) fn finish(
    name: String,
    body: &str,
    body_at: usize,
    expr: Expr,
    operands: Vec<PatternSpec>,
    operand_names: Vec<String>,
    mut captures: Vec<usize>,
    anchor: Option<usize>,
) -> Result<QuerySpec, ParseError> {
    if operands.is_empty() {
        return Err(ParseError::new(body_at, "expression has no patterns"));
    }
    let span = operands[0].len();
    if operands.iter().any(|o| o.len() != span) {
        return Err(ParseError::new(
            body_at,
            format!(
                "operands describe the same window, so they must be the same length:\n{}",
                align_lengths(&operands)
            ),
        ));
    }

    let anchor = anchor.unwrap_or(0);
    if anchor >= span {
        return Err(ParseError::new(body_at, "anchor is past the end of the window"));
    }
    // The anchor is always recorded. Its coordinates are written out anyway, so
    // withholding the base observed there buys nothing, and every real query
    // wants it -- a methylation call is exactly "what was at the anchor".
    if !captures.contains(&anchor) {
        captures.push(anchor);
    }
    captures.retain(|&c| c < span);
    captures.sort_unstable();
    captures.dedup();

    Ok(QuerySpec {
        name,
        text: body.trim().to_string(),
        expr,
        operands,
        operand_names: operand_names.into_boxed_slice(),
        anchor,
        captures: captures.into_boxed_slice(),
        span,
    })
}

/// Check a set of queries assembled from several sources.
///
/// A single query file cannot contain a conflict -- duplicate `pat` and `query`
/// names are rejected as it is parsed -- but a run may combine several files
/// and any number of `-q` strings, and nothing has compared them until now.
///
/// Two rules, both about a name meaning one thing:
///
/// - Query names must be unique, or the output rows cannot be told apart.
/// - A pattern name must denote the same columns everywhere. Re-declaring it
///   identically is fine, and expected when two files share a definition.
pub fn validate_set(specs: &[QuerySpec]) -> Result<(), String> {
    let mut queries: Vec<&str> = Vec::new();
    for q in specs {
        if q.name.is_empty() {
            continue;
        }
        if queries.contains(&q.name.as_str()) {
            return Err(format!(
                "two queries are named '{}'; their rows would be indistinguishable",
                q.name
            ));
        }
        queries.push(&q.name);
    }

    let mut pats: Vec<(&str, &PatternSpec, &str)> = Vec::new();
    for q in specs {
        for (i, op) in q.operands.iter().enumerate() {
            let Some(name) = q.operand_names.get(i).filter(|n| !n.is_empty()) else {
                continue;
            };
            match pats.iter().find(|(n, _, _)| *n == name.as_str()) {
                Some((_, prev, owner)) if prev.columns != op.columns => {
                    return Err(format!(
                        "pattern '{name}' means '{}' in query '{owner}' but '{}' in query \
                         '{}'; a pattern name must denote one thing across the run",
                        prev.text, op.text, q.name
                    ))
                }
                Some(_) => {}
                None => pats.push((name, op, &q.name)),
            }
        }
    }
    Ok(())
}

/// Render every operand with its column count, so a length mismatch shows which
/// one is the odd one out rather than making the reader count `~` characters.
fn align_lengths(ops: &[PatternSpec]) -> String {
    let w = ops.iter().map(|o| o.text.len()).max().unwrap_or(0);
    ops.iter()
        .map(|o| format!("  {:<w$}  {} columns", o.text, o.columns.len(), w = w))
        .collect::<Vec<_>>()
        .join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn every_lowercase_letter_is_an_alias_and_nothing_else_is() {
        for c in 'a'..='z' {
            assert!(is_alias_char(c), "'{c}' should be usable");
            assert!(Aliases::new().insert(c, Seq::T).is_ok(), "'{c}'");
        }
        for c in 'A'..='Z' {
            assert!(!is_alias_char(c), "'{c}' is a base code");
        }
        for c in ['1', '.', '~', '@', '^', '+', '=', '/', '{'] {
            assert!(!is_alias_char(c), "'{c}'");
        }
    }

    #[test]
    fn alias_errors_state_the_rule() {
        let mut a = Aliases::new();
        let e = a.insert('C', Seq::T).unwrap_err();
        assert!(e.contains("uppercase letters are base codes"), "{e}");
        assert!(e.contains("lowercase letters, a through z"), "{e}");

        let e = a.insert('1', Seq::T).unwrap_err();
        assert!(e.contains("only letters"), "{e}");

        a.insert('j', Seq::T).unwrap();
        assert!(a.insert('j', Seq::C).unwrap_err().contains("already an alias"));
    }

    #[test]
    fn mint_order_prefers_letters_that_are_not_base_codes() {
        // A minted `a` would read as a soft-masked A in a grid.
        assert!(ALIAS_MINT_ORDER.starts_with("fijlopq"));
        let mut seen: Vec<char> = ALIAS_MINT_ORDER.chars().collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), 26, "mint order must cover a-z exactly once");
        assert!(seen.iter().all(|c| c.is_ascii_lowercase()));
    }

    #[test]
    fn alias_tables_merge_when_they_agree() {
        let mut a = Aliases::new();
        a.insert('j', Seq::C | Seq::GAP).unwrap();
        let mut b = Aliases::new();
        b.insert('j', Seq::C | Seq::GAP).unwrap();
        b.insert('l', Seq::G | Seq::PAD).unwrap();
        a.merge(&b).unwrap();
        assert_eq!(a.get('l'), Some(Seq::G | Seq::PAD));

        let mut c = Aliases::new();
        c.insert('j', Seq::T).unwrap();
        assert!(a.merge(&c).unwrap_err().contains("in one place"));
    }

    /// A column renders back to the side it was written on, which is what
    /// keeps a relational column from swapping sides on a round trip.
    #[test]
    fn columns_render_back_to_the_side_they_were_written_on() {
        let t = |read: &str, refr: &str| {
            let spec = crate::test_support::query(&format!(
                "[query.q]\nread = \"{read}\"\nrefr = \"{refr}\"\n"
            ));
            let c = spec.operands[0].columns[0];
            let (r, f) = c.to_text();
            format!("{r}@{f}")
        };
        assert_eq!(t("~", "="), "~@=");
        assert_eq!(t("H", "/"), "H@/");
        assert_eq!(t("=", "C"), "=@C");
        assert_eq!(t("/", "G"), "/@G");
        assert_eq!(t("C~", "CG"), "C@C");
    }


    #[test]
    fn rel_semantics() {
        assert!(Rel::Eq.holds(Seq::C, Seq::C));
        assert!(!Rel::Eq.holds(Seq::C, Seq::G));
        assert!(Rel::Ne.holds(Seq::C, Seq::G));
        assert!(!Rel::Ne.holds(Seq::C, Seq::C));
        // Never on gap, pad, or ambiguity — so they are not complements.
        for r in [Rel::Eq, Rel::Ne] {
            assert!(!r.holds(Seq::GAP, Seq::GAP));
            assert!(!r.holds(Seq::PAD, Seq::PAD));
            assert!(!r.holds(Seq::Y, Seq::C));
            assert!(!r.holds(Seq::C, Seq::Y));
        }
    }

    /// A junction is a value like any other: `,` alone, `J` inside a group, and
    /// the two spell the same set.
    /// What a query pins down about the reference base under its anchor,
    /// which is what lets `extract` fill in `refr_base` without a reference.
    #[test]
    fn what_a_query_fixes_about_the_reference_under_its_anchor() {
        use AnchorRefr::*;
        let r = |toml: &str| crate::test_support::query(toml).anchor_refr();
        let one = |read: &str, refr: &str, mark: &str| {
            format!("[query.q]\nread = \"{read}\"\nrefr = \"{refr}\"\nmark = \"{mark}\"\n")
        };
        let two = |r1: &str, f1: &str, r2: &str, f2: &str, w: &str| {
            format!(
                "[pat.a]\nread = \"{r1}\"\nrefr = \"{f1}\"\n[pat.b]\nread = \"{r2}\"\nrefr = \"{f2}\"\n[query.q]\nmark = \"+.\"\nwhere = \"{w}\"\n"
            )
        };
        assert_eq!(r(&one("C~", "CG", "+.")), Fixed(Seq::C));
        assert_eq!(r(&one("C~~", "CHH", "+..")), Fixed(Seq::C));
        assert_eq!(r(&one("~C", "CG", ".+")), Fixed(Seq::G), "the anchor column, not column 0");
        assert_eq!(r(&two("C~", "CG", "~~", "CA", "a and not b")), Fixed(Seq::C));
        assert_eq!(r(&two("~~", "CG", "~~", "YG", "a and b")), Fixed(Seq::C), "C within Y");
        assert_eq!(r(&two("C~", "CG", "T~", "CA", "a or b")), Fixed(Seq::C), "both branches agree");
        assert_eq!(r(&two("C~", "CG", "C~", "TG", "a or b")), Unknown, "branches disagree");
        assert_eq!(r(&one("C~", "YG", "+.")), Unknown, "two possible bases");
        assert_eq!(r(&one("C~", "~G", "+.")), Unknown);
        assert_eq!(r(&two("C~", "CG", "~~", "CA", "not a")), Unknown, "fires with nothing true");
        assert_eq!(r(&one("=~", "~G", "+.")), SameAsRead);
        assert_eq!(r(&one("/~", "CG", "+.")), Fixed(Seq::C));
        assert_eq!(r(&one("/~", "~G", "+.")), Unknown);
    }

}


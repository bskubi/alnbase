//! Recover a read's strand from a rule that the user declares. The program
//! never assumes the strand.
//!
//! The rules ship as query files under `query/strand/`, one per aligner.
//! This module is the engine those files drive: it parses the two
//! `[strand.*]` tables into a [`StrandRule`] and answers, per record, the two
//! things the walk needs.
//!
//! # What a call carries
//!
//! - **the original strand** (`[strand.original]`): the reference strand the
//!   read's sequence came from. For a conversion assay this is the strand that
//!   carried the conversion. It is the direction
//!   [`crate::alignment::walk_alignment`] runs, and it decides which reference
//!   base is informative.
//! - **the aligned direction** (`[strand.aligned]`): whether the read *as
//!   sequenced* runs along the reference forward or backward. This is what
//!   `read_5p` and `read_3p` count from.
//!
//! Both tables take the keys `forward`, `reverse` and `unknown`.
//!
//! The pair gives the strand of origin (OT, CTOT, OB, CTOB) with no further
//! input, because each of the four strands is one combination of the two:
//! see [`StrandCall::origin`]. A rule never states the origin itself.
//!
//! A rule does not have to read FLAG 0x10. That is what makes the aligners
//! whose 0x10 means the converted reference strand (Bismark single-end,
//! BSBolt read 2, BS-Seeker2) ordinary declarations here rather than special
//! cases.
//!
//! # Why no fallback
//!
//! Exactly one key of each table must match each record. No match, or two
//! matches, is an error naming the record -- not a silent default. A rule that
//! does not cover its input is a rule that is wrong about its input, and the
//! failure mode of guessing is calls placed on the wrong strand, which no
//! downstream check would catch. `unknown` is the one escape, and it is
//! written down.

use std::fmt;

use anyhow::{Result, anyhow, bail};
use rust_htslib::bam::record::{Aux, Record};

use crate::record_field::{AuxKind, FieldSpec, RecordField};
use crate::tags::Strand;

// ------------------------------------------------------------ what a call is --

/// A direction along the reference. Both strand tables answer with one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dir {
    Forward,
    Reverse,
}

impl Dir {
    pub fn name(self) -> &'static str {
        match self {
            Dir::Forward => "forward",
            Dir::Reverse => "reverse",
        }
    }

    pub fn is_reverse(self) -> bool {
        matches!(self, Dir::Reverse)
    }
}

/// What a rule decides about one record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StrandCall {
    /// The reference strand the read's sequence came from.
    pub original: Dir,
    /// Whether the read as sequenced runs along the reference forward or
    /// backward. Not FLAG 0x10: for three of the surveyed aligners the bit
    /// carries the converted reference strand instead.
    pub aligned: Dir,
}

impl StrandCall {
    /// The call for a named strand of origin.
    ///
    /// OT and CTOB are copies of the top-strand sequence, so as sequenced they
    /// run forward. OT and CTOT come from the converted top strand, so their
    /// original strand is forward and the informative reference base is C.
    pub fn from_origin(origin: Strand) -> StrandCall {
        let (original, aligned) = match origin {
            Strand::Ot => (Dir::Forward, Dir::Forward),
            Strand::Ctot => (Dir::Forward, Dir::Reverse),
            Strand::Ob => (Dir::Reverse, Dir::Reverse),
            Strand::Ctob => (Dir::Reverse, Dir::Forward),
        };
        StrandCall { original, aligned }
    }

    /// The strand of origin this call is. Each of the four strands is exactly
    /// one combination of the original strand and the aligned direction, so
    /// this is total, and it is the inverse of [`from_origin`](Self::from_origin).
    pub fn origin(self) -> Strand {
        match (self.original, self.aligned) {
            (Dir::Forward, Dir::Forward) => Strand::Ot,
            (Dir::Forward, Dir::Reverse) => Strand::Ctot,
            (Dir::Reverse, Dir::Reverse) => Strand::Ob,
            (Dir::Reverse, Dir::Forward) => Strand::Ctob,
        }
    }

    /// The original strand as the `strand` column writes it: `+` or `-`.
    pub fn sign(self) -> &'static str {
        match self.original {
            Dir::Forward => "+",
            Dir::Reverse => "-",
        }
    }

    /// Whether the walk runs backwards along the reference. The walk emits
    /// columns 5'->3' along the original strand, so a reverse original strand
    /// is walked from the alignment's far end with both sides reverse
    /// complemented.
    pub fn walk_reversed(self) -> bool {
        self.original.is_reverse()
    }

    /// Whether a walk offset counts from the opposite end of the read from the
    /// sequencer.
    ///
    /// SEQ is stored along the reference by every aligner surveyed, so walk
    /// order runs backwards through SEQ exactly when the original strand is
    /// reverse, and sequencing order runs backwards through it exactly when
    /// the read aligns reverse. When those disagree the walk offset has to be
    /// mirrored to give `read_5p`. For a directional library that is precisely
    /// read 2.
    pub fn mirrors_read(self) -> bool {
        self.walk_reversed() != self.aligned.is_reverse()
    }
}

// ------------------------------------------------------- the condition language --

/// One comparison's right-hand side. Equality only: every rule the survey found
/// is expressible with it, and a rule needing more would be a rule whose input
/// should have been repaired upstream.
#[derive(Debug, Clone, PartialEq)]
enum Value {
    Int(i64),
    Str(String),
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Value::Int(n) => write!(f, "{n}"),
            Value::Str(s) => write!(f, "\"{s}\""),
        }
    }
}

/// A condition over one record: `and`, `or`, `not` and parentheses over bare
/// boolean fields and equality comparisons.
#[derive(Debug, Clone, PartialEq)]
enum Cond {
    /// A bare boolean field name: `is_reverse`.
    Bool(RecordField),
    /// `flags == 99`, `XG == "CT"`.
    Eq(RecordField, Value),
    /// `true` or `false`, which name every record or none. `true` is how a rule
    /// says a property holds whatever the record is — an assay with no
    /// conversion walks everything one way — without the tautology
    /// (`is_reverse or not is_reverse`) that saying so otherwise requires.
    Always(bool),
    And(Box<Cond>, Box<Cond>),
    Or(Box<Cond>, Box<Cond>),
    Not(Box<Cond>),
}

impl Cond {
    fn eval(&self, rec: &Record) -> bool {
        match self {
            Cond::Bool(f) => bool_field(*f, rec),
            Cond::Eq(f, want) => eq_field(*f, want, rec),
            Cond::Always(b) => *b,
            Cond::And(a, b) => a.eval(rec) && b.eval(rec),
            Cond::Or(a, b) => a.eval(rec) || b.eval(rec),
            Cond::Not(a) => !a.eval(rec),
        }
    }
}

fn bool_field(f: RecordField, rec: &Record) -> bool {
    use RecordField as F;
    match f {
        F::IsPaired => rec.is_paired(),
        F::IsProperPair => rec.is_proper_pair(),
        F::IsReverse => rec.is_reverse(),
        F::IsMateReverse => rec.is_mate_reverse(),
        F::IsFirstInTemplate => rec.is_first_in_template(),
        F::IsLastInTemplate => rec.is_last_in_template(),
        F::IsSecondary => rec.is_secondary(),
        F::IsSupplementary => rec.is_supplementary(),
        F::IsDuplicate => rec.is_duplicate(),
        F::IsQcFail => rec.is_quality_check_failed(),
        F::IsUnmapped => rec.is_unmapped(),
        F::IsMateUnmapped => rec.is_mate_unmapped(),
        // Rejected at parse time; unreachable at run time.
        _ => false,
    }
}

/// A tag absent from the record makes every comparison on it false, so the
/// record falls through to no match and is reported by name rather than being
/// quietly assigned a strand.
fn eq_field(f: RecordField, want: &Value, rec: &Record) -> bool {
    match f {
        RecordField::Flags => match want {
            Value::Int(n) => i64::from(rec.flags()) == *n,
            Value::Str(_) => false,
        },
        RecordField::Aux(tag, kind) => {
            let Ok(v) = rec.aux(&tag) else { return false };
            match (kind, want) {
                (AuxKind::Str, Value::Str(s)) => aux_str(&v).as_deref() == Some(s.as_str()),
                (AuxKind::Int, Value::Int(n)) => aux_int(&v) == Some(*n),
                (AuxKind::Float, _) => false,
                // A tag declared one way and compared the other is a rule that
                // cannot fire; caught at parse time, so this is belt and braces.
                _ => false,
            }
        }
        _ => false,
    }
}

fn aux_str(v: &Aux<'_>) -> Option<String> {
    match v {
        Aux::String(s) => Some((*s).to_string()),
        Aux::Char(c) => Some((*c as char).to_string()),
        Aux::HexByteArray(s) => Some((*s).to_string()),
        _ => None,
    }
}

fn aux_int(v: &Aux<'_>) -> Option<i64> {
    match v {
        Aux::I8(x) => Some(i64::from(*x)),
        Aux::U8(x) => Some(i64::from(*x)),
        Aux::I16(x) => Some(i64::from(*x)),
        Aux::U16(x) => Some(i64::from(*x)),
        Aux::I32(x) => Some(i64::from(*x)),
        Aux::U32(x) => Some(i64::from(*x)),
        _ => None,
    }
}

// ------------------------------------------------------------------- parsing --

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Ident(String),
    Str(String),
    Int(i64),
    EqEq,
    LParen,
    RParen,
}

fn tokenize(src: &str) -> Result<Vec<Tok>> {
    let b: Vec<char> = src.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c.is_whitespace() {
            i += 1;
        } else if c == '(' {
            out.push(Tok::LParen);
            i += 1;
        } else if c == ')' {
            out.push(Tok::RParen);
            i += 1;
        } else if c == '=' {
            if i + 1 < b.len() && b[i + 1] == '=' {
                out.push(Tok::EqEq);
                i += 2;
            } else {
                bail!("`=` alone is not a comparison in `{src}`: write `==`");
            }
        } else if c == '"' || c == '\'' {
            let quote = c;
            let start = i + 1;
            let mut j = start;
            while j < b.len() && b[j] != quote {
                j += 1;
            }
            if j >= b.len() {
                bail!("unterminated string in `{src}`");
            }
            out.push(Tok::Str(b[start..j].iter().collect()));
            i = j + 1;
        } else if c.is_ascii_digit() || (c == '-' && i + 1 < b.len() && b[i + 1].is_ascii_digit()) {
            let start = i;
            i += 1;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            let text: String = b[start..i].iter().collect();
            let n = text.parse::<i64>().map_err(|_| anyhow!("`{text}` is not a number in `{src}`"))?;
            out.push(Tok::Int(n));
        } else if c.is_ascii_alphanumeric() || c == '_' || c == ':' {
            let start = i;
            while i < b.len() && (b[i].is_ascii_alphanumeric() || b[i] == '_' || b[i] == ':') {
                i += 1;
            }
            out.push(Tok::Ident(b[start..i].iter().collect()));
        } else {
            bail!("unexpected character `{c}` in `{src}`");
        }
    }
    Ok(out)
}

struct Parser<'a> {
    toks: Vec<Tok>,
    at: usize,
    src: &'a str,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.at)
    }

    fn eat_keyword(&mut self, kw: &str) -> bool {
        if let Some(Tok::Ident(s)) = self.peek()
            && s.eq_ignore_ascii_case(kw)
        {
            self.at += 1;
            return true;
        }
        false
    }

    /// `or` binds loosest, then `and`, then `not` — the same precedence `where`
    /// gives them in a query file, so a condition reads the way a `where` reads.
    fn or(&mut self) -> Result<Cond> {
        let mut left = self.and()?;
        while self.eat_keyword("or") {
            let right = self.and()?;
            left = Cond::Or(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn and(&mut self) -> Result<Cond> {
        let mut left = self.unary()?;
        while self.eat_keyword("and") {
            let right = self.unary()?;
            left = Cond::And(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn unary(&mut self) -> Result<Cond> {
        if self.eat_keyword("not") {
            return Ok(Cond::Not(Box::new(self.unary()?)));
        }
        self.primary()
    }

    fn primary(&mut self) -> Result<Cond> {
        match self.peek().cloned() {
            Some(Tok::LParen) => {
                self.at += 1;
                let inner = self.or()?;
                match self.peek() {
                    Some(Tok::RParen) => {
                        self.at += 1;
                        Ok(inner)
                    }
                    _ => bail!("unclosed `(` in `{}`", self.src),
                }
            }
            Some(Tok::Ident(name)) => {
                if name.eq_ignore_ascii_case("and")
                    || name.eq_ignore_ascii_case("or")
                    || name.eq_ignore_ascii_case("not")
                {
                    bail!("`{name}` needs a condition on both sides in `{}`", self.src);
                }
                if name.eq_ignore_ascii_case("true") || name.eq_ignore_ascii_case("false") {
                    self.at += 1;
                    if matches!(self.peek(), Some(Tok::EqEq)) {
                        bail!("`{name}` is a condition already, so it takes no `==` in `{}`", self.src);
                    }
                    return Ok(Cond::Always(name.eq_ignore_ascii_case("true")));
                }
                self.at += 1;
                let field = resolve_field(&name)?;
                if matches!(self.peek(), Some(Tok::EqEq)) {
                    self.at += 1;
                    let value = match self.peek().cloned() {
                        Some(Tok::Str(s)) => {
                            self.at += 1;
                            Value::Str(s)
                        }
                        Some(Tok::Int(n)) => {
                            self.at += 1;
                            Value::Int(n)
                        }
                        _ => bail!("`{name} ==` needs a value in `{}`", self.src),
                    };
                    check_comparable(&name, field, &value)?;
                    Ok(Cond::Eq(field, value))
                } else {
                    check_boolean(&name, field)?;
                    Ok(Cond::Bool(field))
                }
            }
            Some(t) => bail!("unexpected {t:?} in `{}`", self.src),
            None => bail!("`{}` is not a condition", self.src),
        }
    }
}

/// A bare two-character alphanumeric name is an aux tag read as text, which is
/// how every shipped rule spells one: `XG == "CT"`. A type code may be given
/// explicitly (`XG:A`) with the spelling [`RecordField::parse`] already accepts.
fn resolve_field(name: &str) -> Result<RecordField> {
    let b = name.as_bytes();
    let is_bare_tag = b.len() == 2
        && b.iter().all(|c| c.is_ascii_alphanumeric())
        && !name.eq_ignore_ascii_case("or");
    if is_bare_tag && RecordField::parse(name).is_err() {
        return Ok(RecordField::Aux([b[0], b[1]], AuxKind::Str));
    }
    match RecordField::parse(name) {
        Ok(FieldSpec::One(f)) => Ok(f),
        Ok(FieldSpec::Group(_)) => Err(anyhow!(
            "`{name}` is a group of columns, not a condition: name one field, or an aux tag"
        )),
        Err(e) => Err(anyhow!("{e}")),
    }
}

/// The FLAG bits a condition can name on its own. Every other yes/no column
/// (`read_reverse`) is derived from the strand call, so a rule cannot read it.
pub const FLAG_BITS: &[RecordField] = &[
    RecordField::IsPaired,
    RecordField::IsProperPair,
    RecordField::IsUnmapped,
    RecordField::IsMateUnmapped,
    RecordField::IsReverse,
    RecordField::IsMateReverse,
    RecordField::IsFirstInTemplate,
    RecordField::IsLastInTemplate,
    RecordField::IsSecondary,
    RecordField::IsSupplementary,
    RecordField::IsQcFail,
    RecordField::IsDuplicate,
];

fn check_boolean(name: &str, f: RecordField) -> Result<()> {
    if FLAG_BITS.contains(&f) {
        Ok(())
    } else if matches!(f.data_type(), arrow::datatypes::DataType::Boolean) {
        Err(anyhow!(
            "`{name}` comes from the strand rule itself, so a strand rule cannot read it: \
             use a FLAG bit such as `is_reverse`"
        ))
    } else {
        Err(anyhow!(
            "`{name}` is not a yes/no field, so it cannot stand alone: compare it, as in `{name} == ...`"
        ))
    }
}

fn check_comparable(name: &str, f: RecordField, v: &Value) -> Result<()> {
    match (f, v) {
        (RecordField::Flags, Value::Int(_)) => Ok(()),
        (RecordField::Flags, Value::Str(_)) => Err(anyhow!("`flags` compares to a number, not text")),
        (RecordField::Aux(_, AuxKind::Str), Value::Str(_)) => Ok(()),
        (RecordField::Aux(_, AuxKind::Int), Value::Int(_)) => Ok(()),
        (RecordField::Aux(_, AuxKind::Str), Value::Int(_)) => {
            Err(anyhow!("`{name}` is read as text, so compare it to text: `{name} == \"{v}\"`"))
        }
        (RecordField::Aux(_, AuxKind::Int), Value::Str(_)) => {
            Err(anyhow!("`{name}` is read as a number, so compare it to a number"))
        }
        (RecordField::Aux(_, AuxKind::Float), _) => {
            Err(anyhow!("`{name}` is read as a float; strand rules compare only text and integers"))
        }
        _ => Err(anyhow!(
            "`{name}` cannot be compared in a strand rule: only `flags` and aux tags can"
        )),
    }
}

fn parse_cond(src: &str) -> Result<Cond> {
    let toks = tokenize(src)?;
    if toks.is_empty() {
        bail!("empty condition");
    }
    let mut p = Parser { toks, at: 0, src };
    let cond = p.or()?;
    if p.at != p.toks.len() {
        bail!("trailing text in `{src}`");
    }
    Ok(cond)
}

// ---------------------------------------------------------------- the rule --

/// What a key of either strand table names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Key {
    Dir(Dir),
    Unknown,
}

impl Key {
    fn parse(table: &str, key: &str) -> Result<Key> {
        match key {
            "forward" => Ok(Key::Dir(Dir::Forward)),
            "reverse" => Ok(Key::Dir(Dir::Reverse)),
            "unknown" => Ok(Key::Unknown),
            _ => Err(anyhow!(
                "`{key}` is not a key of [strand.{table}]: it takes forward, reverse or unknown"
            )),
        }
    }

    fn name(self) -> String {
        match self {
            Key::Dir(d) => d.name().to_string(),
            Key::Unknown => "unknown".to_string(),
        }
    }
}

/// The compiled `[strand.*]` tables of one run.
///
/// Built by [`StrandRule::from_tables`], which takes each table as the
/// `(key, condition)` pairs the TOML gave, in declaration order. Order is kept
/// so that an error can name the rules in the order the file writes them.
#[derive(Debug, Clone)]
pub struct StrandRule {
    /// The name the rule was loaded from, for error messages.
    source: String,
    original: Vec<(Key, Cond)>,
    aligned: Vec<(Key, Cond)>,
}

impl StrandRule {
    /// Compile the tables of one query file.
    ///
    /// Each table is the pairs as written; an absent table is an empty slice.
    /// Both tables are required. The rule does not know which file it came
    /// from -- the loader does, and says so with [`named`](Self::named), which
    /// only run-time errors need.
    pub fn from_tables(
        original: &[(String, String)],
        aligned: &[(String, String)],
    ) -> Result<StrandRule> {
        match (original.is_empty(), aligned.is_empty()) {
            (true, true) => bail!("a strand rule needs [strand.original] and [strand.aligned]"),
            (true, false) => bail!(
                "a strand rule needs [strand.original] beside [strand.aligned]: \
                 which reference strand a read came from does not follow from its direction"
            ),
            (false, true) => bail!(
                "a strand rule needs [strand.aligned] beside [strand.original]: \
                 that FLAG 0x10 gives the aligned direction is a claim about the input, \
                 not a default"
            ),
            (false, false) => {}
        }

        let compile = |table: &str, pairs: &[(String, String)]| -> Result<Vec<(Key, Cond)>> {
            let compiled: Vec<(Key, Cond)> = pairs
                .iter()
                .map(|(k, v)| Ok((Key::parse(table, k)?, parse_cond(v).map_err(|e| key_err(k, e))?)))
                .collect::<Result<_>>()?;
            reject_identical(
                &format!("strand.{table}"),
                compiled.iter().map(|(k, c)| (k.name(), c)),
            )?;
            Ok(compiled)
        };
        let original = compile("original", original)?;
        let aligned = compile("aligned", aligned)?;

        Ok(StrandRule { source: String::new(), original, aligned })
    }

    /// Name the file this rule was read from, so that a record it cannot cover
    /// says which rule failed to cover it.
    pub fn named(mut self, source: &str) -> StrandRule {
        self.source = source.to_string();
        self
    }

    /// What [`named`](Self::named) was told, or `""` when nothing was. The
    /// loader reads it back to say which file already declared a rule when a
    /// second one turns up.
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Every key the rule declares, as `table.key`, in declaration order:
    /// `[strand.original]`'s, then `[strand.aligned]`'s.
    ///
    /// This is what a rule promises to be able to say. A file that declares a
    /// key it never gets asked about is a claim with nothing behind it, so the
    /// record fixtures under `query/strand/records/` are required to exercise
    /// each of these at least once.
    pub fn declared_keys(&self) -> Vec<String> {
        let original = self.original.iter().map(|(k, _)| format!("original.{}", k.name()));
        let aligned = self.aligned.iter().map(|(k, _)| format!("aligned.{}", k.name()));
        original.chain(aligned).collect()
    }

    /// The call for one record, or `None` when either table says `unknown`.
    ///
    /// `[strand.original]` is asked first. A record it calls `unknown` is not
    /// put to `[strand.aligned]`, so that table need not cover it.
    ///
    /// An error means the rule did not cover the record: either nothing matched
    /// or more than one key did. Both name the record, because the fix is to
    /// look at it.
    pub fn call(&self, rec: &Record) -> Result<Option<StrandCall>> {
        let Key::Dir(original) = pick(&self.original, rec, "strand.original", &self.source)? else {
            return Ok(None);
        };
        let Key::Dir(aligned) = pick(&self.aligned, rec, "strand.aligned", &self.source)? else {
            return Ok(None);
        };
        Ok(Some(StrandCall { original, aligned }))
    }

    /// The keys a record matches, as [`declared_keys`](Self::declared_keys)
    /// writes them, asked in the order [`call`](Self::call) asks them.
    #[cfg(test)]
    fn matched_keys(&self, rec: &Record) -> Result<Vec<String>> {
        let original = pick(&self.original, rec, "strand.original", &self.source)?;
        let mut keys = vec![format!("original.{}", original.name())];
        if original != Key::Unknown {
            let aligned = pick(&self.aligned, rec, "strand.aligned", &self.source)?;
            keys.push(format!("aligned.{}", aligned.name()));
        }
        Ok(keys)
    }
}

fn key_err(key: &str, e: anyhow::Error) -> anyhow::Error {
    anyhow!("in the rule for `{key}`: {e}")
}

/// Exactly one key must match. Naming every key that matched, rather than the
/// first, is what makes an overlapping pair of conditions a five-second fix.
fn pick(table: &[(Key, Cond)], rec: &Record, table_name: &str, source: &str) -> Result<Key> {
    let mut hits = table.iter().filter(|(_, c)| c.eval(rec));
    let first = hits.next();
    let second = hits.next();
    let qname = String::from_utf8_lossy(rec.qname()).into_owned();
    // The rule names itself only when the loader told it which file it is.
    let source = if source.is_empty() { String::new() } else { format!("{source}: ") };
    match (first, second) {
        (Some((k, _)), None) => Ok(*k),
        (None, _) => Err(anyhow!(
            "{source}no rule in [{table_name}] matches read `{qname}` (FLAG {}). \
             Every record must match exactly one; add a rule, or `unknown` if the \
             input cannot say.",
            rec.flags()
        )),
        (Some((a, _)), Some((b, _))) => Err(anyhow!(
            "{source}read `{qname}` (FLAG {}) matches both `{}` and `{}` in [{table_name}]; \
             exactly one must match.",
            rec.flags(),
            a.name(),
            b.name(),
        )),
    }
}

/// Two keys with the same condition can never both be right, and the per-record
/// check would report every record. Catching it at load names the two keys once.
fn reject_identical<'a>(
    table: &str,
    pairs: impl Iterator<Item = (String, &'a Cond)>,
) -> Result<()> {
    let pairs: Vec<(String, &Cond)> = pairs.collect();
    for i in 0..pairs.len() {
        for j in (i + 1)..pairs.len() {
            if pairs[i].1 == pairs[j].1 {
                bail!(
                    "`{}` and `{}` in [{table}] have the same condition, \
                     so no record could match one and not the other",
                    pairs[i].0,
                    pairs[j].0
                );
            }
        }
    }
    Ok(())
}

/// The shipped `query/strand/directional.toml`, parsed once.
///
/// Tests everywhere in the crate need a strand for the walk, and they take it
/// from the file a user would take it from: a fixture written out in Rust
/// could drift from what alnbase ships, and then the tests would be checking a
/// rule nobody runs.
#[cfg(test)]
pub fn directional() -> std::sync::Arc<StrandRule> {
    use std::sync::{Arc, OnceLock};
    static RULE: OnceLock<Arc<StrandRule>> = OnceLock::new();
    RULE.get_or_init(|| {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("query/strand/directional.toml");
        let src = std::fs::read_to_string(&path).expect("the shipped directional rule");
        let (original, aligned) = tables_of(&src);
        let rule = StrandRule::from_tables(&original, &aligned)
            .expect("the shipped directional rule compiles");
        Arc::new(rule.named(&path.display().to_string()))
    })
    .clone()
}

/// One strand table as `(key, condition)` pairs.
#[cfg(test)]
type Pairs = Vec<(String, String)>;

/// Read the two `[strand.*]` tables out of a query file's text, without the
/// query file's own parser.
#[cfg(test)]
fn tables_of(src: &str) -> (Pairs, Pairs) {
    let v: toml::Value = toml::from_str(src).expect("valid TOML");
    let table = |name: &str| -> Pairs {
        v.get("strand")
            .and_then(|s| s.get(name))
            .and_then(|t| t.as_table())
            .map(|t| {
                t.iter()
                    .map(|(k, v)| (k.clone(), v.as_str().expect("a condition is a string").to_string()))
                    .collect()
            })
            .unwrap_or_default()
    };
    (table("original"), table("aligned"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Aux;

    fn pairs(kv: &[(&str, &str)]) -> Vec<(String, String)> {
        kv.iter().map(|(k, v)| ((*k).to_string(), (*v).to_string())).collect()
    }

    /// The aligned direction as a conformant FLAG gives it.
    const FLAG_ALIGNED: &[(&str, &str)] = &[("forward", "not is_reverse"), ("reverse", "is_reverse")];

    /// The rule for a conformant directional library, which is
    /// `query/strand/directional.toml`.
    fn directional() -> StrandRule {
        StrandRule::from_tables(
            &pairs(&[
                (
                    "forward",
                    "(not is_reverse and not is_last_in_template) or (is_reverse and is_last_in_template)",
                ),
                (
                    "reverse",
                    "(is_reverse and not is_last_in_template) or (not is_reverse and is_last_in_template)",
                ),
            ]),
            &pairs(FLAG_ALIGNED),
        )
        .unwrap()
    }

    fn rec(flags: u16, aux: &[(&[u8; 2], &str)]) -> Record {
        let mut r = Record::new();
        r.set(b"read1", None, b"ACGT", &[40, 40, 40, 40]);
        r.set_flags(flags);
        for (tag, value) in aux {
            r.push_aux(*tag, Aux::String(value)).unwrap();
        }
        r
    }

    #[test]
    fn directional_reproduces_the_hardcoded_rule() {
        // (flags, strand) as `Library::Directional` has always derived them:
        // read 1 forward is OT, read 2 reverse is CTOT, and so on.
        let cases = [
            (0x1 | 0x40, Strand::Ot),
            (0x1 | 0x40 | 0x10, Strand::Ob),
            (0x1 | 0x80 | 0x10, Strand::Ctot),
            (0x1 | 0x80, Strand::Ctob),
            // An unpaired read has is_last_in_template false, so it is read 1.
            (0, Strand::Ot),
            (0x10, Strand::Ob),
        ];
        let rule = directional();
        for (flags, want) in cases {
            let call = rule.call(&rec(flags, &[])).unwrap().unwrap();
            assert_eq!(call.origin(), want, "flags {flags:#x}");
        }
    }

    /// Each strand of origin is one combination of the two directions, and
    /// each combination is one strand. Stated strand by strand so that a change
    /// to either half has to be made here too.
    #[test]
    fn the_origin_and_the_two_directions_determine_each_other() {
        let cases = [
            (Strand::Ot, Dir::Forward, Dir::Forward),
            (Strand::Ctot, Dir::Forward, Dir::Reverse),
            (Strand::Ob, Dir::Reverse, Dir::Reverse),
            (Strand::Ctob, Dir::Reverse, Dir::Forward),
        ];
        for (origin, original, aligned) in cases {
            let call = StrandCall::from_origin(origin);
            assert_eq!((call.original, call.aligned), (original, aligned), "{}", origin.name());
            assert_eq!(StrandCall { original, aligned }.origin(), origin, "{}", origin.name());
        }
    }

    #[test]
    fn a_tag_rule_ignores_a_flag_that_lies() {
        // BSBolt sets 0x10 from the converted reference strand on every *_G2A
        // record, so the FLAG says "reverse" for a read that was sequenced
        // forward. Reading both directions from YS places it correctly anyway.
        let rule = StrandRule::from_tables(
            &pairs(&[
                ("forward", r#"YS == "W_C2T" or YS == "W_G2A""#),
                ("reverse", r#"YS == "C_C2T" or YS == "C_G2A""#),
                ("unknown", r#"YS == "WC""#),
            ]),
            &pairs(&[
                ("forward", r#"YS == "W_C2T" or YS == "C_G2A""#),
                ("reverse", r#"YS == "W_G2A" or YS == "C_C2T""#),
            ]),
        )
        .unwrap();

        let r = rec(0x1 | 0x80 | 0x10, &[(b"YS", "W_G2A")]);
        let call = rule.call(&r).unwrap().unwrap();
        assert_eq!(call.origin(), Strand::Ctot);

        // `unknown` is a call of its own, not an error and not a guess.
        assert!(rule.call(&rec(0x4, &[(b"YS", "WC")])).unwrap().is_none());
    }

    #[test]
    fn bismark_two_tags_together_name_all_four() {
        let rule = StrandRule::from_tables(
            &pairs(&[("forward", r#"XG == "CT""#), ("reverse", r#"XG == "GA""#)]),
            &pairs(&[
                ("forward", r#"(XR == "CT" and XG == "CT") or (XR == "GA" and XG == "GA")"#),
                ("reverse", r#"(XR == "GA" and XG == "CT") or (XR == "CT" and XG == "GA")"#),
            ]),
        )
        .unwrap();
        let cases = [
            (("CT", "CT"), Strand::Ot),
            (("GA", "CT"), Strand::Ctot),
            (("CT", "GA"), Strand::Ob),
            (("GA", "GA"), Strand::Ctob),
        ];
        for ((xr, xg), want) in cases {
            let r = rec(0, &[(b"XR", xr), (b"XG", xg)]);
            assert_eq!(rule.call(&r).unwrap().unwrap().origin(), want, "XR={xr} XG={xg}");
        }
    }

    /// bwa-meth's `YD` names only the original strand, and a strand and its
    /// PCR copy share it. The copy runs the other way along the reference,
    /// which the FLAG says, so the pair still names all four strands.
    #[test]
    fn a_tag_for_the_original_strand_and_the_flag_give_all_four() {
        let rule = StrandRule::from_tables(
            &pairs(&[("forward", r#"YD == "f""#), ("reverse", r#"YD == "r""#)]),
            &pairs(FLAG_ALIGNED),
        )
        .unwrap();
        let cases = [
            (0, "f", Strand::Ot),
            (0x10, "f", Strand::Ctot),
            (0x10, "r", Strand::Ob),
            (0, "r", Strand::Ctob),
        ];
        for (flags, yd, want) in cases {
            let call = rule.call(&rec(flags, &[(b"YD", yd)])).unwrap().unwrap();
            assert_eq!(call.origin(), want, "flags {flags:#x} YD {yd}");
        }
    }

    #[test]
    fn whole_flag_comparisons_work() {
        // asTair's own table is written as whole FLAG values.
        let rule = StrandRule::from_tables(
            &pairs(&[
                ("forward", "flags == 99 or flags == 147 or flags == 0"),
                ("reverse", "flags == 83 or flags == 163 or flags == 16"),
            ]),
            &pairs(FLAG_ALIGNED),
        )
        .unwrap();
        assert_eq!(rule.call(&rec(99, &[])).unwrap().unwrap().original, Dir::Forward);
        assert_eq!(rule.call(&rec(163, &[])).unwrap().unwrap().original, Dir::Reverse);
    }

    #[test]
    fn a_record_no_rule_covers_is_named_not_guessed() {
        let rule =
            StrandRule::from_tables(&pairs(&[("forward", "not is_reverse")]), &pairs(FLAG_ALIGNED))
                .unwrap();
        let err = rule.call(&rec(0x10, &[])).unwrap_err().to_string();
        assert!(err.contains("no rule in [strand.original]"), "{err}");
        assert!(err.contains("read1"), "{err}");
    }

    #[test]
    fn two_matching_rules_name_both() {
        let rule = StrandRule::from_tables(
            &pairs(&[("forward", "not is_reverse"), ("reverse", "not is_last_in_template")]),
            &pairs(FLAG_ALIGNED),
        )
        .unwrap();
        let err = rule.call(&rec(0, &[])).unwrap_err().to_string();
        assert!(err.contains("`forward` and `reverse`"), "{err}");
    }

    #[test]
    fn an_absent_tag_matches_nothing() {
        let rule = StrandRule::from_tables(
            &pairs(&[("forward", r#"XG == "CT""#), ("reverse", r#"XG == "GA""#)]),
            &pairs(FLAG_ALIGNED),
        )
        .unwrap();
        // Not "false, therefore reverse": no match at all, which is reported.
        let err = rule.call(&rec(0, &[])).unwrap_err().to_string();
        assert!(err.contains("no rule"), "{err}");
    }

    #[test]
    fn both_tables_are_required() {
        let original = pairs(&[("forward", "true")]);
        let aligned = pairs(FLAG_ALIGNED);
        for (o, a, want) in [
            (&original[..], &[][..], "needs [strand.aligned] beside"),
            (&[][..], &aligned[..], "needs [strand.original] beside"),
            (&[][..], &[][..], "needs [strand.original] and [strand.aligned]"),
        ] {
            let err = StrandRule::from_tables(o, a).unwrap_err().to_string();
            assert!(err.contains(want), "{err}");
        }
    }

    /// `unknown` in either table skips the record. A record the original table
    /// declines is not put to the aligned table, so that table need not cover
    /// it: an unmapped bwa-meth record has no `YD`, and its FLAG says nothing.
    #[test]
    fn unknown_in_either_table_skips_the_record() {
        let rule = StrandRule::from_tables(
            &pairs(&[
                ("forward", r#"YD == "f""#),
                ("reverse", r#"YD == "r""#),
                ("unknown", "not (YD == \"f\" or YD == \"r\")"),
            ]),
            &pairs(&[("forward", "not is_reverse and not is_supplementary"), ("reverse", "is_reverse")]),
        )
        .unwrap();
        // Original unknown, and the aligned table would have matched nothing.
        assert!(rule.call(&rec(0x800, &[])).unwrap().is_none());

        let rule = StrandRule::from_tables(
            &pairs(&[("forward", "true")]),
            &pairs(&[("forward", "not is_reverse"), ("reverse", "is_reverse and not is_secondary"), ("unknown", "is_reverse and is_secondary")]),
        )
        .unwrap();
        assert!(rule.call(&rec(0x10 | 0x100, &[])).unwrap().is_none());
        assert_eq!(rule.call(&rec(0x10, &[])).unwrap().unwrap().origin(), Strand::Ctot);
    }

    #[test]
    fn two_keys_with_one_condition_are_rejected_at_load() {
        let err = StrandRule::from_tables(
            &pairs(&[("forward", "not is_reverse"), ("reverse", "not is_reverse")]),
            &pairs(FLAG_ALIGNED),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("same condition") && err.contains("[strand.original]"), "{err}");
    }

    #[test]
    fn a_key_that_is_not_a_direction_is_rejected() {
        for (o, a, want) in [
            (pairs(&[("OT", "true")]), pairs(FLAG_ALIGNED), "`OT` is not a key of [strand.original]"),
            (pairs(&[("+", "true")]), pairs(FLAG_ALIGNED), "`+` is not a key of [strand.original]"),
            (pairs(&[("forward", "true")]), pairs(&[("fwd", "true")]), "`fwd` is not a key of [strand.aligned]"),
        ] {
            let err = StrandRule::from_tables(&o, &a).unwrap_err().to_string();
            assert!(err.contains(want), "{err}");
        }
    }

    #[test]
    fn precedence_matches_a_where_clause() {
        // `a or b and c` is `a or (b and c)`, and parens override it.
        let loose = parse_cond("is_reverse or is_paired and is_secondary").unwrap();
        let tight = parse_cond("is_reverse or (is_paired and is_secondary)").unwrap();
        assert_eq!(loose, tight);
        let other = parse_cond("(is_reverse or is_paired) and is_secondary").unwrap();
        assert_ne!(loose, other);
    }

    /// Every FLAG bit a condition accepts is one it can evaluate.
    #[test]
    fn every_flag_bit_parses_and_reads_its_bit() {
        for &f in FLAG_BITS {
            let name = f.static_name().unwrap();
            let cond = parse_cond(name).unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(!cond.eval(&rec(0, &[])), "{name} with no flags set");
            assert!(cond.eval(&rec(0xFFF, &[])), "{name} with every flag set");
        }
    }

    /// `true` is what an assay with no conversion needs: one original strand
    /// for every record, said plainly rather than as a tautology over some
    /// field the rule does not otherwise care about.
    #[test]
    fn true_and_false_name_every_record_and_none() {
        let fwd = rec(0, &[]);
        let rev = rec(16, &[]);
        assert!(parse_cond("true").unwrap().eval(&fwd));
        assert!(parse_cond("true").unwrap().eval(&rev));
        assert!(!parse_cond("false").unwrap().eval(&fwd));
        assert!(!parse_cond("false").unwrap().eval(&rev));
        // It composes like any other condition, and is not a field.
        assert!(!parse_cond("true and false").unwrap().eval(&fwd));
        assert!(parse_cond("not false").unwrap().eval(&fwd));
        assert_eq!(parse_cond("TRUE").unwrap(), parse_cond("true").unwrap());
        assert!(parse_cond("true == 1").unwrap_err().to_string().contains("takes no `==`"));
    }

    #[test]
    fn unparseable_conditions_say_what_is_wrong() {
        let cases = [
            ("is_reverse ==", "needs a value"),
            ("is_reverse and", "is not a condition"),
            ("(is_reverse", "unclosed"),
            ("mapq", "not a yes/no field"),
            ("flags == \"99\"", "not text"),
            (r#"XG == 3"#, "read as text"),
            ("is_reverse = true", "write `==`"),
            ("no_such_field", "unknown field"),
            // Derived from the strand call, so it would read as false on every
            // record and send every record to one key.
            ("read_reverse", "comes from the strand rule itself"),
        ];
        for (src, want) in cases {
            let err = parse_cond(src).unwrap_err().to_string();
            assert!(err.contains(want), "`{src}` gave `{err}`, wanted `{want}`");
        }
    }

    /// Every rule alnbase ships compiles. The files are the interface, so a
    /// file that the engine cannot load is a broken release, not a broken test.
    #[test]
    fn every_shipped_rule_compiles() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("query/strand");
        let mut seen = 0;
        for entry in std::fs::read_dir(&dir).expect("query/strand exists") {
            let path = entry.unwrap().path();
            if path.extension().and_then(|e| e.to_str()) != Some("toml") {
                continue;
            }
            let name = path.file_name().unwrap().to_string_lossy().into_owned();
            let src = std::fs::read_to_string(&path).unwrap();
            let (original, aligned) = tables_of(&src);
            StrandRule::from_tables(&original, &aligned).unwrap_or_else(|e| panic!("{name}: {e}"));
            seen += 1;
        }
        assert_eq!(seen, 12, "query/strand should hold twelve rules");
    }

    /// `directional.toml` is what `--library directional` did, and the point of
    /// shipping it is that replacing the hardcoded rule changes nothing. This
    /// checks the shipped file itself, not the copy written out in this module.
    #[test]
    fn the_shipped_directional_file_matches_the_hardcoded_rule() {
        let rule = super::directional();

        // Every combination of the two bits the old rule read, plus the
        // unpaired forms, compared against `Library::Directional` itself.
        for flags in [0u16, 0x10, 0x1 | 0x40, 0x1 | 0x40 | 0x10, 0x1 | 0x80, 0x1 | 0x80 | 0x10] {
            let r = rec(flags, &[]);
            let want = crate::tags::Library::Directional.strand(&r);
            let got = rule.call(&r).unwrap().unwrap();
            assert_eq!(got.origin(), want, "flags {flags:#x}");
            // And the two values the walk actually uses agree with what the
            // walk derives from the FLAG today.
            let bottom_today = r.is_last_in_template() != r.is_reverse();
            assert_eq!(got.walk_reversed(), bottom_today, "walk direction disagrees at flags {flags:#x}");
            assert_eq!(got.aligned.is_reverse(), r.is_reverse(), "flags {flags:#x}");
            // And the value the written offsets use: mirroring a walk offset
            // into an as-sequenced one was "the record is read 2", which for
            // this rule is what `mirrors_read` comes to.
            assert_eq!(
                got.mirrors_read(),
                r.is_last_in_template(),
                "offset mirroring disagrees at flags {flags:#x}"
            );
        }
    }

    /// `mirrors_read` is the original strand against the aligned direction,
    /// and it is the one thing standing between a walk offset and `read_5p`.
    #[test]
    fn the_offsets_mirror_when_the_walk_runs_against_the_sequencer() {
        for (origin, mirror) in [
            (Strand::Ot, false),   // walked forward, sequenced forward
            (Strand::Ob, false),   // walked backward, sequenced backward
            (Strand::Ctot, true),  // walked forward, sequenced backward
            (Strand::Ctob, true),  // walked backward, sequenced forward
        ] {
            assert_eq!(StrandCall::from_origin(origin).mirrors_read(), mirror, "{origin:?}");
        }
    }

    // ------------------------------------------------ the shipped rules --

    /// The call a fixture record's QNAME asks for. `None` for `unknown`, which
    /// is a rule declining a record rather than a strand.
    ///
    /// The grammar: `unknown`, or `<original>:<aligned>` with each side
    /// `forward` or `reverse`, optionally followed by `/` and a label for the
    /// reader.
    fn expectation(qname: &str) -> Option<StrandCall> {
        let want = qname.split('/').next().unwrap();
        if want == "unknown" {
            return None;
        }
        let dir = |s: &str| match s {
            "forward" => Dir::Forward,
            "reverse" => Dir::Reverse,
            _ => panic!("`{s}` in `{qname}` is not forward or reverse"),
        };
        let (original, aligned) = want
            .split_once(':')
            .unwrap_or_else(|| panic!("`{want}` is not `unknown` or `original:aligned`"));
        Some(StrandCall { original: dir(original), aligned: dir(aligned) })
    }

    /// Every shipped rule is run against records, not only described in prose.
    ///
    /// A file that says Bismark's `XR:Z:GA` with `XG:Z:CT` is CTOT has a record
    /// making that claim, so an edit that breaks it fails here rather than
    /// going on producing calls on the wrong strand. The fixtures are written
    /// from each aligner's source, the same reading the rule itself came from,
    /// so they catch a rule that stops matching its own documentation -- not an
    /// aligner that changes what it writes. Only that aligner's own BAM can
    /// catch the second, which is why the validation demos remain the place for
    /// it.
    #[test]
    fn every_shipped_rule_calls_its_records_as_their_names_say() {
        use rust_htslib::bam::Read as _;

        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("query/strand");
        let mut seen = 0;
        for entry in std::fs::read_dir(&root).unwrap() {
            let path = entry.unwrap().path();
            if path.extension().and_then(|e| e.to_str()) != Some("toml") {
                continue;
            }
            let stem = path.file_stem().unwrap().to_str().unwrap().to_string();
            let text = std::fs::read_to_string(&path).unwrap();
            let rule = crate::query_toml::parse_file(&text)
                .unwrap_or_else(|e| panic!("{stem}.toml: {e}"))
                .strand
                .unwrap_or_else(|| panic!("{stem}.toml declares no [strand.*] table"))
                .named(&format!("{stem}.toml"));

            let sam = root.join("records").join(format!("{stem}.sam"));
            assert!(sam.exists(), "{stem}.toml ships no records at {}", sam.display());
            let mut reader = rust_htslib::bam::Reader::from_path(&sam).unwrap();
            let mut covered: Vec<String> = Vec::new();
            let mut records = 0;
            for rec in reader.records() {
                let rec = rec.unwrap();
                let qname = String::from_utf8(rec.qname().to_vec()).unwrap();
                let got = rule.call(&rec).unwrap_or_else(|e| panic!("{stem}.sam `{qname}`: {e}"));
                assert_eq!(got, expectation(&qname), "{stem}.sam `{qname}`");
                covered.extend(rule.matched_keys(&rec).unwrap());
                records += 1;
            }
            assert!(records > 0, "{stem}.sam has no records");

            // A declared key no record reaches is a claim with nothing behind
            // it, which is the failure this fixture exists to prevent.
            let declared = rule.declared_keys();
            let missing: Vec<&String> = declared.iter().filter(|k| !covered.contains(k)).collect();
            assert!(missing.is_empty(), "{stem}.toml declares {missing:?}, which no record exercises");
            seen += 1;
        }
        assert_eq!(seen, 12, "a rule was added or removed without its records");
    }

    /// The other direction: a records file with no rule is a fixture for
    /// something that no longer ships, and would otherwise sit there passing.
    #[test]
    fn no_records_file_outlives_its_rule() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("query/strand");
        for entry in std::fs::read_dir(root.join("records")).unwrap() {
            let path = entry.unwrap().path();
            if path.extension().and_then(|e| e.to_str()) != Some("sam") {
                continue;
            }
            let stem = path.file_stem().unwrap().to_str().unwrap();
            assert!(
                root.join(format!("{stem}.toml")).exists(),
                "{stem}.sam has no {stem}.toml to check"
            );
        }
    }
}

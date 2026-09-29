//! Recover a read's original strand from a rule that the user declares. The
//! program never assumes the strand.
//!
//! The design is in `docs/design/strand-rules.md`; the rules themselves ship as
//! query files under `queries/strand/`, one per aligner. This module is the
//! engine those files drive: it parses the `[strand.*]` tables into a
//! [`StrandRule`] and answers, per record, the three things the walk needs.
//!
//! # What a call carries
//!
//! - **the conversion strand** ([`ConvStrand`]): which reference strand carried
//!   the conversion. This is the direction [`crate::alignment::walk_alignment`]
//!   runs and which reference base is informative.
//! - **the sequenced direction** ([`SeqDir`]): whether the read *as sequenced*
//!   runs along the reference forward or backward, which is what `off_5p` and
//!   `off_3p` count from.
//! - **the strand of origin** ([`Strand`]), when the input distinguishes all
//!   four. `None` when it does not. bwa-meth, BISCUIT, HISAT-3N and dnmtools
//!   say only which strand was converted, and not which of the two reads on
//!   that strand this record is.
//!
//! The first two follow from the third by biology, not by convention, so a rule
//! that names the origin needs no FLAG at all. That is what makes the aligners
//! whose 0x10 means the converted reference strand (Bismark single-end, BSBolt
//! read 2, BS-Seeker2) ordinary declarations here rather than special cases.
//!
//! # Why no fallback
//!
//! Exactly one key of a table must match each record. No match, or two matches,
//! is an error naming the record — not a silent default. A rule that does not
//! cover its input is a rule that is wrong about its input, and the failure
//! mode of guessing is calls placed on the wrong strand, which no downstream
//! check would catch. `unknown` is the one escape, and it is written down.

use std::fmt;

use anyhow::{Result, anyhow, bail};
use rust_htslib::bam::record::{Aux, Record};

use crate::record_field::{AuxKind, FieldSpec, RecordField};
use crate::tags::Strand;

// ------------------------------------------------------------ what a call is --

/// The reference strand that carried the conversion: the `strand` column's
/// `+`/`-`, and the direction the walk runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvStrand {
    Plus,
    Minus,
}

impl ConvStrand {
    pub fn name(self) -> &'static str {
        match self {
            ConvStrand::Plus => "+",
            ConvStrand::Minus => "-",
        }
    }
}

/// Whether the read as sequenced runs along the reference forward or backward.
///
/// Not FLAG 0x10: for three of the surveyed aligners the bit carries the
/// converted reference strand instead, and this is the value that says so.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqDir {
    Forward,
    Reverse,
}

impl SeqDir {
    pub fn name(self) -> &'static str {
        match self {
            SeqDir::Forward => "forward",
            SeqDir::Reverse => "reverse",
        }
    }

    /// Whether the read runs opposite the reference — the `read_reverse`
    /// column.
    pub fn is_reverse(self) -> bool {
        matches!(self, SeqDir::Reverse)
    }
}

/// What a rule decides about one record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StrandCall {
    pub conversion: ConvStrand,
    pub sequenced: SeqDir,
    /// `None` when the input distinguishes only two of the four strands.
    pub origin: Option<Strand>,
}

impl StrandCall {
    /// The call a named strand of origin implies. See the table in
    /// `docs/design/strand-rules.md`: OT and CTOB are copies of the top-strand
    /// sequence and so read forward; OT and CTOT are the converted top strand
    /// and its complement, so the informative reference base is C.
    pub fn from_origin(origin: Strand) -> StrandCall {
        let (conversion, sequenced) = match origin {
            Strand::Ot => (ConvStrand::Plus, SeqDir::Forward),
            Strand::Ctot => (ConvStrand::Plus, SeqDir::Reverse),
            Strand::Ob => (ConvStrand::Minus, SeqDir::Reverse),
            Strand::Ctob => (ConvStrand::Minus, SeqDir::Forward),
        };
        StrandCall { conversion, sequenced, origin: Some(origin) }
    }

    /// Whether the walk runs backwards along the reference: a, and nothing
    /// else. The walk emits columns 5'->3' along the conversion strand, so a
    /// minus-strand conversion is walked from the alignment's far end with
    /// both sides reverse complemented.
    pub fn walk_reversed(self) -> bool {
        matches!(self.conversion, ConvStrand::Minus)
    }

    /// Whether a walk offset counts from the opposite end of the read from the
    /// sequencer: a against b.
    ///
    /// SEQ is stored along the reference by every aligner surveyed, so walk
    /// order runs backwards through SEQ exactly when the conversion strand is
    /// minus, and sequencing order runs backwards through it exactly when the
    /// read is reverse. When those disagree the walk offset has to be mirrored
    /// to give `off_5p`. For a directional library that is precisely read 2,
    /// which is the rule this replaces.
    pub fn mirrors_read(self) -> bool {
        self.walk_reversed() != self.sequenced.is_reverse()
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

fn check_boolean(name: &str, f: RecordField) -> Result<()> {
    if matches!(f.data_type(), arrow::datatypes::DataType::Boolean) {
        Ok(())
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

/// What a key of `[strand.origin]` names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OriginKey {
    /// One of the four strands.
    Origin(Strand),
    /// `"+"`/`"-"`: the conversion strand is known, the origin is not. The
    /// `[strand.sequenced]` table then applies to exactly these records.
    Conversion(ConvStrand),
    Unknown,
}

impl OriginKey {
    fn parse(key: &str) -> Result<OriginKey> {
        if let Some(s) = Strand::from_name(key) {
            return Ok(OriginKey::Origin(s));
        }
        match key {
            "+" => Ok(OriginKey::Conversion(ConvStrand::Plus)),
            "-" => Ok(OriginKey::Conversion(ConvStrand::Minus)),
            "unknown" => Ok(OriginKey::Unknown),
            _ => Err(anyhow!(
                "`{key}` is not a strand: [strand.origin] takes OT, CTOT, OB, CTOB, \
                 \"+\", \"-\" or unknown"
            )),
        }
    }

    fn name(self) -> String {
        match self {
            OriginKey::Origin(s) => s.name().to_string(),
            OriginKey::Conversion(c) => c.name().to_string(),
            OriginKey::Unknown => "unknown".to_string(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConvKey {
    Conversion(ConvStrand),
    Unknown,
}

impl ConvKey {
    fn parse(key: &str) -> Result<ConvKey> {
        match key {
            "+" => Ok(ConvKey::Conversion(ConvStrand::Plus)),
            "-" => Ok(ConvKey::Conversion(ConvStrand::Minus)),
            "unknown" => Ok(ConvKey::Unknown),
            _ => Err(anyhow!(
                "`{key}` is not a conversion strand: [strand.conversion] takes \"+\", \"-\" or unknown"
            )),
        }
    }

    fn name(self) -> String {
        match self {
            ConvKey::Conversion(c) => c.name().to_string(),
            ConvKey::Unknown => "unknown".to_string(),
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
    origin: Vec<(OriginKey, Cond)>,
    conversion: Vec<(ConvKey, Cond)>,
    sequenced: Vec<(SeqDir, Cond)>,
}

impl StrandRule {
    /// Compile the tables of one query file.
    ///
    /// Each table is the pairs as written; an absent table is an empty slice.
    /// The rule does not know which file it came from — the loader does, and
    /// says so with [`named`](Self::named), which only run-time errors need.
    pub fn from_tables(
        origin: &[(String, String)],
        conversion: &[(String, String)],
        sequenced: &[(String, String)],
    ) -> Result<StrandRule> {
        if origin.is_empty() && conversion.is_empty() {
            bail!("a strand rule needs [strand.origin] or [strand.conversion]");
        }
        if !origin.is_empty() && !conversion.is_empty() {
            bail!(
                "[strand.origin] and [strand.conversion] are alternatives — \
                 the conversion strand follows from the origin, and a second declaration \
                 of it can only disagree"
            );
        }

        let origin: Vec<(OriginKey, Cond)> = origin
            .iter()
            .map(|(k, v)| Ok((OriginKey::parse(k)?, parse_cond(v).map_err(|e| key_err(k, e))?)))
            .collect::<Result<_>>()?;
        let conversion: Vec<(ConvKey, Cond)> = conversion
            .iter()
            .map(|(k, v)| Ok((ConvKey::parse(k)?, parse_cond(v).map_err(|e| key_err(k, e))?)))
            .collect::<Result<_>>()?;
        let sequenced: Vec<(SeqDir, Cond)> = sequenced
            .iter()
            .map(|(k, v)| {
                let dir = match k.as_str() {
                    "forward" => SeqDir::Forward,
                    "reverse" => SeqDir::Reverse,
                    _ => bail!(
                        "`{k}` is not a direction: [strand.sequenced] takes forward or reverse"
                    ),
                };
                Ok((dir, parse_cond(v).map_err(|e| key_err(k, e))?))
            })
            .collect::<Result<_>>()?;

        // A four-way origin table derives the sequenced direction, so declaring
        // it as well is declaring a theorem: rejected rather than checked.
        let origin_needs_sequenced =
            origin.iter().any(|(k, _)| matches!(k, OriginKey::Conversion(_)));
        if !conversion.is_empty() && sequenced.is_empty() {
            bail!(
                "[strand.conversion] needs [strand.sequenced] beside it — \
                 that the FLAG's 0x10 gives the sequenced direction is a claim about \
                 the input, not a default"
            );
        }
        if origin_needs_sequenced && sequenced.is_empty() {
            bail!(
                "[strand.origin] claims records with \"+\" or \"-\", whose \
                 sequenced direction does not follow — add [strand.sequenced]"
            );
        }
        if !origin.is_empty() && !origin_needs_sequenced && !sequenced.is_empty() {
            bail!(
                "[strand.sequenced] is already implied by [strand.origin]; \
                 remove it, or name the records it applies to with \"+\"/\"-\" keys"
            );
        }

        reject_identical("strand.origin", origin.iter().map(|(k, c)| (k.name(), c)))?;
        reject_identical("strand.conversion", conversion.iter().map(|(k, c)| (k.name(), c)))?;
        reject_identical(
            "strand.sequenced",
            sequenced.iter().map(|(k, c)| (k.name().to_string(), c)),
        )?;

        Ok(StrandRule { source: String::new(), origin, conversion, sequenced })
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

    /// Whether every record this rule classifies gets a strand of origin.
    ///
    /// False for a rule built on `[strand.conversion]`, and for a
    /// `[strand.origin]` table with a `"+"` or `"-"` key: both reach the
    /// conversion strand and stop, which does not say whether a `+` record is
    /// OT or CTOT. A `strand` tag writes one value per strand of origin, so
    /// this is what lets that pairing be refused before the run opens the BAM
    /// rather than on its first record.
    ///
    /// A rule's `unknown` key does not make this false: a record the rule
    /// declines is skipped and counted, never written.
    pub fn names_origin(&self) -> bool {
        !self.origin.is_empty()
            && !self.origin.iter().any(|(k, _)| matches!(k, OriginKey::Conversion(_)))
    }

    /// Every key the rule declares, as the file writes it, in declaration
    /// order: `[strand.origin]`'s or `[strand.conversion]`'s, then
    /// `[strand.sequenced]`'s.
    ///
    /// This is what a rule promises to be able to say. A file that declares a
    /// strand it never gets asked about is a claim with nothing behind it, so
    /// the record fixtures under `queries/strand/records/` are required to
    /// exercise each of these at least once.
    pub fn declared_keys(&self) -> Vec<String> {
        let mut keys: Vec<String> = if self.origin.is_empty() {
            self.conversion.iter().map(|(k, _)| k.name()).collect()
        } else {
            self.origin.iter().map(|(k, _)| k.name()).collect()
        };
        keys.extend(self.sequenced.iter().map(|(d, _)| d.name().to_string()));
        keys
    }

    /// The call for one record, or `None` when the rule says `unknown`.
    ///
    /// An error means the rule did not cover the record: either nothing matched
    /// or more than one key did. Both name the record, because the fix is to
    /// look at it.
    pub fn call(&self, rec: &Record) -> Result<Option<StrandCall>> {
        if !self.origin.is_empty() {
            let key = pick(&self.origin, rec, |k| k.name(), "strand.origin", &self.source)?;
            return match key {
                OriginKey::Unknown => Ok(None),
                OriginKey::Origin(s) => Ok(Some(StrandCall::from_origin(s))),
                OriginKey::Conversion(c) => {
                    let dir = self.sequenced_of(rec)?;
                    Ok(Some(StrandCall { conversion: c, sequenced: dir, origin: None }))
                }
            };
        }

        let key = pick(&self.conversion, rec, |k| k.name(), "strand.conversion", &self.source)?;
        match key {
            ConvKey::Unknown => Ok(None),
            ConvKey::Conversion(c) => {
                let dir = self.sequenced_of(rec)?;
                Ok(Some(StrandCall { conversion: c, sequenced: dir, origin: None }))
            }
        }
    }

    fn sequenced_of(&self, rec: &Record) -> Result<SeqDir> {
        pick(&self.sequenced, rec, |d| d.name().to_string(), "strand.sequenced", &self.source)
    }
}

fn key_err(key: &str, e: anyhow::Error) -> anyhow::Error {
    anyhow!("in the rule for `{key}`: {e}")
}

/// Exactly one key must match. Naming every key that matched, rather than the
/// first, is what makes an overlapping pair of conditions a five-second fix.
fn pick<K: Copy>(
    table: &[(K, Cond)],
    rec: &Record,
    name: impl Fn(K) -> String,
    table_name: &str,
    source: &str,
) -> Result<K> {
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
            name(*a),
            name(*b),
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

/// The shipped `queries/strand/directional.toml`, parsed once.
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
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("queries/strand/directional.toml");
        let src = std::fs::read_to_string(&path).expect("the shipped directional rule");
        let v: toml::Value = toml::from_str(&src).expect("valid TOML");
        let table = |name: &str| -> Vec<(String, String)> {
            v.get("strand")
                .and_then(|s| s.get(name))
                .and_then(|t| t.as_table())
                .map(|t| {
                    t.iter()
                        .map(|(k, v)| {
                            (k.clone(), v.as_str().expect("a condition is a string").to_string())
                        })
                        .collect()
                })
                .unwrap_or_default()
        };
        let rule = StrandRule::from_tables(&table("origin"), &table("conversion"), &table("sequenced"))
            .expect("the shipped directional rule compiles");
        Arc::new(rule.named(&path.display().to_string()))
    })
    .clone()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Aux;

    fn pairs(kv: &[(&str, &str)]) -> Vec<(String, String)> {
        kv.iter().map(|(k, v)| ((*k).to_string(), (*v).to_string())).collect()
    }

    /// The four-way table every conformant directional aligner gets, which is
    /// `queries/strand/directional.toml`.
    fn directional() -> StrandRule {
        StrandRule::from_tables(
            &pairs(&[
                ("OT", "not is_reverse and not is_last_in_template"),
                ("CTOT", "is_reverse and is_last_in_template"),
                ("OB", "is_reverse and not is_last_in_template"),
                ("CTOB", "not is_reverse and is_last_in_template"),
            ]),
            &[],
            &[],
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
            assert_eq!(call.origin, Some(want), "flags {flags:#x}");
        }
    }

    #[test]
    fn the_walk_direction_and_sequenced_direction_follow_from_the_origin() {
        // The table in the design doc, which is biology and not convention.
        let cases = [
            (Strand::Ot, ConvStrand::Plus, SeqDir::Forward),
            (Strand::Ctot, ConvStrand::Plus, SeqDir::Reverse),
            (Strand::Ob, ConvStrand::Minus, SeqDir::Reverse),
            (Strand::Ctob, ConvStrand::Minus, SeqDir::Forward),
        ];
        for (origin, conv, dir) in cases {
            let call = StrandCall::from_origin(origin);
            assert_eq!(call.conversion, conv, "{}", origin.name());
            assert_eq!(call.sequenced, dir, "{}", origin.name());
        }
    }

    #[test]
    fn a_tag_rule_ignores_a_flag_that_lies() {
        // BSBolt sets 0x10 from the converted reference strand on every *_G2A
        // record, so the FLAG says "reverse" for a read that was sequenced
        // forward. Naming the origin from YS places it correctly anyway.
        let rule = StrandRule::from_tables(
            &pairs(&[
                ("OT", r#"YS == "W_C2T""#),
                ("CTOT", r#"YS == "W_G2A""#),
                ("OB", r#"YS == "C_C2T""#),
                ("CTOB", r#"YS == "C_G2A""#),
                ("unknown", r#"YS == "WC""#),
            ]),
            &[],
            &[],
        )
        .unwrap();

        let r = rec(0x1 | 0x80 | 0x10, &[(b"YS", "W_G2A")]);
        let call = rule.call(&r).unwrap().unwrap();
        assert_eq!(call.origin, Some(Strand::Ctot));
        assert_eq!(call.conversion, ConvStrand::Plus);
        assert_eq!(call.sequenced, SeqDir::Reverse);

        // `unknown` is a call of its own, not an error and not a guess.
        assert!(rule.call(&rec(0x4, &[(b"YS", "WC")])).unwrap().is_none());
    }

    #[test]
    fn bismark_two_tags_together_name_all_four() {
        let rule = StrandRule::from_tables(
            &pairs(&[
                ("OT", r#"XR == "CT" and XG == "CT""#),
                ("CTOT", r#"XR == "GA" and XG == "CT""#),
                ("OB", r#"XR == "CT" and XG == "GA""#),
                ("CTOB", r#"XR == "GA" and XG == "GA""#),
            ]),
            &[],
            &[],
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
            assert_eq!(rule.call(&r).unwrap().unwrap().origin, Some(want), "XR={xr} XG={xg}");
        }
    }

    #[test]
    fn a_two_way_rule_reports_the_origin_as_unknown() {
        let rule = StrandRule::from_tables(
            &[],
            &pairs(&[("+", r#"YD == "f""#), ("-", r#"YD == "r""#)]),
            &pairs(&[("forward", "not is_reverse"), ("reverse", "is_reverse")]),
        )
        .unwrap();
        let call = rule.call(&rec(0x10, &[(b"YD", "r")])).unwrap().unwrap();
        assert_eq!(call.conversion, ConvStrand::Minus);
        assert_eq!(call.sequenced, SeqDir::Reverse);
        assert_eq!(call.origin, None, "YD says which strand was converted, not which read this is");
    }

    #[test]
    fn whole_flag_comparisons_work() {
        // asTair's own table is written as whole FLAG values.
        let rule = StrandRule::from_tables(
            &[],
            &pairs(&[
                ("+", "flags == 99 or flags == 147 or flags == 0"),
                ("-", "flags == 83 or flags == 163 or flags == 16"),
            ]),
            &pairs(&[("forward", "not is_reverse"), ("reverse", "is_reverse")]),
        )
        .unwrap();
        assert_eq!(rule.call(&rec(99, &[])).unwrap().unwrap().conversion, ConvStrand::Plus);
        assert_eq!(rule.call(&rec(163, &[])).unwrap().unwrap().conversion, ConvStrand::Minus);
    }

    #[test]
    fn a_record_no_rule_covers_is_named_not_guessed() {
        let rule = StrandRule::from_tables(
            &pairs(&[("OT", "not is_reverse")]),
            &[],
            &[],
        )
        .unwrap();
        let err = rule.call(&rec(0x10, &[])).unwrap_err().to_string();
        assert!(err.contains("no rule"), "{err}");
        assert!(err.contains("read1"), "{err}");
    }

    #[test]
    fn two_matching_rules_name_both() {
        let rule = StrandRule::from_tables(
            &pairs(&[("OT", "not is_reverse"), ("OB", "not is_last_in_template")]),
            &[],
            &[],
        )
        .unwrap();
        let err = rule.call(&rec(0, &[])).unwrap_err().to_string();
        assert!(err.contains("OT") && err.contains("OB"), "{err}");
    }

    #[test]
    fn an_absent_tag_matches_nothing() {
        let rule = StrandRule::from_tables(
            &pairs(&[("OT", r#"XG == "CT""#), ("OB", r#"XG == "GA""#)]),
            &[],
            &[],
        )
        .unwrap();
        // Not "false, therefore OB": no match at all, which is reported.
        let err = rule.call(&rec(0, &[])).unwrap_err().to_string();
        assert!(err.contains("no rule"), "{err}");
    }

    #[test]
    fn the_two_forms_are_alternatives() {
        let err = StrandRule::from_tables(
            &pairs(&[("OT", "not is_reverse")]),
            &pairs(&[("+", "not is_reverse")]),
            &[],
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("alternatives"), "{err}");
    }

    #[test]
    fn conversion_without_sequenced_is_an_error() {
        let err = StrandRule::from_tables(
            &[],
            &pairs(&[("+", r#"YD == "f""#), ("-", r#"YD == "r""#)]),
            &[],
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("[strand.sequenced]"), "{err}");
    }

    #[test]
    fn a_four_way_origin_table_may_not_also_declare_the_sequenced_direction() {
        let err = StrandRule::from_tables(
            &pairs(&[
                ("OT", "not is_reverse and not is_last_in_template"),
                ("CTOT", "is_reverse and is_last_in_template"),
                ("OB", "is_reverse and not is_last_in_template"),
                ("CTOB", "not is_reverse and is_last_in_template"),
            ]),
            &[],
            &pairs(&[("forward", "not is_reverse"), ("reverse", "is_reverse")]),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("already implied"), "{err}");
    }

    #[test]
    fn a_partly_known_origin_table_needs_the_sequenced_direction() {
        let tables = pairs(&[("OT", r#"XG == "CT""#), ("+", r#"XG == "CT?""#)]);
        let err = StrandRule::from_tables(&tables, &[], &[]).unwrap_err().to_string();
        assert!(err.contains("[strand.sequenced]"), "{err}");

        // With the table present it is accepted, and the "+" records take their
        // direction from it.
        let rule = StrandRule::from_tables(
            &pairs(&[("OT", r#"XG == "CT""#), ("+", r#"XG == "CU""#)]),
            &[],
            &pairs(&[("forward", "not is_reverse"), ("reverse", "is_reverse")]),
        )
        .unwrap();
        let call = rule.call(&rec(0x10, &[(b"XG", "CU")])).unwrap().unwrap();
        assert_eq!(call.conversion, ConvStrand::Plus);
        assert_eq!(call.sequenced, SeqDir::Reverse);
        assert_eq!(call.origin, None);
    }

    #[test]
    fn two_keys_with_one_condition_are_rejected_at_load() {
        let err = StrandRule::from_tables(
            &pairs(&[("OT", "not is_reverse"), ("OB", "not is_reverse")]),
            &[],
            &[],
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("same condition"), "{err}");
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

    /// `true` is what an assay with no conversion needs: one walk direction for
    /// every record, said plainly rather than as a tautology over some field the
    /// rule does not otherwise care about.
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
        ];
        for (src, want) in cases {
            let err = parse_cond(src).unwrap_err().to_string();
            assert!(err.contains(want), "`{src}` gave `{err}`, wanted `{want}`");
        }
    }

    /// Read the `[strand.*]` tables out of a shipped file, without the query
    /// file's own parser, which does not know about them yet.
    fn tables_of(src: &str) -> (Vec<(String, String)>, Vec<(String, String)>, Vec<(String, String)>) {
        let v: toml::Value = toml::from_str(src).expect("valid TOML");
        let table = |name: &str| -> Vec<(String, String)> {
            v.get("strand")
                .and_then(|s| s.get(name))
                .and_then(|t| t.as_table())
                .map(|t| {
                    t.iter()
                        .map(|(k, v)| {
                            (k.clone(), v.as_str().expect("a condition is a string").to_string())
                        })
                        .collect()
                })
                .unwrap_or_default()
        };
        (table("origin"), table("conversion"), table("sequenced"))
    }

    /// Every rule alnbase ships compiles. The files are the interface, so a
    /// file that the engine cannot load is a broken release, not a broken test.
    #[test]
    fn every_shipped_rule_compiles() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("queries/strand");
        let mut seen = 0;
        for entry in std::fs::read_dir(&dir).expect("queries/strand exists") {
            let path = entry.unwrap().path();
            if path.extension().and_then(|e| e.to_str()) != Some("toml") {
                continue;
            }
            let name = path.file_name().unwrap().to_string_lossy().into_owned();
            let src = std::fs::read_to_string(&path).unwrap();
            let (origin, conversion, sequenced) = tables_of(&src);
            StrandRule::from_tables(&origin, &conversion, &sequenced)
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            seen += 1;
        }
        assert_eq!(seen, 12, "queries/strand should hold twelve rules");
    }

    /// `directional.toml` is what `--library directional` does, and the point of
    /// shipping it is that replacing the hardcoded rule changes nothing. This
    /// checks the shipped file itself, not the copy written out in this module.
    #[test]
    fn the_shipped_directional_file_matches_the_hardcoded_rule() {
        let src = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("queries/strand/directional.toml"),
        )
        .unwrap();
        let (origin, conversion, sequenced) = tables_of(&src);
        let rule = StrandRule::from_tables(&origin, &conversion, &sequenced).unwrap();

        // Every combination of the two bits the old rule read, plus the
        // unpaired forms, compared against `Library::Directional` itself.
        for flags in [0u16, 0x10, 0x1 | 0x40, 0x1 | 0x40 | 0x10, 0x1 | 0x80, 0x1 | 0x80 | 0x10] {
            let r = rec(flags, &[]);
            let want = crate::tags::Library::Directional.strand(&r);
            let got = rule.call(&r).unwrap().unwrap();
            assert_eq!(got.origin, Some(want), "flags {flags:#x}");
            // And the two values the walk actually uses agree with what the
            // walk derives from the FLAG today.
            let bottom_today = r.is_last_in_template() != r.is_reverse();
            assert_eq!(
                got.conversion == ConvStrand::Minus,
                bottom_today,
                "walk direction disagrees at flags {flags:#x}"
            );
            assert_eq!(got.sequenced.is_reverse(), r.is_reverse(), "flags {flags:#x}");
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

    /// `mirrors_read` is a against b, and it is the one thing standing between
    /// a walk offset and `off_5p`. Stated strand by strand so that a change to
    /// either half of the table has to be made here too.
    #[test]
    fn the_offsets_mirror_when_the_walk_runs_against_the_sequencer() {
        use crate::tags::Strand;
        for (origin, mirror) in [
            (Strand::Ot, false),   // walked forward, sequenced forward
            (Strand::Ob, false),   // walked backward, sequenced backward
            (Strand::Ctot, true),  // walked forward, sequenced backward
            (Strand::Ctob, true),  // walked backward, sequenced forward
        ] {
            assert_eq!(StrandCall::from_origin(origin).mirrors_read(), mirror, "{origin:?}");
        }

        // A rule naming only the conversion strand still answers, because both
        // halves were declared; this is the bwa-meth `YD` shape.
        let call = StrandCall {
            conversion: ConvStrand::Plus,
            sequenced: SeqDir::Reverse,
            origin: None,
        };
        assert!(call.mirrors_read());
        assert!(!call.walk_reversed());
    }

    #[test]
    fn a_strand_name_that_is_not_one_is_rejected() {
        let err = StrandRule::from_tables(&pairs(&[("TOP", "is_reverse")]), &[], &[])
            .unwrap_err()
            .to_string();
        assert!(err.contains("is not a strand"), "{err}");
    }

    // ------------------------------------------------ the shipped rules --

    /// The call a fixture record's QNAME asks for, and the rule keys that
    /// asking it exercises. `None` for `unknown`, which is a rule declining a
    /// record rather than a strand.
    ///
    /// The grammar is the one `queries/strand/records/README.md` states:
    /// `OT`, `CTOT`, `OB`, `CTOB`, `unknown`, or `<conversion>:<sequenced>`
    /// for a rule that reaches only the conversion strand, optionally followed
    /// by `/` and a label for the reader.
    fn expectation(qname: &str) -> (Option<StrandCall>, Vec<String>) {
        let want = qname.split('/').next().unwrap();
        if want == "unknown" {
            return (None, vec!["unknown".to_string()]);
        }
        if let Some(origin) = Strand::from_name(want) {
            return (Some(StrandCall::from_origin(origin)), vec![origin.name().to_string()]);
        }
        let (conv, seq) = want
            .split_once(':')
            .unwrap_or_else(|| panic!("`{want}` is not a strand, `unknown`, or `conv:sequenced`"));
        let conversion = match conv {
            "+" => ConvStrand::Plus,
            "-" => ConvStrand::Minus,
            _ => panic!("`{conv}` is not a conversion strand"),
        };
        let sequenced = match seq {
            "forward" => SeqDir::Forward,
            "reverse" => SeqDir::Reverse,
            _ => panic!("`{seq}` is not a sequenced direction"),
        };
        (
            Some(StrandCall { conversion, sequenced, origin: None }),
            vec![conv.to_string(), seq.to_string()],
        )
    }

    /// Every shipped rule is run against records, not only described in prose.
    ///
    /// This is the other half of open question 3 of
    /// `docs/design/strand-rules.md`: a file that says Bismark's `XR:Z:GA` with
    /// `XG:Z:CT` is CTOT now has a record making that claim, so an edit that
    /// breaks it fails here rather than going on producing calls on the wrong
    /// strand. The fixtures are written from each aligner's source, the same
    /// reading the rule itself came from, so they catch a rule that stops
    /// matching its own documentation — not an aligner that changes what it
    /// writes. Only that aligner's own BAM can catch the second, which is why
    /// the validation demos remain the place for it.
    #[test]
    fn every_shipped_rule_calls_its_records_as_their_names_say() {
        use rust_htslib::bam::Read as _;

        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("queries/strand");
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
                let (want, keys) = expectation(&qname);
                let got = rule.call(&rec).unwrap_or_else(|e| panic!("{stem}.sam `{qname}`: {e}"));
                assert_eq!(got, want, "{stem}.sam `{qname}`");
                covered.extend(keys);
                records += 1;
            }
            assert!(records > 0, "{stem}.sam has no records");

            // A declared strand no record asks for is a claim with nothing
            // behind it, which is the failure this fixture exists to prevent.
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
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("queries/strand");
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

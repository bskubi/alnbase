//! Declare which BAM tags a scan writes. The user makes the declaration in a
//! TOML query file.
//!
//! ```toml
//! [tag.XM.bases]              # one character per position, Bismark-style
//! fill = "."
//! z = "TG"
//! Z = "CG"
//!
//! [tag.XR.strand]             # one value per read, by conversion strand
//! CT = ["OT", "OB"]
//! GA = ["CTOT", "CTOB"]
//! ```
//!
//! This module is the validated configuration and nothing else: it knows no
//! TOML, reads no BAM and writes no tag. [`crate::query_toml`] builds it from a
//! file; [`crate::bam_out`] writes the tags it declares, and
//! [`crate::tag_extract`] reads a `bases` tag back into hit rows.
//!
//! # One kind per tag, named by its table
//!
//! A tag's table says what kind of value it holds -- `bases`, `strand` -- and
//! the contents are that kind's settings. New kinds arrive as new table names
//! rather than as flags on existing ones, so a file written today keeps meaning
//! what it meant.
//!
//! # Strands are an input, not a policy
//!
//! A `strand` tag says what to write for a read from each strand of origin. It
//! does not say how that strand is determined, which is the run's strand rule's
//! job ([`crate::strand_rule`]). All four strands exist in every protocol -- a
//! directional paired-end library puts read 1 on OT/OB and read 2 on
//! CTOT/CTOB, a non-directional one puts either read anywhere -- so a mapping
//! must cover all four whatever protocol it is later used with.

use std::fmt;

/// A validated two-character SAM tag name: a letter, then a letter or digit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TagName([u8; 2]);

impl TagName {
    pub fn parse(s: &str) -> Result<Self, String> {
        let b = s.as_bytes();
        let ok = b.len() == 2 && b[0].is_ascii_alphabetic() && b[1].is_ascii_alphanumeric();
        if !ok {
            return Err(format!(
                "'{s}' is not a BAM tag name: a tag is two characters, a letter then a \
                 letter or digit, e.g. XM"
            ));
        }
        Ok(TagName([b[0], b[1]]))
    }

    pub fn bytes(self) -> [u8; 2] {
        self.0
    }
}

impl fmt::Display for TagName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.0[0] as char, self.0[1] as char)
    }
}

/// A bisulfite conversion strand.
///
/// OT and OB are the original top and bottom strands; CTOT and CTOB are the
/// strands complementary to them, which PCR copies of the originals produce.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Ot,
    Ctot,
    Ob,
    Ctob,
}

impl Strand {
    /// Every strand, in the order values are indexed by.
    pub const ALL: [Strand; 4] = [Strand::Ot, Strand::Ctot, Strand::Ob, Strand::Ctob];

    pub fn name(self) -> &'static str {
        match self {
            Strand::Ot => "OT",
            Strand::Ctot => "CTOT",
            Strand::Ob => "OB",
            Strand::Ctob => "CTOB",
        }
    }

    /// Exact, case-sensitive: the names are conventional spellings, and
    /// accepting `ot` would invite `Ot`, `oT` and a lookup table of guesses.
    pub fn from_name(s: &str) -> Option<Strand> {
        Strand::ALL.into_iter().find(|st| st.name() == s)
    }

    fn index(self) -> usize {
        self as usize
    }

    /// Whether a read from this strand reports the bottom strand of the
    /// reference: OB, and CTOB, its complement's copy.
    ///
    /// This is the orientation [`crate::alignment::walk_alignment`] walks a
    /// record in, so a per-base tag can map a walk offset back onto SEQ.
    pub fn is_bottom(self) -> bool {
        matches!(self, Strand::Ob | Strand::Ctob)
    }
}

/// The strand rule alnbase used to compile in, kept only as a golden value to
/// check the shipped `queries/strand/directional.toml` against.
///
/// A run gets its strand from a query file's `[strand.*]` tables and from
/// nowhere else (`docs/design/strand-rules.md`), so this is not a route a user
/// can take. It survives because the evidence that retiring it changed nothing
/// is a test comparing the shipped file to it, record by record, and that test
/// needs something to compare against.
///
/// # Directional
///
/// Read 1 is the converted original strand, read in its own 5'->3' direction,
/// and read 2 its complement: so read 1 aligning forward came from OT and in
/// reverse from OB, and read 2 aligning in reverse came from CTOT and forward
/// from CTOB. A read with no mate is read 1. Nothing but the FLAG is consulted.
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Library {
    Directional,
}

#[cfg(test)]
impl Library {
    /// Everything a caller needs to know about `record`'s strand: which
    /// reference strand carried the conversion, which way the read was
    /// sequenced, and which of the four strands it came from.
    pub fn call(self, record: &rust_htslib::bam::Record) -> crate::strand_rule::StrandCall {
        crate::strand_rule::StrandCall::from_origin(self.strand(record))
    }

    /// The strand `record` is reported against.
    pub fn strand(self, record: &rust_htslib::bam::Record) -> Strand {
        match self {
            Library::Directional => {
                directional_strand(record.is_last_in_template(), record.is_reverse())
            }
        }
    }
}

#[cfg(test)]
fn directional_strand(is_read2: bool, is_reverse: bool) -> Strand {
    match (is_read2, is_reverse) {
        (false, false) => Strand::Ot,
        (false, true) => Strand::Ob,
        (true, true) => Strand::Ctot,
        (true, false) => Strand::Ctob,
    }
}

/// A value that may be written as one string or as a list of them.
///
/// `list` records which, so that a setting taking exactly one value can reject
/// `fill = ["."]` rather than quietly accepting a one-element list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Names {
    pub items: Vec<String>,
    pub list: bool,
}

/// Every tag one run writes, in declaration order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TagConfig {
    pub tags: Vec<TagSpec>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TagSpec {
    pub name: TagName,
    pub kind: TagKind,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagKind {
    Bases(BasesTag),
    Strand(StrandTag),
}

/// One character per position: a code where a query hits, `fill` elsewhere.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BasesTag {
    /// Written at every position no code claims.
    pub fill: u8,
    /// In declaration order.
    pub codes: Vec<BasesCode>,
}

/// A character, and the query whose hits write it.
///
/// One query per character: a character that could stand for either of two
/// queries would say a base was marked without saying how, and nothing reading
/// the tag could tell which.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BasesCode {
    pub code: u8,
    /// A query name, resolved against the file that declared the tag.
    pub query: String,
}

/// One value per read, chosen by the read's conversion strand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StrandTag {
    values: [String; 4],
}

impl StrandTag {
    pub fn value(&self, s: Strand) -> &str {
        &self.values[s.index()]
    }
}

impl TagConfig {
    pub fn is_empty(&self) -> bool {
        self.tags.is_empty()
    }

    /// Fold another file's tags in. One tag written from two declarations
    /// would silently keep whichever came last, so a repeat is an error even
    /// when the two agree.
    pub fn merge(&mut self, other: TagConfig) -> Result<(), String> {
        for t in other.tags {
            if self.tags.iter().any(|mine| mine.name == t.name) {
                return Err(format!("tag {} is declared in more than one query file", t.name));
            }
            self.tags.push(t);
        }
        Ok(())
    }
}

/// Check that every query a bases tag uses has a read base under its anchor
/// whenever it fires; see [`crate::dsl::QuerySpec::anchor_can_lack_read_base`]. `queries`
/// must include every query the tag names.
pub fn check_anchors(tag: &BasesTag, queries: &[crate::dsl::QuerySpec]) -> Result<(), String> {
    for code in &tag.codes {
        let spec = queries.iter().find(|q| q.name == code.query);
        if spec.is_some_and(|q| q.anchor_can_lack_read_base()) {
            return Err(format!(
                "code '{}': query '{}' can fire with its anchor on a deletion or intron, where \
                 the read has no base to mark. Make the anchor column's read side exclude gaps \
                 -- e.g. N rather than ~ -- or AND the query with a pattern that does.",
                code.code as char, code.query
            ));
        }
    }
    Ok(())
}

/// Printable ASCII other than space: what a one-character code in a `Z` tag
/// can be without being invisible in a SAM dump.
fn is_code_char(c: char) -> bool {
    c.is_ascii_graphic()
}

/// Build a `bases` tag from its table's entries.
///
/// One-character keys are codes; the only longer key is `fill`. A code can
/// never be spelled like a setting, so nothing about the table needs to be
/// remembered -- a misspelled setting is reported as unknown.
///
/// `queries` is every query name the file declares.
pub fn bases_tag(entries: &[(String, Names)], queries: &[&str]) -> Result<BasesTag, String> {
    let mut fill: Option<u8> = None;
    let mut codes: Vec<BasesCode> = Vec::new();

    for (key, names) in entries {
        let mut chars = key.chars();
        let single = match (chars.next(), chars.next()) {
            (Some(c), None) => Some(c),
            _ => None,
        };

        let Some(code) = single else {
            if key != "fill" {
                return Err(format!(
                    "'{key}' is neither a one-character code nor a setting; a bases table \
                     takes one-character codes and 'fill'"
                ));
            }
            let value = match (names.list, names.items.as_slice()) {
                (false, [v]) => v,
                _ => return Err("'fill' takes one character, not a list".to_string()),
            };
            let mut vc = value.chars();
            fill = match (vc.next(), vc.next()) {
                (Some(c), None) if is_code_char(c) => Some(c as u8),
                _ => {
                    return Err(format!(
                        "'fill' must be one printable character other than a space, \
                         got '{value}'"
                    ))
                }
            };
            continue;
        };

        if !is_code_char(code) {
            return Err(format!(
                "code '{code}' must be a printable character other than a space"
            ));
        }
        let q = match (names.list, names.items.as_slice()) {
            (false, [q]) => q,
            _ => {
                return Err(format!(
                    "code '{code}' takes one query name, not a list: a character stands for \
                     exactly one query"
                ))
            }
        };
        if !queries.contains(&q.as_str()) {
            return Err(format!(
                "code '{code}': no query named '{q}' (defined: {})",
                if queries.is_empty() { "none".to_string() } else { queries.join(", ") }
            ));
        }
        // A query under two codes would write two characters at the same
        // position. It is the one overlap that is visible without comparing
        // what the queries match.
        if let Some(prev) = codes.iter().find(|c| &c.query == q) {
            return Err(format!(
                "query '{q}' is used by code '{}' and again by code '{code}'",
                prev.code as char
            ));
        }
        codes.push(BasesCode { code: code as u8, query: q.clone() });
    }

    let fill = fill.ok_or_else(|| {
        "no 'fill': say which character to write where no code applies, e.g. fill = \".\""
            .to_string()
    })?;
    if codes.is_empty() {
        return Err("declares no codes".to_string());
    }
    if let Some(c) = codes.iter().find(|c| c.code == fill) {
        return Err(format!(
            "'{}' is both the fill and a code, so the tag could not tell them apart",
            c.code as char
        ));
    }
    Ok(BasesTag { fill, codes })
}

/// Build a `strand` tag from its table's entries: each key is a value to
/// write, each list the strands that write it. Every strand must get exactly
/// one value.
pub fn strand_tag(entries: &[(String, Names)]) -> Result<StrandTag, String> {
    let mut values: [Option<String>; 4] = Default::default();

    for (value, names) in entries {
        if value.is_empty() || !value.chars().all(|c| c == ' ' || c.is_ascii_graphic()) {
            return Err(format!(
                "'{value}' cannot be a tag value: it must be non-empty printable ASCII"
            ));
        }
        if names.items.is_empty() {
            return Err(format!("value '{value}' lists no strands"));
        }
        for n in &names.items {
            let Some(st) = Strand::from_name(n) else {
                let hint = Strand::from_name(&n.to_ascii_uppercase())
                    .map(|st| format!(" (strand names are uppercase: {})", st.name()))
                    .unwrap_or_default();
                return Err(format!(
                    "value '{value}': '{n}' is not a strand; expected OT, CTOT, OB or CTOB{hint}"
                ));
            };
            if let Some(prev) = &values[st.index()] {
                return Err(format!(
                    "strand {} is given value '{prev}' and also '{value}'",
                    st.name()
                ));
            }
            values[st.index()] = Some(value.clone());
        }
    }

    let missing: Vec<&str> = Strand::ALL
        .into_iter()
        .filter(|st| values[st.index()].is_none())
        .map(Strand::name)
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "no value for {}; every strand needs one, since which strands occur depends \
             on the library, not on this file",
            missing.join(", ")
        ));
    }
    Ok(StrandTag { values: values.map(|v| v.expect("checked above")) })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one(s: &str) -> Names {
        Names { items: vec![s.to_string()], list: false }
    }
    fn many(s: &[&str]) -> Names {
        Names { items: s.iter().map(|x| x.to_string()).collect(), list: true }
    }
    fn e(k: &str, v: Names) -> (String, Names) {
        (k.to_string(), v)
    }

    #[test]
    fn tag_names() {
        for good in ["XM", "xm", "X1", "zz"] {
            assert!(TagName::parse(good).is_ok(), "{good}");
        }
        for bad in ["X", "XMM", "1X", "X-", "", "Xé"] {
            assert!(TagName::parse(bad).is_err(), "{bad}");
        }
        assert_eq!(TagName::parse("XM").unwrap().to_string(), "XM");
    }

    #[test]
    fn a_directional_library_reads_the_strand_from_the_flag() {
        assert_eq!(directional_strand(false, false), Strand::Ot);
        assert_eq!(directional_strand(false, true), Strand::Ob);
        assert_eq!(directional_strand(true, true), Strand::Ctot);
        assert_eq!(directional_strand(true, false), Strand::Ctob);
        let bottom: Vec<_> = Strand::ALL.into_iter().filter(|s| s.is_bottom()).collect();
        assert_eq!(bottom, [Strand::Ob, Strand::Ctob]);
    }

    #[test]
    fn strand_names_round_trip_and_are_exact() {
        for st in Strand::ALL {
            assert_eq!(Strand::from_name(st.name()), Some(st));
        }
        assert_eq!(Strand::from_name("ot"), None);
    }

    #[test]
    fn a_chars_tag_keeps_its_declaration() {
        let t = bases_tag(
            &[e("fill", one(".")), e("z", one("TG")), e("Z", one("CG"))],
            &["TG", "CG"],
        )
        .unwrap();
        assert_eq!(t.fill, b'.');
        assert_eq!(t.codes.len(), 2);
        assert_eq!((t.codes[0].code, t.codes[0].query.as_str()), (b'z', "TG"));
        assert_eq!((t.codes[1].code, t.codes[1].query.as_str()), (b'Z', "CG"));
    }

    #[test]
    fn chars_tag_errors() {
        let q = &["a", "b"][..];
        let bad = |entries: &[(String, Names)]| bases_tag(entries, q).unwrap_err();
        assert!(bad(&[e("z", one("a"))]).contains("no 'fill'"));
        assert!(bad(&[e("fill", one("."))]).contains("no codes"));
        assert!(bad(&[e("fill", one("..")), e("z", one("a"))]).contains("one printable"));
        assert!(bad(&[e("fill", one(" ")), e("z", one("a"))]).contains("one printable"));
        assert!(bad(&[e("fill", many(&["."])), e("z", one("a"))]).contains("not a list"));
        assert!(bad(&[e("fil", one(".")), e("z", one("a"))]).contains("neither"));
        assert!(bad(&[e("fill", one(".")), e(".", one("a"))]).contains("both the fill"));
        assert!(bad(&[e("fill", one(".")), e(" ", one("a"))]).contains("printable"));
        assert!(bad(&[e("fill", one(".")), e("z", one("nope"))]).contains("no query named"));
        assert!(bad(&[e("fill", one(".")), e("z", one("a")), e("Z", one("a"))])
            .contains("used by code 'z' and again by code 'Z'"));
        for list in [many(&["a", "b"]), many(&["a"]), many(&[])] {
            assert!(
                bad(&[e("fill", one(".")), e("z", list.clone())]).contains("one query name, not a list"),
                "{list:?}"
            );
        }
    }

    #[test]
    fn bismark_xr_and_xg() {
        let xr = strand_tag(&[e("CT", many(&["OT", "CTOB"])), e("GA", many(&["OB", "CTOT"]))])
            .unwrap();
        let xg = strand_tag(&[e("CT", many(&["OT", "CTOT"])), e("GA", many(&["OB", "CTOB"]))])
            .unwrap();
        let got: Vec<_> = Strand::ALL.iter().map(|&s| (xr.value(s), xg.value(s))).collect();
        // OT, CTOT, OB, CTOB: Bismark's read and genome conversions.
        assert_eq!(got, [("CT", "CT"), ("GA", "CT"), ("GA", "GA"), ("CT", "GA")]);
    }

    #[test]
    fn strand_tag_errors() {
        let bad = |entries: &[(String, Names)]| strand_tag(entries).unwrap_err();
        let msg = bad(&[e("CT", many(&["OT", "CTOB"]))]);
        assert!(msg.contains("no value for CTOT, OB"), "{msg}");
        let msg = bad(&[e("CT", many(&["OT", "ot"]))]);
        assert!(msg.contains("uppercase: OT"), "{msg}");
        let msg = bad(&[e("CT", many(&["OT", "CTOT", "OB", "CTOB"])), e("GA", one("OB"))]);
        assert!(msg.contains("strand OB is given value 'CT' and also 'GA'"), "{msg}");
        assert!(bad(&[e("", many(&["OT", "CTOT", "OB", "CTOB"]))]).contains("non-empty"));
        assert!(bad(&[e("CT", many(&[]))]).contains("lists no strands"));
        // One value for every strand is legitimate, if unusual.
        assert!(strand_tag(&[e("x", many(&["OT", "CTOT", "OB", "CTOB"]))]).is_ok());
    }

    #[test]
    fn merging_rejects_a_tag_declared_twice() {
        let xr = || TagSpec {
            name: TagName::parse("XR").unwrap(),
            kind: TagKind::Strand(
                strand_tag(&[e("CT", many(&["OT", "CTOT", "OB", "CTOB"]))]).unwrap(),
            ),
        };
        let mut a = TagConfig { tags: vec![xr()] };
        let msg = a.merge(TagConfig { tags: vec![xr()] }).unwrap_err();
        assert!(msg.contains("tag XR is declared in more than one"), "{msg}");
    }
}

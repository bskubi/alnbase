//! Render a query in English.
//!
//! The grid shows what each operand constrains and the tree shows how the
//! operands combine, but both still require reading the DSL. This says it in
//! words, which is what you want when checking a query against a protocol
//! description, pasting it into a methods section, or handing it to someone who
//! has never seen the syntax.
//!
//! Two rules keep it readable rather than exhaustive:
//!
//! - Fully unconstrained columns (`~@~`) are not mentioned. In a padded window
//!   most columns are free, and listing them buries the two that matter.
//! - A side constrained at a single column of a run is named by that column
//!   ("the read is C or T at column 3") rather than spelled as a sequence with
//!   "anything" padding it out.

use crate::dsl::{ColumnSpec, PatternSpec, QuerySpec, Rel};
use crate::predicate::Expr;
use crate::seq::Seq;

/// "A", "A or C", "A, C or T".
fn join_or(parts: &[String]) -> String {
    match parts.len() {
        0 => String::new(),
        1 => parts[0].clone(),
        2 => format!("{} or {}", parts[0], parts[1]),
        _ => format!("{} or {}", parts[..parts.len() - 1].join(", "), parts[parts.len() - 1]),
    }
}

/// A base set in words. `~` is "anything"; `N` is "any base"; anything else is
/// spelled out, because "H" means nothing to a reader who is not already
/// fluent in IUPAC.
pub fn describe(s: Seq) -> String {
    if s == Seq::ANY {
        return "anything".to_string();
    }
    // What `~` meant before junctions were their own value, and still what the
    // relational codes constrain.
    if s == Seq::N_GAP_PAD {
        return "anything the read spans".to_string();
    }
    if s.is_empty() {
        return "nothing".to_string();
    }
    let mut parts: Vec<String> = Vec::new();
    let bases = s.0 & Seq::BASES;
    if bases == Seq::N.0 {
        parts.push("any base".to_string());
    } else if bases != 0 {
        let letters: Vec<String> = Seq(bases)
            .individual_bits()
            .iter()
            .map(|b| b.to_char().to_string())
            .collect();
        parts.push(join_or(&letters));
    }
    if s.0 & Seq::GAP.0 != 0 {
        parts.push("a gap".to_string());
    }
    if s.0 & Seq::PAD.0 != 0 {
        parts.push("a pad".to_string());
    }
    if s.0 & Seq::CLIP.0 != 0 {
        parts.push("a clipped base".to_string());
    }
    if s.0 & Seq::SKIP.0 != 0 {
        parts.push("a junction".to_string());
    }
    join_or(&parts)
}

/// Unconstrained: says nothing about the column. Both spellings count, since a
/// relational column carries `N_GAP_PAD` and constrains through the relation
/// rather than through its base set.
fn free(s: Seq) -> bool {
    s == Seq::ANY || s == Seq::N_GAP_PAD
}

fn constrained(c: &ColumnSpec) -> bool {
    c.rel.is_some() || !free(c.read) || !free(c.refr)
}

/// "columns 3-4" or "column 3".
fn cols(a: usize, b: usize) -> String {
    if a == b {
        format!("column {a}")
    } else {
        format!("columns {a}-{b}")
    }
}

/// Spell one side across a run: `CG` when every column is one unambiguous base,
/// otherwise a "then" chain.
fn spell(seqs: &[Seq]) -> String {
    let plain = seqs
        .iter()
        .all(|s| s.0 & Seq::FLAGS == 0 && (s.0 & Seq::BASES).count_ones() == 1);
    if plain {
        return seqs.iter().map(|s| s.to_char()).collect();
    }
    seqs.iter().map(|s| describe(*s)).collect::<Vec<_>>().join(" then ")
}

/// One side's contribution to a run, or `None` if it constrains nothing there.
fn side_phrase(name: &str, seqs: &[Seq], start: usize) -> Option<String> {
    let live: Vec<usize> = (0..seqs.len()).filter(|&i| !free(seqs[i])).collect();
    match live.len() {
        0 => None,
        // Naming the one column it touches beats spelling a sequence that is
        // mostly "anything" -- but a one-column run already said which column
        // in its own prefix, so do not say it twice.
        1 => {
            let i = live[0];
            let d = describe(seqs[i]);
            Some(if seqs.len() == 1 {
                format!("the {name} is {d}")
            } else {
                format!("the {name} is {d} at column {}", start + i)
            })
        }
        _ => Some(format!("the {name} is {}", spell(seqs))),
    }
}

fn rel_phrase(c: &ColumnSpec, start: usize, end: usize) -> String {
    let verb = match c.rel {
        Some(Rel::Eq) => "agrees with",
        _ => "differs from",
    };
    // The relational code always leaves its own side unconstrained, so at most
    // one side carries a base set here.
    let extra = if !free(c.read) {
        format!("the read is {} and ", describe(c.read))
    } else if !free(c.refr) {
        format!("the reference is {} and ", describe(c.refr))
    } else {
        String::new()
    };
    // The "both unambiguous" caveat is stated once at the end rather than on
    // every clause; repeated three times across a mismatch run it drowns out
    // the sentence.
    format!("at {} {extra}the read {verb} the reference", cols(start, end))
}

/// One operand in words.
pub fn pattern_prose(p: &PatternSpec) -> String {
    let n = p.columns.len();
    let mut clauses: Vec<String> = Vec::new();
    let mut i = 0usize;

    while i < n {
        if !constrained(&p.columns[i]) {
            i += 1;
            continue;
        }
        // A relational column couples the two sides, so it cannot join a run
        // that describes them independently. Identical adjacent ones do collapse
        // -- "at columns 0-2 the read differs" beats saying it three times.
        if p.columns[i].rel.is_some() {
            let start = i;
            while i < n && p.columns[i] == p.columns[start] {
                i += 1;
            }
            clauses.push(rel_phrase(&p.columns[start], start, i - 1));
            continue;
        }
        let start = i;
        while i < n && constrained(&p.columns[i]) && p.columns[i].rel.is_none() {
            i += 1;
        }
        let end = i - 1;
        let reads: Vec<Seq> = p.columns[start..=end].iter().map(|c| c.read).collect();
        let refrs: Vec<Seq> = p.columns[start..=end].iter().map(|c| c.refr).collect();

        // Read first, then reference: the same order as the `read@reference`
        // syntax, so the sentence and the pattern can be read side by side.
        let mut parts: Vec<String> = Vec::new();
        if let Some(s) = side_phrase("read", &reads, start) {
            parts.push(s);
        }
        if let Some(s) = side_phrase("reference", &refrs, start) {
            parts.push(s);
        }
        clauses.push(format!("at {} {}", cols(start, end), parts.join(", and ")));
    }

    let width = if n == 1 { "1 column".to_string() } else { format!("{n} columns") };
    if clauses.is_empty() {
        return format!("{width}, unconstrained");
    }
    format!("{width}: {}", clauses.join("; "))
}

/// The expression in words, with operands referred to by name where they have
/// one and by number where they do not.
fn logic_prose(e: &Expr, spec: &QuerySpec, parent_is_or: bool) -> String {
    match e {
        Expr::Const(b) => if *b { "always" } else { "never" }.to_string(),
        Expr::Pattern(i) => format!("{} matches", spec.operand_label(*i)),
        Expr::Not(x) => match &**x {
            Expr::Pattern(i) => format!("{} does not match", spec.operand_label(*i)),
            other => format!("it is not the case that ({})", logic_prose(other, spec, false)),
        },
        Expr::And(xs) => {
            let inner = xs
                .iter()
                .map(|x| logic_prose(x, spec, false))
                .collect::<Vec<_>>()
                .join(" and ");
            if parent_is_or {
                format!("({inner})")
            } else {
                inner
            }
        }
        Expr::Or(xs) => {
            let inner = xs
                .iter()
                .map(|x| logic_prose(x, spec, true))
                .collect::<Vec<_>>()
                .join(" or ");
            format!("({inner})")
        }
        Expr::Count { lo, hi, patterns } => format!(
            "between {lo} and {hi} of {} match",
            patterns
                .iter()
                .map(|p| spec.operand_label(*p))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// The whole query in words.
pub fn prose(spec: &QuerySpec) -> String {
    let mut out = String::new();

    out.push_str(&format!(
        "  Fires at a position where {}.\n",
        logic_prose(&spec.expr, spec, false)
    ));

    for (i, op) in spec.operands.iter().enumerate() {
        out.push_str(&format!("    {} {}\n", spec.operand_label(i), pattern_prose(op)));
    }

    if spec.operands.iter().any(|o| o.columns.iter().any(|c| c.rel.is_some())) {
        out.push_str(
            "    Note: 'agrees with' and 'differs from' both require an unambiguous\n\
             \x20         A/C/G/T on each side, so neither holds at a gap, a pad, or an\n\
             \x20         ambiguity code.\n",
        );
    }

    out.push_str(&format!(
        "  Reports the coordinates of column {}",
        spec.anchor
    ));
    match spec.captures.len() {
        0 => out.push_str(".\n"),
        1 => out.push_str(&format!(
            ", and records what was observed at column {}.\n",
            spec.captures[0]
        )),
        _ => out.push_str(&format!(
            ", and records what was observed at columns {}.\n",
            spec.captures
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )),
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(toml: &str) -> String {
        prose(&crate::test_support::query(toml))
    }

    /// The operand descriptions only, excluding the footnote -- which mentions
    /// the same verbs and would otherwise inflate any count.
    fn operands(s: &str) -> String {
        p(s).lines()
            .filter(|l| l.trim_start().starts_with("p1") || l.trim_start().starts_with("p2"))
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn base_sets_are_spelled_out() {
        assert_eq!(describe(Seq::C), "C");
        assert_eq!(describe(Seq::Y), "C or T");
        assert_eq!(describe(Seq::H), "A, C or T");
        assert_eq!(describe(Seq::N), "any base");
        assert_eq!(describe(Seq::ANY), "anything");
        assert_eq!(describe(Seq::GAP), "a gap");
        assert_eq!(describe(Seq::C | Seq::GAP), "C or a gap");
    }

    #[test]
    fn the_methyl_hic_query_reads_as_english() {
        let t = p("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n[pat.p2]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.mCpG_hic]\nmark = \"...+....\"\nwhere = \"p1 and not p2\"\n");
        println!("{t}");
        assert!(t.contains("p1 matches and p2 does not match"), "{t}");
        // The free columns either side are not mentioned.
        // Read side first, matching `read@reference`.
        assert!(
            t.contains("at columns 3-4 the read is C or T at column 3, and the reference is CG"),
            "{t}"
        );
        assert!(t.contains("at columns 0-7 the reference is GATCGATC"), "{t}");
        assert!(t.contains("coordinates of column 3"), "{t}");
        assert!(t.contains("observed at column 3"), "{t}");
    }

    #[test]
    fn a_single_constrained_column_is_named_not_spelled() {
        // The read is free at column 4, so spelling it would give
        // "C-or-T then anything", which is noise.
        let t = p("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n\n[query.x]\nmark = \"+.......\"\nwhere = \"p1\"\n");
        assert!(t.contains("the read is C or T at column 3"), "{t}");
        assert!(!t.contains("then anything"), "{t}");
    }

    #[test]
    fn relational_columns_get_their_own_clause() {
        let t = p("[pat.p1]\nread = \"~\"\nrefr = \"=\"\n\n[query.m]\nmark = \"+\"\nwhere = \"p1\"\n");
        assert!(t.contains("agrees with the reference"), "{t}");
        // The caveat is a footnote, stated once, not on every clause.
        assert_eq!(t.matches("unambiguous").count(), 1, "{t}");
        assert!(!operands("[pat.p1]\nread = \"~\"\nrefr = \"=\"\n\n[query.m]\nmark = \"+\"\nwhere = \"p1\"\n").contains("unambiguous"), "{t}");
        let t = p("[pat.p1]\nread = \"~\"\nrefr = \"/\"\n\n[query.x]\nmark = \"+\"\nwhere = \"p1\"\n");
        assert!(t.contains("differs from the reference"), "{t}");
        let t = p("[pat.p1]\nread = \"H\"\nrefr = \"/\"\n\n[query.y]\nmark = \"+\"\nwhere = \"p1\"\n");
        assert!(t.contains("the read is A, C or T and the read differs"), "{t}");
    }

    #[test]
    fn precedence_survives_into_prose() {
        let t = p("[pat.p1]\nread = \"A\"\nrefr = \"A\"\n[pat.p2]\nread = \"B\"\nrefr = \"B\"\n[pat.p3]\nread = \"C\"\nrefr = \"C\"\n\n[query.x]\nmark = \"+\"\nwhere = \"p1 or p2 and p3\"\n");
        assert!(t.contains("(p1 matches or (p2 matches and p3 matches))"), "{t}");
        let t = p("[pat.p1]\nread = \"A\"\nrefr = \"A\"\n[pat.p2]\nread = \"B\"\nrefr = \"B\"\n[pat.p3]\nread = \"C\"\nrefr = \"C\"\n\n[query.y]\nmark = \"+\"\nwhere = \"(p1 or p2) and p3\"\n");
        assert!(t.contains("(p1 matches or p2 matches) and p3 matches"), "{t}");
    }

    #[test]
    fn negated_groups_are_spelled_out() {
        let t = p("[pat.p1]\nread = \"A\"\nrefr = \"A\"\n[pat.p2]\nread = \"C\"\nrefr = \"C\"\n[pat.p3]\nread = \"G\"\nrefr = \"G\"\n\n[query.x]\nmark = \"+\"\nwhere = \"p1 and not (p2 and p3)\"\n");
        assert!(t.contains("it is not the case that"), "{t}");
    }

    #[test]
    fn identical_relational_columns_collapse() {
        let t = operands("[pat.p1]\nread = \"~~~\"\nrefr = \"///\"\n\n[query.mm3]\nmark = \"+..\"\nwhere = \"p1\"\n");
        assert!(t.contains("at columns 0-2 the read differs from the reference"), "{t}");
        assert_eq!(t.matches("differs from").count(), 1, "{t}");
    }

    #[test]
    fn differing_relational_columns_stay_separate() {
        let t = operands("[pat.p1]\nread = \"~~~\"\nrefr = \"=/=\"\n\n[query.x]\nmark = \"+..\"\nwhere = \"p1\"\n");
        assert_eq!(t.matches("agrees with").count(), 2, "{t}");
        assert_eq!(t.matches("differs from").count(), 1, "{t}");
    }

    #[test]
    fn captures_are_listed() {
        let t = p("[pat.p1]\nread = \"C~C\"\nrefr = \"CGC\"\n\n[query.x]\nmark = \"+.^\"\nwhere = \"p1\"\n");
        assert!(t.contains("observed at columns 0, 2"), "{t}");
        // The anchor is recorded even when nothing was marked, so this reads
        // as one capture rather than none.
        let t = p("[pat.p1]\nread = \"C~\"\nrefr = \"CG\"\n\n[query.y]\nmark = \"+.\"\nwhere = \"p1\"\n");
        assert!(t.contains("coordinates of column 0, and records what was observed at column 0."), "{t}");
    }

    #[test]
    fn gaps_and_pads_are_named() {
        let t = p("[pat.p1]\nread = \"N\"\nrefr = \".\"\n\n[query.ins]\nmark = \"+\"\nwhere = \"p1\"\n");
        assert!(t.contains("a gap"), "{t}");
        let t = p("[pat.p1]\nread = \"_N\"\nrefr = \"NN\"\n\n[query.edge]\nmark = \"+.\"\nwhere = \"p1\"\n");
        assert!(t.contains("a pad"), "{t}");
    }
}
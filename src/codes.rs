//! Print the code table for `--list-codes`.
//!
//! Every character a pattern can contain, in one place the user can actually
//! reach. An error message that names the free alias characters only helps
//! after you have already guessed wrong; this is the answer to "what can I
//! type", which is the question people have first.
//!
//! Built from [`Seq`] and [`Aliases`] rather than written out, so it cannot
//! drift from what the parser accepts.

use crate::dsl::Aliases;
use crate::seq::Seq;

/// Codes with a one-character form, which is what a matrix grid can hold.
const IUPAC: &[(char, &str)] = &[
    ('A', "A"),
    ('C', "C"),
    ('G', "G"),
    ('T', "T"),
    ('R', "A or G"),
    ('Y', "C or T"),
    ('S', "C or G"),
    ('W', "A or T"),
    ('K', "G or T"),
    ('M', "A or C"),
    ('B', "C, G or T"),
    ('D', "A, G or T"),
    ('H', "A, C or T"),
    ('V', "A, C or G"),
    ('N', "any base"),
];

fn row(code: &str, meaning: &str) -> String {
    format!("  {code:<8}{meaning}\n")
}

/// The full table. `aliases` is whatever the current invocation declared.
pub fn table(aliases: &Aliases) -> String {
    let mut out = String::new();

    out.push_str("BASE CODES (uppercase)\n");
    for (c, m) in IUPAC {
        out.push_str(&row(&c.to_string(), m));
    }

    out.push_str("\nEXTENSIONS\n");
    out.push_str(&row(". - Z", "gap: a deletion, or an insertion on the other side"));
    out.push_str(&row("_ X", "pad: past the end of the read, or off the contig"));
    out.push_str(&row(": L", "clip: a soft-clipped base, in the flank beside the aligned part"));
    out.push_str(&row(", J", "junction: reference skipped by a CIGAR N, i.e. an intron"));
    out.push_str(&row("~", "anything at all — a base, a gap, a pad, a clip, or a junction"));
    out.push_str(&row("{ACG}", "union of any codes; {C.} is C-or-a-gap, {NJ} base-or-junction"));
    out.push_str(
        "  Write a junction ',' on its own and 'J' inside a group: a comma in a\n  \
         group is rejected, because '{A,C,G}' is the natural guess at a union and\n  \
         would quietly mean A, C, G or a junction.\n",
    );

    out.push_str("\nJUNCTIONS\n");
    out.push_str(
        "  An intron is not emitted base by base -- it would be tens of thousands\n  \
         of columns the read never covered. It emits k reference bases at each\n  \
         end against ',' on the read side, where k is the widest query, and one\n  \
         ',@,' marker for everything between. A pattern that does not name a\n  \
         junction cannot span one, which is the point: two exonic bases either\n  \
         side of an intron are not adjacent in the genome.\n\n  \
         An intron short enough to fit in the two windows is shown whole and has\n  \
         no marker, so ',@,' means \"a junction with something elided\".\n",
    );

    out.push_str("\nRELATIONAL (one side only; the other gives the base set)\n");
    out.push_str(&row("=", "both sides the same unambiguous A/C/G/T"));
    out.push_str(&row("/", "both sides unambiguous A/C/G/T and different"));
    out.push_str(
        "  Neither fires on a gap, a pad, a junction, or an ambiguity code, so they\n  \
         are not complements: '~@=' and '~@/' together mean \"both sides are real\".\n",
    );

    out.push_str("\nMARK ROW (the `mark` key of a query, as wide as its patterns)\n");
    out.push_str(&row(".", "nothing"));
    out.push_str(&row("^", "also record what was observed at this column"));
    out.push_str(&row("+", "report this column's coordinates (default: column 0)"));
    out.push_str(
        "  Markers go only in `mark`, never in a `read` or `refr` row, and rows take\n  \
         no repeat counts: write every column out. The anchor is always recorded,\n  \
         so '+' implies '^'.\n",
    );

    out.push_str("\nALIASES\n");
    out.push_str(&format!(
        "  Uppercase letters are base codes. Lowercase letters are yours.\n\n  \
         A row holds one character per column, but only {} of the {} base\n  \
         sets have a one-character form. An alias names one of the rest, in the\n  \
         [alias] table of a --query-file:\n\n    \
         [alias]\n    \
         j = \"{{C._}}\"\n\n  \
         Any of a-z will do. Prefer f i j l o p q, the letters that are not a\n  \
         base code in either case -- a lowercase 'a' reads as a soft-masked A.\n",
        single_char_count(),
        Seq::COUNT - 1
    ));

    if aliases.is_empty() {
        out.push_str("  None declared.\n");
    } else {
        out.push_str("\n  Declared here:\n");
        for (c, s) in aliases.iter() {
            out.push_str(&format!("    {c} = {}\n", s.name()));
        }
    }

    out
}

/// How many `Seq` values render as a single character.
///
/// Counted rather than stated: it is the reason aliases exist, and a wrong
/// number in the help text would be a small lie about the tool's own limits.
fn single_char_count() -> usize {
    Seq::all()
        .filter(|s| !s.is_empty() && s.name().chars().count() == 1)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_iupac_code_round_trips_through_the_parser() {
        // The table is a promise about what the parser accepts.
        for (c, _) in IUPAC {
            assert!(
                Seq::from_name(&c.to_string()).is_some(),
                "table lists '{c}' but the parser rejects it"
            );
        }
        for c in [".", "-", "Z", "_", "X", ":", "L", ",", "J", "~"] {
            assert!(Seq::from_name(c).is_some(), "table lists '{c}'");
        }
    }

    #[test]
    fn the_alias_count_matches_reality() {
        let n = single_char_count();
        // Every code the table lists, plus the ones with several spellings.
        assert!(n >= IUPAC.len(), "{n}");
        let t = table(&Aliases::new());
        assert!(t.contains(&format!("only {n} of")), "{t}");
    }

    #[test]
    fn the_rule_is_stated_before_the_details() {
        let t = table(&Aliases::new());
        assert!(
            t.contains("Uppercase letters are base codes. Lowercase letters are yours."),
            "{t}"
        );
        assert!(t.contains("Any of a-z"), "{t}");
        assert!(t.contains("None declared."), "{t}");
    }

    #[test]
    fn declared_aliases_are_shown() {
        let mut a = Aliases::new();
        a.insert('j', Seq::C | Seq::GAP).unwrap();
        let t = table(&a);
        assert!(t.contains("j = {C.}"), "{t}");
    }

    #[test]
    fn the_relational_caveat_is_stated() {
        // The single most misread part of the language.
        let t = table(&Aliases::new());
        assert!(t.contains("not"), "{t}");
        assert!(t.contains("complements"), "{t}");
    }
}
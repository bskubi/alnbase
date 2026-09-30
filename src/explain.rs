//! Show what a query says, column by column. This is what `--explain` prints.
//!
//! The DSL is dense. In the worked example of
//! [`crate::query_toml::EXAMPLE`], `mCpG` is twelve columns wide, and the whole
//! point of its `junc_end` operand is that the junction's final C lands on the
//! same column as the `cpg` operand's cytosine — which is invisible unless the
//! operands are printed stacked and aligned.
//!
//! This is also the on-ramp to the language. People learn a DSL by reading a
//! working example and modifying it, not from grammar documentation, so an
//! example plus `--explain` teaches more than a syntax section does.
//!
//! Pure rendering: it takes a parsed query and, optionally, the shape the
//! predicate compiled to. No BAM, no automaton.

use crate::dsl::{Aliases, QuerySpec};
use crate::predicate::Expr;
use crate::prose;

/// How each operand enters the expression.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Polarity {
    /// Only ever positive: a match makes the query more likely to fire.
    Positive,
    /// Only ever under `not`: a match suppresses the query.
    Negative,
    /// Both, which is legal but almost always a mistake worth seeing.
    Mixed,
}

impl Polarity {
    fn label(self) -> &'static str {
        match self {
            Polarity::Positive => "   ",
            Polarity::Negative => "not",
            Polarity::Mixed => " +-",
        }
    }
}

/// Walk the expression recording, for each operand, whether it is reached under
/// an even or odd number of `not`s.
fn polarities(e: &Expr, neg: bool, seen: &mut [(bool, bool)]) {
    match e {
        Expr::Const(_) => {}
        Expr::Pattern(p) => {
            let slot = &mut seen[*p];
            if neg {
                slot.1 = true;
            } else {
                slot.0 = true;
            }
        }
        Expr::Not(x) => polarities(x, !neg, seen),
        Expr::And(xs) | Expr::Or(xs) => {
            for x in xs {
                polarities(x, neg, seen);
            }
        }
        Expr::Count { patterns, .. } => {
            for p in patterns {
                let slot = &mut seen[*p];
                if neg {
                    slot.1 = true;
                } else {
                    slot.0 = true;
                }
            }
        }
    }
}

fn operand_polarities(spec: &QuerySpec) -> Vec<Polarity> {
    let mut seen = vec![(false, false); spec.operands.len()];
    polarities(&spec.expr, false, &mut seen);
    seen.iter()
        .map(|&(pos, negv)| match (pos, negv) {
            (true, true) => Polarity::Mixed,
            (_, true) => Polarity::Negative,
            _ => Polarity::Positive,
        })
        .collect()
}

/// A node's label and the children to draw beneath it.
///
/// A `not` applied directly to a pattern goes into the label. It does not
/// become a level of its own. The tree therefore shows `not [2]`, and not a
/// `not` node with one child below it. A negated operand is the common case,
/// so this keeps the tree short.
fn node<'e>(e: &'e Expr, spec: &QuerySpec) -> (String, Vec<&'e Expr>) {
    let ops = &spec.operands;
    let leaf = |i: usize, neg: bool| {
        let text = ops.get(i).map(|o| o.text.as_str()).unwrap_or("?");
        (
            format!(
                "{}{}  {}",
                if neg { "not " } else { "" },
                spec.operand_label(i),
                text
            ),
            Vec::new(),
        )
    };
    match e {
        Expr::Const(b) => (b.to_string(), Vec::new()),
        Expr::Pattern(i) => leaf(*i, false),
        Expr::Not(x) => match &**x {
            Expr::Pattern(i) => leaf(*i, true),
            other => ("not".to_string(), vec![other]),
        },
        Expr::And(xs) => ("and".to_string(), xs.iter().collect()),
        Expr::Or(xs) => ("or".to_string(), xs.iter().collect()),
        Expr::Count { lo, hi, patterns } => (
            format!(
                "{}..{} of {}",
                lo,
                hi,
                patterns
                    .iter()
                    .map(|p| spec.operand_label(*p))
                    .collect::<Vec<_>>()
                    .join(" ")
            ),
            Vec::new(),
        ),
    }
}

fn draw(e: &Expr, spec: &QuerySpec, prefix: &str, connector: &str, out: &mut String) {
    let (label, kids) = node(e, spec);
    out.push_str(&format!("  {prefix}{connector}{label}\n"));
    // A child of the last branch has nothing below it to connect to, so its
    // continuation is blank rather than a trailing spine.
    let child_prefix = match connector {
        "" => prefix.to_string(),
        c if c.starts_with('`') => format!("{prefix}   "),
        _ => format!("{prefix}|  "),
    };
    for (i, k) in kids.iter().enumerate() {
        let last = i + 1 == kids.len();
        draw(k, spec, &child_prefix, if last { "`- " } else { "|- " }, out);
    }
}

/// The expression as a tree, so precedence is visible rather than inferred.
///
/// A user can reasonably write `A@A or B@B and C@C`. A reader can also easily
/// misread it. The tree shows which way the expression bound.
fn structure(spec: &QuerySpec) -> String {
    let mut out = String::new();
    draw(&spec.expr, spec, "", "", &mut out);
    out
}

/// Which kind of line a rendered row is.
enum Row {
    /// The `0 1 2 ...` header.
    Index,
    /// `^` under each captured column.
    Capture,
    /// `+` under the anchor column.
    Anchor,
    /// `(operand index, 0 for read | 1 for reference)`.
    Side(usize, usize),
}


/// Render one query for `--explain`, showing declared aliases by name.
///
/// A grid showing `J` where the user wrote `J` is the point of having aliases;
/// expanding them back to `{C._}` would both widen the column and hide the
/// abbreviation they chose.
pub fn explain_with(
    spec: &QuerySpec,
    shape: Option<(&str, usize)>,
    aliases: &Aliases,
) -> String {
    let pol = operand_polarities(spec);
    let n = spec.span;

    let cells: Vec<Vec<(String, String)>> = spec
        .operands
        .iter()
        .map(|op| {
            op.columns
                .iter()
                .map(|c| c.to_text_with(aliases))
                .collect()
        })
        .collect();

    // Every cell in a column shares a width, so the grid lines up regardless of
    // how wide a `{ACG}` group or a two-digit index renders.
    let mut widths = vec![1usize; n];
    for (i, w) in widths.iter_mut().enumerate() {
        *w = (*w).max(i.to_string().len());
        for row in &cells {
            *w = (*w).max(row[i].0.len()).max(row[i].1.len());
        }
    }

    // Build every row's label first, then size the gutter from the labels
    // themselves. Computing that width by hand is how the columns end up one
    // space out of true.
    // Anchor and capture get a row each rather than sharing a cell. Stacking
    // them would make any column carrying both two characters wide, which pads
    // every base in that column and breaks the eye's read down the grid --
    // exactly on the column the reader most wants to follow.
    let mut rows: Vec<(String, Row)> = vec![("col".to_string(), Row::Index)];
    if !spec.captures.is_empty() {
        rows.push(("capture".to_string(), Row::Capture));
    }
    rows.push(("anchor".to_string(), Row::Anchor));
    for oi in 0..cells.len() {
        rows.push((
            format!("{} {}  read", pol[oi].label(), spec.operand_label(oi)),
            Row::Side(oi, 0),
        ));
        rows.push(("refr".to_string(), Row::Side(oi, 1)));
    }
    let gutter = rows.iter().map(|(p, _)| p.len()).max().unwrap_or(0);

    let mut out = String::new();
    out.push_str(if spec.name.is_empty() { "(unnamed)" } else { &spec.name });
    out.push('\n');
    out.push_str("  ");
    out.push_str(&spec.text);
    out.push_str("\n\n");

    // A lone pattern has no structure worth drawing; the echoed text above says
    // everything the tree would.
    if !matches!(spec.expr, Expr::Pattern(_)) {
        out.push_str(&structure(spec));
        out.push('\n');
    }

    // Words before the grid: the reader who needs the grid has already been
    // told what to look for, and the reader who does not can stop here.
    out.push_str(&prose::prose(spec));
    out.push('\n');

    for (label, row) in &rows {
        out.push_str(&format!("{label:>gutter$}  "));
        for i in 0..n {
            let w = widths[i];
            let cell = match row {
                Row::Index => i.to_string(),
                Row::Capture => {
                    if spec.captures.contains(&i) { "^".into() } else { String::new() }
                }
                Row::Anchor => {
                    if spec.anchor == i { "+".into() } else { String::new() }
                }
                Row::Side(oi, 0) => cells[*oi][i].0.clone(),
                Row::Side(oi, _) => cells[*oi][i].1.clone(),
            };
            out.push_str(&format!("{cell:>w$} "));
        }
        out.push('\n');
    }

    let captures = if spec.captures.is_empty() {
        "none".to_string()
    } else {
        spec.captures.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", ")
    };
    out.push_str(&format!(
        "\n  anchor   column {}\n  captures {}\n  span     {} columns\n",
        spec.anchor, captures, n
    ));
    if let Some((form, groups)) = shape {
        out.push_str(&format!("  compiled {form}, {groups} groups\n"));
    }
    if names_a_pad(spec) {
        out.push_str(READ2_PAD_NOTE);
    }
    out
}

/// Whether some operand asks for a flank column (a pad or a clip) on the read
/// side, rather than merely admitting one the way `~` admits everything.
fn names_a_pad(spec: &QuerySpec) -> bool {
    use crate::seq::Seq;
    spec.operands.iter().flat_map(|op| op.columns.iter()).any(|c| {
        c.read.0 & (Seq::PAD.0 | Seq::CLIP.0) != 0 && c.read != Seq::ANY && c.read != Seq::N_GAP_PAD
    })
}

/// Printed under a query that places a pad. This is the one case where the
/// direction of a pattern and the direction of an offset differ.
const READ2_PAD_NOTE: &str = "\n  note     patterns run along the conversion strand, but off_5p/off_3p count from\n           \
     the ends as sequenced. For read 1 a pad left of the read is past its 5' end;\n           \
     for read 2 (CTOT, CTOB) it is past its 3' end, and a pad on the right is past\n           \
     its 5' end.\n";

#[cfg(test)]
mod tests {
    use super::*;

    fn ex(toml: &str) -> String {
        explain_with(&crate::test_support::query(toml), Some(("cnf", 2)), &Aliases::new())
    }

    /// A query that places a pad is told how its direction relates to read 2's
    /// offsets; one that merely admits pads through `~` is not.
    #[test]
    fn a_pad_in_a_pattern_gets_the_read_2_note() {
        assert!(ex("[query.edge]\nread = \"_Y\"\nrefr = \"~C\"\n").contains("for read 2 (CTOT, CTOB)"));
        assert!(ex("[alias]\ng = \"{._}\"\n[query.gp]\nread = \"gY\"\nrefr = \"~C\"\n").contains("note"));
        assert!(!ex("[query.cg]\nread = \"Y~\"\nrefr = \"CG\"\n").contains("note"));
    }

    /// Every grid line must have exactly the same length. If two lines differ,
    /// the columns do not line up, and the grid shows nothing useful.
    fn assert_aligned(text: &str) {
        let grid: Vec<&str> = text
            .lines()
            .skip_while(|l| !l.trim_start().starts_with("col "))
            .take_while(|l| !l.is_empty())
            .collect();
        assert!(grid.len() >= 4, "too few grid lines in:\n{text}");
        for l in &grid {
            assert_eq!(l.len(), grid[0].len(), "misaligned:\n{text}");
        }
    }

    #[test]
    fn the_methyl_hic_query_lines_up() {
        let t = ex("[pattern.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n[pattern.p2]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.mCpG_hic]\nmark = \"...+....\"\nwhere = \"p1 and not p2\"\n");
        assert_aligned(&t);
        assert!(t.contains("mCpG_hic"), "{t}");
        assert!(t.contains("not"), "{t}");
        assert!(t.contains("anchor   column 3"), "{t}");
        assert!(t.contains("captures 3"), "{t}");
        assert!(t.contains("span     8 columns"), "{t}");
        assert!(t.contains("compiled cnf, 2 groups"), "{t}");
        println!("{t}");
    }

    #[test]
    fn both_markers_on_one_column_do_not_widen_it() {
        // The anchored-and-captured column must stay the same width as its
        // neighbours; that is the whole reason the markers get a row each.
        let t = ex("[pattern.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n[pattern.p2]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.mCpG_hic]\nmark = \"...+....\"\nwhere = \"p1 and not p2\"\n");
        assert_aligned(&t);
        let read = t
            .lines()
            .skip_while(|l| !l.trim_start().starts_with("col "))
            .find(|l| l.contains("  read"))
            .unwrap();
        let cells: Vec<&str> = read.split_whitespace().skip(2).collect();
        assert_eq!(cells, ["~", "~", "~", "Y", "~", "~", "~", "~"], "{t}");
        // One row each, not one cell holding both. Counted over the grid only:
        // the echoed expression above it contains the marker characters too.
        let grid: Vec<&str> = t
            .lines()
            .skip_while(|l| !l.trim_start().starts_with("col "))
            .take_while(|l| !l.is_empty())
            .collect();
        assert_eq!(grid.iter().filter(|l| l.contains('^')).count(), 1, "{t}");
        assert_eq!(grid.iter().filter(|l| l.contains('+')).count(), 1, "{t}");
        assert!(grid[1].trim_start().starts_with("capture"), "{t}");
        assert!(grid[2].trim_start().starts_with("anchor"), "{t}");
    }

    #[test]
    fn the_capture_row_shows_the_implicit_anchor_capture() {
        // The anchor is always recorded, so the capture row is always present
        // and always marks at least the anchor column.
        let t = ex("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n[pattern.p2]\nread = \"C\"\nrefr = \"C\"\n\n[query.x]\nmark = \"+\"\nwhere = \"p1 and p2\"\n");
        let grid: Vec<&str> = t
            .lines()
            .skip_while(|l| !l.trim_start().starts_with("col "))
            .take_while(|l| !l.is_empty())
            .collect();
        assert!(grid[1].trim_start().starts_with("capture"), "{t}");
        assert!(grid[2].trim_start().starts_with("anchor"), "{t}");
        assert_aligned(&t);
    }

    #[test]
    fn wide_groups_do_not_break_alignment() {
        assert_aligned(&ex("[alias]\nf = \"{C._}\"\n\n[pattern.p1]\nread = \"V~f\"\nrefr = \"CGA\"\n[pattern.p2]\nread = \"~~~\"\nrefr = \"KCA\"\n\n[query.x]\nmark = \"+..\"\nwhere = \"p1 and p2\"\n"));
    }

    #[test]
    fn double_digit_columns_do_not_break_alignment() {
        assert_aligned(&ex("[pattern.p1]\nread = \"~~~~~~~~~~~~~~~C\"\nrefr = \"~~~~~~~~~~~~~~~C\"\n\n[query.x]\nmark = \"+...............\"\nwhere = \"p1\"\n"));
    }

    #[test]
    fn relational_columns_render_on_the_written_side() {
        let t = ex("[pattern.p1]\nread = \"~\"\nrefr = \"=\"\n\n[query.m]\nmark = \"+\"\nwhere = \"p1\"\n");
        assert!(t.contains("="), "{t}");
        let t = ex("[pattern.p1]\nread = \"=\"\nrefr = \"C\"\n\n[query.m]\nmark = \"+\"\nwhere = \"p1\"\n");
        assert!(t.contains("="), "{t}");
    }

    #[test]
    fn polarity_is_reported_per_operand() {
        let t = ex("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n[pattern.p2]\nread = \"C\"\nrefr = \"C\"\n\n[query.x]\nmark = \"+\"\nwhere = \"p1 and not p2\"\n");
        let not_lines: Vec<&str> = t.lines().filter(|l| l.starts_with("not")).collect();
        assert_eq!(not_lines.len(), 1, "{t}");
        // The negated operand is the second one.
        assert!(not_lines[0].contains("p2"), "{t}");
    }

    #[test]
    fn an_operand_used_both_ways_is_flagged() {
        let t = ex("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n[pattern.p2]\nread = \"C\"\nrefr = \"C\"\n\n[query.x]\nmark = \"+\"\nwhere = \"p1 and (not p1 or p2)\"\n");
        assert!(t.contains("+-"), "expected a mixed-polarity marker in:\n{t}");
    }

    /// The grid region only; the tree and the echoed expression sit above it.
    fn tree_of(t: &str) -> Vec<String> {
        t.lines()
            .skip(3)
            .take_while(|l| !l.trim().is_empty())
            .map(|l| l.trim_end().to_string())
            .collect()
    }

    #[test]
    fn precedence_is_visible() {
        // `or` binds loosest, so it must be the root and `and` the subtree.
        // This is the misreading the tree exists to prevent.
        let t = tree_of(&ex("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n[pattern.p2]\nread = \"B\"\nrefr = \"B\"\n[pattern.p3]\nread = \"C\"\nrefr = \"C\"\n\n[query.p]\nmark = \"+\"\nwhere = \"p1 or p2 and p3\"\n"));
        assert_eq!(
            t,
            [
                "  or",
                "  |- p1  A@A",
                "  `- and",
                "     |- p2  B@B",
                "     `- p3  C@C",
            ]
        );
    }

    #[test]
    fn parentheses_change_the_root() {
        let t = tree_of(&ex("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n[pattern.p2]\nread = \"B\"\nrefr = \"B\"\n[pattern.p3]\nread = \"C\"\nrefr = \"C\"\n\n[query.p]\nmark = \"+\"\nwhere = \"(p1 or p2) and p3\"\n"));
        assert_eq!(t[0], "  and");
        assert_eq!(t[1], "  |- or");
    }

    #[test]
    fn negated_patterns_fold_into_one_line() {
        // `not` on a bare pattern is the common case and should not cost a
        // level of nesting.
        let t = tree_of(&ex("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n[pattern.p2]\nread = \"C\"\nrefr = \"C\"\n\n[query.p]\nmark = \"+\"\nwhere = \"p1 and not p2\"\n"));
        assert_eq!(t, ["  and", "  |- p1  A@A", "  `- not p2  C@C"]);
    }

    #[test]
    fn negated_groups_keep_their_own_node() {
        let t = tree_of(&ex("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n[pattern.p2]\nread = \"C\"\nrefr = \"C\"\n[pattern.p3]\nread = \"G\"\nrefr = \"G\"\n\n[query.p]\nmark = \"+\"\nwhere = \"p1 and not (p2 and p3)\"\n"));
        assert!(t.iter().any(|l| l.trim_end().ends_with("not")), "{t:?}");
        assert!(t.iter().any(|l| l.contains("p3  G@G")), "{t:?}");
    }

    #[test]
    fn a_lone_pattern_has_no_tree() {
        let t = ex("[pattern.p1]\nread = \"C~\"\nrefr = \"CG\"\n\n[query.p]\nmark = \"+.\"\nwhere = \"p1\"\n");
        assert!(!t.contains("|-"), "{t}");
        assert!(!t.contains("`-"), "{t}");
    }

    #[test]
    fn tree_operand_numbers_match_the_grid() {
        let t = ex("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n[pattern.p2]\nread = \"C\"\nrefr = \"C\"\n\n[query.p]\nmark = \"+\"\nwhere = \"p1 and not p2\"\n");
        // [2] is the negated one in both renderings.
        assert!(t.contains("`- not p2  C@C"), "{t}");
        assert!(t.lines().any(|l| l.starts_with("not p2  read")), "{t}");
    }

    #[test]
    fn the_grid_uses_declared_aliases() {
        let mut a = Aliases::new();
        a.insert('j', crate::seq::Seq::C | crate::seq::Seq::GAP).unwrap();
        let spec = crate::test_support::query("[alias]\nj = \"{C.}\"\n\n[pattern.p1]\nread = \"j~\"\nrefr = \"CG\"\n\n[query.x]\nmark = \"+.\"\nwhere = \"p1\"\n");

        let with = explain_with(&spec, None, &a);
        assert!(with.contains(" j "), "{with}");
        assert!(!with.contains("{C.}"), "{with}");
        assert_aligned(&with);

        // Without the table the set is spelled out, which is wider but still
        // has to line up.
        let without = explain_with(&spec, None, &Aliases::new());
        assert!(without.contains("{C.}"), "{without}");
        assert_aligned(&without);
    }

    #[test]
    fn named_operands_replace_the_numbers_everywhere() {
        let spec = crate::test_support::query(
            "[pattern.cpg]\nread = \"Y~\"\nrefr = \"CG\"\n[pattern.junc]\nread = \"~~\"\nrefr = \"GA\"\n\
             [query.mCpG]\nwhere = \"cpg and not junc\"\n",
        );
        let t = explain_with(&spec, Some(("cnf", 2)), &Aliases::new());
        // Tree, grid gutter and prose all use the name.
        assert!(t.contains("|- cpg  "), "{t}");
        assert!(t.contains("`- not junc  "), "{t}");
        assert!(t.contains("cpg matches and junc does not match"), "{t}");
        assert!(t.lines().any(|l| l.starts_with("not junc  read")), "{t}");
        assert!(!t.contains("p1  "), "{t}");
        assert_aligned(&t);
    }



    #[test]
    fn works_without_a_compiled_shape() {
        let spec = crate::test_support::query("[pattern.p1]\nread = \"A\"\nrefr = \"A\"\n\n[query.unnamed]\nmark = \"+\"\nwhere = \"p1\"\n");
        let t = explain_with(&spec, None, &Aliases::new());
        assert!(!t.contains("compiled"), "{t}");
    }
}
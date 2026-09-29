//! Run queries over a single record and show each step of the match.
//!
//! Two questions bring people here, and the second is the hard one:
//!
//! - *Why did this fire?* Answered by evaluating each operand at the position
//!   and annotating the expression tree with what was true.
//! - *Why did it not fire?* Answered by the automaton's own state. Each pattern's
//!   active traces say how many of its columns currently match, so a pattern
//!   that never completes can be reported with the furthest it ever got and the
//!   column that broke it. That is the thing you cannot work out by staring at
//!   the query, and it is why this reads `State` rather than just `Hits`.
//!
//! Input is either a synthetic column pair typed on the command line, or real
//! columns from a BAM. Both arrive here as `&[Column]`, so the two paths cannot
//! disagree about what a trace means.

use std::fmt::Write as _;

use crate::column::Column;
use crate::byg_search::Bitmap;
use crate::dsl::QuerySpec;
use crate::predicate::{Expr, Predicate};
use crate::query::QuerySet;
use crate::seq::Seq;

/// Build columns from a `READ@REFR` pair, for tracing without a BAM.
///
/// Coordinates are the column index: synthetic input has no alignment, so
/// pretending otherwise would put made-up genomic positions in the output.
/// Quality is absent for the same reason.
pub fn columns_from_pair(s: &str) -> Result<Vec<Column>, String> {
    let (read, refr) = s
        .split_once('@')
        .ok_or_else(|| format!("trace input needs a '@' between read and reference: '{s}'"))?;
    let (read, refr) = (read.trim(), refr.trim());
    if read.chars().count() != refr.chars().count() {
        return Err(format!(
            "read is {} columns, reference is {}; a trace is a column-by-column \
             alignment, so the two must line up",
            read.chars().count(),
            refr.chars().count()
        ));
    }
    if read.is_empty() {
        return Err("trace input is empty".to_string());
    }

    // Trace input is *observed* data, not a pattern. `Seq::from_name` is the
    // right decoder because it understands `_` and `Z` as well as IUPAC, but
    // `~` has to be rejected: an observation is a concrete value, and feeding
    // "anything" in would match only patterns that are themselves `~`, which
    // looks like a bug in the query rather than in the input.
    let decode = |c: char, side: &str, i: usize| -> Result<Seq, String> {
        if c == '~' {
            return Err(format!(
                "{side} column {i}: '~' means \"anything\" in a pattern, but a trace takes \
                 observed data. Use a concrete base, 'N', '.' for a gap, '_' for a pad, ':' \
                 for a clip or ',' for a junction."
            ));
        }
        Seq::from_name(&c.to_string()).filter(|s| !s.is_empty()).ok_or_else(|| {
            format!("{side} column {i}: '{c}' is not a base, a gap, a pad, a clip or a junction")
        })
    };

    read.chars()
        .zip(refr.chars())
        .enumerate()
        .map(|(i, (r, f))| {
            Ok(Column {
                read: decode(r, "read", i)?,
                refr: decode(f, "reference", i)?,
                refr_pos: i as i64,
                read_off: i as i64,
                qual: -1,
            })
        })
        .collect()
}

/// How far pattern `p` has matched, given the automaton's live traces.
///
/// A pattern occupies `len` adjacent bits ending at `end_bits[p]`, with higher
/// bits meaning longer prefixes, so the highest live bit is the current match
/// length. Zero means the pattern is not started.
fn progress(state: &Bitmap, end_bit: usize, len: usize) -> usize {
    let start = end_bit + 1 - len;
    (start..=end_bit)
        .rev()
        .find(|&b| state.get(b) != 0)
        .map_or(0, |b| b - start + 1)
}

/// One pattern's history over the whole record.
struct PatternRun {
    label: String,
    /// Match length at each column.
    progress: Vec<usize>,
    /// Columns where the pattern completed.
    hits: Vec<usize>,
}

impl PatternRun {
    /// Furthest this pattern ever got, and the first column where it did.
    ///
    /// First rather than last: when a pattern reaches the same depth several
    /// times, the earliest is where it first ran out of road, which is what you
    /// want to look at. `max_by_key` would hand back the last one.
    ///
    /// Written as a fold rather than an iterator chain because a closure taking
    /// `&(usize, &usize)` cannot destructure the inner reference -- explicitly
    /// dereferencing inside an already-borrowing pattern is rejected under the
    /// 2024 binding-mode rules.
    fn best(&self) -> (usize, usize) {
        let mut best = (0usize, 0usize);
        for (i, &n) in self.progress.iter().enumerate() {
            if n > best.0 {
                best = (n, i);
            }
        }
        best
    }
}

/// Render an expression with each leaf's truth value at one position.
fn annotate(e: &Expr, spec: &QuerySpec, truth: &[bool], indent: usize) -> String {
    let pad = " ".repeat(indent);
    match e {
        Expr::Const(b) => format!("{pad}{b}\n"),
        Expr::Pattern(i) => {
            format!("{pad}{:<24} {}\n", spec.operand_label(*i), truth[*i])
        }
        Expr::Not(x) => match &**x {
            Expr::Pattern(i) => format!(
                "{pad}{:<24} {}\n",
                format!("not {}", spec.operand_label(*i)),
                !truth[*i]
            ),
            other => {
                let mut s = format!("{pad}not\n");
                s.push_str(&annotate(other, spec, truth, indent + 2));
                s
            }
        },
        Expr::And(xs) | Expr::Or(xs) => {
            let op = if matches!(e, Expr::And(_)) { "and" } else { "or" };
            let mut s = format!("{pad}{op}\n");
            for x in xs {
                s.push_str(&annotate(x, spec, truth, indent + 2));
            }
            s
        }
        Expr::Count { lo, hi, patterns } => {
            let n = patterns.iter().filter(|&&p| truth[p]).count();
            format!("{pad}{lo}..{hi} of {} matched\n", n)
        }
    }
}

/// What a column is, for the header grid.
fn cell(s: Seq) -> char {
    s.to_char()
}

/// Trace `specs` (compiled as `qs`) over one record's columns.
///
/// `label` names the record in the output. `verbose` prints the per-column
/// progress grid, which is the part that tells you where a match broke.
pub fn trace(
    qs: &QuerySet,
    specs: &[QuerySpec],
    cols: &[Column],
    label: &str,
    verbose: bool,
) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "{label}  ({} columns)\n", cols.len());

    if cols.is_empty() {
        out.push_str("  no columns\n");
        return out;
    }

    let (runs, fired) = replay(qs, specs, cols);

    if verbose {
        write_grid(&mut out, qs, cols, &runs);
    }

    write_queries(&mut out, qs, specs, cols, &runs, &fired);
    out
}

/// Per query, per column: where it fired and which operands were true there.
type Fired = Vec<Vec<(usize, Vec<bool>)>>;

/// Run the automaton over `cols`, keeping what the report needs from every
/// position rather than only the final state.
///
/// Tracing therefore reads `State` and not only `Hits`. A pattern that never
/// completes still has live traces, and how far they got is what says why the
/// pattern did not fire.
fn replay(qs: &QuerySet, specs: &[QuerySpec], cols: &[Column]) -> (Vec<PatternRun>, Fired) {
    let mut state = qs.plan.new_state();
    let mut hits = qs.plan.new_hits();

    // Label an interned pattern with the name some query gave it, falling back
    // to its text. Patterns are shared, so the name may come from any query
    // that uses it -- and within a run a name denotes one pattern, which
    // `dsl::validate_set` guarantees.
    let mut runs: Vec<PatternRun> = qs
        .patterns
        .iter()
        .enumerate()
        .map(|(pi, p)| PatternRun {
            label: specs
                .iter()
                .find_map(|s| {
                    (0..s.operands.len())
                        .find(|&i| pattern_id(qs, s, i) == Some(pi))
                        .and_then(|i| s.operand_names.get(i).filter(|n| !n.is_empty()).cloned())
                })
                .unwrap_or_else(|| p.text.clone()),
            progress: Vec::with_capacity(cols.len()),
            hits: Vec::new(),
        })
        .collect();
    let mut fired: Fired = vec![Vec::new(); specs.len()];

    for (ci, c) in cols.iter().enumerate() {
        state.step(QuerySet::compat_idx(c.read, c.refr), &mut hits);

        for (pi, run) in runs.iter_mut().enumerate() {
            let end = qs.plan.end_bits[pi];
            run.progress.push(progress(&state.inner, end, qs.plan.pattern_lengths[pi]));
            if hits.inner.get(end) != 0 {
                run.hits.push(ci);
            }
        }

        for (qi, (spec, q)) in specs.iter().zip(qs.queries.iter()).enumerate() {
            if q.span > ci + 1 {
                continue; // the window has not been filled yet
            }
            if q.predicate.eval(hits.inner.words()) {
                let truth = operand_truth(qs, spec, &hits.inner);
                fired[qi].push((ci, truth));
            }
        }
    }
    (runs, fired)
}

/// The column grid and each pattern's progress across it.
fn write_grid(out: &mut String, qs: &QuerySet, cols: &[Column], runs: &[PatternRun]) {
    {
        let w = cols.len().saturating_sub(1).to_string().len().max(1);
        let idx: String = (0..cols.len())
            .map(|i| format!("{i:>w$} "))
            .collect();
        let rd: String = cols.iter().map(|c| format!("{:>w$} ", cell(c.read))).collect();
        let rf: String = cols.iter().map(|c| format!("{:>w$} ", cell(c.refr))).collect();
        let _ = writeln!(out, "  {:<12}{idx}", "col");
        let _ = writeln!(out, "  {:<12}{rd}", "read");
        let _ = writeln!(out, "  {:<12}{rf}", "refr");

        out.push_str("\n  pattern progress: columns matched so far, * where complete\n");
        let namew = runs.iter().map(|r| r.label.len()).max().unwrap_or(0).min(28);
        for (pi, r) in runs.iter().enumerate() {
            let row: String = r
                .progress
                .iter()
                .map(|&n| {
                    let done = n == qs.plan.pattern_lengths[pi];
                    let s = if n == 0 {
                        ".".to_string()
                    } else if done {
                        "*".to_string()
                    } else {
                        n.to_string()
                    };
                    format!("{s:>w$} ")
                })
                .collect();
            let short: String = r.label.chars().take(namew).collect();
            let _ = writeln!(out, "  {short:<namew$}  {row}");
        }
        out.push('\n');
    }
}

/// What fired where, and for anything that did not, which operand fell short.
fn write_queries(
    out: &mut String,
    qs: &QuerySet,
    specs: &[QuerySpec],
    cols: &[Column],
    runs: &[PatternRun],
    fired: &Fired,
) {
    for (qi, spec) in specs.iter().enumerate() {
        let name = if spec.name.is_empty() { "(unnamed)" } else { &spec.name };
        if fired[qi].is_empty() {
            let _ = writeln!(out, "  {name}: no match");
            // The useful part: which operand never completed, and where it got
            // closest. A query fails because one of its patterns did.
            for i in 0..spec.operands.len() {
                let Some(pi) = pattern_id(qs, spec, i) else { continue };
                let r = &runs[pi];
                let len = qs.plan.pattern_lengths[pi];
                if r.hits.is_empty() {
                    let (best, at) = r.best();
                    let _ = if best == 0 {
                        // "best at column 7" would be noise: it never started,
                        // so no column is more informative than any other.
                        writeln!(
                            out,
                            "    {:<24} never started; no column matches its first",
                            spec.operand_label(i)
                        )
                    } else {
                        writeln!(
                            out,
                            "    {:<24} never completed; got {best} of {len} columns, \
                             furthest at column {at}",
                            spec.operand_label(i)
                        )
                    };
                } else {
                    let _ = writeln!(
                        out,
                        "    {:<24} completed at {}",
                        spec.operand_label(i),
                        join(&r.hits)
                    );
                }
            }
        } else {
            let _ = writeln!(
                out,
                "  {name}: matches ending at column {}",
                join(&fired[qi].iter().map(|(c, _)| *c).collect::<Vec<_>>())
            );
            for (ci, truth) in &fired[qi] {
                let anchor = *ci + 1 - specs[qi].span + specs[qi].anchor;
                let a = &cols[anchor];
                let _ = writeln!(
                    out,
                    "    column {ci}, anchor column {anchor}: read {}, refr {}, \
                     walk offset {}, refr_pos {} (0-based)",
                    cell(a.read),
                    cell(a.refr),
                    if a.read == Seq::PAD { "none".to_string() } else { a.read_off.to_string() },
                    a.refr_pos
                );
                out.push_str(&annotate(&spec.expr, spec, truth, 6));
            }
        }
    }
}

/// Which of a query's operands are satisfied right now.
fn operand_truth(qs: &QuerySet, spec: &QuerySpec, hits: &Bitmap) -> Vec<bool> {
    (0..spec.operands.len())
        .map(|i| {
            pattern_id(qs, spec, i).is_some_and(|pi| hits.get(qs.plan.end_bits[pi]) != 0)
        })
        .collect()
}

/// The pattern id for a query's operand.
///
/// `QuerySet` keeps one copy of a pattern that several queries share. The
/// mapping is therefore by columns, and not by position. The search is linear.
/// That costs little, because a trace runs once for each record, and not once
/// for each column.
fn pattern_id(qs: &QuerySet, spec: &QuerySpec, operand: usize) -> Option<usize> {
    let want = &spec.operands.get(operand)?.columns;
    qs.patterns.iter().position(|p| &p.columns == want)
}

fn join(v: &[usize]) -> String {
    v.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compile(toml: &str) -> (QuerySet, Vec<QuerySpec>) {
        let specs = crate::test_support::queries(toml);
        (QuerySet::compile(&specs).unwrap(), specs)
    }

    fn run(toml: &str, pair: &str) -> String {
        let (qs, specs) = compile(toml);
        let cols = columns_from_pair(pair).unwrap();
        trace(&qs, &specs, &cols, "test", true)
    }

    #[test]
    fn synthetic_columns_parse() {
        let c = columns_from_pair("TTC@TTG").unwrap();
        assert_eq!(c.len(), 3);
        assert_eq!(c[2].read, Seq::C);
        assert_eq!(c[2].refr, Seq::G);
        assert_eq!(c[2].read_off, 2);
        assert!(columns_from_pair("TTC@TT").unwrap_err().contains("3 columns"));
        assert!(columns_from_pair("TTC").unwrap_err().contains("'@'"));
        assert!(columns_from_pair("TT?@TTG").unwrap_err().contains("column 2"));
    }

    #[test]
    fn gaps_and_pads_are_accepted() {
        let c = columns_from_pair("_A.@CAG").unwrap();
        assert_eq!(c[0].read, Seq::PAD);
        assert_eq!(c[2].read, Seq::GAP);
    }

    #[test]
    fn a_match_reports_the_position_and_the_reasoning() {
        let t = run("[pat.p1]\nread = \"Y~\"\nrefr = \"CG\"\n\n[query.mCpG]\nmark = \"+.\"\nwhere = \"p1\"\n", "TTCGTT@TTCGTT");
        // The CpG completes on column 3.
        assert!(t.contains("mCpG: matches ending at column 3"), "{t}");
        assert!(t.contains("anchor column 2"), "{t}");
        assert!(t.contains("p1"), "{t}");
        assert!(t.contains("true"), "{t}");
    }

    #[test]
    fn the_progress_grid_shows_where_a_match_broke() {
        // GATCGATC against a reference that diverges at column 4.
        let t = run("[pat.p1]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.j]\nmark = \"+.......\"\nwhere = \"p1\"\n", "GATCTATC@GATCTATC");
        assert!(t.contains("no match"), "{t}");
        // It got four columns in before the mismatch.
        assert!(t.contains("got 4 of 8 columns"), "{t}");
        assert!(t.contains("furthest at column 3"), "{t}");
    }

    #[test]
    fn an_exclusion_that_fires_is_shown_as_the_reason() {
        let (qs, specs) = compile("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n[pat.p2]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.m]\nmark = \"...+....\"\nwhere = \"p1 and not p2\"\n");

        // A genuine CpG: fires, with the exclusion false.
        let good = columns_from_pair("TTTCGTTT@TTTCGTTT").unwrap();
        let t = trace(&qs, &specs, &good, "good", false);
        assert!(t.contains("m: matches ending at column 7"), "{t}");
        assert!(t.contains("not p2") && t.contains("true"), "{t}");

        // The same CpG inside a junction: suppressed, and the trace says which
        // operand did it.
        let bad = columns_from_pair("GATCGATC@GATCGATC").unwrap();
        let t = trace(&qs, &specs, &bad, "junction", false);
        assert!(t.contains("m: no match"), "{t}");
        assert!(t.contains("p1") && t.contains("completed at 7"), "{t}");
        assert!(t.contains("p2") && t.contains("completed at 7"), "{t}");
    }

    #[test]
    fn named_operands_are_used_when_present() {
        let specs = crate::test_support::queries(
            "[pat.cpg]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n\
             [pat.junc]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\
             [query.m]\nmark = \"...+....\"\nwhere = \"cpg and not junc\"\n",
        );
        let qs = QuerySet::compile(&specs).unwrap();
        let cols = columns_from_pair("GATCGATC@GATCGATC").unwrap();
        let t = trace(&qs, &specs, &cols, "r", false);
        assert!(t.contains("cpg"), "{t}");
        assert!(t.contains("junc"), "{t}");
        assert!(!t.contains("p1"), "{t}");
        // The progress grid uses the names too, not the raw pattern text.
        let t = trace(&qs, &specs, &cols, "r", true);
        assert!(t.lines().any(|l| l.trim_start().starts_with("cpg  ")), "{t}");
        assert!(!t.contains("~~~Y~~~~@"), "{t}");
    }

    #[test]
    fn a_wildcard_in_trace_input_is_rejected() {
        let e = columns_from_pair("~C@GC").unwrap_err();
        assert!(e.contains("observed data"), "{e}");
        assert!(e.contains("read column 0"), "{e}");
    }

    #[test]
    fn the_furthest_column_is_the_first_one_that_got_there() {
        // A pattern that stalls at the same depth repeatedly should point at
        // where it first stalled, not the last time it was still stuck.
        let r = PatternRun {
            label: "p".into(),
            progress: vec![1, 2, 3, 2, 3, 1],
            hits: Vec::new(),
        };
        assert_eq!(r.best(), (3, 2));

        // All zeros: nothing to point at.
        let r = PatternRun { label: "p".into(), progress: vec![0; 5], hits: Vec::new() };
        assert_eq!(r.best(), (0, 0));
    }

    #[test]
    fn a_pattern_that_never_starts_is_reported_as_zero() {
        let t = run("[query.x]\nread = \"AAAA\"\nrefr = \"AAAA\"\n", "TTTT@TTTT");
        assert!(t.contains("never started"), "{t}");
        assert!(!t.contains("furthest at"), "{t}");
    }

    #[test]
    fn several_queries_are_traced_together() {
        let t = run("[pat.p1]\nread = \"C~\"\nrefr = \"CG\"\n[pat.p2]\nread = \"~\"\nrefr = \"/\"\n\n[query.a]\nmark = \"+.\"\nwhere = \"p1\"\n\n[query.b]\nmark = \"+\"\nwhere = \"p2\"\n", "CG@CA");
        assert!(t.contains("a: no match"), "{t}");
        assert!(t.contains("b: matches ending at column 1"), "{t}");
    }

    #[test]
    fn the_window_must_be_full_before_a_query_can_fire() {
        // A negated-only operand would otherwise "match" at column 0 of an
        // 8-wide query, before eight columns have even been seen.
        let t = run("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n[pat.p2]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.m]\nmark = \"...+....\"\nwhere = \"p1 and not p2\"\n", "TTTCGTTT@TTTCGTTT");
        assert!(t.contains("matches ending at column 7"), "{t}");
        assert!(!t.contains("matches ending at column 0"), "{t}");
    }
}
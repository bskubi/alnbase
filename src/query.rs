//! Compile parsed queries into an automaton and one predicate per query.
//!
//! Know these three properties before you read the rest.
//!
//! **The run keeps one copy of each pattern.** Two queries that name the same
//! pattern share one set of bits in the automaton. An overlap between two
//! queries therefore costs the boolean layer nothing extra. Queries do overlap
//! in practice. An exclusion clause such as `not ~8@GATCGATC` is usually
//! repeated in every query of a run. The key is the set of lowered columns, and
//! not the text. Two spellings of the same pattern therefore share one set of
//! bits as well. One spelling can come through an alias, and the other can be
//! written out in full.
//!
//! **The compat table understands relational codes.** To build the compat set of
//! a column, the code lists the observed values that are subsets of each side,
//! and then filters that list by the relation. `=` and `/` therefore cost
//! nothing during the scan, because this module resolves them once. The test is
//! subset and not equality, because an observed value is itself a set. An IUPAC
//! ambiguity code in the reference is several nucleotide bits at the same time.
//! It satisfies a column when every value that it could be is a value that the
//! column allows.
//!
//! **A query is a predicate, not a pattern.** The two layers meet at `end_bits`,
//! which is the `bit_of` table that [`crate::predicate`] needs, together with
//! the width of the bitmap. Nothing about columns or bases crosses that line.
//!
//! This module also refuses a query that can never do useful work, instead of
//! letting it disappoint the user during the scan. A column that no observed
//! value satisfies is a `DeadColumn`. A query that an empty window already
//! satisfies is an `AlwaysOnEmpty`. A `not X` with no positive term is an
//! example of the second, and it would fire almost everywhere. That check is
//! semantic and not syntactic: `not X` alone is refused, but `~@= and not X` is
//! accepted. To tell the two apart, the code asks the compiled predicate what it
//! says about an empty bitmap. Both checks need the compiled automaton, which is
//! why they are here and not in [`crate::dsl`].

use std::collections::HashMap;
use std::fmt;

use crate::byg_search::Plan;
use crate::dsl::{ColumnSpec, PatternSpec, QuerySpec};
use crate::predicate::{self, AnyPredicate, Expr, Predicate};
use crate::seq::Seq;

/// The limit on DNF and CNF expansion. Above this limit the tree-walking
/// fallback takes over. alnbase logs the size for each query at startup, so that
/// real usage can show whether the limit is set well.
const NORMAL_FORM_BUDGET: usize = 4096;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CompileError {
    /// No observation can satisfy this column. Its operand can therefore never
    /// match, and the query would never fire and never say why.
    DeadColumn { query: String, pattern: String, column: usize, why: String },
    /// A window with no pattern hits at all satisfies the query. Such a query
    /// is true at every padded position past the end of every read.
    AlwaysOnEmpty { query: String },
    Predicate { query: String, source: predicate::CompileError },
}

impl fmt::Display for CompileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CompileError::DeadColumn { query, pattern, column, why } => write!(
                f,
                "in query '{query}', pattern '{pattern}' column {column}: {why}. \
                 No observation can satisfy this column, so the pattern can never match."
            ),
            CompileError::AlwaysOnEmpty { query } => write!(
                f,
                "query '{query}' is satisfied by a window containing no matches at all, so it \
                 would fire at essentially every position, including past the ends of reads. \
                 Combine the negation with at least one positive pattern."
            ),
            CompileError::Predicate { query, source } => {
                write!(f, "in query '{query}': {source}")
            }
        }
    }
}

impl std::error::Error for CompileError {}

/// A column that parsed correctly but adds nothing to the query. This is not a
/// fatal error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Warning {
    pub query: String,
    pub pattern: String,
    pub column: usize,
    pub message: String,
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "warning: in query '{}', pattern '{}' column {}: {}",
            self.query, self.pattern, self.column, self.message
        )
    }
}

/// One compiled query: what makes it fire, where it reports, and what it
/// records.
#[derive(Debug)]
pub struct Query {
    pub name: String,
    /// The expression for the output column, with every operand written out in
    /// the terse base-level form.
    ///
    /// This is not the text that the user wrote. A query in matrix form says
    /// `cpg and not junc_end`, and those names live in the file that declared
    /// them. You could read a hits row that carried those names only if you also
    /// had that file. The code therefore resolves the names one time, here. The
    /// substitution is the same for every row that a query ever writes, and a
    /// query can write very many rows.
    pub predicate: AnyPredicate,
    /// The column whose coordinates the output reports.
    pub anchor: usize,
    /// The columns to record, in ascending order.
    pub captures: Box<[usize]>,
    pub span: usize,
    /// The `ColumnRing::get` argument for the anchor column.
    ///
    /// The ring increments its index after a store, so `get(1)` is the column
    /// that the ring stored last. That column is the *last* column of the span,
    /// and column `i` is `span - i` positions back. An error of one here moves
    /// every reported coordinate and reports nothing. The code therefore
    /// resolves the value once, here, and the scan loop does no arithmetic.
    pub anchor_back: usize,
    /// The `ColumnRing::get` arguments for the captures, in the same order as
    /// `captures`.
    pub capture_back: Box<[usize]>,
}

/// Everything that the scanner needs: one automaton and one predicate per
/// query.
#[derive(Debug)]
pub struct QuerySet {
    pub queries: Box<[Query]>,
    /// The patterns, with one copy of each. The index of a pattern is the
    /// pattern id that the automaton uses.
    pub patterns: Box<[PatternSpec]>,
    pub plan: Plan,
    /// The length of the longest pattern. It sets the flank of the walk and the
    /// depth of the column ring. Both read it from here, so one place derives
    /// it.
    pub max_span: usize,
    /// True when every query records its anchor column and nothing else. This
    /// is the common case, because the anchor is always captured.
    ///
    /// The output can then be flat. The coordinates and the quality of the
    /// capture are those of the anchor, which the row already holds as scalars,
    /// so only the observed bases are new. Lists would have length one, and four
    /// of their five columns would repeat what the row already says.
    pub flat_captures: bool,
    pub warnings: Vec<Warning>,
}

impl QuerySet {
    /// One line per query: the normal form that its predicate compiled to, the
    /// size of that form, and the span of the query.
    ///
    /// `--explain` prints these lines after the detail for each query. The
    /// detail says what one query means. These lines say what the whole *set*
    /// cost to compile, which is worth a look before a long run. A query that
    /// fell back to the tree walk, or a span that is much wider than you
    /// intended, appears here as one unusual line among many ordinary ones.
    pub fn shapes(&self) -> Vec<String> {
        self.queries
            .iter()
            .map(|q| {
                let (form, n) = q.predicate.shape();
                format!("{}: {form}, {n} groups, span {}", q.name, q.span)
            })
            .collect()
    }

    /// The compat index of an observed symbol. It must match the layout that
    /// `fill_compat` writes.
    #[inline]
    pub fn compat_idx(read: Seq, refr: Seq) -> usize {
        read.0 as usize * Seq::COUNT + refr.0 as usize
    }

    pub fn compile(specs: &[QuerySpec]) -> Result<Self, CompileError> {
        let mut patterns: Vec<PatternSpec> = Vec::new();
        let mut interned: HashMap<Box<[ColumnSpec]>, usize> = HashMap::new();
        // Per query: local operand index -> global pattern id.
        let mut mapping: Vec<Vec<usize>> = Vec::with_capacity(specs.len());

        for spec in specs {
            let mut ids = Vec::with_capacity(spec.operands.len());
            for op in &spec.operands {
                let id = *interned.entry(op.columns.clone()).or_insert_with(|| {
                    patterns.push(op.clone());
                    patterns.len() - 1
                });
                ids.push(id);
            }
            mapping.push(ids);
        }

        let lengths: Vec<usize> = patterns.iter().map(|p| p.len()).collect();
        let mut plan = Plan::from_lengths(&lengths, Seq::COUNT * Seq::COUNT);
        let max_span = lengths.iter().copied().max().unwrap_or(0);

        let warnings = fill_compat(&mut plan, &patterns, specs, &mapping)?;

        // --- one predicate per query ---------------------------------------
        let mut queries = Vec::with_capacity(specs.len());
        for (spec, ids) in specs.iter().zip(mapping.iter()) {
            let expr = remap(&spec.expr, ids);
            let predicate =
                predicate::compile(&expr, &plan.end_bits, plan.word_count, NORMAL_FORM_BUDGET)
                    .map_err(|source| CompileError::Predicate {
                        query: spec.text.clone(),
                        source,
                    })?;

            // The empty-window check is semantic, not syntactic: `not X` alone
            // is rejected, `~@= and not X` is fine, and the engine already
            // knows the answer.
            if predicate.zero_value() {
                return Err(CompileError::AlwaysOnEmpty { query: spec.text.clone() });
            }

            queries.push(Query {
                name: spec.name.clone(),
                predicate,
                anchor: spec.anchor,
                capture_back: spec.captures.iter().map(|&i| spec.span - i).collect(),
                captures: spec.captures.clone(),
                span: spec.span,
                anchor_back: spec.span - spec.anchor,
            });
        }

        let flat_captures = queries
            .iter()
            .all(|q| q.captures.len() == 1 && q.captures[0] == q.anchor);

        Ok(QuerySet {
            flat_captures,
            queries: queries.into_boxed_slice(),
            patterns: patterns.into_boxed_slice(),
            plan,
            max_span,
            warnings,
        })
    }

}

/// Rewrite an expression that uses query-local pattern ids so that it uses the
/// global pattern ids.
fn remap(e: &Expr, ids: &[usize]) -> Expr {
    match e {
        Expr::Const(b) => Expr::Const(*b),
        Expr::Pattern(p) => Expr::Pattern(ids[*p]),
        Expr::Not(x) => Expr::not(remap(x, ids)),
        Expr::And(xs) => Expr::And(xs.iter().map(|x| remap(x, ids)).collect()),
        Expr::Or(xs) => Expr::Or(xs.iter().map(|x| remap(x, ids)).collect()),
        Expr::Count { lo, hi, patterns } => Expr::Count {
            lo: *lo,
            hi: *hi,
            patterns: patterns.iter().map(|p| ids[*p]).collect(),
        },
    }
}

/// Set the compat bits of every pattern column. This function also checks that
/// each column is satisfiable, because it already enumerates the values that the
/// check needs.
fn fill_compat(
    plan: &mut Plan,
    patterns: &[PatternSpec],
    specs: &[QuerySpec],
    mapping: &[Vec<usize>],
) -> Result<Vec<Warning>, CompileError> {
    // For error messages: which query first mentioned each pattern.
    let mut owner: HashMap<usize, &str> = HashMap::new();
    for (spec, ids) in specs.iter().zip(mapping.iter()) {
        for &id in ids {
            owner.entry(id).or_insert(&spec.text);
        }
    }

    let mut warnings = Vec::new();
    let mut bit = 0usize;

    for (pid, pattern) in patterns.iter().enumerate() {
        for (ci, col) in pattern.columns.iter().enumerate() {
            let mut n_set = 0usize;
            // Which declared read/refr bits actually appear in a live pair.
            let mut live_read = 0u8;
            let mut live_refr = 0u8;

            for read_obs in col.read.bit_combinations() {
                for refr_obs in col.refr.bit_combinations() {
                    if let Some(rel) = col.rel {
                        if !rel.holds(read_obs, refr_obs) {
                            continue;
                        }
                    }
                    plan.set_compat(QuerySet::compat_idx(read_obs, refr_obs), bit);
                    live_read |= read_obs.0;
                    live_refr |= refr_obs.0;
                    n_set += 1;
                }
            }

            let query = owner.get(&pid).copied().unwrap_or("");
            if n_set == 0 {
                return Err(CompileError::DeadColumn {
                    query: query.to_string(),
                    pattern: pattern.text.clone(),
                    column: ci,
                    why: match col.rel {
                        Some(r) => format!(
                            "'{}' requires unambiguous A/C/G/T on both sides, which no value \
                             here can provide",
                            r.to_char()
                        ),
                        // Unreachable in practice: without a relation the two
                        // sides constrain independently, and an empty set is
                        // already rejected at parse time. Kept as a guard.
                        None => "the two sides have no compatible values".to_string(),
                    },
                });
            }

            for (declared, live, side) in [
                (col.read.0, live_read, "read"),
                (col.refr.0, live_refr, "reference"),
            ] {
                let dead = declared & !live;
                if dead != 0 {
                    warnings.push(Warning {
                        query: query.to_string(),
                        pattern: pattern.text.clone(),
                        column: ci,
                        message: format!(
                            "'{}' on the {side} side is unreachable; the column reduces to '{}'",
                            Seq(dead).name(),
                            Seq(live).name()
                        ),
                    });
                }
            }

            bit += 1;
        }
    }

    Ok(warnings)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(toml: &str) -> QuerySet {
        QuerySet::compile(&crate::test_support::queries(toml)).unwrap_or_else(|e| panic!("{e}"))
    }

    fn err(toml: &str) -> CompileError {
        QuerySet::compile(&crate::test_support::queries(toml)).unwrap_err()
    }

    /// Feed a sequence of columns through the automaton, and report which
    /// queries fire at each position.
    fn run(qs: &QuerySet, cols: &[(Seq, Seq)]) -> Vec<Vec<usize>> {
        let mut state = qs.plan.new_state();
        let mut hits = qs.plan.new_hits();
        let mut out = Vec::new();
        for &(r, f) in cols {
            state.step(QuerySet::compat_idx(r, f), &mut hits);
            let words = hits.inner.words();
            out.push(
                qs.queries
                    .iter()
                    .enumerate()
                    .filter(|(_, q)| q.predicate.eval(words))
                    .map(|(i, _)| i)
                    .collect(),
            );
        }
        out
    }

    fn cols(read: &str, refr: &str) -> Vec<(Seq, Seq)> {
        read.chars()
            .zip(refr.chars())
            .map(|(a, b)| (Seq::FROM_FASTA[a as usize].unwrap(), Seq::FROM_FASTA[b as usize].unwrap()))
            .collect()
    }

    #[test]
    fn a_plain_pattern_fires_at_the_last_column() {
        let qs = set("[pat.p1]\nread = \"Y~\"\nrefr = \"CG\"\n\n[query.unnamed]\nmark = \"+.\"\nwhere = \"p1\"\n");
        //                    0    1    2    3
        let fired = run(&qs, &cols("TTCA", "TTCG"));
        assert_eq!(fired[0], Vec::<usize>::new());
        assert_eq!(fired[1], Vec::<usize>::new());
        assert_eq!(fired[2], Vec::<usize>::new());
        // The CG pattern completes on column 3.
        assert_eq!(fired[3], vec![0]);
    }

    #[test]
    fn patterns_are_interned_across_queries() {
        let qs = set("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n[pat.p2]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.a]\nmark = \"+.......\"\nwhere = \"p1 and not p2\"\n\n[query.b]\nmark = \"+.......\"\nwhere = \"p2 or p1\"\n");
        // Four operand mentions, two distinct patterns.
        assert_eq!(qs.patterns.len(), 2);
        assert_eq!(qs.max_span, 8);
    }

    #[test]
    fn negation_needs_a_positive_operand() {
        assert_eq!(
            err("[pat.p1]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.q]\nmark = \"+.......\"\nwhere = \"not p1\"\n"),
            CompileError::AlwaysOnEmpty { query: "not p1".into() }
        );
        // With a positive operand it compiles.
        set("[pat.p1]\nread = \"~~~~~~~~\"\nrefr = \"========\"\n[pat.p2]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.unnamed]\nmark = \"+.......\"\nwhere = \"p1 and not p2\"\n");
    }

    #[test]
    fn boolean_exclusion_works_end_to_end() {
        // A CpG at columns 3-4 of an 8-column window, excluded when the
        // reference reads GATCGATC — whose own CG sits at exactly 3-4.
        let qs = set("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n[pat.p2]\nread = \"~~~~~~~~\"\nrefr = \"GATCGATC\"\n\n[query.mCpG]\nmark = \"...+....\"\nwhere = \"p1 and not p2\"\n");
        assert_eq!(qs.queries[0].span, 8);
        assert_eq!(qs.queries[0].anchor, 3);
        assert_eq!(&*qs.queries[0].captures, &[3]);

        // A genuine CpG in a non-junction context: fires on the last column.
        let good = run(&qs, &cols("TTTCGTTT", "TTTCGTTT"));
        assert_eq!(good[7], vec![0]);

        // The same CpG inside a ligation junction: excluded.
        let junction = run(&qs, &cols("GATCGATC", "GATCGATC"));
        assert_eq!(junction[7], Vec::<usize>::new());
    }

    #[test]
    fn relational_columns_match_and_mismatch() {
        let qs = set("[pat.p1]\nread = \"~\"\nrefr = \"=\"\n[pat.p2]\nread = \"~\"\nrefr = \"/\"\n\n[query.m]\nmark = \"+\"\nwhere = \"p1\"\n\n[query.x]\nmark = \"+\"\nwhere = \"p2\"\n");
        let fired = run(&qs, &cols("ACGT", "ACTT"));
        assert_eq!(fired[0], vec![0]); // A vs A
        assert_eq!(fired[1], vec![0]); // C vs C
        assert_eq!(fired[2], vec![1]); // G vs T
        assert_eq!(fired[3], vec![0]); // T vs T
    }

    #[test]
    fn relational_codes_ignore_gaps_pads_and_ambiguity() {
        let qs = set("[pat.p1]\nread = \"~\"\nrefr = \"=\"\n[pat.p2]\nread = \"~\"\nrefr = \"/\"\n\n[query.m]\nmark = \"+\"\nwhere = \"p1\"\n\n[query.x]\nmark = \"+\"\nwhere = \"p2\"\n");
        let none: Vec<usize> = vec![];
        // Neither fires: they are not complements.
        assert_eq!(run(&qs, &[(Seq::GAP, Seq::GAP)])[0], none);
        assert_eq!(run(&qs, &[(Seq::PAD, Seq::PAD)])[0], none);
        assert_eq!(run(&qs, &[(Seq::N, Seq::C)])[0], none);
        assert_eq!(run(&qs, &[(Seq::C, Seq::N)])[0], none);
    }

    #[test]
    fn dead_columns_are_rejected() {
        let e = err("[query.q]\nread = \".\"\nrefr = \"=\"\n");
        match e {
            CompileError::DeadColumn { column, .. } => assert_eq!(column, 0),
            other => panic!("{other:?}"),
        }
        assert!(matches!(err("[alias]\ng = \"{._}\"\n\n[query.q]\nread = \"g\"\nrefr = \"/\"\n"), CompileError::DeadColumn { .. }));
        // Without a relation the sides constrain independently, so a column
        // whose sides disagree is a mismatch, not a contradiction.
        set("[pat.p1]\nread = \"M\"\nrefr = \"G\"\n\n[query.unnamed]\nmark = \"+\"\nwhere = \"p1\"\n");
    }

    #[test]
    fn unreachable_components_warn_without_failing() {
        // The gap can never satisfy '=', so the column reduces to (C, C).
        let qs = set("[alias]\nj = \"{C.}\"\n\n[pat.p1]\nread = \"j\"\nrefr = \"=\"\n\n[query.unnamed]\nmark = \"+\"\nwhere = \"p1\"\n");
        assert_eq!(qs.warnings.len(), 2, "{:?}", qs.warnings);
        assert!(qs.warnings[0].message.contains("unreachable"), "{:?}", qs.warnings[0]);
    }

    #[test]
    fn backsteps_are_resolved_at_compile_time() {
        // The ring post-increments, so get(1) is the last column of the span
        // and column i is span - i back. Column 3 of an 8-wide span is 5.
        let qs = set("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n\n[query.unnamed]\nmark = \"...+....\"\nwhere = \"p1\"\n");
        let q = &qs.queries[0];
        assert_eq!(&*q.capture_back, &[5]);
        assert_eq!(q.anchor_back, 5);
        // Default anchor is column 0: the far end, matching the old
        // `ring.get(k_pair.len())`.
        let qs = set("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n\n[query.unnamed]\nmark = \"+.......\"\nwhere = \"p1\"\n");
        assert_eq!(qs.queries[0].anchor_back, 8);
        assert_eq!(qs.queries[0].span, 8);
    }

    #[test]
    fn flat_captures_tracks_the_common_case() {
        // Anchor only, explicit or defaulted: nothing to put in a list.
        assert!(set("[pat.p1]\nread = \"Y~\"\nrefr = \"CG\"\n\n[query.unnamed]\nmark = \"+.\"\nwhere = \"p1\"\n").flat_captures);
        assert!(set("[pat.p1]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n\n[query.unnamed]\nmark = \"...+....\"\nwhere = \"p1\"\n").flat_captures);
        assert!(set("[pat.p1]\nread = \"Y~\"\nrefr = \"CG\"\n[pat.p2]\nread = \"~~~Y~~~~\"\nrefr = \"~~~CG~~~\"\n\n[query.a]\nmark = \"+.\"\nwhere = \"p1\"\n\n[query.b]\nmark = \"...+....\"\nwhere = \"p2\"\n").flat_captures);
        // One extra capture anywhere turns the whole file into the list form,
        // since a shard has one schema.
        assert!(!set("[pat.p1]\nread = \"Y~\"\nrefr = \"CG\"\n\n[query.unnamed]\nmark = \"+^\"\nwhere = \"p1\"\n").flat_captures);
        assert!(!set("[pat.p1]\nread = \"Y~\"\nrefr = \"CG\"\n\n[query.a]\nmark = \"+.\"\nwhere = \"p1\"\n\n[query.b]\nmark = \"+^\"\nwhere = \"p1\"\n").flat_captures);
    }


    #[test]
    fn shapes_are_reportable() {
        let qs = set("[pat.p1]\nread = \"Y~\"\nrefr = \"CG\"\n[pat.p2]\nread = \"~~\"\nrefr = \"GA\"\n\n[query.mCpG]\nmark = \"+.\"\nwhere = \"p1 and not p2\"\n");
        let s = &qs.shapes()[0];
        assert!(s.starts_with("mCpG: "), "{s}");
        assert!(s.contains("span 2"), "{s}");
    }
}
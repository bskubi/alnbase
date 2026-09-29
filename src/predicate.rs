//! Evaluate a boolean predicate over a bitmap of pattern hits.
//!
//! This module knows nothing about BYG search, BAM records or CIGARs. It has
//! only two points of contact with the search layer:
//!
//!   * a slice `bit_of: &[usize]`, which maps a `PatternId` to a bit index. For
//!     the BYG search this slice is `Plan::end_bits`.
//!   * the `&[u64]` words of a hit bitmap, which the caller gives to
//!     `Predicate::eval`.
//!
//! To change the evaluation strategy you write a new `Predicate` impl. No code
//! upstream or downstream has to change.
//!
//! # Normal forms
//!
//! [`compile`] turns an expression into sum-of-products form (DNF) or
//! product-of-sums form (CNF). Both forms lower to the *same* machine shape: a
//! masked word compare, plus optional popcount guards. The two evaluators
//! therefore differ only in the sense of the comparison and in the direction of
//! the early exit:
//!
//! ```text
//! DNF group (a conjunction): satisfied if  words & mask == cmp   ... OR  over the groups
//! CNF group (a disjunction): satisfied if  words & mask != cmp   ... AND over the groups
//! ```
//!
//! Neither form is small in the general case, and each form grows very large on
//! the inputs that suit the other. `(a|b)&(c|d)&(e|f)` is 3 CNF groups and 8 DNF
//! groups, and its dual is the reverse. [`compile`] therefore expands both forms
//! and keeps the smaller one. If both forms are over the budget, it falls back
//! to a tree walk. Expressions that behave like parity are exponential in both
//! forms and always take the fallback.
//!
//! # Counts
//!
//! The popcount guards come from [`Expr::Count`], which means "between `lo` and
//! `hi` of these patterns hit". It lowers to a masked popcount, and not to the
//! `C(n, k)` conjunctions that the same statement needs in plain booleans.
//!
//! Nothing that a user can write reaches this code yet. [`crate::query_toml`]
//! has no syntax for a count, so the only callers of [`Expr::at_least`] and the
//! related functions are the tests in this module. The representation and the
//! evaluator are present and tested. It is the query language that has not
//! caught up. Know this before you look for the syntax, and before you decide
//! that the code is dead.

use std::collections::BTreeMap;
use std::fmt;

/// Index into the pattern list, *not* a bit index. Callers speak in pattern ids;
/// lowering to bit positions happens once, at compile time.
pub type PatternId = usize;

// ---------------------------------------------------------------------------
// Specification layer
// ---------------------------------------------------------------------------

/// A boolean expression over pattern hits. This is the stable interface a CLI
/// parser (or anything else) targets.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Expr {
    Const(bool),
    Pattern(PatternId),
    Not(Box<Expr>),
    And(Vec<Expr>),
    Or(Vec<Expr>),
    /// True when the number of listed patterns that hit falls in `lo..=hi`.
    ///
    /// A closed range rather than a bare threshold, so that negation stays
    /// inside the representation: the complement of a range is a union of at
    /// most two ranges, which the compiler handles without expanding to
    /// `C(n, k)` conjunctions. "At most 2 mismatches" is therefore just as cheap
    /// as "at least 3".
    Count {
        lo: usize,
        hi: usize,
        patterns: Vec<PatternId>,
    },
}

impl Expr {
    pub fn pat(p: PatternId) -> Expr {
        Expr::Pattern(p)
    }
    pub fn not(e: Expr) -> Expr {
        Expr::Not(Box::new(e))
    }
    pub fn all<I: IntoIterator<Item = Expr>>(xs: I) -> Expr {
        Expr::And(xs.into_iter().collect())
    }
    pub fn any<I: IntoIterator<Item = Expr>>(xs: I) -> Expr {
        Expr::Or(xs.into_iter().collect())
    }

    /// At least `k` of `patterns` hit.
    pub fn at_least(k: usize, patterns: Vec<PatternId>) -> Expr {
        let hi = patterns.len();
        Expr::Count { lo: k, hi, patterns }
    }
    /// At most `k` of `patterns` hit.
    pub fn at_most(k: usize, patterns: Vec<PatternId>) -> Expr {
        Expr::Count { lo: 0, hi: k, patterns }
    }
    /// Exactly `k` of `patterns` hit.
    pub fn exactly(k: usize, patterns: Vec<PatternId>) -> Expr {
        Expr::Count { lo: k, hi: k, patterns }
    }

    /// Every pattern the expression mentions, in ascending order.
    pub fn patterns(&self) -> Vec<PatternId> {
        let mut v = Vec::new();
        self.collect_patterns(&mut v);
        v.sort_unstable();
        v.dedup();
        v
    }

    fn collect_patterns(&self, out: &mut Vec<PatternId>) {
        match self {
            Expr::Const(_) => {}
            Expr::Pattern(p) => out.push(*p),
            Expr::Not(x) => x.collect_patterns(out),
            Expr::And(xs) | Expr::Or(xs) => {
                for x in xs {
                    x.collect_patterns(out)
                }
            }
            Expr::Count { patterns, .. } => out.extend_from_slice(patterns),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CompileError {
    /// A pattern id had no entry in the supplied `bit_of` table.
    UnknownPattern(PatternId),
    /// `bit_of` pointed past the end of the bitmap.
    BitOutOfRange { pattern: PatternId, bit: usize },
    /// Normalization exceeded the group budget. Soft: [`compile`] falls back.
    TooComplex { budget: usize },
}

impl fmt::Display for CompileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CompileError::UnknownPattern(p) => write!(f, "unknown pattern id {p}"),
            CompileError::BitOutOfRange { pattern, bit } => {
                write!(f, "pattern {pattern} maps to bit {bit}, past end of bitmap")
            }
            CompileError::TooComplex { budget } => {
                write!(f, "expression exceeded normalization budget of {budget} groups")
            }
        }
    }
}

impl std::error::Error for CompileError {}

// ---------------------------------------------------------------------------
// Evaluation interface
// ---------------------------------------------------------------------------

/// A compiled predicate. Object safe, but prefer monomorphising the record loop
/// over `P: Predicate` so dispatch is hoisted out of the per-symbol work.
pub trait Predicate {
    /// `words` is the raw hit bitmap. Must be at least `word_count` long.
    fn eval(&self, words: &[u64]) -> bool;

    /// Value when the hit bitmap is all zeros. Known at compile time; lets the
    /// caller skip evaluation entirely at symbols that produced no hits.
    fn zero_value(&self) -> bool;
}

impl<P: Predicate + ?Sized> Predicate for Box<P> {
    #[inline]
    fn eval(&self, words: &[u64]) -> bool {
        (**self).eval(words)
    }
    #[inline]
    fn zero_value(&self) -> bool {
        (**self).zero_value()
    }
}

// ---------------------------------------------------------------------------
// Normalization: a single expansion shared by both normal forms
// ---------------------------------------------------------------------------

/// A count constraint with a normalized, deduplicated pattern list.
#[derive(Clone, Debug, PartialEq, Eq)]
struct CountSpec {
    lo: usize,
    hi: usize,
    patterns: Vec<PatternId>,
}

enum CountOutcome {
    AlwaysTrue,
    AlwaysFalse,
    Spec(CountSpec),
}

fn norm_count(lo: usize, hi: usize, patterns: &[PatternId]) -> CountOutcome {
    let mut ps = patterns.to_vec();
    ps.sort_unstable();
    ps.dedup();
    let n = ps.len();
    let hi = hi.min(n);
    if lo > hi {
        return CountOutcome::AlwaysFalse;
    }
    if lo == 0 && hi == n {
        return CountOutcome::AlwaysTrue;
    }
    CountOutcome::Spec(CountSpec { lo, hi, patterns: ps })
}

/// The complement of a non-trivial range: `[0, lo-1]` and `[hi+1, n]`, whichever
/// are non-empty. Never `AlwaysTrue`, since the input is never `AlwaysFalse`.
fn complement(spec: &CountSpec) -> Vec<CountSpec> {
    let n = spec.patterns.len();
    let mut out = Vec::with_capacity(2);
    if spec.lo > 0 {
        out.push(CountSpec {
            lo: 0,
            hi: spec.lo - 1,
            patterns: spec.patterns.clone(),
        });
    }
    if spec.hi < n {
        out.push(CountSpec {
            lo: spec.hi + 1,
            hi: n,
            patterns: spec.patterns.clone(),
        });
    }
    out
}

/// A conjunction of literals and count constraints.
///
/// In DNF this is read as written (everything ANDed). After dualization the very
/// same structure is read as a CNF clause (everything ORed) — see [`dualize`].
#[derive(Clone, Debug)]
struct Group {
    lits: BTreeMap<PatternId, bool>,
    counts: Vec<CountSpec>,
}

impl Group {
    /// The empty conjunction: TRUE in DNF, FALSE (empty clause) in CNF.
    fn empty() -> Group {
        Group {
            lits: BTreeMap::new(),
            counts: Vec::new(),
        }
    }

    fn size(&self) -> usize {
        self.lits.len() + self.counts.len()
    }

    /// Conjoin two groups. `None` if they contradict on a literal.
    fn merge(&self, other: &Group) -> Option<Group> {
        let mut lits = self.lits.clone();
        for (&p, &pol) in &other.lits {
            if let Some(prev) = lits.insert(p, pol) {
                if prev != pol {
                    return None; // p AND NOT p
                }
            }
        }
        let mut counts = self.counts.clone();
        for c in &other.counts {
            if !counts.contains(c) {
                counts.push(c.clone());
            }
        }
        Some(Group { lits, counts })
    }

    /// Does `self` subsume `other`? `(A) OR (A AND B) == A`.
    fn subsumes(&self, other: &Group) -> bool {
        self.lits.iter().all(|(p, pol)| other.lits.get(p) == Some(pol))
            && self.counts.iter().all(|c| other.counts.contains(c))
    }
}

/// Expand to disjunctive normal form, pushing negations down on the fly.
/// `neg` is the accumulated polarity of the enclosing context.
fn dnf(e: &Expr, neg: bool, budget: &mut usize) -> Result<Vec<Group>, CompileError> {
    Ok(match e {
        Expr::Const(b) => {
            if *b ^ neg {
                vec![Group::empty()] // TRUE
            } else {
                vec![] // FALSE: the empty disjunction
            }
        }
        Expr::Pattern(p) => {
            let mut g = Group::empty();
            g.lits.insert(*p, !neg);
            vec![g]
        }
        Expr::Not(x) => dnf(x, !neg, budget)?,

        // Conjunctive context: And in positive position, or (De Morgan) Or in
        // negative position. Cartesian product of the children's expansions.
        Expr::And(xs) if !neg => conjoin(xs, neg, budget)?,
        Expr::Or(xs) if neg => conjoin(xs, neg, budget)?,

        // Disjunctive context: concatenate.
        Expr::Or(xs) | Expr::And(xs) => {
            let mut acc = Vec::new();
            for x in xs {
                acc.extend(dnf(x, neg, budget)?);
            }
            spend(&acc, budget)?;
            acc
        }

        Expr::Count { lo, hi, patterns } => {
            let outcome = norm_count(*lo, *hi, patterns);
            match (outcome, neg) {
                (CountOutcome::AlwaysTrue, false) | (CountOutcome::AlwaysFalse, true) => {
                    vec![Group::empty()]
                }
                (CountOutcome::AlwaysFalse, false) | (CountOutcome::AlwaysTrue, true) => vec![],
                (CountOutcome::Spec(spec), false) => {
                    let mut g = Group::empty();
                    g.counts.push(spec);
                    vec![g]
                }
                (CountOutcome::Spec(spec), true) => {
                    // Complement is a union of ranges, so each becomes its own
                    // disjunct. No C(n, k) expansion.
                    complement(&spec)
                        .into_iter()
                        .map(|c| Group {
                            lits: BTreeMap::new(),
                            counts: vec![c],
                        })
                        .collect()
                }
            }
        }
    })
}

fn conjoin(xs: &[Expr], neg: bool, budget: &mut usize) -> Result<Vec<Group>, CompileError> {
    let mut acc = vec![Group::empty()];
    for x in xs {
        let d = dnf(x, neg, budget)?;
        let mut out = Vec::new();
        for a in &acc {
            for b in &d {
                if let Some(g) = a.merge(b) {
                    out.push(g);
                }
            }
        }
        spend(&out, budget)?;
        acc = out;
        if acc.is_empty() {
            break; // contradiction; nothing downstream can revive it
        }
    }
    Ok(acc)
}

fn spend(acc: &[Group], budget: &usize) -> Result<(), CompileError> {
    if acc.len() > *budget {
        Err(CompileError::TooComplex { budget: *budget })
    } else {
        Ok(())
    }
}

/// Drop any group subsumed by a more general one. Cheap at compile time, and
/// every group removed is work saved at every symbol. Valid in both forms:
/// `a | (a & b) == a` dualizes to `a & (a | b) == a`.
fn absorb(mut gs: Vec<Group>) -> Vec<Group> {
    gs.sort_by_key(|g| g.size());
    let mut kept: Vec<Group> = Vec::with_capacity(gs.len());
    'candidate: for g in gs {
        for k in &kept {
            if k.subsumes(&g) {
                continue 'candidate;
            }
        }
        kept.push(g);
    }
    kept
}

/// CNF by duality: if `DNF(!e) == OR_i AND_j l_ij` then `e == AND_i OR_j !l_ij`.
///
/// The literal negation is free at the mask level. A DNF group's `cmp` word
/// holds its *positive* literals; the dual clause's `cmp` must hold the bits
/// whose set-ness *fails* the clause, which are exactly the clause's negative
/// literals — i.e. the original group's positive literals. Same word, different
/// reading. Only the counts need real work, since `!(C1 & C2)` is
/// `!C1 | !C2` and each complement is a union of ranges.
fn dualize(groups: Vec<Group>) -> Vec<Group> {
    groups
        .into_iter()
        .map(|g| {
            let mut counts = Vec::new();
            for c in &g.counts {
                for d in complement(c) {
                    if !counts.contains(&d) {
                        counts.push(d);
                    }
                }
            }
            Group { lits: g.lits, counts }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Lowering to word masks
// ---------------------------------------------------------------------------

/// Matches iff `words[lo + i] & mask[i] == cmp[i]` for all i.
///
/// In DNF, `cmp` is the required polarity pattern and a match means the
/// conjunction of literals holds. In CNF, `cmp` is the all-literals-wrong
/// pattern and a match means the clause's literals are all unsatisfied.
#[derive(Clone, Debug)]
struct WordMask {
    lo: usize,
    mask: Box<[u64]>,
    cmp: Box<[u64]>,
}

impl WordMask {
    #[inline(always)]
    fn matches(&self, words: &[u64]) -> bool {
        // Single-word is the overwhelmingly common case: patterns are laid out
        // contiguously, so a group over related patterns rarely straddles a word
        // boundary. Special-cased to avoid the loop.
        if self.mask.len() == 1 {
            return words[self.lo] & self.mask[0] == self.cmp[0];
        }
        for i in 0..self.mask.len() {
            if words[self.lo + i] & self.mask[i] != self.cmp[i] {
                return false;
            }
        }
        true
    }
}

/// A count constraint lowered to a masked popcount.
#[derive(Clone, Debug)]
struct CountMask {
    lo: usize,
    mask: Box<[u64]>,
    min: u32,
    max: u32,
}

impl CountMask {
    #[inline(always)]
    fn matches(&self, words: &[u64]) -> bool {
        let mut n = 0u32;
        if self.mask.len() == 1 {
            n = (words[self.lo] & self.mask[0]).count_ones();
        } else {
            for i in 0..self.mask.len() {
                n += (words[self.lo + i] & self.mask[i]).count_ones();
            }
        }
        n >= self.min && n <= self.max
    }
}

#[derive(Clone, Debug)]
struct GroupMask {
    lits: WordMask,
    counts: Box<[CountMask]>,
}

impl GroupMask {
    /// DNF reading: everything ANDed.
    #[inline(always)]
    fn all_hold(&self, words: &[u64]) -> bool {
        self.lits.matches(words) && self.counts.iter().all(|c| c.matches(words))
    }

    /// CNF reading: everything ORed. The literal part is satisfied when the
    /// all-wrong pattern does *not* match.
    #[inline(always)]
    fn any_holds(&self, words: &[u64]) -> bool {
        !self.lits.matches(words) || self.counts.iter().any(|c| c.matches(words))
    }
}

// ---------------------------------------------------------------------------
// Backend 1a: sum of products
// ---------------------------------------------------------------------------

/// OR over conjunctive groups. Early-exits on the first group that holds.
#[derive(Clone, Debug)]
pub struct DnfPredicate {
    groups: Box<[GroupMask]>,
    zero: bool,
}

impl DnfPredicate {
    /// Number of product terms. Useful for logging what users' expressions cost.
    pub fn group_count(&self) -> usize {
        self.groups.len()
    }
}

impl Predicate for DnfPredicate {
    #[inline]
    fn eval(&self, words: &[u64]) -> bool {
        for g in self.groups.iter() {
            if g.all_hold(words) {
                return true;
            }
        }
        false
    }
    #[inline]
    fn zero_value(&self) -> bool {
        self.zero
    }
}

// ---------------------------------------------------------------------------
// Backend 1b: product of sums
// ---------------------------------------------------------------------------

/// AND over disjunctive clauses. Early-exits on the first clause that fails,
/// which is the common case once a symbol has any hits at all.
#[derive(Clone, Debug)]
pub struct CnfPredicate {
    groups: Box<[GroupMask]>,
    zero: bool,
}

impl CnfPredicate {
    /// Number of clauses.
    pub fn group_count(&self) -> usize {
        self.groups.len()
    }
}

impl Predicate for CnfPredicate {
    #[inline]
    fn eval(&self, words: &[u64]) -> bool {
        for g in self.groups.iter() {
            if !g.any_holds(words) {
                return false;
            }
        }
        true
    }
    #[inline]
    fn zero_value(&self) -> bool {
        self.zero
    }
}

// ---------------------------------------------------------------------------
// Compilation entry points
// ---------------------------------------------------------------------------

/// Check every pattern id up front so that normalization can only ever fail
/// with the soft `TooComplex`, which keeps the fallback logic honest.
fn validate(expr: &Expr, bit_of: &[usize], word_count: usize) -> Result<(), CompileError> {
    for p in expr.patterns() {
        let bit = *bit_of.get(p).ok_or(CompileError::UnknownPattern(p))?;
        if bit / 64 >= word_count {
            return Err(CompileError::BitOutOfRange { pattern: p, bit });
        }
    }
    Ok(())
}

fn lower_group(g: &Group, bit_of: &[usize]) -> GroupMask {
    let lits = if g.lits.is_empty() {
        WordMask {
            lo: 0,
            mask: Box::from(&[][..]),
            cmp: Box::from(&[][..]),
        }
    } else {
        let bits: Vec<usize> = g.lits.keys().map(|&p| bit_of[p]).collect();
        let lo = bits.iter().map(|&b| b / 64).min().unwrap();
        let hi = bits.iter().map(|&b| b / 64).max().unwrap();
        let mut mask = vec![0u64; hi - lo + 1];
        let mut cmp = vec![0u64; hi - lo + 1];
        for ((_, &pol), &bit) in g.lits.iter().zip(bits.iter()) {
            let w = bit / 64 - lo;
            let m = 1u64 << (bit % 64);
            mask[w] |= m;
            if pol {
                cmp[w] |= m;
            }
        }
        WordMask {
            lo,
            mask: mask.into_boxed_slice(),
            cmp: cmp.into_boxed_slice(),
        }
    };

    let counts = g
        .counts
        .iter()
        .map(|c| {
            let bits: Vec<usize> = c.patterns.iter().map(|&p| bit_of[p]).collect();
            let lo = bits.iter().map(|&b| b / 64).min().unwrap();
            let hi = bits.iter().map(|&b| b / 64).max().unwrap();
            let mut mask = vec![0u64; hi - lo + 1];
            for &bit in &bits {
                mask[bit / 64 - lo] |= 1u64 << (bit % 64);
            }
            CountMask {
                lo,
                mask: mask.into_boxed_slice(),
                min: c.lo as u32,
                max: c.hi as u32,
            }
        })
        .collect::<Vec<_>>()
        .into_boxed_slice();

    GroupMask { lits, counts }
}

/// Compile to sum-of-products.
pub fn compile_dnf(
    expr: &Expr,
    bit_of: &[usize],
    word_count: usize,
    budget: usize,
) -> Result<DnfPredicate, CompileError> {
    validate(expr, bit_of, word_count)?;
    let mut remaining = budget;
    let groups = absorb(dnf(expr, false, &mut remaining)?);
    let groups: Box<[GroupMask]> = groups
        .iter()
        .map(|g| lower_group(g, bit_of))
        .collect();
    let mut p = DnfPredicate { groups, zero: false };
    p.zero = p.eval(&vec![0u64; word_count]);
    Ok(p)
}

/// Compile to product-of-sums.
pub fn compile_cnf(
    expr: &Expr,
    bit_of: &[usize],
    word_count: usize,
    budget: usize,
) -> Result<CnfPredicate, CompileError> {
    validate(expr, bit_of, word_count)?;
    let mut remaining = budget;
    // DNF of the negation, then read every group as a clause.
    let groups = absorb(dualize(dnf(expr, true, &mut remaining)?));
    let groups: Box<[GroupMask]> = groups
        .iter()
        .map(|g| lower_group(g, bit_of))
        .collect();
    let mut p = CnfPredicate { groups, zero: false };
    p.zero = p.eval(&vec![0u64; word_count]);
    Ok(p)
}

// ---------------------------------------------------------------------------
// Backend 2: tree walk (fallback, and the reference for differential tests)
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum Node {
    Const(bool),
    Bit { word: usize, mask: u64 },
    Not(Box<Node>),
    And(Box<[Node]>),
    Or(Box<[Node]>),
    Count(CountMask),
}

/// Walks a lowered expression tree. Slower than the normal-form backends but
/// never exceeds the budget, so it is the fallback for a pathological
/// expression.
#[derive(Clone, Debug)]
pub struct TreePredicate {
    root: Node,
    zero: bool,
}

pub fn compile_tree(
    expr: &Expr,
    bit_of: &[usize],
    word_count: usize,
) -> Result<TreePredicate, CompileError> {
    validate(expr, bit_of, word_count)?;
    let root = lower_node(expr, bit_of);
    let zero = eval_node(&root, &vec![0u64; word_count]);
    Ok(TreePredicate { root, zero })
}

fn lower_node(e: &Expr, bit_of: &[usize]) -> Node {
    match e {
        Expr::Const(b) => Node::Const(*b),
        Expr::Pattern(p) => {
            let bit = bit_of[*p];
            Node::Bit {
                word: bit / 64,
                mask: 1u64 << (bit % 64),
            }
        }
        Expr::Not(x) => Node::Not(Box::new(lower_node(x, bit_of))),
        Expr::And(xs) => Node::And(xs.iter().map(|x| lower_node(x, bit_of)).collect()),
        Expr::Or(xs) => Node::Or(xs.iter().map(|x| lower_node(x, bit_of)).collect()),
        Expr::Count { lo, hi, patterns } => match norm_count(*lo, *hi, patterns) {
            CountOutcome::AlwaysTrue => Node::Const(true),
            CountOutcome::AlwaysFalse => Node::Const(false),
            CountOutcome::Spec(spec) => {
                let g = Group {
                    lits: BTreeMap::new(),
                    counts: vec![spec],
                };
                let gm = lower_group(&g, bit_of);
                Node::Count(gm.counts[0].clone())
            }
        },
    }
}

fn eval_node(n: &Node, words: &[u64]) -> bool {
    match n {
        Node::Const(b) => *b,
        Node::Bit { word, mask } => words[*word] & mask != 0,
        Node::Not(x) => !eval_node(x, words),
        Node::And(xs) => xs.iter().all(|x| eval_node(x, words)),
        Node::Or(xs) => xs.iter().any(|x| eval_node(x, words)),
        Node::Count(c) => c.matches(words),
    }
}

impl Predicate for TreePredicate {
    #[inline]
    fn eval(&self, words: &[u64]) -> bool {
        eval_node(&self.root, words)
    }
    #[inline]
    fn zero_value(&self) -> bool {
        self.zero
    }
}

// ---------------------------------------------------------------------------
// Backend selection + multi-expression driver
// ---------------------------------------------------------------------------

/// Static dispatch over the available backends. Add a variant to add a backend
/// without touching call sites.
#[derive(Clone, Debug)]
pub enum AnyPredicate {
    Dnf(DnfPredicate),
    Cnf(CnfPredicate),
    Tree(TreePredicate),
}

impl AnyPredicate {
    /// Which form was chosen, and how many groups it needed. Worth logging once
    /// per run: it tells you whether the budget is set sensibly and whether the
    /// fallback ever fires on real user expressions.
    pub fn shape(&self) -> (&'static str, usize) {
        match self {
            AnyPredicate::Dnf(p) => ("dnf", p.group_count()),
            AnyPredicate::Cnf(p) => ("cnf", p.group_count()),
            AnyPredicate::Tree(_) => ("tree", 0),
        }
    }
}

impl Predicate for AnyPredicate {
    #[inline]
    fn eval(&self, words: &[u64]) -> bool {
        match self {
            AnyPredicate::Dnf(p) => p.eval(words),
            AnyPredicate::Cnf(p) => p.eval(words),
            AnyPredicate::Tree(p) => p.eval(words),
        }
    }
    #[inline]
    fn zero_value(&self) -> bool {
        match self {
            AnyPredicate::Dnf(p) => p.zero_value(),
            AnyPredicate::Cnf(p) => p.zero_value(),
            AnyPredicate::Tree(p) => p.zero_value(),
        }
    }
}

/// Expand both normal forms, keep the smaller, fall back to a tree walk if both
/// blow the budget. Ties go to CNF: its early exit fires on the *false* result,
/// which is the common one at symbols that have hits but do not satisfy this
/// particular expression.
pub fn compile(
    expr: &Expr,
    bit_of: &[usize],
    word_count: usize,
    budget: usize,
) -> Result<AnyPredicate, CompileError> {
    validate(expr, bit_of, word_count)?;

    let d = match compile_dnf(expr, bit_of, word_count, budget) {
        Ok(p) => Some(p),
        Err(CompileError::TooComplex { .. }) => None,
        Err(e) => return Err(e),
    };
    let c = match compile_cnf(expr, bit_of, word_count, budget) {
        Ok(p) => Some(p),
        Err(CompileError::TooComplex { .. }) => None,
        Err(e) => return Err(e),
    };

    Ok(match (d, c) {
        (Some(d), Some(c)) => {
            if d.group_count() < c.group_count() {
                AnyPredicate::Dnf(d)
            } else {
                AnyPredicate::Cnf(c)
            }
        }
        (Some(d), None) => AnyPredicate::Dnf(d),
        (None, Some(c)) => AnyPredicate::Cnf(c),
        (None, None) => AnyPredicate::Tree(compile_tree(expr, bit_of, word_count)?),
    })
}


// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Deliberately straddles a word boundary: patterns 2..5 live in word 1.
    const BITS: [usize; 6] = [3, 10, 70, 71, 120, 127];
    const WORDS: usize = 2;
    const N: usize = BITS.len();

    fn words_from(set: &[bool]) -> Vec<u64> {
        let mut w = vec![0u64; WORDS];
        for (i, &on) in set.iter().enumerate() {
            if on {
                w[BITS[i] / 64] |= 1u64 << (BITS[i] % 64);
            }
        }
        w
    }

    /// Reference semantics, evaluated directly on a bool vector.
    fn reference(e: &Expr, set: &[bool]) -> bool {
        match e {
            Expr::Const(b) => *b,
            Expr::Pattern(p) => set[*p],
            Expr::Not(x) => !reference(x, set),
            Expr::And(xs) => xs.iter().all(|x| reference(x, set)),
            Expr::Or(xs) => xs.iter().any(|x| reference(x, set)),
            Expr::Count { lo, hi, patterns } => {
                let mut ps = patterns.clone();
                ps.sort_unstable();
                ps.dedup();
                let n = ps.iter().filter(|&&p| set[p]).count();
                n >= *lo && n <= *hi
            }
        }
    }

    /// Differential test: DNF, CNF and the tree walker must agree with the
    /// reference on every one of the 2^N possible hit bitmaps.
    fn exhaustive(e: &Expr) {
        let d = compile_dnf(e, &BITS, WORDS, 100_000).expect("dnf");
        let c = compile_cnf(e, &BITS, WORDS, 100_000).expect("cnf");
        let t = compile_tree(e, &BITS, WORDS).expect("tree");
        for m in 0..(1u32 << N) {
            let set: Vec<bool> = (0..N).map(|i| m >> i & 1 == 1).collect();
            let w = words_from(&set);
            let want = reference(e, &set);
            assert_eq!(d.eval(&w), want, "dnf mismatch on {set:?} for {e:?}");
            assert_eq!(c.eval(&w), want, "cnf mismatch on {set:?} for {e:?}");
            assert_eq!(t.eval(&w), want, "tree mismatch on {set:?} for {e:?}");
        }
        let zeros = vec![0u64; WORDS];
        assert_eq!(d.zero_value(), d.eval(&zeros));
        assert_eq!(c.zero_value(), c.eval(&zeros));
        assert_eq!(t.zero_value(), t.eval(&zeros));
    }

    #[test]
    fn simple_conjunction() {
        exhaustive(&Expr::all([Expr::pat(0), Expr::pat(1)]));
    }

    #[test]
    fn negation_and_word_straddle() {
        exhaustive(&Expr::any([
            Expr::all([Expr::pat(0), Expr::not(Expr::pat(1))]),
            Expr::pat(2),
        ]));
        exhaustive(&Expr::all([Expr::pat(1), Expr::pat(2), Expr::pat(4)]));
    }

    #[test]
    fn de_morgan() {
        exhaustive(&Expr::not(Expr::any([
            Expr::all([Expr::pat(0), Expr::pat(3)]),
            Expr::not(Expr::pat(2)),
        ])));
    }

    #[test]
    fn xor_like() {
        exhaustive(&Expr::all([
            Expr::any([Expr::pat(1), Expr::pat(2)]),
            Expr::not(Expr::all([Expr::pat(1), Expr::pat(2)])),
        ]));
    }

    #[test]
    fn constants_and_contradictions() {
        exhaustive(&Expr::all([Expr::pat(0), Expr::not(Expr::pat(0))]));
        exhaustive(&Expr::any([Expr::pat(0), Expr::not(Expr::pat(0))]));
        exhaustive(&Expr::all([Expr::Const(true), Expr::pat(3)]));
        exhaustive(&Expr::all([Expr::Const(false), Expr::pat(3)]));
        exhaustive(&Expr::Const(true));
        exhaustive(&Expr::Const(false));
    }

    #[test]
    fn count_ranges() {
        let all: Vec<PatternId> = (0..N).collect();
        exhaustive(&Expr::at_least(2, all.clone()));
        exhaustive(&Expr::at_most(2, all.clone()));
        exhaustive(&Expr::exactly(3, all.clone()));
        exhaustive(&Expr::Count { lo: 1, hi: 3, patterns: all.clone() });
        exhaustive(&Expr::at_least(0, all.clone())); // trivially true
        exhaustive(&Expr::at_least(N + 1, all.clone())); // trivially false
        exhaustive(&Expr::at_least(2, vec![1, 2, 2, 1])); // duplicates
    }

    #[test]
    fn negated_counts_do_not_expand() {
        let all: Vec<PatternId> = (0..N).collect();
        exhaustive(&Expr::not(Expr::at_least(2, all.clone())));
        exhaustive(&Expr::not(Expr::exactly(2, all.clone())));
        exhaustive(&Expr::not(Expr::Count { lo: 1, hi: 3, patterns: all.clone() }));

        // "At most 2" used to fall back via C(n, k) expansion. It must now be a
        // single group in both forms.
        let e = Expr::at_most(2, all);
        assert_eq!(compile_dnf(&e, &BITS, WORDS, 4).unwrap().group_count(), 1);
        assert_eq!(compile_cnf(&e, &BITS, WORDS, 4).unwrap().group_count(), 1);
    }

    #[test]
    fn counts_inside_conjunctions() {
        // The case that previously forced a fallback: a count ANDed with
        // structural conditions. "Heavily converted island AND motif present."
        let e = Expr::all([
            Expr::at_least(3, vec![0, 1, 2, 3]),
            Expr::pat(4),
            Expr::not(Expr::pat(0)),
        ]);
        exhaustive(&e);
        assert!(matches!(
            compile(&e, &BITS, WORDS, 1024).unwrap(),
            AnyPredicate::Dnf(_) | AnyPredicate::Cnf(_)
        ));

        exhaustive(&Expr::any([
            Expr::at_most(1, vec![0, 1, 2]),
            Expr::all([Expr::pat(3), Expr::pat(4)]),
        ]));
        exhaustive(&Expr::all([
            Expr::at_least(2, vec![0, 1, 2]),
            Expr::at_most(1, vec![2, 3, 4]),
        ]));
    }

    #[test]
    fn absorption_shrinks_groups() {
        let e = Expr::any([Expr::pat(0), Expr::all([Expr::pat(0), Expr::pat(1)])]);
        assert_eq!(compile_dnf(&e, &BITS, WORDS, 4096).unwrap().group_count(), 1);
        exhaustive(&e);
    }

    #[test]
    fn form_selection_picks_the_smaller() {
        let pairs = [[0usize, 1], [2, 3], [4, 5]];

        // Product of sums over disjoint pairs: 3 clauses, but 8 product terms.
        let pos = Expr::all(pairs.iter().map(|p| Expr::any([Expr::pat(p[0]), Expr::pat(p[1])])));
        assert_eq!(compile_dnf(&pos, &BITS, WORDS, 4096).unwrap().group_count(), 8);
        assert_eq!(compile_cnf(&pos, &BITS, WORDS, 4096).unwrap().group_count(), 3);
        assert_eq!(compile(&pos, &BITS, WORDS, 4096).unwrap().shape(), ("cnf", 3));
        exhaustive(&pos);

        // Sum of products: the exact mirror image.
        let sop = Expr::any(pairs.iter().map(|p| Expr::all([Expr::pat(p[0]), Expr::pat(p[1])])));
        assert_eq!(compile_dnf(&sop, &BITS, WORDS, 4096).unwrap().group_count(), 3);
        assert_eq!(compile_cnf(&sop, &BITS, WORDS, 4096).unwrap().group_count(), 8);
        assert_eq!(compile(&sop, &BITS, WORDS, 4096).unwrap().shape(), ("dnf", 3));
        exhaustive(&sop);

        // A budget that only one form fits still compiles, via that form.
        assert_eq!(compile(&pos, &BITS, WORDS, 4).unwrap().shape(), ("cnf", 3));
        assert_eq!(compile(&sop, &BITS, WORDS, 4).unwrap().shape(), ("dnf", 3));
    }

    #[test]
    fn parity_blows_up_in_both_forms_and_falls_back() {
        // XOR over all N variables: 2^(N-1) groups either way.
        let parity = (0..N).fold(Expr::Const(false), |acc, i| {
            Expr::any([
                Expr::all([acc.clone(), Expr::not(Expr::pat(i))]),
                Expr::all([Expr::not(acc), Expr::pat(i)]),
            ])
        });
        let p = compile(&parity, &BITS, WORDS, 8).unwrap();
        assert_eq!(p.shape().0, "tree");
        exhaustive(&parity);

        // With a generous budget it compiles, but stays large in both forms.
        let half = 1 << (N - 1);
        assert_eq!(compile_dnf(&parity, &BITS, WORDS, 100_000).unwrap().group_count(), half);
        assert_eq!(compile_cnf(&parity, &BITS, WORDS, 100_000).unwrap().group_count(), half);
    }

    #[test]
    fn errors_are_hard_not_soft() {
        assert_eq!(
            compile(&Expr::pat(99), &BITS, WORDS, 16).unwrap_err(),
            CompileError::UnknownPattern(99)
        );
        // bit 120 needs word 1; a one-word bitmap cannot hold it.
        assert_eq!(
            compile(&Expr::pat(4), &BITS, 1, 16).unwrap_err(),
            CompileError::BitOutOfRange { pattern: 4, bit: 120 }
        );
    }


    /// Randomized differential test over generated expressions.
    #[test]
    fn fuzz_all_backends_agree() {
        let mut seed = 0x243f_6a88_85a3_08d3u64;
        let mut rand = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        fn build(depth: usize, r: &mut impl FnMut() -> u64) -> Expr {
            if depth == 0 {
                return match r() % 8 {
                    0 => Expr::Const(r() % 2 == 0),
                    1..=2 => {
                        let k = (r() % 4) as usize;
                        let hi = k + (r() % 3) as usize;
                        let ps: Vec<usize> =
                            (0..N).filter(|_| r() % 2 == 0).collect();
                        if ps.is_empty() {
                            Expr::pat((r() % N as u64) as usize)
                        } else {
                            Expr::Count { lo: k, hi, patterns: ps }
                        }
                    }
                    _ => Expr::pat((r() % N as u64) as usize),
                };
            }
            match r() % 4 {
                0 => Expr::not(build(depth - 1, r)),
                1 => Expr::all((0..2 + r() % 2).map(|_| build(depth - 1, r)).collect::<Vec<_>>()),
                2 => Expr::any((0..2 + r() % 2).map(|_| build(depth - 1, r)).collect::<Vec<_>>()),
                _ => build(depth - 1, r),
            }
        }
        for _ in 0..300 {
            let e = build(3, &mut rand);
            exhaustive(&e);
        }
    }
}
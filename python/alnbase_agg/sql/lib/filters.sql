-- alnbase-agg: measures, relations, and the filter library.
--
--     .read lib/filters.sql
--
-- Read this after your script has defined `src`, `measure` and `relation`,
-- and before it defines `resolved`. DuckDB binds a macro's tables when a macro is
-- created, so reading it too early names a view you have not written yet rather
-- than failing later on a filter that silently matched nothing.
--
-- WHAT YOUR SCRIPT OWES THIS FILE
--
--   src        one row per hit, with at least:
--                read      what a read-level filter judges. Any expression
--                          that identifies one molecule; {'q': qname, 'm':
--                          mate} is the usual one. Never record_id: a chimeric
--                          Methyl-HiC read is several alignment records.
--                block     the partition over reads that the outputs are
--                          split by: a cell barcode, a spatial spot, or
--                          'bulk' when there is nothing to split by. Not
--                          called 'cell': the same column is a cell in one
--                          assay and a spot or a library in another, and
--                          every level below works the same way whichever
--                          you meant. It need not be one column -- a
--                          condition and a barcode joined together is one
--                          block -- but it must be a single value, because
--                          it is written out as the barcode of a cellInfo
--                          file and as a dataset name inside an H5.
--                name      the alnbase query that fired on this row.
--                contig    reference sequence name.
--                refr_pos  position, in the base the hit table declares.
--                strand    the strand the call sits on: '+', '-' or '.'.
--                qname     the read name. Needed only by read_fails_listed.
--
--   measure    (measure, name, category)
--   relation   (relation, measure, main, other)
--
-- MEASURES
--
-- A measure sorts alnbase query names into named buckets called categories.
-- It is the only place the biology lives.
--
--     ('mCG', 'CG',  'methylated')    a methylated CG
--     ('mCG', 'TG',  'unmethylated')  the same position, seen converted
--
--     ('gt',  'A>A', 'A')             a genotype measure: four categories and
--     ('gt',  'A>C', 'C')             no privileged one. A query name encodes
--     ('gt',  'A>G', 'G')             the reference base as well as the read
--     ('gt',  'A>T', 'T')             base, so "is this the reference allele"
--                                     needs no per-site lookup.
--
-- A query name belongs to at most one category *of a given measure*, and to no
-- category of the measures it has nothing to do with. Declaring two measures
-- over disjoint sets of query names is normal and is how an assay with two
-- independent readouts is handled: a NOMe-seq run declares an accessibility
-- measure over its GC queries and a methylation measure over its CG queries,
-- and no query name appears in both.
--
-- Within one measure the categories partition the hits that join it, which is
-- what lets a denominator be a plain sum with no double counting.
--
-- RELATIONS
--
-- A category count is a number on its own. A relation puts two sets of
-- categories against each other, and everything a filter or a file format
-- reads is an interpretation of one:
--
--     ('mCG', 'mCG', ['methylated'], ['unmethylated'])
--     ('G_vs_A', 'snp', ['G'], ['A'])
--     ('G_vs_rest', 'snp', ['G'], ['A','C','T'])
--
-- `main` is the side in focus and `other` is what it is set against.
-- The two are disjoint and neither repeats the other, exactly as the
-- categories they are built from do not overlap. What gets read off a relation
-- is derived:
--
--     main    hits on the side in focus
--     other   hits on the side it is set against
--     total   both together, which is what a format writes as `den`
--     ratio   main / total
--
-- So `ratio` is a way of reading a relation rather than a thing you declare,
-- which is why a run can ask for a count and a total off the same declaration
-- without writing it twice. A measure with exactly two categories has one
-- possible relation and normally declares it under the measure's own name,
-- which is why the common case never mentions relations at all.
--
-- WHAT THIS FILE GIVES BACK
--
--   read_category, block_category   one row per category: how many hits
--   read_relation, block_relation   num and den, per declared relation
--   block_extent                    reads and distinct positions, per block
--   read_fails*(...)                the reads that fail a test
--   block_fails*(...)               the blocks that fail a test
--
-- Every filter returns what FAILS it, so filters compose by union and a run
-- reads as one `WHERE ... NOT IN (...)`. The alternative -- each filter
-- returning what passes -- composes by intersection, which reads worse and
-- makes "rejected by nothing" indistinguishable from "judged by nothing".


-- A mistyped measure, category, relation or comparison is the one mistake here
-- that is otherwise silent: a filter naming something nothing declares matches
-- no rows, rejects nobody and says nothing. It cannot be caught inside a query,
-- because DuckDB's optimizer prunes an unreferenced `error()` call out of the
-- plan. Deriving a type from the declarations themselves makes it a bind-time
-- error instead, with no hand-written validator anywhere:
-- `read_fails('mCG_typo', 'ratio', '>', 0.5)` fails with
--     Conversion Error: Could not convert string 'mCG_typo' to UINT8
DROP TYPE IF EXISTS measure_name;
CREATE TYPE measure_name AS ENUM (SELECT DISTINCT measure FROM measure);
DROP TYPE IF EXISTS category_name;
CREATE TYPE category_name AS ENUM (SELECT DISTINCT category FROM measure);
DROP TYPE IF EXISTS relation_name;
CREATE TYPE relation_name AS ENUM (SELECT DISTINCT relation FROM relation);
DROP TYPE IF EXISTS comparison;
CREATE TYPE comparison AS ENUM ('<', '<=', '>', '>=', '=', '<>');
DROP TYPE IF EXISTS quantity;
CREATE TYPE quantity AS ENUM ('main', 'other', 'total', 'ratio');


-- Every filter with a threshold takes the comparison as an argument, so the
-- call says which way the test points and whether the boundary is included.
--
-- Spelling it out is not pedantry. `> 0.4` and `>= 0.4` are different filters,
-- and premethyst's read cut uses the second; with the direction baked into a
-- macro name there was no way to ask for it, and a run needing it had to write
-- the SQL out by hand. One macro per quantity, with the comparison visible at
-- the call site, is both shorter to document and more expressive.
--
-- What the operators mean lives here and nowhere else, so no filter below can
-- disagree with another about a boundary.
CREATE OR REPLACE MACRO compare(x, op, t) AS
  CASE op::comparison::VARCHAR
    WHEN '<'  THEN x <  t
    WHEN '<=' THEN x <= t
    WHEN '>'  THEN x >  t
    WHEN '>=' THEN x >= t
    WHEN '='  THEN x =  t
    ELSE           x <> t
  END;


-- ---------------------------------------------------------------------------
-- Counting: one row per category, then one row per relation.
-- ---------------------------------------------------------------------------

-- How many hits each read contributed to each category of each measure.
CREATE OR REPLACE VIEW read_category AS
SELECT s.read, m.measure, m.category, count(*) AS n
FROM src s JOIN measure m USING (name)
GROUP BY ALL;

-- The same one level up. Note that this reads `src` and not `resolved`: a
-- block is judged on everything that was sequenced from it, including reads a
-- read filter is about to throw away. That is the right order for a QC cut on
-- depth -- a block is thin because little was sequenced, not because a filter
-- was strict -- and the wrong order if you meant "blocks with too few
-- surviving reads". For that, define a second view over `resolved` in your own
-- script; it is the same four lines.
CREATE OR REPLACE VIEW block_category AS
SELECT s.block, m.measure, m.category, count(*) AS n
FROM src s JOIN measure m USING (name)
GROUP BY ALL;

-- Each relation, read off per read: `num` is its first side and `den` the two
-- sides together, which is the pair every file format writes and every
-- proportion is taken over.
--
-- `coalesce(..., 0)` is load bearing: a read that covered ten CH sites and
-- methylated none of them has a methylation ratio of 0, not no ratio, and a
-- sum over an empty FILTER is NULL rather than zero.
CREATE OR REPLACE VIEW read_relation AS
SELECT c.read,
       r.relation,
       coalesce(sum(c.n) FILTER (WHERE list_contains(r.main, c.category)), 0) AS num,
       coalesce(sum(c.n) FILTER (WHERE list_contains(r.main, c.category)
                                    OR list_contains(r.other, c.category)), 0) AS den
FROM read_category c JOIN relation r USING (measure)
GROUP BY ALL;

CREATE OR REPLACE VIEW block_relation AS
SELECT c.block,
       r.relation,
       coalesce(sum(c.n) FILTER (WHERE list_contains(r.main, c.category)), 0) AS num,
       coalesce(sum(c.n) FILTER (WHERE list_contains(r.main, c.category)
                                    OR list_contains(r.other, c.category)), 0) AS den
FROM block_category c JOIN relation r USING (measure)
GROUP BY ALL;



-- What a relation yields, in one place. A filter is a grouping, a quantity and
-- a fail condition, and this is the middle one.
-- The first two are named for the relation's own two sides, so there is
-- nothing to learn: you declared `main` and `other`, and `main` and `other`
-- are how many hits fell on each. `total` is the two together, and `ratio` is
-- the side in focus over that total.
--
-- Note the change of vocabulary between layers, which is deliberate. A
-- relation is declared and filtered in relation words -- main, other -- while
-- a file format writes a ratio, and `num` and `den` are the right words there.
-- The parameters below are the counts the views produce, not the declaration.
--
-- Every quantity is named, including the last. Falling through to the ratio
-- would mean that adding a fifth quantity to the ENUM above and forgetting a
-- branch here silently computed a ratio instead -- which is the shape of
-- failure this library exists to refuse, so it raises rather than guesses.
CREATE OR REPLACE MACRO quantity_value(q, num, den) AS
  CASE q::quantity::VARCHAR
    WHEN 'main'  THEN num::DOUBLE
    WHEN 'other' THEN (den - num)::DOUBLE
    WHEN 'total' THEN den::DOUBLE
    WHEN 'ratio' THEN CASE WHEN den > 0 THEN num::DOUBLE / den END
    ELSE error('quantity_value: no branch for quantity ' || q)
  END;


-- How much of a relation each block saw, counted in reads and in distinct
-- positions rather than in hits.
--
-- Its own view because neither number can be derived from `num` and `den`: a
-- block with 2000 hits may be one read seen 2000 times or 2000 reads seen
-- once, and every single-cell pipeline cuts on the second. premethyst's
-- `unique_reads >= 10000` and its `calls-filter -G` minimum CG *sites* are the
-- two commonest block cuts in this ecosystem, and neither is a hit count.
CREATE OR REPLACE VIEW block_extent AS
SELECT s.block,
       r.relation,
       count(DISTINCT s.read) AS reads,
       count(DISTINCT (s.contig, s.refr_pos, s.strand)) AS sites
FROM src s
JOIN measure m USING (name)
JOIN relation r ON r.measure = m.measure
WHERE list_contains(r.main, m.category) OR list_contains(r.other, m.category)
GROUP BY ALL;


-- ---------------------------------------------------------------------------
-- The filters.
--
-- A filter is three choices, and each is written where it belongs:
--
--     what it groups on     the macro: read_fails, block_fails, site_fails
--     what it computes      the quantity: 'main', 'other', 'total', 'ratio'
--     when it fails         the comparison and the threshold
--
-- so "blocks with fewer than 1000 hits" is `block_fails('mCG', 'total', '<',
-- 1000)` and needs no macro of its own. Every filter returns what FAILS it, so
-- filters compose by union and a run reads as one `WHERE ... NOT IN (...)`.
--
-- Below these are the quantities that are not a plain aggregate over one
-- group's own rows -- a run along the read, a count of reads or of covered
-- positions rather than of hits, the reads a file names -- which keep a macro
-- each because they cannot be written as one. The site-level members of that
-- set are in lib/count_sites.sql, because they need `site_all` to exist.
-- ---------------------------------------------------------------------------

CREATE OR REPLACE MACRO read_fails(rel, q, op, x) AS TABLE
  SELECT read FROM read_relation
  WHERE relation = rel::relation_name::VARCHAR
    AND compare(quantity_value(q, num, den), op, x);

CREATE OR REPLACE MACRO block_fails(rel, q, op, x) AS TABLE
  SELECT block FROM block_relation
  WHERE relation = rel::relation_name::VARCHAR
    AND compare(quantity_value(q, num, den), op, x);

-- How many reads the block contributed, and how many distinct positions it
-- covered. Their own macros because neither is an aggregate over the block's
-- own num and den; see `block_extent` above.
CREATE OR REPLACE MACRO block_fails_reads(rel, op, n) AS TABLE
  SELECT block FROM block_extent
  WHERE relation = rel::relation_name::VARCHAR AND compare(reads, op, n);

CREATE OR REPLACE MACRO block_fails_sites(rel, op, n) AS TABLE
  SELECT block FROM block_extent
  WHERE relation = rel::relation_name::VARCHAR AND compare(sites, op, n);


-- The longest run of consecutive main-side calls along the read, in reference
-- order, counted among the positions the relation covers.
--
-- Its own macro because it is a window over positions rather than an aggregate
-- over one read's totals. It is also the test a ratio provably cannot express:
-- a read with four methylated CH in a row and a read with the same four
-- scattered have identical ratios, and only the first looks like a patch of
-- failed conversion. The shape is gaps-and-islands -- the difference between a
-- row number over all calls and a row number within each run of like calls is
-- constant exactly within a run -- so it costs two window functions and no
-- procedural code.
CREATE OR REPLACE MACRO read_fails_consecutive(rel, op, n) AS TABLE
WITH spec AS (
  SELECT measure, main, other FROM relation
  WHERE relation = rel::relation_name::VARCHAR
), calls AS (
  SELECT DISTINCT
         s.read,
         s.refr_pos,
         list_contains((SELECT main FROM spec), m.category) AS in_main
  FROM src s JOIN measure m USING (name)
  WHERE m.measure = (SELECT measure FROM spec)
    AND (list_contains((SELECT main FROM spec), m.category)
      OR list_contains((SELECT other FROM spec), m.category))
), islands AS (
  SELECT read,
         in_main,
         row_number() OVER (PARTITION BY read ORDER BY refr_pos)
       - row_number() OVER (PARTITION BY read, in_main ORDER BY refr_pos) AS island
  FROM calls
)
SELECT DISTINCT read FROM (
  SELECT read, in_main, island, count(*) AS run
  FROM islands
  GROUP BY read, in_main, island
) WHERE in_main AND compare(run, op, n);

-- No observations of the measure at all, which no threshold can reach.
--
-- The counting views above are sparse: a read gets a row because it had a hit,
-- so a read covering no CH site is ABSENT from read_relation rather than
-- present with a total of 0. `read_fails(rel, 'total', '<', 1)` therefore
-- binds, runs and matches nothing, ever. This macro reads `src` instead, where
-- every read is present whatever it covered.
--
-- Whether such reads should go is a real choice -- premethyst discards exactly
-- them -- so it is spelled, not assumed. The one read filter about a measure
-- rather than a relation.
CREATE OR REPLACE MACRO read_fails_missing(m) AS TABLE
  SELECT DISTINCT read FROM src
  WHERE read NOT IN (SELECT read FROM read_category
                     WHERE measure = m::measure_name::VARCHAR);

-- A list of read names in a file: a blacklist, a duplicate set, or the output
-- of another tool entirely. The file is read at query time, so regenerating it
-- and re-running changes the answer without touching the script.
CREATE OR REPLACE MACRO read_fails_listed(path) AS TABLE
  SELECT DISTINCT read FROM src
  WHERE qname IN (SELECT qname FROM read_csv(path, header = false,
                                             columns = {'qname': 'VARCHAR'}));

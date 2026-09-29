-- alnbase-agg: count the surviving hits into sites.
--
--     .read lib/count_sites.sql
--
-- REQUIRES resolved
--
-- Read this after your script has defined `resolved`, the view that says which
-- observations survived. Everything under `formats/` reads what this defines.
--
-- WHAT THIS FILE GIVES BACK
--
--   site_stranded   one row per (block, contig, position, strand, measure,
--                   category), with a count. One base, on its own strand.
--   site_category   what merging redefines. The same rows, until a script
--                   reads lib/merge_strands.sql.
--   site_all        the same rows pivoted against the declared relations: num
--                   and den per (block, position, strand, relation).
--   site_relation   what the formats read, and what dropping sites redefines.
--   block_summary   one row per (block, relation): num, den, and the count of
--                   distinct positions covered.
--   pos1, pos0      the coordinate convention, in one place.
--   site_key        what identifies one site to a filter.
--   site_fails*     the site-level filters, which need site_all to exist.
--
-- CHANGING THE DIMENSIONS
--
-- `GROUP BY ALL` means the grouping is exactly the non-aggregate columns in
-- the SELECT, so adding or dropping a dimension is adding or dropping a line.
-- A run that wants to count by the read's own strand as well, or by a capture
-- column the alnbase query wrote, copies this file next to its script and adds
-- the column; `.read lib/count_sites.sql` then finds the copy.


-- The coordinate convention, in one place.
--
-- A position off by one produces a file that writes cleanly, counts
-- plausibly, plots with the right shape and is wrong at every base, and
-- nothing downstream can catch it. So no format file below writes `+ 1`: they
-- call these, and the one arithmetic statement in the whole toolchain is here.
--
-- alnbase records `coordinate_base = 0` in every hit table's footer, which is
-- the default, so a run normally sets nothing. `alnbase-agg check` is what
-- confirms the positions really are where they claim to be, by looking the
-- reference base up in a FASTA.
CREATE OR REPLACE MACRO pos1(p) AS
  p + (1 - coalesce(getvariable('coordinate_base'), 0));   -- 1-based inclusive
CREATE OR REPLACE MACRO pos0(p) AS
  p - coalesce(getvariable('coordinate_base'), 0);         -- 0-based inclusive


-- One row per base per block per category. Long rather than pivoted, because
-- a measure may have any number of categories; the pivot into num and den
-- happens once, in `site_all` below.
--
-- A position whose only reads a filter rejected is absent here entirely,
-- rather than present with zero counts. An uncovered cytosine is not an
-- unmethylated cytosine, and every downstream format in this ecosystem is
-- read that way.
CREATE OR REPLACE VIEW site_stranded AS
SELECT r.block,
       r.contig,
       r.refr_pos,
       r.strand,
       m.measure,
       m.category,
       count(*) AS n
FROM resolved r JOIN measure m USING (name)
GROUP BY ALL;

-- The indirection is what lets strand merging be an optional include. A view
-- cannot be redefined in terms of itself, so `merge_strands.sql` needs
-- something underneath to redefine `site_category` from.
CREATE OR REPLACE VIEW site_category AS SELECT * FROM site_stranded;


-- What every format reads: one row per position per declared relation, with the
-- numerator and the denominator the ratio names. `den - num` is the negative
-- count every bed-like format wants.
--
-- `coalesce(..., 0)` because a position covered only by the denominator's
-- other category has no numerator rows at all, and a sum over an empty FILTER
-- is NULL. An unmethylated cytosine has num 0, not num NULL.
CREATE OR REPLACE VIEW site_all AS
SELECT c.block,
       c.contig,
       c.refr_pos,
       c.strand,
       r.relation,
       coalesce(sum(c.n) FILTER (WHERE list_contains(r.main, c.category)), 0)::BIGINT AS num,
       coalesce(sum(c.n) FILTER (WHERE list_contains(r.main, c.category)
                                    OR list_contains(r.other, c.category)), 0)::BIGINT AS den
FROM site_category c JOIN relation r USING (measure)
GROUP BY ALL;

-- The second indirection, and the same trick as `site_category` above: a run
-- that drops sites redefines `site_relation` in terms of `site_all`, and every
-- format follows without being passed anything.
CREATE OR REPLACE VIEW site_relation AS SELECT * FROM site_all;

-- What identifies one site, for a filter that rejects some. The struct is
-- built the same way everywhere, so a site filter composes with another
-- exactly as two read filters do.
CREATE OR REPLACE MACRO site_key(b, c, p, s) AS
  {'block': b, 'contig': c, 'pos': p, 'strand': s};


-- Per block, per relation, after filtering. `sites` counts rows of `site_relation`,
-- which is one row per position, so it is the distinct-position count a
-- cellInfo file means by "coverage" -- a different number from `den`, which
-- counts calls. A block covering one CG twenty times has den 20 and sites 1.
CREATE OR REPLACE VIEW block_summary AS
SELECT block,
       relation,
       sum(num)::BIGINT AS num,
       sum(den)::BIGINT AS den,
       count(*)::BIGINT AS sites
FROM site_relation
WHERE den > 0
GROUP BY ALL;


-- ---------------------------------------------------------------------------
-- Site-level filters. Each returns the sites that FAIL.
--
-- They live here rather than in lib/filters.sql because they read `site_all`,
-- which does not exist until this file has run; DuckDB binds a macro's tables
-- when the macro is created.
--
-- A site is one position in one block, so these judge a position per block:
-- "this cell did not see this CG often enough" rather than "nobody did". A
-- filter meant to apply everywhere returns every block's row for the position,
-- which is what four of the five below do.
--
-- This is not the `min_depth` a format takes. `min_depth` is what one file
-- shows and is set per output; these are what the run keeps, are recorded by
-- the run, and appear in its report beside every other filter.
-- ---------------------------------------------------------------------------

-- The same three choices as every other filter: the grouping is the site, the
-- quantity is 'main', 'other', 'total' or 'ratio', and the condition is a
-- comparison. `site_fails('mCG', 'total', '<', 5)` is the coverage cut.
CREATE OR REPLACE MACRO site_fails(rel, q, op, x) AS TABLE
  SELECT site_key(block, contig, refr_pos, strand) AS site
  FROM site_all
  WHERE relation = rel::relation_name::VARCHAR
    AND compare(quantity_value(q, num, den), op, x);

-- How many blocks covered the position at all. Its own macro because it is a
-- count across blocks rather than an aggregate over one site's own rows: it
-- judges the position once and then rejects it in every block.
CREATE OR REPLACE MACRO site_fails_blocks(rel, op, n) AS TABLE
  SELECT site_key(s.block, s.contig, s.refr_pos, s.strand) AS site
  FROM site_all s
  JOIN (SELECT contig, refr_pos, strand, count(DISTINCT block) AS blocks
        FROM site_all
        WHERE relation = rel::relation_name::VARCHAR AND den > 0
        GROUP BY contig, refr_pos, strand) b
    USING (contig, refr_pos, strand)
  WHERE s.relation = rel::relation_name::VARCHAR
    AND compare(b.blocks, op, n);

-- Positions inside any region of a BED file, and positions outside every one.
--
-- Both directions exist because both are asked for and they are not each
-- other's complement in any way a run can spell: `in_regions` is an ENCODE
-- blacklist or a list of known SNPs, `outside_regions` is a capture panel or
-- premethyst's `context-extract`, which keeps only the positions a site file
-- names. Without the second, restricting to a panel meant inverting the BED by
-- hand.
--
-- Three columns are read and the rest of the line ignored, and the coordinates
-- are BED's own -- 0-based, half open -- which is why `pos0` appears rather
-- than a bare comparison.
CREATE OR REPLACE MACRO bed_regions(path) AS TABLE
  SELECT * FROM read_csv(path, header = false, delim = E'\t',
                         columns = {'contig': 'VARCHAR', 'start': 'BIGINT',
                                    'end': 'BIGINT'},
                         null_padding = true);

CREATE OR REPLACE MACRO site_fails_in_regions(path) AS TABLE
  SELECT site_key(s.block, s.contig, s.refr_pos, s.strand) AS site
  FROM site_all s
  WHERE EXISTS (SELECT 1 FROM bed_regions(path) r
                WHERE r.contig = s.contig
                  AND pos0(s.refr_pos) >= r.start
                  AND pos0(s.refr_pos) <  r."end");

CREATE OR REPLACE MACRO site_fails_outside_regions(path) AS TABLE
  SELECT site_key(s.block, s.contig, s.refr_pos, s.strand) AS site
  FROM site_all s
  WHERE NOT EXISTS (SELECT 1 FROM bed_regions(path) r
                    WHERE r.contig = s.contig
                      AND pos0(s.refr_pos) >= r.start
                      AND pos0(s.refr_pos) <  r."end");


-- Every site at a position that fails a test, pooling the evidence across
-- blocks and strands and rejecting the position in all of them.
--
-- This is the SNP veto. A C read as a T may be bisulfite conversion or it may
-- be a C-to-T variant, and the two are indistinguishable on that strand. The
-- evidence is on the other one: an alnbase query anchored on the base opposite
-- the cytosine counts reference-matching and variant bases there, unaffected
-- by conversion. That evidence is a relation of its own at the same position
-- and the other strand, in whichever blocks happened to cover it -- so a
-- filter keyed to one block's own row cannot reach it, and this one is keyed
-- to the position instead.
--
-- It generalises past the SNP case: any "reject this position everywhere
-- because of what some other relation says about it" is this macro.
-- The evidence is POOLED before it is compared, not tested block by block.
-- Testing per block and rejecting if any one fails would veto a position
-- because a single cell saw a single variant read, which is the opposite of
-- what evidence across cells is for.
CREATE OR REPLACE MACRO site_fails_position(rel, q, op, x) AS TABLE
  SELECT site_key(s.block, s.contig, s.refr_pos, s.strand) AS site
  FROM site_all s
  JOIN (SELECT contig, refr_pos, sum(num) AS num, sum(den) AS den
        FROM site_all
        WHERE relation = rel::relation_name::VARCHAR
        GROUP BY contig, refr_pos) v
    USING (contig, refr_pos)
  WHERE compare(quantity_value(q, v.num, v.den), op, x);

-- A single-cell methylation run, end to end.
--
--     alnbase-agg run single_cell_methylation.sql \
--         -s hits='calls_*_*.parquet' \
--         -s fai=genome.fa.fai \
--         -s sample=pbmc \
--         -s exclude=exclude.txt
--
-- The hit tables come from an alnbase run that selected the record fields this
-- script names:
--
--     alnbase query --parquet --query-file queries.toml \
--         --query-file query/strand/bsbolt.toml \
--         -F qname,is_first_in_template,ref_name,strand,mapq,XB:Z \
--         reads.bam ref.aref calls
--
-- Read this top to bottom. It is the whole run: there is no second file that
-- says what these declarations mean, and nothing is filled in elsewhere.


-- 1. THE HITS, AND THE HIT-LEVEL FILTER ------------------------------------
--
-- `src` is the adapter between what alnbase wrote and what the library
-- expects. It is the only place that knows your column names, so a run against
-- a differently-selected hit table changes here and nowhere else.
--
-- The WHERE is the hit-level filter, and it runs before anything is counted.
-- This is where per-position artefacts go: end repair sits a fixed distance
-- from the 3' end and random-priming bias the same distance from the 5' end,
-- and both offsets count from the read as it was sequenced, including
-- hard-clipped bases.

CREATE OR REPLACE VIEW src AS
SELECT *,
       -- `read` is the read, not the alignment record: a chimeric Methyl-HiC
       -- read is several records, and judging them separately would filter
       -- half a molecule.
       ref_name AS contig,
       {'qname': qname, 'mate': is_first_in_template} AS read,
       XB AS block
FROM read_parquet(getvariable('hits'))
WHERE qual >= 20
  AND mapq >= 30
  AND read_5p >= 10
  AND read_3p >= 2;


-- 2. WHAT THE QUERY NAMES MEAN ---------------------------------------------
--
-- The only place in the run where biology appears. `name` is the alnbase query
-- that fired on a hit row; this sorts those names into disjoint categories.
-- Here `methylated` is the call and `unmethylated` the same position seen not
-- called, so the two together are the denominator and no name is written twice.
--
-- A measure is not limited to two categories -- a SNP measure has one per
-- base, and the filters take the category by name. What a file format writes
-- is two of them set against each other, because every format here puts one
-- ratio on a line.
--
-- Everything downstream is arithmetic over these three columns, so an assay
-- this toolchain has never heard of is a different set of rows here and no
-- change anywhere else.

CREATE OR REPLACE VIEW measure (measure, name, category) AS VALUES
  ('mCG', 'CG',  'methylated'), ('mCG', 'TG',  'unmethylated'),
  ('mCH', 'CHG', 'methylated'), ('mCH', 'CHH', 'methylated'),
  ('mCH', 'THG', 'unmethylated'), ('mCH', 'THH', 'unmethylated');


-- And what is put against what. Every filter and every file format reads an
-- interpretation of a relation: `main` is the side in focus, `other` what it
-- is set against, `total` is both together and `ratio` is `main` over `total`.
-- The sides are disjoint and neither repeats the other. A measure with two categories has one possible relation
-- and normally declares it under the measure's own name.

CREATE OR REPLACE VIEW relation (relation, measure, main, other) AS VALUES
  ('mCG', 'mCG', ['methylated'], ['unmethylated']),
  ('mCH', 'mCH', ['methylated'], ['unmethylated']);


.read lib/filters.sql


-- 3. WHICH OBSERVATIONS SURVIVE --------------------------------------------
--
-- Every level of filtering composes in one WHERE, because every filter returns
-- the things that fail it. Read it as a list of reasons to throw something
-- away.

CREATE OR REPLACE VIEW resolved AS
SELECT * FROM src
WHERE read NOT IN (
        -- Incomplete bisulfite conversion: a read with too much CH methylation
        -- did not convert, so its CG calls cannot be trusted either.
        SELECT read FROM read_fails('mCH', 'ratio', '>', 0.4)

        -- A patch of adjacent CH calls, which a rate cannot see: four in a row
        -- and four scattered give the same rate, and only the first looks like
        -- a conversion failure.
        UNION ALL SELECT read FROM read_fails_consecutive('mCH', '>', 3)

        -- Reads named in a file. Regenerating the file and re-running changes
        -- the answer without touching this script.
        UNION ALL SELECT read FROM read_fails_listed(getvariable('exclude'))

        -- A threshold that depends on how much evidence there is, which is
        -- how premethyst judges a read. The library ships no macro for it:
        -- here it is the CASE it is, and the first branch that matches
        -- decides. The chained front end says the same thing as two calls to
        -- `both(...)`, unioned like any other pair of reasons.
        UNION ALL
        SELECT read FROM read_relation
        WHERE relation = 'mCH'
          AND CASE WHEN den < 3  THEN false          -- too few sites to judge
                   WHEN den <= 5 THEN num > 1        -- one call is tolerable
                   ELSE num::DOUBLE / den >= 0.4 END)

  -- And one level up: blocks with too little data to be worth keeping.
  AND block NOT IN (SELECT block FROM block_fails('mCG', 'total', '<', 1000));


-- 4. WHAT THE RUN THREW AWAY -----------------------------------------------
--
-- A `.sql` file cannot be asked what filters it applied, so a run that wants
-- an account of itself prints one. These are ordinary queries; `alnbase-agg
-- run` shows the result of any statement that returns rows.

SELECT 'reads seen'     AS what, count(DISTINCT read) AS n FROM src
UNION ALL
SELECT 'reads kept',    count(DISTINCT read) FROM resolved
UNION ALL
SELECT 'high mCH',      count(*) FROM read_fails('mCH', 'ratio', '>', 0.4)
UNION ALL
SELECT 'mCH in a row',  count(*) FROM read_fails_consecutive('mCH', '>', 3)
UNION ALL
SELECT 'listed',        count(*) FROM read_fails_listed(getvariable('exclude'))
UNION ALL
SELECT 'thin blocks', count(*) FROM block_fails('mCG', 'total', '<', 1000);


-- 5. COUNTING ---------------------------------------------------------------

.read lib/count_sites.sql

-- Add the two cytosines of each CG together, onto the plus-strand one. A
-- choice, not a fact -- a run studying hemimethylation must not do it -- so it
-- is written down rather than assumed.
SET VARIABLE merge = 'mCG';
.read lib/merge_strands.sql

-- Contig order, for anything that sorts across contigs. `ORDER BY contig` is
-- lexical and would put chr10 between chr1 and chr2.
.read lib/load_contigs.sql


-- 6. THE FILES --------------------------------------------------------------
--
-- A format file is a parameterised include: set what it reads, then read it.
-- The variables are late bound, so the same file written once produces both
-- the CG and the CH output.

SET VARIABLE relation = 'mCG';
SET VARIABLE out = getvariable('sample') || '.CG.cov';
.read formats/bismark_cov.sql

SET VARIABLE relation = 'mCH';
SET VARIABLE out = getvariable('sample') || '.CH.cov';
.read formats/bismark_cov.sql

-- One line per block: the metadata table Amethyst builds its object from.
SET VARIABLE cg_relation = 'mCG';
SET VARIABLE ch_relation = 'mCH';
SET VARIABLE out = getvariable('sample') || '.cellInfo.txt';
.read formats/cellinfo.sql

-- And the per-block HDF5. The only output in this run that needs Python, and
-- only because the file is a dataset per block rather than one table. Read
-- twice, once per context, because an Amethyst file holds a group for each.
SET VARIABLE out = getvariable('sample') || '.h5';

SET VARIABLE relation = 'mCG';
SET VARIABLE context = 'CG';
.read formats/amethyst_h5.sql
.python formats/amethyst_h5.py

SET VARIABLE relation = 'mCH';
SET VARIABLE context = 'CH';
.read formats/amethyst_h5.sql
.python formats/amethyst_h5.py

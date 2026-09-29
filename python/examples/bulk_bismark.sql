-- The smallest useful run: bulk CG methylation, one coverage file.
--
--     alnbase-agg run bulk_bismark.sql -s hits='calls_*_*.parquet' -s out=bulk.CG.cov
--
-- Four declarations and four includes. Compare it with
-- `single_cell_methylation.sql`, which is the same shape with more in it:
-- nothing here is a cut-down mode, it is the same machinery with fewer rows in
-- `measure` and a shorter predicate in `resolved`.

-- The adapter. `block` is a constant because bulk data has nothing to split
-- by, and every level below still works: the block filters have one block to
-- judge.
CREATE OR REPLACE VIEW src AS
SELECT *,
       ref_name AS contig,
       {'qname': qname, 'mate': is_first_in_template} AS read,
       'bulk' AS block
FROM read_parquet(getvariable('hits'))
WHERE qual >= 20;

CREATE OR REPLACE VIEW measure (measure, name, category) AS VALUES
  ('mCG', 'CG', 'methylated'), ('mCG', 'TG', 'unmethylated');

CREATE OR REPLACE VIEW relation (relation, measure, main, other) AS VALUES
  ('mCG', 'mCG', ['methylated'], ['unmethylated']);

.read lib/filters.sql

-- No read filter at all, which is a position worth taking explicitly: `src`
-- and `resolved` differ by nothing, and a reader can see that at a glance.
CREATE OR REPLACE VIEW resolved AS SELECT * FROM src;

.read lib/count_sites.sql

SET VARIABLE merge = 'mCG';
.read lib/merge_strands.sql

SET VARIABLE relation = 'mCG';
.read formats/bismark_cov.sql

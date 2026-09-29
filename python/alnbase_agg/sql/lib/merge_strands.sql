-- alnbase-agg: fold a palindromic context onto one coordinate.
--
--     SET VARIABLE merge = 'mCG';
--     .read lib/count_sites.sql
--     .read lib/merge_strands.sql
--
-- REQUIRES lib/count_sites.sql
--
-- A palindromic context is two bases, one on each strand, at two different
-- coordinates. A CG is the common case: a cytosine on the plus strand at p,
-- and another on the minus strand at p+1, opposite the G. alnbase reports both,
-- because they are two observations of two different bases and throwing that
-- away is not a walk's business. Most downstream formats want them added
-- together, and that addition is this file.
--
-- It is a separate include because it is a choice, not a fact. A run studying
-- strand-specific methylation, or hemimethylation, must not do it; a run
-- writing a `.cov` for methylKit almost certainly must. Reading this file is
-- how the choice gets written down.
--
-- Everything downstream follows, because DuckDB resolves a view by name each
-- time it is queried: redefining `site_category` here changes what
-- `site_relation`, `block_summary` and every file under `formats/` see, with no
-- argument passed anywhere. It folds the raw category counts, so a measure
-- with four categories merges exactly as a two-category one does.
--
-- VARIABLES
--   merge   the measure to fold, e.g. 'mCG'. Every other measure is left
--           exactly as it was, because a context that is not palindromic has
--           no partner base and folding it would add two unrelated positions
--           together.
--   span    how far apart the two bases are. Default 1, which is a CG.
--
-- SPAN
-- The number is the distance between the two bases of the palindrome, and it
-- is the only thing here that knows anything about a context:
--
--   CG        span 1    C at p, C at p+1 on the minus strand
--   CAG/CTG   span 2    the symmetric part of plant CHG methylation; the two
--                       cytosines sit two apart. CCG is not symmetric this way
--                       -- its partner is a CG -- so a run folding CHG should
--                       declare CCG as its own measure and leave it alone.
--   GATC      span 1    the two adenines of a Dam site, for 6mA
--
-- COORDINATES
-- The minus-strand base moves back onto its partner, which is correct whatever
-- base the hit table counts from: the shift is a difference between two
-- positions, not an absolute. The merged row is reported on the plus strand,
-- which is what every format below expects to see.

CREATE OR REPLACE VIEW site_category AS
SELECT block,
       contig,
       CASE WHEN measure = getvariable('merge') AND strand = '-'
            THEN refr_pos - coalesce(getvariable('span'), 1)
            ELSE refr_pos END AS refr_pos,
       CASE WHEN measure = getvariable('merge')
            THEN '+' ELSE strand END AS strand,
       measure,
       category,
       sum(n)::BIGINT AS n
FROM site_stranded
GROUP BY ALL;

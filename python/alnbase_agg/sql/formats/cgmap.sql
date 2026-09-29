-- BS-Seeker2 / cgmaptools CGmap, bulk.
--
--     SET VARIABLE out = 'sample.CGmap';
--     .read formats/cgmap.sql
--
-- REQUIRES lib/count_sites.sql
--
-- Eight columns, tab separated, no header, one line per cytosine, both
-- strands, every context in the one file:
--
--     <chr> <nucleotide on the + strand> <pos> <class> <dinucleotide>
--     <level> <mC> <total>
--
-- This is the one format here that does not take a `relation`, and the reason
-- is its fourth and fifth columns: they name the context on the line itself,
-- so the file holds all of them and the context has to be a value rather than
-- a choice made once per file.
--
-- WHAT IT NEEDS FROM YOUR MEASURES
-- Both columns are read off the cytosine's trinucleotide, 5'->3' on its own
-- strand -- the class is CG, CHG or CHH, the dinucleotide is the cytosine plus
-- its 3' neighbour -- and `CHG` does not say which base H was. So this file
-- reports the relations whose names are trinucleotides, CAA through CTT, which
-- is what an alnbase query file enumerating the sixteen `C**` patterns puts in
-- the hit table's `name` column. Other relations the run declares, such as an
-- mCH used only by a read filter, are left out.
--
-- VARIABLES
--   out         the path to write.
--   min_depth   drop positions seen fewer than this many times. Default 1.
--
-- COLUMN 2 IS STRAND RELATIVE AND NOTHING ELSE IS
-- It is the base on the Watson strand: `C` for a cytosine on the plus strand
-- and `G` for one on the minus strand. The position is the cytosine's own
-- either way.
--
-- THE LEVEL IS A FRACTION, FORMATTED AS PYTHON PRINTS IT
-- BS-Seeker2 writes `str(round(level, 2))`. Two things about that are not
-- `printf('%.2f')`: Python's round breaks a tie to the even digit, so 5/8 is
-- 0.62 and not 0.63 -- at eight-read coverage, not an exotic case -- and `str`
-- keeps one decimal on a whole number, so a fully methylated site is `1.0` and
-- not `1`. `round_even` gives the first and the equality test gives the second.

-- An empty CGmap is what you get if no relation is named for a trinucleotide,
-- and an empty file is the failure mode worth catching: the run looks like it
-- worked. DuckDB's optimizer prunes `error()` out of a predicate or an unused
-- projection, but evaluates it here, where it is the value of a scalar
-- aggregate that the statement returns.
SELECT CASE WHEN count(*) > 0 THEN 'ok'
            ELSE error('formats/cgmap.sql found no relation named for a '
                       || 'trinucleotide. Its context columns are read off '
                       || 'names like CAA, CAC ... CTT, one per pattern, so a '
                       || 'run writing CGmap declares those as its relations.')
       END
FROM relation
WHERE regexp_matches(relation, '^C[ACGT]{2}$');

COPY (
  SELECT contig,
         CASE strand WHEN '+' THEN 'C' WHEN '-' THEN 'G' END  AS nucleotide,
         pos1(refr_pos)                                       AS pos,
         CASE WHEN substr(relation, 2, 1) = 'G' THEN 'CG'
              WHEN substr(relation, 3, 1) = 'G' THEN 'CHG'
              ELSE 'CHH' END                                  AS class,
         substr(relation, 1, 2)                                AS dinucleotide,
         CASE WHEN round_even(sum(num) * 1.0 / sum(den), 2)
                 = floor(round_even(sum(num) * 1.0 / sum(den), 2))
              THEN printf('%.1f',  round_even(sum(num) * 1.0 / sum(den), 2))
              ELSE printf('%.15g', round_even(sum(num) * 1.0 / sum(den), 2))
         END                                                  AS level,
         sum(num)                                             AS mC,
         sum(den)                                             AS total
  FROM site_relation
  WHERE regexp_matches(relation, '^C[ACGT]{2}$')
  GROUP BY contig, refr_pos, strand, relation
  HAVING sum(den) >= coalesce(getvariable('min_depth'), 1)
  ORDER BY contig, refr_pos
) TO (getvariable('out')) (FORMAT csv, DELIMITER E'\t', HEADER false);
-- No ELSE on the nucleotide CASE: alnbase writes '+' or '-' and nothing else,
-- so an unexpected strand comes out empty and visible rather than silently
-- counted as a minus strand.

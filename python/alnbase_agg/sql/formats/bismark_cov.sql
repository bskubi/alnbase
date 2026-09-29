-- Bismark coverage (`.cov`), bulk.
--
--     SET VARIABLE relation = 'mCG';
--     SET VARIABLE out = 'sample.CG.cov';
--     .read formats/bismark_cov.sql
--
-- REQUIRES lib/count_sites.sql
--
-- Six columns, tab separated, no header:
--
--     <chr> <start> <end> <percent methylated> <methylated> <unmethylated>
--
-- `start` and `end` are the same 1-based position. This is the format most of
-- the ecosystem accepts -- `coverage2cytosine`, methylKit, DSS and bsseq all
-- read it -- so it is the one to write when in doubt.
--
-- VARIABLES
--   relation  which relation this file reports. One number per line, so one.
--   out         the path to write.
--   min_depth   drop positions seen fewer than this many times. Default 1,
--               which drops nothing, because a position with any coverage at
--               all is a real observation.
--
-- THE PERCENTAGE
-- Bismark's comes out of Perl stringifying `100 * c / (c + t)`, which is
-- fifteen significant digits: `66.6666666666667`, not the seventeen DuckDB
-- would print. `printf('%.15g', ...)` matches it exactly, so a diff against
-- Bismark's own output is a diff about methylation rather than about floating
-- point.
--
-- BLOCKS
-- Summed away. A `.cov` file is bulk; `formats/amethyst_h5.sql` is the
-- per-block counterpart over the same `site_relation`.

COPY (
  SELECT contig,
         pos1(refr_pos)                                         AS start,
         pos1(refr_pos)                                         AS "end",
         printf('%.15g', 100.0 * sum(num) / sum(den))           AS percent,
         sum(num)                                               AS methylated,
         sum(den) - sum(num)                                    AS unmethylated
  FROM site_relation
  WHERE relation = getvariable('relation')
  GROUP BY contig, refr_pos
  HAVING sum(den) >= coalesce(getvariable('min_depth'), 1)
  ORDER BY contig, refr_pos
) TO (getvariable('out')) (FORMAT csv, DELIMITER E'\t', HEADER false);
-- The parentheses around the target are load bearing: `TO getvariable('out')`
-- is a parser error, because COPY's target is a literal unless it is wrapped
-- into an expression.

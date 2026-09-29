-- MethylDackel's default output, bulk.
--
--     SET VARIABLE relation = 'mCG';
--     SET VARIABLE out = 'sample_CpG.bedGraph';
--     .read formats/methyldackel_bedgraph.sql
--
-- REQUIRES lib/count_sites.sql
--
-- Six columns, tab separated, no header:
--
--     <chr> <start> <end> <percent methylated> <methylated> <unmethylated>
--
-- The same six fields as Bismark's coverage file, and neither of the other two
-- bed-like formats here: 0-based half-open where `bismark_cov.sql` is 1-based,
-- six columns where `bedgraph.sql` has four. It is its own file rather than an
-- option on either, because a coordinate convention chosen by a flag is the
-- kind of thing that produces output which is wrong everywhere and looks right
-- anywhere.
--
-- VARIABLES
--   relation  which relation this file reports.
--   out         the path to write.
--   min_depth   drop positions seen fewer than this many times. Default 1.
--
-- THE PERCENTAGE IS TRUNCATED, NOT ROUNDED
-- MethylDackel writes `(int)(100.0 * m / (m + u))` (`extract.c:46-51`), a C
-- cast, so two thirds methylated is `66` and a position that is 99.6%
-- methylated is `99`. Bismark, on the same numbers, writes `66.6666666666667`.
-- `floor` and not `CAST(... AS INTEGER)`: DuckDB's cast rounds, so it would
-- write 67 where MethylDackel writes 66. The counts are non-negative, so
-- flooring is exactly C's truncation toward zero.
--
-- MethylDackel also writes a `track type="bedGraph" description=...` line
-- first, which is a display directive rather than data and cannot be part of a
-- COPY. A diff against MethylDackel's own output skips its first line.

COPY (
  SELECT contig,
         pos0(refr_pos)                                       AS start,
         pos0(refr_pos) + 1                                   AS "end",
         floor(100.0 * sum(num) / sum(den))::INTEGER          AS percent,
         sum(num)                                             AS methylated,
         sum(den) - sum(num)                                  AS unmethylated
  FROM site_relation
  WHERE relation = getvariable('relation')
  GROUP BY contig, refr_pos
  HAVING sum(den) >= coalesce(getvariable('min_depth'), 1)
  ORDER BY contig, refr_pos
) TO (getvariable('out')) (FORMAT csv, DELIMITER E'\t', HEADER false);

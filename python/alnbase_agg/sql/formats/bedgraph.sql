-- Bismark bedGraph, bulk.
--
--     SET VARIABLE relation = 'mCG';
--     SET VARIABLE out = 'sample.CG.bedGraph';
--     .read formats/bedgraph.sql
--
-- REQUIRES lib/count_sites.sql
--
-- Four columns, tab separated, no header:
--
--     <chr> <start> <end> <percent methylated>
--
-- 0-based half-open, which is where it differs from `bismark_cov.sql` as much
-- as in its column count. A genome browser's idea of one number per interval,
-- and nothing about the format is specific to methylation: point `relation` at
-- an editing rate or a mismatch rate and the file is still a bedGraph.
--
-- VARIABLES
--   relation  which relation this file reports.
--   out         the path to write.
--   min_depth   drop positions seen fewer than this many times. Default 1.
--
-- The `track type=bedGraph` line Bismark emits first is not written. It is a
-- UCSC display directive rather than data, every reader tolerates its absence,
-- and emitting it would mean the file could not simply be a COPY.

COPY (
  SELECT contig,
         pos0(refr_pos)                                 AS start,
         pos0(refr_pos) + 1                             AS "end",
         printf('%.15g', 100.0 * sum(num) / sum(den))   AS percent
  FROM site_relation
  WHERE relation = getvariable('relation')
  GROUP BY contig, refr_pos
  HAVING sum(den) >= coalesce(getvariable('min_depth'), 1)
  ORDER BY contig, refr_pos
) TO (getvariable('out')) (FORMAT csv, DELIMITER E'\t', HEADER false);

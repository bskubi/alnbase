-- ENCODE bedMethyl, bulk.
--
--     SET VARIABLE relation = 'mCG';
--     SET VARIABLE out = 'sample.CG.bedmethyl';
--     .read formats/bedmethyl.sql
--
-- REQUIRES lib/count_sites.sql
--
-- Eleven columns, tab separated, no header: the nine BED columns, then
-- coverage and percent.
--
--     <chr> <start> <end> <name> <score> <strand>
--     <thickStart> <thickEnd> <itemRgb> <coverage> <percent methylated>
--
-- This is what modkit writes and what the ENCODE portal distributes, so it is
-- the format to reach for when the consumer is a browser or a pipeline built
-- around long-read modification calls rather than around Bismark.
--
-- VARIABLES
--   relation  which relation this file reports. Also the `name` column, which
--               is what makes a CG file and a CHH file distinguishable once
--               they have been concatenated.
--   out         the path to write.
--   min_depth   drop positions seen fewer than this many times. Default 1.
--
-- STRAND
-- Each cytosine keeps its own strand, which is the point of the column. A run
-- that read `lib/merge_strands.sql` has already collapsed each CG onto its plus
-- strand cytosine, so this file will say `+` throughout -- honest, but a
-- strand-merged bedMethyl conventionally says `.`, so merge or write this,
-- not both.
--
-- `score` is the coverage capped at 1000, which is the spec's rule and the
-- reason a browser can shade by it. `thickStart`/`thickEnd`/`itemRgb` carry no
-- information we have; the spec wants all nine BED columns present, so they
-- are filled rather than omitted.

COPY (
  SELECT contig,
         pos0(refr_pos)                                    AS start,
         pos0(refr_pos) + 1                                AS "end",
         getvariable('relation')                            AS name,
         least(sum(den), 1000)                             AS score,
         strand,
         pos0(refr_pos)                                    AS "thickStart",
         pos0(refr_pos) + 1                                AS "thickEnd",
         '0,0,0'                                           AS "itemRgb",
         sum(den)                                          AS coverage,
         printf('%.2f', 100.0 * sum(num) / sum(den))       AS percent
  FROM site_relation
  WHERE relation = getvariable('relation')
  GROUP BY contig, refr_pos, strand
  HAVING sum(den) >= coalesce(getvariable('min_depth'), 1)
  ORDER BY contig, refr_pos
) TO (getvariable('out')) (FORMAT csv, DELIMITER E'\t', HEADER false);

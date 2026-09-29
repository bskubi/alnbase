-- premethyst / Amethyst cellInfo, one line per block.
--
--     SET VARIABLE cg_relation = 'mCG';
--     SET VARIABLE ch_relation = 'mCH';
--     SET VARIABLE out = 'sample.cellInfo.txt';
--     .read formats/cellinfo.sql
--
-- REQUIRES lib/count_sites.sql
--
-- Six columns, tab separated, no header:
--
--     <barcode> <all calls> <CG sites> <CG percent> <CH sites> <CH percent>
--
-- This is the shape Amethyst's own vignette ships
-- (`pbmc_vignette_cellInfo.txt`) and what `obj@metadata` is built from: the
-- `mcg_pct` and `mch_pct` that amethyst-facet reads for its score and ratio
-- metrics are columns 4 and 6.
--
-- VARIABLES
--   cg_relation, ch_relation   which relations fill the two pairs of columns.
--                            Nothing here is specific to CG and CH beyond the
--                            names: a NOMe run points them at GCH and HCG.
--   out                      the path to write.
--   min_sites                drop blocks covering fewer than this many CG
--                            positions. Default 0, which drops nothing.
--
-- WHAT "COV" MEANS IN EACH COLUMN, WHICH IS NOT THE SAME THING TWICE
-- Column 2 counts *calls* across every relation the run declares -- how much
-- data the block produced. Columns 3 and 5 count *distinct positions* covered,
-- which is what premethyst and sciMETv2 mean by `CG_Cov`. A block covering one
-- CG twenty times contributes 20 to the first and 1 to the second.
--
-- premethyst's own file has ten columns: four more counting fragments, pairs,
-- singletons and an excluded-mCH percentage. Those are read accounting rather
-- than site accounting, and a run that wants them adds them here from
-- `read_relation`, which still holds one row per read.

COPY (
  SELECT block,
         sum(den)                                   AS all_calls,
         sum(sites) FILTER (WHERE relation = getvariable('cg_relation'))  AS cg_cov,
         printf('%.2f', 100.0
                * sum(num) FILTER (WHERE relation = getvariable('cg_relation'))
                / nullif(sum(den) FILTER (WHERE relation = getvariable('cg_relation')), 0))
                                                    AS cg_pct,
         sum(sites) FILTER (WHERE relation = getvariable('ch_relation'))  AS ch_cov,
         printf('%.2f', 100.0
                * sum(num) FILTER (WHERE relation = getvariable('ch_relation'))
                / nullif(sum(den) FILTER (WHERE relation = getvariable('ch_relation')), 0))
                                                    AS ch_pct
  FROM block_summary
  GROUP BY block
  HAVING coalesce(sum(sites) FILTER (WHERE relation = getvariable('cg_relation')), 0)
         >= coalesce(getvariable('min_sites'), 0)
  ORDER BY block
) TO (getvariable('out')) (FORMAT csv, DELIMITER E'\t', HEADER false);
-- `nullif(..., 0)` rather than a guard in the WHERE: a block with CG coverage
-- and no CH at all is a real block and belongs in the file, with an empty
-- percentage rather than a division error or a fabricated zero.

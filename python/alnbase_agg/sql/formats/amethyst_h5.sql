-- Amethyst HDF5 v2, one dataset per block. The SQL half.
--
--     SET VARIABLE fai = 'genome.fa.fai';
--     .read lib/load_contigs.sql
--     SET VARIABLE relation = 'mCG';
--     SET VARIABLE context = 'CG';
--     SET VARIABLE out = 'sample.h5';
--     .read formats/amethyst_h5.sql
--     .python formats/amethyst_h5.py     -- the chain adds this line itself
--
-- REQUIRES lib/count_sites.sql, lib/load_contigs.sql
--
-- This is the only format here that SQL cannot finish on its own, because the
-- file is not one table: it holds a dataset per block, at `/<context>/<block>/1`.
-- So the work is split at the natural seam. This file aggregates and sorts;
-- `amethyst_h5.py` walks the sorted stream and writes one group at a time,
-- holding one block in memory rather than a genome.
--
-- VARIABLES
--   relation  which relation fills the file.
--   context   the HDF5 group name, e.g. 'CG'. Amethyst puts no restriction on
--             it: GCH, HCG and CAC are all in use.
--   out       the path to write.
--   (fai)     read by lib/load_contigs.sql, which this file needs.
--
-- THE SORT IS THE CONTRACT
-- Amethyst's index stores only the first row number and a count per contig
-- (`R/index.R:58-78`), so every contig's rows must be contiguous in the
-- dataset. The ORDER BY here is what makes that true, and `check_sorted` in
-- the Python half refuses to write if this file is ever edited into something
-- that is not sorted by block.
--
-- Contigs are ordered by their position in the reference rather than by name,
-- because `ORDER BY chr` is lexical and would put chr10 between chr1 and chr2.
-- Either satisfies Amethyst, which only needs contiguity; reference order is
-- chosen because it is also what every coordinate-sorted format wants, so
-- there is one contig order in the toolchain and not two.
--
-- v2 DROPPED pct
-- The v2 dtype is `(chr S10, pos i8, c i8, t i8)` and no Amethyst reader uses
-- a percentage: `R/diff.R:297` explicitly drops the column when it finds one.
-- So it is not computed here.

CREATE OR REPLACE VIEW amethyst AS
SELECT s.block,
       s.contig          AS chr,
       pos1(s.refr_pos)    AS pos,
       s.num::BIGINT       AS c,
       (s.den - s.num)::BIGINT AS t
FROM site_relation s JOIN contigs k ON k.name = s.contig
WHERE s.relation = getvariable('relation')
  AND s.den > 0
ORDER BY s.block, k.i, s.refr_pos;
-- The join is an inner join on purpose. A contig in the data that the
-- reference index does not name is a sign that the hit tables and the FASTA
-- are not the same assembly, and dropping it silently would be the worst
-- outcome -- so the Python half counts what the join removed and refuses to
-- write a file that lost anything.

-- alnbase-agg: load the contig order from a FASTA index.
--
--     SET VARIABLE fai = 'genome.fa.fai';
--     .read lib/load_contigs.sql
--
-- Defines `contigs(name, i)`, where `i` is the contig's position in the
-- reference. Any format that has to sort by coordinate across contigs needs
-- it, because `ORDER BY contig` is lexical and puts chr10 between chr1 and
-- chr2 -- which an Amethyst H5 or a tabix index will not forgive.
--
-- A `.fai` is the convenient source because samtools has already written one
-- beside every reference anyone is using, and its first column is the contig
-- name in reference order. A `.genome`/`.chromSizes` file works identically.
-- A run with a handful of contigs can skip this file and write them out:
--
--     CREATE OR REPLACE VIEW contigs (name, i) AS
--       VALUES ('chr1', 1), ('chr2', 2);

CREATE OR REPLACE VIEW contigs AS
SELECT name, row_number() OVER () AS i
FROM read_csv(getvariable('fai'),
              header = false,
              delim = E'\t',
              columns = {'name': 'VARCHAR', 'length': 'BIGINT',
                         'offset': 'BIGINT', 'linebases': 'BIGINT',
                         'linewidth': 'BIGINT'},
              -- A .genome file has two columns and a .fai has five. Reading
              -- the first two and ignoring the rest covers both, and a file
              -- with a different shape fails here rather than silently
              -- ordering contigs by something that is not their order.
              null_padding = true);

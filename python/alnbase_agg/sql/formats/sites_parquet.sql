-- The site table itself, as parquet. The escape hatch.
--
--     SET VARIABLE out = 'sample.sites.parquet';
--     .read formats/sites_parquet.sql
--
-- REQUIRES lib/count_sites.sql
--
-- Every other file under `formats/` drops something: a block key, a strand, a
-- relation, a count that the target format has no column for. This one drops
-- nothing. It is the format to write when the aggregation is an intermediate
-- rather than an end -- another DuckDB query, a converter, an R session --
-- because parquet keeps the types, so whatever reads it back does not have to
-- be told what the columns were.
--
-- It is also the answer when a format this toolchain does not ship has a
-- converter that does: write this, convert from it, and nobody has to
-- implement the format twice.
--
-- VARIABLES
--   out           the path to write.
--   compression   passed through to DuckDB. Default is DuckDB's own (snappy).

COPY (
  SELECT * FROM site_relation
  WHERE den > 0
  ORDER BY block, contig, refr_pos, relation
) TO (getvariable('out'))
  (FORMAT parquet, COMPRESSION (coalesce(getvariable('compression'), 'snappy')));

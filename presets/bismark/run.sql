-- Bismark's methylation extractor output, as an alnbase-agg run.
--
--     python generate.py > queries.toml
--     alnbase query -q queries.toml --insertions emit \
--         -F qname,is_first_in_template,ref_name,strand reads.bam calls
--     alnbase-agg run run.sql -s hits='calls_*_*.parquet' -s sample=sample
--
-- Pair it with queries.toml, which says what counts as a call; this file says
-- how the calls are grouped and written. README.md has each of Bismark's
-- choices and where in its source it is made.
--
-- This is a specification, not a comparison harness. It is the tool's implicit
-- extraction query written out explicitly, so that reading it is reading the
-- list of decisions Bismark makes without saying so.

CREATE OR REPLACE VIEW src AS
SELECT *,
       {'qname': qname, 'mate': is_first_in_template} AS unit,
       'bulk' AS cell
FROM read_parquet(getvariable('hits'));
-- No WHERE. bismark_methylation_extractor applies no MAPQ or base-quality
-- threshold of its own -- filtering happens before it, in the aligner and in
-- deduplicate_bismark -- so a run that added one would not be reproducing it.

-- The measures are named after Bismark's own XM letters, because that is what
-- queries.toml puts in the hit table's `name` column: uppercase is methylated,
-- lowercase unconverted.
CREATE OR REPLACE VIEW measure (measure, name, role) AS VALUES
  ('CG',  'Z', 'num'), ('CG',  'Z', 'den'), ('CG',  'z', 'den'),
  ('CHG', 'X', 'num'), ('CHG', 'X', 'den'), ('CHG', 'x', 'den'),
  ('CHH', 'H', 'num'), ('CHH', 'H', 'den'), ('CHH', 'h', 'den');

-- `U` and `u` are Bismark's "context unknown": a cytosine whose next two
-- read-projected reference characters include an N or an inserted base.
-- Bismark writes them into the XM string and its extractor then discards them,
-- so they reach no .cov file. They are declared in queries.toml so the tag
-- round-trips, and are deliberately not a measure here. To count them, add:
--
--   ('unknown', 'U', 'num'), ('unknown', 'U', 'den'), ('unknown', 'u', 'den')

.read lib/filters.sql

-- No read filter either, for the same reason.
CREATE OR REPLACE VIEW resolved AS SELECT * FROM src;

.read lib/count_sites.sql

-- Deliberately no lib/merge_strands.sql. bismark2bedGraph does not merge the two
-- cytosines of a CG; `coverage2cytosine --merge_CpG` is the separate step that
-- does, and a preset that merged here would be reproducing the wrong command.

SET VARIABLE measure = 'CG';
SET VARIABLE out = getvariable('sample') || '.CG.bismark.cov';
.read formats/bismark_cov.sql

-- Bismark writes CG only unless bismark_methylation_extractor is given --CX.
-- Uncommenting these is that flag.
--
-- SET VARIABLE measure = 'CHG';
-- SET VARIABLE out = getvariable('sample') || '.CHG.bismark.cov';
-- .read formats/bismark_cov.sql
--
-- SET VARIABLE measure = 'CHH';
-- SET VARIABLE out = getvariable('sample') || '.CHH.bismark.cov';
-- .read formats/bismark_cov.sql

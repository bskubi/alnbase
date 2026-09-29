-- MethylDackel extract's default output, as an alnbase-agg run.
--
--     alnbase query -q queries.toml --insertions skip \
--         -F qname,is_first_in_template,ref_name,strand,mapq,flags,\
--is_paired,is_proper_pair,is_mate_unmapped,qual \
--         reads.bam calls
--     alnbase-agg run run.sql -s hits='calls_*_*.parquet' -s sample=sample
--
-- Pair it with queries.toml, which is the other half: that file says what
-- counts as a call, this one says how the calls are grouped and written.
-- README.md has what each of MethylDackel's choices is and where it comes from.

-- MethylDackel's default read and base filters (`common.c:415-430`,
-- `common.c:127`). They are options rather than hardcoded behaviour, which is
-- why they are here and not in queries.toml: a query file says what a
-- methylation call *is*, and these say which alignments this run chose to
-- believe. To compare against MethylDackel run with different thresholds,
-- change these to match -- or drop them and turn them off on both sides.
--
-- Not applied: --ignoreNH, because Bismark writes no NH tag.
CREATE OR REPLACE VIEW src AS
SELECT *,
       {'qname': qname, 'mate': is_first_in_template} AS unit,
       'bulk' AS cell
FROM read_parquet(getvariable('hits'))
WHERE mapq >= 10                                 -- -q 10
  AND flags & 3840 = 0                           -- -F 0xF00: secondary, QC fail,
                                                 --    duplicate, supplementary
  AND NOT (is_paired AND is_mate_unmapped)       -- per --keepSingleton
  AND NOT (is_paired AND NOT is_proper_pair)     -- per --keepDiscordant
  AND qual >= 5;                                 -- -p 5, on the cytosine's own base

CREATE OR REPLACE VIEW measure (measure, name, role) AS VALUES
  ('CG',  'CpG_methylated', 'num'), ('CG',  'CpG_methylated', 'den'),
  ('CG',  'CpG_unmethylated', 'den'),
  ('CHG', 'CHG_methylated', 'num'), ('CHG', 'CHG_methylated', 'den'),
  ('CHG', 'CHG_unmethylated', 'den'),
  ('CHH', 'CHH_methylated', 'num'), ('CHH', 'CHH_methylated', 'den'),
  ('CHH', 'CHH_unmethylated', 'den');

.read lib/filters.sql

-- MethylDackel judges rows, not reads: there is no per-read filter in it.
CREATE OR REPLACE VIEW resolved AS SELECT * FROM src;

.read lib/count_sites.sql

-- No merge. MethylDackel pileups per position and does not merge the two
-- cytosines of a CG unless --mergeContext is given. Both strands are already
-- here: alnbase's walk orients each read along its conversion strand, so a
-- cytosine on either strand is reported at its own coordinate.

SET VARIABLE measure = 'CG';
SET VARIABLE out = getvariable('sample') || '_CpG.bedGraph';
.read formats/methyldackel_bedgraph.sql

-- MethylDackel writes CpG only unless --CHG or --CHH is given. Uncommenting a
-- pair of lines is the whole difference from `MethylDackel extract --CHG`.
--
-- SET VARIABLE measure = 'CHG';
-- SET VARIABLE out = getvariable('sample') || '_CHG.bedGraph';
-- .read formats/methyldackel_bedgraph.sql
--
-- SET VARIABLE measure = 'CHH';
-- SET VARIABLE out = getvariable('sample') || '_CHH.bedGraph';
-- .read formats/methyldackel_bedgraph.sql

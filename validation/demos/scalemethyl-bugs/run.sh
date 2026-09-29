#!/usr/bin/env bash
# Reproduce the ScaleMethyl bug candidates in README.md.
#
#   SCALEMETHYL=/path/to/ScaleMethyl SCALEMETHYL_ENV=/path/to/scaleMethylTools-env \
#   ALLCOOLS_ENV=/path/to/allcools-env METHSCAN_ENV=/path/to/methscan-env \
#   PY=/path/to/python-with-duckdb ALNBASE=/path/to/alnbase ./run.sh [workdir]
#
# SCALEMETHYL is a checkout of github.com/ScaleBio/ScaleMethyl at 2a4973c. SCALEMETHYL_ENV is
# built from that checkout's envs/scaleMethylTools.conda.yml (../../envs/scalemethyl.yaml is
# a copy); it provides the pinned Python stack, bwameth.py, bwa-mem2, samtools, bgzip and
# tabix. ALLCOOLS_ENV is ../../envs/allcools.yaml; METHSCAN_ENV is ../../envs/methscan.yaml.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${SCALEMETHYL:?set SCALEMETHYL}" "${SCALEMETHYL_ENV:?set SCALEMETHYL_ENV}"
: "${ALLCOOLS_ENV:?set ALLCOOLS_ENV}" "${METHSCAN_ENV:?set METHSCAN_ENV}"
: "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$SCALEMETHYL_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }
promise() { printf -- '-- promised: %s\n' "$*"; }
shows() { printf -- '-- %s:\n' "$*"; }

# Reference, 0-based positions (ScaleMethyl's parquet and alnbase's hits are 0-based):
#   chrZ  3000 random bases (seed 1), for reads aligned by bwa-meth
#   chrC  two copies of a 60 bp unit: 5 A, then CAAA x 10 (CHH Cs at 5, 9, ..., 41),
#         then ACGA (a CG at 46-47), then 11 A. The second copy starts at 60.
"$PY" - <<'PYEOF'
import random
random.seed(1)
unit = "A" * 5 + "CAAA" * 10 + "ACGA" + "A" * 11
contigs = {
    "chrZ": "".join(random.choice("ACGT") for _ in range(3000)),
    "chrC": unit * 2,
}
with open("ref.fa", "w") as f:
    for name, seq in contigs.items():
        f.write(f">{name}\n{seq}\n")
PYEOF
samtools faidx ref.fa
"$ALNBASE" index ref.fa ref.aref >/dev/null 2>&1

# met_extract DIR BAM [OPTIONS]: run met_extract.py in DIR as the pipeline's Extract step
# does (bwa-meth aligner, one subprocess); prints cellInfo.txt with its rows sorted, as
# their order varies between runs
met_extract() {
  local dir="$1" bam="$2"; shift 2
  mkdir -p "$dir"
  (cd "$dir" && python "$SCALEMETHYL/bin/met_extract.py" "$WORK/$bam" --sample s \
     --subprocesses 1 --aligner bwa-meth --ref "$WORK/ref.fa" "$@" >/dev/null 2>&1 &&
   { head -n 1 s.cellInfo.txt; tail -n +2 s.cellInfo.txt | sort; } | column -t)
}
# calls DIR: ScaleMethyl's per-cell totals by context and strand from DIR/s.met_*.parquet
calls() {
  duck "select barcode, context, strand, sum(methylated) as methylated,
               sum(unmethylated) as unmethylated
        from (select *, 'CG' as context from '$1/s.met_CG.parquet'
              union all select *, 'CH' as context from '$1/s.met_CH.parquet')
        group by all order by all"
}
# alnbase_hits NAME: CG/CH methylated/unmethylated hits for NAME.bam as NAME_*.parquet,
# with each record's name and FLAG
alnbase_hits() {
  "$ALNBASE" query --query-file "$HERE/cg_ch.toml" --parquet -F qname -F flags "$1.bam" \
    ref.aref "$1.parquet" >/dev/null 2>&1
}
duck() { "$PY" -c 'import duckdb, sys; duckdb.sql(sys.argv[1]).show(max_width=200)' "$1"; }

say "SM-1: bwa-meth path takes strand from the alignment orientation (met_extract.py)"
# Paired-end reads from a directional library, aligned by bwa-meth as the pipeline does.
# 50 fragments of 300 bp across chrZ, alternately from the original top (OT) and original
# bottom (OB) strand; mates are the 100 bp at each end. Two cells read the same fragments:
#   ACGT+ACGT+AAAA  every CG C retained, every other C converted (CG 100%, CH 0%)
#   ACGT+ACGT+CCCC  every C converted (CG 0%, CH 0%)
"$PY" - <<'PYEOF'
ref = open("ref.fa").read().split()[1]
complement = str.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s.translate(complement)[::-1]
def bisulfite(s, keep_cg):
    return "".join("T" if b == "C" and not (keep_cg and s[i + 1:i + 2] == "G") else b
                   for i, b in enumerate(s))
with open("r1.fq", "w") as r1, open("r2.fq", "w") as r2:
    for cell, keep_cg in [("ACGT+ACGT+AAAA", True), ("ACGT+ACGT+CCCC", False)]:
        # Convert each whole strand, then cut fragments, so a CG split by a fragment end
        # keeps its state. The bottom strand is converted in its own 5'->3' orientation.
        top, bottom = bisulfite(ref, keep_cg), bisulfite(revcomp(ref), keep_cg)
        for i in range(50):
            strand = "OT" if i % 2 == 0 else "OB"
            start, end = 100 + 50 * i, 400 + 50 * i
            converted = top[start:end] if strand == "OT" else bottom[len(ref) - end:len(ref) - start]
            name = f"{strand}{i}:{cell}"
            r1.write(f"@{name}\n{converted[:100]}\n+\n{'I' * 100}\n")
            r2.write(f"@{name}\n{revcomp(converted)[:100]}\n+\n{'I' * 100}\n")
PYEOF
bwameth.py index-mem2 ref.fa >/dev/null 2>&1
bwameth.py --reference ref.fa r1.fq r2.fq 2>/dev/null | samtools sort -o sm1.bam - 2>/dev/null
samtools index sm1.bam
promise "in a directional library both mates of a fragment carry the same conversion; for OT,"
promise "C->T at reference C, whichever way the mate is aligned. Cell AAAA: CG 100%, CH 0%."
promise "Cell CCCC: CG 0%, CH 0%. Every call on + for OT fragments, on - for OB fragments"
shows "bwa-meth output: FLAG, YD tag and fragment strand (first two letters of the name)"
samtools view sm1.bam | awk '{split($NF, yd, ":"); print substr($1, 1, 2), $2, yd[3]}' |
  sort | uniq -c
shows "ScaleMethyl met_extract.py --threshold 1.0 (CH read filter off): cellInfo.txt"
met_extract sm1_t1 sm1.bam --threshold 1.0
shows "ScaleMethyl calls by cell, context and strand"
calls sm1_t1
shows "ScaleMethyl met_extract.py --threshold 0.5 (the pipeline default): cellInfo.txt"
met_extract sm1_t05 sm1.bam --threshold 0.5
shows "alnbase: calls by cell and context"
alnbase_hits sm1
duck "select regexp_extract(qname, ':(.*)', 1) as barcode, split_part(name, '_', 1) as context,
             count(*) filter (where name like '%\_methylated' escape '\') as methylated,
             count(*) filter (where name like '%unmethylated') as unmethylated
      from 'sm1_*.parquet' group by all order by all"
shows "alnbase on first-in-pair records only (FLAG 99, 83), the mates ScaleMethyl keeps at 0.5, vs ScaleMethyl's --threshold 0.5 parquet: positions whose per-cell counts differ"
duck "with alnbase as (
        select regexp_extract(qname, ':(.*)', 1) as barcode, split_part(name, '_', 1) as context,
               refr_pos as pos,
               count(*) filter (where name like '%\_methylated' escape '\') as methylated,
               count(*) filter (where name like '%unmethylated') as unmethylated
        from 'sm1_*.parquet' where flags in (99, 83) group by all),
      scalemethyl as (
        select barcode, 'CG' as context, pos, methylated, unmethylated from 'sm1_t05/s.met_CG.parquet'
        union all
        select barcode, 'CH', pos, methylated, unmethylated from 'sm1_t05/s.met_CH.parquet')
      select count(*) filter (where a.methylated is distinct from s.methylated
                                 or a.unmethylated is distinct from s.unmethylated) as differing,
             count(*) as positions
      from alnbase a full join scalemethyl s using (barcode, context, pos)"

say "SM-2: CH read filter tests the running ratio (met_extract.py)"
# Two forward single-end reads in one cell, 44M over one copy each of the chrC unit. Each
# covers 10 CHH Cs and retains exactly one of them, plus the CG C: first retains its
# leftmost CHH C (read at 0-based 5, C at 5), last its rightmost (read at 65, C at 101).
"$PY" - <<'PYEOF'
ref = open("ref.fa").read().split()[3]
def read(start, retained_chh):
    seq = "".join("C" if b == "C" and (i == retained_chh or ref[i + 1] == "G") else
                  "T" if b == "C" else b for i, b in enumerate(ref[start:start + 44], start))
    return seq
header = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrZ\tLN:3000\n@SQ\tSN:chrC\tLN:120\n"
with open("sm2.sam", "w") as f:
    f.write(header)
    for name, start, retained in [("first", 5, 5), ("last", 65, 101)]:
        f.write("\t".join([f"{name}:ACGT+ACGT+GGGG", "0", "chrC", str(start + 1), "60", "44M",
                           "*", "0", "0", read(start, retained), "I" * 44, "YD:Z:f"]) + "\n")
PYEOF
samtools sort -o sm2.bam sm2.sam 2>/dev/null && samtools index sm2.bam
promise "docs: 'Reads with greater than this percentage CH methylation will be discarded'."
promise "Both reads have 1 of 10 CH methylated (10%), so both are kept at --threshold 0.5,"
promise "with CG calls at 46 and 106 and CH_high 0"
shows "ScaleMethyl met_extract.py --threshold 0.5: cellInfo.txt"
met_extract sm2 sm2.bam --threshold 0.5
shows "ScaleMethyl CG rows"
duck "select chr, pos, strand, methylated, unmethylated from 'sm2/s.met_CG.parquet'"
shows "alnbase: calls per read"
alnbase_hits sm2
duck "select qname, count(*) filter (where name = 'CH_methylated') as CH_methylated,
             count(*) filter (where name like 'CH%') as CH_calls,
             list(refr_pos) filter (where name like 'CG%') as CG_positions
      from 'sm2_*.parquet' group by qname order by qname"

say "SM-3, SM-4, SM-14: per-cell ALLC and Bismark .cov exports (write_allc.py, write_bismark.py)"
# The pipeline's MATRIX_GENERATION steps on the SM-2 calls: one passing cell. Its one CG
# call is the C at 0-based 106, 1-based 107, in reference context CGA.
cp sm2/s.met_CG.parquet export/ 2>/dev/null || { mkdir -p export && cp sm2/s.met_CG.parquet export/; }
printf 'cell_id,sample,pass\nACGT+ACGT+GGGG,s,pass\n' > export/allCells.csv
(cd export &&
 python "$SCALEMETHYL/bin/write_allc.py" --met_calls s.met_CG.parquet --all_cells allCells.csv \
   --sample s --context CG >/dev/null 2>&1 &&
 tabix -b2 -e2 -s1 CG/ACGT+ACGT+GGGG.allc.tsv.gz &&
 python "$SCALEMETHYL/bin/write_bismark.py" --met_calls s.met_CG.parquet \
   --all_cells allCells.csv --sample s --context CG >/dev/null 2>&1)
ALLC=export/CG/ACGT+ACGT+GGGG.allc.tsv.gz
COV=export/CG/ACGT+ACGT+GGGG.CG.cov.gz
promise "docs: ALLC columns 'as described' by ALLCools, whose spec gives 1-based positions and a"
promise "3-base context; the C is at chrC 107 (reference bases 107-109: $(samtools faidx ref.fa chrC:107-109 | tail -1))"
promise "so: chrC  107  +  CGA  1  1  1"
shows "ScaleMethyl ALLC"
zcat "$ALLC"
shows "reference bases at the written position, 106-108"
samtools faidx ref.fa chrC:106-108 | tail -1
shows "tabix (as the pipeline indexes it) for the C, chrC:107"
tabix "$ALLC" chrC:107-107 | wc -l
shows "ALLCools extract-allc --mc_contexts CGN: rows written"
printf 'chrZ\t3000\nchrC\t120\n' > chrom.sizes
"$ALLCOOLS_ENV/bin/allcools" extract-allc --allc_path "$ALLC" --output_prefix export/extract \
  --mc_contexts CGN --chrom_size_path chrom.sizes --strandness both >/dev/null 2>&1 &&
  zcat export/extract.CGN-Both.allc.tsv.gz | wc -l
promise "docs: 'bismark .cov format ... columns for the chr, pos, percent_methylated,"
promise "methylated_count and unmethylated_count'. Bismark's .cov (which MethSCAn, linked from"
promise "write_bismark.py, reads by default): chrC  107  107  100  1  0"
shows "ScaleMethyl Bismark .cov"
zcat "$COV"
shows "MethSCAn prepare (default --input-format bismark)"
mkdir -p export/methscan_in && cp "$COV" export/methscan_in/
"$METHSCAN_ENV/bin/methscan" prepare export/methscan_in/*.cov.gz export/methscan 2>&1 |
  grep -v '^\s*$' | tail -3 || true

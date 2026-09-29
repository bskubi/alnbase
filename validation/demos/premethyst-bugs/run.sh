#!/usr/bin/env bash
# Reproduce the premethyst bug candidates in README.md.
#
#   PREMETHYST=/path/to/premethyst PREMETHYST_ENV=/path/to/premethyst-env \
#   BISMARK_ENV=/path/to/bismark-env PY=/path/to/python-with-duckdb \
#   ALNBASE=/path/to/alnbase ./run.sh [workdir]
#
# PREMETHYST is a checkout of github.com/adeylab/premethyst at fa5a47b. PREMETHYST_ENV is
# ../../envs/premethyst.yaml with BSBolt built into it by ../../envs/build_bsbolt.sh; it
# provides perl, samtools and bsbolt. BISMARK_ENV is ../../envs/bismark.yaml.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${PREMETHYST:?set PREMETHYST}" "${PREMETHYST_ENV:?set PREMETHYST_ENV}"
: "${BISMARK_ENV:?set BISMARK_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
mkdir -p "$WORK" && cd "$WORK"

# premethyst is committed without the executable bit, and bam-extract starts its per-cell
# workers by running `premethyst` from PATH, so run an executable copy of the checkout.
rm -rf premethyst && cp -r "$PREMETHYST" premethyst && chmod +x premethyst/premethyst
export PATH="$WORK/premethyst:$PREMETHYST_ENV/bin:$PATH"

say() { printf '\n== %s\n' "$*"; }
promise() { printf -- '-- promised: %s\n' "$*"; }
shows() { printf -- '-- %s:\n' "$*"; }

# Reference, 1-based positions (premethyst's .cov files are 1-based; alnbase's refr_pos
# is 0-based and is shown plus 1):
#   chrP  3000 random bases (seed 1)
#   chrQ  300 random bases over A, G and T with a CG after every 12, so its top strand
#         has Cs only in CG context
"$PY" - <<'PYEOF'
import random
random.seed(1)
contigs = {
    "chrP": "".join(random.choice("ACGT") for _ in range(3000)),
    "chrQ": "".join("".join(random.choice("AGT") for _ in range(12)) + "CG" for _ in range(22))[:300],
}
with open("ref.fa", "w") as f:
    for name, seq in contigs.items():
        f.write(f">{name}\n{seq}\n")
PYEOF
samtools faidx ref.fa
"$ALNBASE" index ref.fa ref.aref >/dev/null 2>&1

# Paired-end reads from a directional library: one 300 bp fragment per cell, with CG Cs
# retained and every other C converted. Read names follow premethyst's
# CELL:FRAGMENT#0 layout. The directional read (the one that reads the converted strand
# 5'->3') is named r1 and is given to the aligners as the first read.
#   cell    contig  fragment  strand  edit in the directional read, 50 bases from its 5' end
#   plain   chrP    200-500   OT      none
#   otins   chrP    700-1000  OT      2 inserted bases (AA)
#   otdel   chrP    1200-1500 OT      2 deleted reference bases
#   obins   chrP    1700-2000 OB      2 inserted bases (AA)
#   obdel   chrP    2200-2500 OB      2 deleted reference bases
#   noch    chrQ    0-300     OT      none (no CH C on the strand read)
"$PY" - <<'PYEOF'
contigs = dict(zip(*[iter(open("ref.fa").read().split())] * 2))
complement = str.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s.translate(complement)[::-1]
def bisulfite(s):
    return "".join("T" if b == "C" and s[i + 1:i + 2] != "G" else b for i, b in enumerate(s))
with open("r1.fq", "w") as r1, open("r2.fq", "w") as r2, open("r1_half.csv", "w") as frags:
    frags.write("cell,chr,first,last\n")
    for cell, contig, start, strand, edit in [
            ("plain", ">chrP", 200, "OT", None), ("otins", ">chrP", 700, "OT", "ins"),
            ("otdel", ">chrP", 1200, "OT", "del"), ("obins", ">chrP", 1700, "OB", "ins"),
            ("obdel", ">chrP", 2200, "OB", "del"), ("noch", ">chrQ", 0, "OT", None)]:
        ref = contigs[contig]
        end = start + 300
        converted = (bisulfite(ref)[start:end] if strand == "OT"
                     else bisulfite(revcomp(ref))[len(ref) - end:len(ref) - start])
        first = {None: converted[:100],
                 "ins": converted[:50] + "AA" + converted[50:98],
                 "del": converted[:50] + converted[52:102]}[edit]
        second = revcomp(converted)[:100]
        name = f"{cell}:1#0"
        r1.write(f"@{name}\n{first}\n+\n{'I' * len(first)}\n")
        r2.write(f"@{name}\n{second}\n+\n{'I' * len(second)}\n")
        # The half of the fragment read by r1 (1-based, inclusive); r2 reads the other half
        half = (start + 1, start + 150) if strand == "OT" else (start + 151, end)
        frags.write(f"{cell},{contig[1:]},{half[0]},{half[1]}\n")
PYEOF

# extract NAME BAM [OPTIONS]: premethyst bam-extract on BAM, name-sorted first as
# premethyst's own fastq-align step leaves it; prints cellInfo.txt with a header row
# (columns from bam_extract.pm:301) and its rows sorted, as their order varies between runs
extract() {
  local name="$1" bam="$2"; shift 2
  rm -rf "$name" "$name".*
  samtools sort -n -o "$name.nsrt.bam" "$bam" 2>/dev/null
  premethyst bam-extract "$@" -O "$name" "$name.nsrt.bam" >/dev/null 2>&1
  { printf 'cell\tall_calls\tCG_sites\tCG_pct\tCH_sites\tCH_pct\tfragments\tpairs\tsingles\tno_CH\n'
    sort "$name.cellInfo.txt"; } | awk -F '\t' -v OFS='\t' '{for (i = 1; i <= NF; i++) if ($i == "") $i = "\"\""; print}' |
    column -t -s $'\t'
}
# compare NAME BAM: premethyst's calls in NAME/*.cov against alnbase's hits for BAM's
# first-in-pair records, per cell, context and 1-based position, within the half of each
# fragment that r1 reads (r1_half.csv). The mates do not overlap, so premethyst's calls there
# come from r1 alone.
compare() {
  "$ALNBASE" query --query-file "$HERE/cg_ch.toml" --parquet -F qname -F flags "$2" ref.aref \
    "$1.alnbase.parquet" >/dev/null 2>&1
  duck "with alnbase as (
          select split_part(qname, ':', 1) as cell, split_part(name, '_', 1) as context,
                 refr_pos + 1 as pos,
                 count(*) filter (where name like '%\_methylated' escape '\') as methylated,
                 count(*) filter (where name like '%unmethylated') as unmethylated
          from '$1.alnbase_*.parquet' where flags & 64 = 64 group by all),
        premethyst as (
          select regexp_extract(filename, '([^/]+)\.(C[GH])\.cov$', 1) as cell,
                 regexp_extract(filename, '([^/]+)\.(C[GH])\.cov$', 2) as context,
                 column1 as pos, column4 as methylated, column3 as unmethylated
          from read_csv('$1/*.cov', header = false, delim = '\t', filename = true))
        select cell,
               count(*) filter (where a.pos is not null and p.pos is not null) as same_position,
               count(*) filter (where p.pos is null) as alnbase_only,
               count(*) filter (where a.pos is null) as premethyst_only
        from alnbase a full join premethyst p using (cell, context, pos)
          join 'r1_half.csv' h using (cell)
        where pos between h.first and h.last
        group by cell order by cell"
}
duck() { "$PY" -c 'import duckdb, sys; duckdb.sql(sys.argv[1]).show(max_width=200)' "$1"; }

# positions NAME CELL: CG positions in CELL's r1 half, from alnbase's first-in-pair hits and
# from premethyst's NAME/CELL.CG.cov
positions() {
  duck "select 'alnbase' as source, list(refr_pos + 1 order by refr_pos) as CG_positions
        from '$1.alnbase_*.parquet', 'r1_half.csv' h
        where qname like '$2:%' and name like 'CG%' and flags & 64 = 64 and h.cell = '$2'
          and refr_pos + 1 between h.first and h.last
        union all
        select 'premethyst', list(column1 order by column1)
        from read_csv('$1/$2.CG.cov', header = false, delim = '\t'), 'r1_half.csv' h
        where h.cell = '$2' and column1 between h.first and h.last"
}

say "PM-1, PM-3, PM-4: premethyst bam-extract on BSBolt alignments (XB tag, the default)"
bsbolt Index -G ref.fa -DB bsbolt_db >/dev/null 2>&1
bsbolt Align -DB bsbolt_db -F1 r1.fq -F2 r2.fq -O bsbolt >/dev/null 2>&1
shows "BSBolt alignments of r1: name, FLAG, POS, CIGAR"
samtools view -f 64 bsbolt.bam | cut -f 1,2,4,6
promise "every call at the reference position of its base, so the same positions as alnbase's"
promise "hits for the same reads. PM-3: the noch fragment's CG calls are reported. PM-4 (usage:"
promise "'-m Minimum chromosome size to retain (def = 10000000)'): no rows on chrP (3 kb) or chrQ"
shows "premethyst bam-extract: cellInfo.txt"
extract pm_bsbolt bsbolt.bam
shows "calls in each r1 half: positions in both, only in alnbase, only in premethyst"
compare pm_bsbolt bsbolt.bam
for cell in otins otdel obins obdel; do
  shows "$cell: CG positions in the r1 half"
  positions pm_bsbolt "$cell"
done

say "PM-2: premethyst bam-extract -B on Bismark alignments (XM tag)"
# The same reads, plus a cell whose r1 starts with 10 bases that are not in the reference, so
# that Bismark's --local mode soft-clips them; --local is used for every cell.
#   clip    chrP    2600-2900 OT      10 extra bases before r1 (soft clip), none in r2
"$PY" - <<'PYEOF'
contigs = dict(zip(*[iter(open("ref.fa").read().split())] * 2))
ref = contigs[">chrP"]
complement = str.maketrans("ACGT", "TGCA")
bisulfite = lambda s: "".join("T" if b == "C" and s[i + 1:i + 2] != "G" else b for i, b in enumerate(s))
converted = bisulfite(ref)[2600:2900]
first = "GATTAGATTA" + converted[:90]
second = converted.translate(complement)[::-1][:100]
with open("r1.fq", "a") as r1, open("r2.fq", "a") as r2, open("r1_half.csv", "a") as frags:
    r1.write(f"@clip:1#0\n{first}\n+\n{'I' * 100}\n")
    r2.write(f"@clip:1#0\n{second}\n+\n{'I' * 100}\n")
    frags.write("clip,chrP,2601,2750\n")
PYEOF
mkdir -p bismark_genome && cp ref.fa bismark_genome/
"$BISMARK_ENV/bin/bismark_genome_preparation" --path_to_aligner "$BISMARK_ENV/bin" bismark_genome >/dev/null 2>&1
PATH="$BISMARK_ENV/bin:$PATH" bismark --local --genome bismark_genome -1 r1.fq -2 r2.fq \
  -o bismark --basename bismark >/dev/null 2>&1
shows "Bismark alignments of r1: name, FLAG, POS, CIGAR"
samtools view -f 64 bismark/bismark_pe.bam | cut -f 1,2,4,6
promise "the same positions as alnbase's hits (usage: '-B Methylation call field is in Bismark"
promise "format (XM:Z:) Eg. Bismark or UA-Meth as aligner')"
shows "premethyst bam-extract -B: cellInfo.txt"
extract pm_bismark bismark/bismark_pe.bam -B
shows "calls in each r1 half: positions in both, only in alnbase, only in premethyst"
compare pm_bismark bismark/bismark_pe.bam
for cell in clip otins otdel; do
  shows "$cell: CG positions in the r1 half"
  positions pm_bismark "$cell"
done

say "PM-11: the last record of a cell is dropped when it has no mate (bam_extract.pm)"
# The BSBolt BAM with r2 of the plain fragment removed, as bam-rmdup's per-record MAPQ
# filter (samtools view -q 10) removes a low-quality mate; and single-end alignments of r1.
samtools view -b -e '!(qname =~ "^plain:" && flag.read2)' -o no_plain_r2.bam bsbolt.bam
promise "a fragment with one record is extracted as a single (cellInfo counts 'singles');"
promise "usage: fastq-align takes -2 as optional, so single-end input is accepted"
shows "premethyst bam-extract, plain without r2: cellInfo.txt"
extract pm_no_plain_r2 no_plain_r2.bam
bsbolt Align -DB bsbolt_db -F1 r1.fq -O bsbolt_se >/dev/null 2>&1
shows "premethyst bam-extract, single-end r1 alignments: cellInfo.txt"
extract pm_se bsbolt_se.bam

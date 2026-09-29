#!/usr/bin/env bash
# Reproduce the MethylDackel bug candidates in README.md on hand-written reads.
#
#   METHYLDACKEL=/path/to/MethylDackel SAMTOOLS=/path/to/samtools \
#   PY=/path/to/python-with-duckdb ALNBASE=/path/to/alnbase \
#   ./run.sh [workdir]
#
# METHYLDACKEL is the source build from ../../envs/build_methyldackel.sh (commit 3c77bda).
# Every case prints what MethylDackel's help or usage text promises, then what it does.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
MIRROR="$HERE/../../../presets/methyldackel/queries.toml"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${METHYLDACKEL:?set METHYLDACKEL}" "${SAMTOOLS:?set SAMTOOLS}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }
promise() { printf -- '-- promised: %s\n' "$*"; }
shows() { printf -- '-- %s:\n' "$*"; }

# Reference, 0-based positions:
#   chrA  CGCGCGCGCG + 30 A        C of a CG at 0, 2, 4, 6, 8
#   chrB  10 A + CG + 28 A         C of a CG at 10
#   chrC  1 Mb of A                C of a CG at 100 and at 900000
"$PY" - <<'EOF'
def contig(name, length, cgs):
    s = bytearray(b"A" * length)
    for c in cgs:
        s[c:c + 2] = b"CG"
    return f">{name}\n" + "\n".join(s[i:i + 60].decode() for i in range(0, length, 60)) + "\n"
with open("ref.fa", "w") as f:
    f.write(contig("chrA", 40, [0, 2, 4, 6, 8]))
    f.write(contig("chrB", 40, [10]))
    f.write(contig("chrC", 1_000_000, [100, 900_000]))
EOF
"$SAMTOOLS" faidx ref.fa
"$ALNBASE" index ref.fa ref.aref >/dev/null 2>&1

# bam NAME: SAM records on stdin (tab-separated) -> sorted, indexed NAME.bam
bam() {
  { printf '@HD\tVN:1.6\tSO:unsorted\n'
    awk '{printf "@SQ\tSN:%s\tLN:%s\n", $1, $2}' ref.fa.fai
    tr ' ' '\t'; } | "$SAMTOOLS" sort -o "$1.bam" - && "$SAMTOOLS" index "$1.bam"
}

# alnbase hit rows from the MethylDackel mirror query, as OUT_PREFIX_*.parquet.
alnbase_hits() {  # BAM OUT_PREFIX
  "$ALNBASE" query --query-file "$MIRROR" --parquet -F qname "$1" ref.aref "$2.parquet" >/dev/null 2>&1
}
duck() { "$PY" -c 'import duckdb, sys; duckdb.sql(sys.argv[1]).show()' "$1"; }

say "MD-8: --OT left bound (extract)"
bam md8 <<'EOF'
r1 0 chrA 1 60 20M * 0 0 CGCGCGCGCGAAAAAAAAAA IIIIIIIIIIIIIIIIIIII
EOF
promise "--help: '--OT A,B,C,D' includes calls at 1-based read positions A through B;"
promise "'--OT 5,0,0,0 would include all but the first 4 bases'."
promise "Read positions of the CG Cs are 1, 3, 5, 7, 9 (0-based 0, 2, 4, 6, 8)."
for ot in 0,0,0,0 3,0,0,0 0,5,0,0; do
  "$METHYLDACKEL" extract -q 0 --OT "$ot" -o "md8_$ot" ref.fa md8.bam 2>/dev/null
  printf 'MethylDackel --OT %s: starts ' "$ot"
  tail -n +2 "md8_${ot}_CpG.bedGraph" | cut -f2 | paste -sd' '
done
echo "expected: 0,0,0,0 -> 0 2 4 6 8; 3,0,0,0 -> 2 4 6 8; 0,5,0,0 -> 0 2 4"
alnbase_hits md8.bam md8
shows "alnbase, 1-based read position A = off_5p A-1, so --OT 3,0,0,0 is off_5p >= 2"
duck "select refr_pos, off_5p from 'md8_*.parquet' where name like 'CpG%' and off_5p >= 2 order by refr_pos"

say "MD-14: bedGraph percentage (extract)"
bam md14 <<'EOF'
m1 0 chrA 1 60 10M * 0 0 CGCGCGCGCG IIIIIIIIII
m2 0 chrA 1 60 10M * 0 0 CGCGCGCGCG IIIIIIIIII
u1 0 chrA 1 60 10M * 0 0 TGCGCGCGCG IIIIIIIIII
EOF
promise "README: column 4 is 'the methylation percentage rounded to an integer'; 2 of 3 is 66.67, so 67."
"$METHYLDACKEL" extract -q 0 -o md14 ref.fa md14.bam 2>/dev/null
shows "MethylDackel, site 0"
awk 'NR > 1 && $2 == 0' md14_CpG.bedGraph

say "PR-3: NH filter (perRead)"
bam pr3 <<'EOF'
r3 0 chrB 1 60 20M * 0 0 AAAAAAAAAACGAAAAAAAA IIIIIIIIIIIIIIIIIIII NH:i:2
r3 256 chrB 1 60 20M * 0 0 AAAAAAAAAACGAAAAAAAA IIIIIIIIIIIIIIIIIIII NH:i:2
EOF
promise "usage: 'if an NH tag is present and its value is >1 then an entry is ignored', so no rows;"
promise "and --ignoreNH is an option."
shows "MethylDackel perRead"
"$METHYLDACKEL" perRead ref.fa pr3.bam 2>/dev/null
shows "MethylDackel perRead --ignoreNH"
"$METHYLDACKEL" perRead --ignoreNH ref.fa pr3.bam 2>&1 | head -n 1 || true

say "PR-4: minimum base quality (perRead)"
# r4a: bases 9 and 10 (0-based) have Phred 2; base 10 is the C of the CG.
# r4b: base 9, the last aligned base, has Phred 2; base 10 (the C) is soft-clipped.
bam pr4 <<'EOF'
r4a 0 chrB 1 60 20M * 0 0 AAAAAAAAAACGAAAAAAAA IIIIIIIII##IIIIIIIII
r4b 0 chrB 1 60 10M5S * 0 0 AAAAAAAAAACAAAA IIIIIIIII#IIIII
EOF
promise "usage: '-p Minimum Phred threshold to include a base' (default 5); soft clips are not aligned."
promise "The CG C is either below Q5 (r4a) or soft-clipped (r4b), so both rows are 0.0 with 0 bases."
shows "MethylDackel perRead"
"$METHYLDACKEL" perRead ref.fa pr4.bam 2>/dev/null
alnbase_hits pr4.bam pr4
shows "alnbase, CG calls at Phred >= 5 per read (a soft-clipped base is a clip column, never a C)"
duck "select qname, count(*) filter (where name like 'CpG%' and qual >= 5) as cg_calls
      from 'pr4_*.parquet' group by qname order by qname"

say "PR-5: -l BED (perRead)"
bam pr5 <<'EOF'
in 0 chrC 101 60 20M * 0 0 CGAAAAAAAAAAAAAAAAAA IIIIIIIIIIIIIIIIIIII
out 0 chrC 900001 60 20M * 0 0 CGAAAAAAAAAAAAAAAAAA IIIIIIIIIIIIIIIIIIII
EOF
printf 'chrC\t50\t150\n' > pr5.bed
promise "usage: '-l FILE A BED file listing regions for inclusion', so only read 'in' (POS 100 0-based)."
shows "MethylDackel perRead -l pr5.bed"
"$METHYLDACKEL" perRead -l pr5.bed ref.fa pr5.bam 2>/dev/null

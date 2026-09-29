# premethyst bugs, reproduced

Five premethyst `bam-extract` behaviours that put calls at the wrong position, drop calls, or
ignore a documented option. Four come from reading the source
(`docs/design/research/methylation-callers-hardcoded.md`, IDs PM-1, PM-2, PM-3, PM-4); PM-11
was found while building this demo. This directory confirms each one by running the tool.
sciMETv2's `sciMET_BSBolt2cellCalls.pl` has the same extraction logic.

Tool: premethyst GitHub master `fa5a47b`, run from a copy of the checkout with perl and
samtools 1.24 (`../../envs/premethyst.yaml`). Inputs come from the aligners premethyst
reads: BSBolt GitHub master `ea4870e` (1.6.0; built by `../../envs/build_bsbolt.sh`) for the
default `XB` tag, and Bismark 0.25.1 with bowtie2 (`../../envs/bismark.yaml`) for `-B`.
Checked on 2026-09-17 with alnbase 0.1.10.

```
PREMETHYST=... PREMETHYST_ENV=... BISMARK_ENV=... PY=python-with-duckdb ALNBASE=... \
  ./run.sh [workdir]
```

`expected.txt` is the output of that run (reruns are identical). The reads are simulated
paired-end fragments from a directional library, one 300 bp fragment per cell, with CG Cs
retained and all other Cs converted, aligned by BSBolt and by Bismark (`--local`). The
directional read (r1) of some fragments carries a 2 bp insertion or deletion, or 10 bases
that Bismark soft-clips. alnbase's calls for the same BAMs come from `cg_ch.toml` in this
directory. The comparison uses first-in-pair records and the half of each fragment that r1
reads; the mates do not overlap, so premethyst's calls there come from r1 alone (see the
BSBolt FLAG note below for why r2 is left out).

| ID | Command | Expected | Observed | Cause |
|---|---|---|---|---|
| PM-1 | `bam-extract` on BSBolt BAMs | each call at its base's reference position: the CG positions alnbase reports for the same record | the fragment without indels matches at all 22 positions. After a 2 bp insertion every call is 2 positions to the right (`otins`: 768, 785, 795 written as 770, 787, 797); after a 2 bp deletion, 2 to the left (`otdel`: 1257 as 1255); the same on original-bottom reads | the position starts at POS and advances by 1 per call letter and by XB's digit counts (`bam_extract.pm:352-381`). BSBolt's digits count query bases, including inserted bases and not deleted ones (BSBolt `bwa.c:285-316`), and premethyst never reads the CIGAR |
| PM-2 | `bam-extract -B` on Bismark BAMs | usage: "-B Methylation call field is in Bismark format (XM:Z:) Eg. Bismark or UA-Meth as aligner"; calls at the same positions as alnbase's | `10S90M`: every call 10 positions to the right (2603 written as 2613). Insertions and deletions shift calls as in PM-1 | the position advances by 1 per XM character, including soft-clipped and inserted bases, and not over deleted ones (`:334-351`) |
| PM-3 | `bam-extract` | usage documents only "-M Max allowed fraction mCH sites methylated"; a fragment on chrQ whose strand has Cs only in CG context gives its CG calls (7 in r1's half) | no rows; cellInfo counts the fragment in its last column | a fragment is used only if it has at least one CH call (`:255-273`); the zero-CH case is counted as excluded |
| PM-4 | `bam-extract` | usage: "-m Minimum chromosome size to retain (def = 10000000) Used to exclude random and other small contigs"; no rows on a 3 kb or 300 bp contig | rows on both | `$minSize` is set (`:13`) and printed in the usage (`:36`) but never read |
| PM-11 | `bam-extract` on a BAM with one mate removed, and on single-end BSBolt BAMs | a fragment with one record is extracted as a single (cellInfo has a `singles` column); `fastq-align` takes `-2` as optional, so single-end input is expected | the cell whose r2 was removed has 0 fragments and no calls. With single-end input every cell has 0 fragments | a record is held until the next record shows whether it is a mate (`:183-197`). When the barcode changes, the held record is discarded without being written (`:164-177`), and so is the last one in the file (`:204-206`) |

## What this means for comparisons

- **PM-1 and PM-2** shift every call downstream of an indel along the read by the indel
  length, and PM-2 shifts every call on a soft-clipped Bismark record by the clip length.
  The shifted calls land on other positions, often not Cs, and are mixed into the per-cell
  `.cov` files with correct calls from other fragments. Indels are a few per cent of reads in
  typical data; soft clips are rare with Bismark's default end-to-end mode but common with
  `--local` or other `-B` aligners. The error is silent. A comparison against premethyst
  should expect position differences on every record with an I, D or S in its CIGAR; alnbase
  places each call by the CIGAR.
- **PM-3** drops whole fragments that have CG calls but no CH call, so CG-dense, CH-poor
  fragments (CG islands read on a short fragment) are under-represented. Fragments above the
  `-M` threshold are also dropped, without being counted (follows from the code; not in this
  demo). alnbase leaves per-read filters to queries over its hits.
- **PM-4** means small contigs (unplaced scaffolds, alt contigs, spike-ins, chrM) stay in the
  output even though the usage says they are excluded by default.
- **PM-11** loses one fragment per cell whenever the cell's last record in name order has no
  mate, as when `bam-rmdup`'s per-record MAPQ filter removes one mate. With single-end input
  it loses every read, so single-end premethyst output is empty.
- In the Bismark run, `otdel` also differs at 1250: Bismark itself calls that C as CHH (XM `H`),
  because the deletion removed the G after it. premethyst writes that call at the right
  position; the difference is Bismark's, not premethyst's.
- **BSBolt FLAG (not a premethyst bug).** BSBolt sets FLAG 0x10 from
  the converted reference strand the read matched, not from the read's orientation
  (`bwamem.c:849-869`, `:890`), while SEQ is still reverse-complemented by orientation
  (`:937-950`). Single-end records and first-in-pair records of a directional library come
  out consistent. Second-in-pair records of original-top fragments are written with SEQ
  reverse-complemented but FLAG 0x10 unset (`plain` r2: FLAG 129, `YS:Z:W_G2A`), and the
  reverse for original-bottom fragments. premethyst reads positions from XB and is not
  affected. This was an alnbase limitation when this demo was written, because alnbase then
  took the strand from the FLAG and so called those records on the wrong strand and counted
  `off_5p` from the wrong end. It no longer is: pass `queries/strand/bsbolt.toml`, which
  takes the strand from `YS` and never reads 0x10, and BSBolt BAMs need no FLAG repair. The
  bit and the rule were both measured on BSBolt's own output in `../bsbolt-strand/`.

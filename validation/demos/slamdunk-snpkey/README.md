# SlamDunk's SNP mask collides across chromosomes (SD-2), reproduced

SlamDunk masks T→C conversions that fall on known SNPs, so that germline variation is not
counted as metabolic labelling. The mask is a Python dict keyed by

```python
key = snp[0] + snp[1]        # SNPtools.py:33 — CHROM concatenated with POS, no separator
```

so `"chr11" + "234"` and `"chr1" + "1234"` are the same key. A T→C SNP on one chromosome
masks genuine conversions at the colliding position on the other.

Full entry: `docs/design/research/non-methylation-tool-findings.md` (SD-1..SD-6).

Tool: `t-neumann/slamdunk` at GitHub `14270e5` (v0.4.3, 2024-04-25),
`slamdunk/utils/SNPtools.py`. Checked 2026-09-18.

```
SLAMDUNK=<path to a slamdunk checkout> python probe.py
```

Needs `pybedtools`, which `SNPtools` imports at module scope. `expected.txt` is the output.
The probe exercises the real `SNPDictionary` class rather than a reimplementation, so it
does not require NextGenMap or a BAM.

## 1. One SNP, two chromosomes

```
declared: one T->C SNP at chr11:234 (1-based)

  isTCSnp(chr11 ,    233)  1-based     234  -> True
  isTCSnp(chr1  ,   1233)  1-based    1234  -> True     <-- wrong chromosome
  isTCSnp(chr12 ,    233)  1-based     234  -> False
```

The error is one-directional: a conversion at the colliding position is discarded as
germline. It never creates a conversion, so the effect is a deflated conversion rate and a
missing new-RNA read, not a false positive.

## 2. Footprint on GRCh38

Two names can collide only when one is a prefix of the other **and** the concatenated
digits form a valid position string. `chr1`/`chr10` is therefore safe — `str(pos)` never
carries a leading zero — but `chr1`/`chr11` is not.

| | colliding positions | share |
|---|---|---|
| chr1 ↔ chr11 | 99,999,999 | 40.2% of chr1 |
| chr1 ↔ chr12 | 48,956,422 | 19.7% of chr1 |
| chr2 ↔ chr21 | 46,709,983 | 19.3% of chr2 |
| chr2 ↔ chr22 | 42,193,529 | 17.4% of chr2 |
| chr1 ↔ chr13..chr19 | 9,999,999 each | 4.0% of chr1 each |

**87.9% of chr1 and 36.7% of chr2** lie in a shared key space: 308 Mb, 10.0% of the primary
assembly. An Ensembl-style reference using bare `1`, `11`, `21` collides identically.

How much of that footprint turns into an actual mask depends on T→C SNP density in the
partner region, so the practical damage is far smaller than the footprint. SlamDunk's own
`snp` dunk calls variants from the sample with VarScan and emits only variant sites, which
keeps the VCF small; a dbSNP-scale VCF passed to `alleyoop -v/--vcf` does not.

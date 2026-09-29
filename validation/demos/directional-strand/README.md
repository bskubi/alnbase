# The one rule with no tag to read: `directional.toml` against Bismark's tags

[`queries/strand/directional.toml`](../../../queries/strand/directional.toml) is the eleventh
of the twelve strand rules to be run against real output, and the odd one out. Every other
file describes some aligner's tag. This one reads no tag at all: it is a claim about the
**library** — that read 1 was sequenced from a converted original and read 2 from that
original's PCR copy — so it applies to any conformant aligner's output and rests on none of
the evidence any of them writes.

That is what makes it hard to check. There is nothing in the BAM to compare it against, and
it is also alnbase's historical `--library directional` behaviour written out, so every row
alnbase has ever produced traces back to it.

```
BISMARK_ENV=<env> PY=<python with pyarrow> ALNBASE=<alnbase> ./run.sh [workdir]
```

`BISMARK_ENV` is the environment from [`../../envs/bismark.yaml`](../../envs/bismark.yaml):
Bismark 0.25.1, bowtie2 and samtools.

## Three answers per record

Bismark is the way in. Its `XR`/`XG` pair names all four strands per record — the read's own
conversion and the index that won — and [`bismark.toml`](../../../queries/strand/bismark.toml)
was checked against real output in [`../bismark-nondirectional-pe/`](../bismark-nondirectional-pe/).
So on a directional Bismark BAM every record has two independent four-way calls, one from the
FLAG alone and one from the tags, and a third answer the fixture knows: the strand the read
was actually built from.

The fixture is the 1200-base reference every strand demo uses, as a directional library — each
fragment giving an OT/CTOT duplex or an OB/CTOB duplex, read 1 always the converted original.
Two probe pairs make a readable table; a sweep of 200 pairs, 60 bp from each end of a 200 bp
fragment every 10 bases along both strands, turns the same question into a count.

| arm | records | `directional.toml` vs. truth | `bismark.toml` vs. truth | the two rules vs. each other |
|---|---|---|---|---|
| paired-end | 402 | 0 wrong | 0 wrong | 0 disagreements |
| single-end | 202 | 0 wrong | 0 wrong | 0 disagreements |

alnbase then walks the same records under each rule in turn and asks one CG query. A hit's
anchor is a cytosine in the read, and which reference base it lands on is decided by the
strand the rule declared — the C of a forward CG on a top-strand record, the G of that same CG
on a bottom-strand one — so the positions check the walk rather than the label:

```
-- pe
   directional.toml   2601 CG hits over  402 records; anchors on the wrong base of the CG ...: 0
   bismark.toml       2601 CG hits over  402 records; anchors on the wrong base of the CG ...: 0
   hits where the two rules disagree: 0 (of 2601 and 2601)
```

with 1320 apiece and 0 disagreements on the single-end arm.

## What the single-end arm settles

`bismark.toml`'s comment says that on single-end output Bismark sets 0x10 from the **genome
conversion** rather than from SEQ's orientation, which is why that file has to read the tags
at all. `directional.toml` reads 0x10. Whether it is therefore wrong on single-end Bismark
output is a question its own comment did not answer — it listed "Bismark 0.25.1 directional
paired-end" and stopped.

It is not wrong, and the run says why rather than asserting it:

```
records whose 0x10 contradicts how SEQ was actually stored: 0
```

On a directional library the two ways of setting the bit coincide. `XG == "GA"` means the read
matched the bottom-strand index, which for a directional library means it was sequenced from
the bottom strand, which means SEQ had to be reverse-complemented to store it
reference-forward. The non-conformant bit and the conformant one agree on every record. They
part company exactly on the CTOT and CTOB reads a directional library does not contain — which
is the same boundary the rule's own comment already draws, and the reason `bismark.toml`
remains the file to reach for when the library is not directional.

So the rule's comment now names Bismark single-end as well, on this evidence.

## What this round does not establish

The two rules agreeing is a strong check here, because they read disjoint evidence — one the
FLAG, the other two tags Bismark writes from its index choice — and both were compared against
a truth the fixture holds independently of either. That is not the case everywhere: in
[`../bwameth-strand/`](../bwameth-strand/) `directional.toml` and `bwameth.toml` also agreed
hit for hit, but both trace back to the contig bwa-meth chose, so they would agree even on a
record placed on the wrong strand.

- **Non-directional input is where this rule breaks, by design.** Its comment says so: a CTOT
  read arriving as read 1 is called OT. That is not measured here, because the fixture is a
  directional library; [`../bismark-nondirectional-pe/`](../bismark-nondirectional-pe/) has the
  records a non-directional library actually produces.
- **PBAT** was the twelfth rule, the same rule with the two conversion strands exchanged.
  [`../pbat-strand/`](../pbat-strand/) ran it and it is no longer shipped: Bismark `--pbat`
  and methylpy `--pbat` both hand back a BAM in the directional arrangement, so this file
  covers Bismark `--pbat` paired-end and `bismark.toml` covers both its layouts.
- The other aligners `directional.toml`'s comment names — BISCUIT `-b 1`, HISAT-3N
  `--directional-mapping`, methylpy's output, BSMAP `-n 0` single-end — are not run here. Each
  has its own demo directory; what this round adds is the four-way check against independent
  evidence, which only a four-way tag can give.
- No indels, no soft or hard clipping, no repeats, no multimappers, no unmapped records.
  Mapping efficiency is 100% on both arms and every record is stored reference-forward.

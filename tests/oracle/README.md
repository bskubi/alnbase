# Walk oracle

`walk_oracle.py` checks every column alnbase emits (reference position, read and
reference bases, quality, `off_5p`, `off_3p`) against an independent derivation from
pysam's `get_aligned_pairs(with_cigar=True)`, the CIGAR and pyfaidx. The random
records cover all four flag cases, M/=/X/I/D/N/S/H, contig edges and missing
qualities. The column types checked are aligned bases, deletions, emitted insertions,
pads past each read end, intron context and the `,@,` marker. Each is checked at
three `--end-context` / `--splice-context` settings, one of them zero.

Run it after `cargo test` whenever the walk, offsets, pads or intron handling change,
and before tagging a version:

```bash
cargo build --release
python tests/oracle/walk_oracle.py --alnbase "$CARGO_TARGET_DIR/release/alnbase" --records 3000 --seed 1
```

Needs Python with `pysam`, `pyfaidx` and `duckdb` (bioconda/conda-forge, or pip).
Exits 0 when every column agrees at every setting. CI runs the same command.

It is deliberately outside `cargo test`: the point is that none of its arithmetic is
alnbase's. Checked to fail against:

- alnbase 0.1.0 (read 2 offsets mirrored) and 0.1.3 (hard-clipped bases not counted in
  offsets);
- the oracle itself with any one of these rules broken: a gap takes the offset of the
  aligned base before it in walk order, pads continue from the soft clips, the marker
  sits at the first elided base in walk order.

Until it modelled pads and intron columns it failed on 0.1.8, which gives even a
width-1 query a pad at each end by default; the walk itself was not at fault.

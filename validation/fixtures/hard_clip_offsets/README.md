# A hard clip is still a sequenced base

`chrH` is 40 bp of A with a CG at 10-11 and another at 30-31 (1-based). Each read is 40
bases long as sequenced and is split in two, so each half carries a 20-base hard clip for
the half the other record holds. `fwd` is a forward split (`20M20H` at 1, `20H20M` at 21);
`rev` is the same split on the bottom strand, where SEQ and CIGAR are stored in reference
orientation, so the *trailing* `H` is the one at the sequenced 5' end.

Hard-clipped bases are absent from SEQ but were sequenced, so `off_5p` counts them
(`docs/reference/03-walk-and-matching.md` §3.2):

```
five_in_seq = (flags & 0x80) ? seq_len - 1 - read_off : read_off
off_5p      = five_in_seq + hard_clip_5p
off_3p      = (seq_len - 1 - five_in_seq) + hard_clip_3p
```

Both cytosines sit at SEQ index 9 of their own record, and the second one is nonetheless
`off_5p` 29 -- the 30th cycle, which is what it was. `off_5p + off_3p` is 39 on every row,
the sequenced length minus one, not the 19 that SEQ alone would give. A caller that indexes
into SEQ instead reports cycle 9 for both, which is BISCUIT's BC-2 (the hard clip shifts the
SEQ index in `pileup`) and BC-10 (`-5`/`-3` trimming measured along SEQ); BC-3 aborts on a
hard clip outright. It also matters for any M-bias plot or read-end trim, since those are
`off_5p` cutoffs and a split read's second half is all late cycles.

The two strands are the check on which `H` is which: the forward read's leading `H` and the
reverse read's trailing `H` are both the clip at the sequenced 5' end, and the table is
symmetric because of it.

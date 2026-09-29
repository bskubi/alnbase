"""Simulate hyper-edited RNA-seq reads with per-site ground truth.

    simulate.py OUTDIR [--seed N] [--reads N]

Writes into OUTDIR:

    genome.fa          one contig, chrS
    reads.fastq        single-end 100 bp reads
    truth_reads.tsv    one row per read: class, origin, CIGAR against the genome
    truth_sites.tsv    one row per edited base: the genomic coordinate it is on

The point of the simulation is that every editing event has a known genomic
coordinate, so a pipeline's output can be scored for recall (did it find the
site) and for accuracy (did it put the site in the right place). Those are
different failures and real data cannot tell them apart.

Layout. One synthetic chromosome holds a four-exon gene. Inside the last exon
sits an inverted pair of A-rich elements standing in for the inverted Alus that
ADAR actually binds; reads are drawn from there, so editing clusters where a
real experiment would put it. The elements are synthetic rather than a real Alu
consensus: the pipeline has no model of secondary structure, so nothing
downstream depends on the sequence being a genuine repeat, and a synthetic one
keeps base composition inside the filters the pipeline applies.

Read classes. Reads are drawn from the spliced transcript, so a read crossing
an exon boundary is junction-spanning by construction. Four classes carry
editing and differ only in how they sit on the genome -- contiguous, across a
junction, with a short indel, or with a non-genomic tail -- which is exactly
the axis the pipeline's flat-window assumption is sensitive to. Two unedited
classes are controls: a plain background read, and a read with a few scattered
mismatches that is not hyper-editing and must not be called as such.

Strand. ADAR deaminates adenosine on the transcript. Half the reads are emitted
reverse-complemented, where the same event reads out as T->C, which is what the
pipeline's twelve transform combinations exist to catch.
"""
import argparse
import os
import random

READ_LEN = 100
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


# --- genome ------------------------------------------------------------------

def random_seq(rng, n, gc=0.42):
    at, gcp = (1 - gc) / 2, gc / 2
    return "".join(rng.choices("ACGT", weights=[at, gcp, gcp, at], k=n))


def a_rich_element(rng, n=280):
    """Stand-in for an Alu: A-rich, but balanced enough to clear the pipeline's
    composition filter both before and after heavy A->G editing."""
    return "".join(rng.choices("ACGT", weights=[0.34, 0.20, 0.20, 0.26], k=n))


def build_genome(rng, length=300_000):
    seq = list(random_seq(rng, length))

    # Four exons. The gene sits in the middle of the contig with plenty of
    # flanking sequence, so reads never run off the end.
    exons = [(100_000, 100_900), (120_000, 120_600), (150_000, 150_500),
             (180_000, 181_400)]

    # An inverted pair of A-rich elements inside the last exon: the substrate.
    # The second copy carries 12% divergence from the first, which is roughly
    # how far apart two Alu copies in one hairpin actually are. Identical copies
    # would make every read from the pair a perfect two-locus mapper, and the
    # pipeline's rule for choosing between a read's alignments would then be
    # refusing an unanswerable question rather than being tested on a real one.
    elem = a_rich_element(rng)
    partner = list(revcomp(elem))
    for i in range(len(partner)):
        if rng.random() < 0.12:
            partner[i] = rng.choice([b for b in "ACGT" if b != partner[i]])
    left = 180_300
    right = 180_300 + len(elem) + 120
    seq[left:left + len(elem)] = list(elem)
    seq[right:right + len(elem)] = partner

    # Canonical GT..AG at every intron boundary. Real introns have them, and a
    # splice-aware aligner looks for them -- but note that none of the six
    # three-letter transforms leaves the motif intact (a2g turns AG into GG,
    # t2c turns GT into GC, and so on). Splice-aware alignment against a
    # transformed genome therefore cannot find junctions de novo and has to be
    # given an annotation. That is a property of the method, not of this
    # simulation.
    for (s_prev, e_prev), (s_next, _) in zip(exons, exons[1:]):
        seq[e_prev:e_prev + 2] = list("GT")
        seq[s_next - 2:s_next] = list("AG")

    return "".join(seq), exons, (left, right, len(elem))


def transcript_of(genome, exons):
    """Spliced sequence, plus a map from transcript offset to genomic offset."""
    parts, tmap = [], []
    for start, end in exons:
        parts.append(genome[start:end])
        tmap.extend(range(start, end))
    return "".join(parts), tmap


# --- composition filter ------------------------------------------------------

def passes_composition(read):
    """Mirror of the pipeline's per-read filter, closely enough that simulated
    reads are not silently discarded before the part under test.

    filter_sam.pl rejects a read when one base exceeds 60% of it, when any base
    is under 10%, or when it contains a simple repeat of ten or more units.
    """
    n = len(read)
    counts = {b: read.count(b) for b in "ACGT"}
    if max(counts.values()) / n > 0.6:
        return False
    if min(counts.values()) / n < 0.1:
        return False
    for unit in 1, 2, 3:
        for i in range(n - unit * 10 + 1):
            block = read[i:i + unit * 10]
            if block == block[:unit] * 10:
                return False
    return True


# --- reads -------------------------------------------------------------------

def blocks_from_tmap(tmap, t_start, t_len):
    """Genomic blocks covered by transcript[t_start : t_start+t_len]."""
    coords = tmap[t_start:t_start + t_len]
    blocks, run_start, prev = [], coords[0], coords[0]
    for c in coords[1:]:
        if c != prev + 1:
            blocks.append((run_start, prev + 1))
            run_start = c
        prev = c
    blocks.append((run_start, prev + 1))
    return blocks


def cigar_from_blocks(blocks):
    ops = []
    for i, (s, e) in enumerate(blocks):
        if i:
            ops.append(f"{s - blocks[i - 1][1]}N")
        ops.append(f"{e - s}M")
    return "".join(ops)


def apply_edits(rng, seq, rate):
    """A->G at a fraction of the adenosines, returning the edited sequence and
    the offsets that changed. Offsets are into `seq`, which is transcript
    orientation; the caller converts to genomic and to read orientation."""
    a_pos = [i for i, b in enumerate(seq) if b == "A"]
    k = max(1, int(round(len(a_pos) * rate)))
    chosen = sorted(rng.sample(a_pos, min(k, len(a_pos))))
    out = list(seq)
    for i in chosen:
        out[i] = "G"
    return "".join(out), chosen


def sequencing_errors(rng, seq, rate):
    out, hit = list(seq), []
    for i, b in enumerate(out):
        if rng.random() < rate:
            out[i] = rng.choice([x for x in "ACGT" if x != b])
            hit.append(i)
    return "".join(out), hit


CLASSES = ["flat_edited", "junction_edited", "del_edited", "ins_edited",
           "clip_edited", "flat_unedited", "snp_read"]


def draw_read(rng, genome, transcript, tmap, g_to_t, elem_span, cls, rid, err_rate):
    """Return (record, sites) or None if the draw has to be retried.

    `record` describes the read; `sites` is one entry per edited base, already
    in genomic coordinates and in read-as-sequenced offsets.
    """
    left, right, elem_len = elem_span

    # Draw from the A-rich elements for edited classes, anywhere for controls.
    if cls in ("flat_unedited", "snp_read"):
        g_lo, g_hi = 100_000, 181_400 - READ_LEN
    else:
        g_lo, g_hi = left - 40, right + elem_len - READ_LEN + 40

    if cls == "junction_edited":
        # Straddle one of the three junctions, off centre by a random amount.
        exon_end_t = rng.choice(junction_offsets(tmap))
        t_start = exon_end_t - rng.randint(25, READ_LEN - 25)
    else:
        g_start = rng.randint(g_lo, g_hi)
        if g_start not in g_to_t:
            return None
        t_start = g_to_t[g_start]

    span = READ_LEN
    if cls == "del_edited":
        del_len = rng.randint(2, 5)
        span = READ_LEN + del_len
    if t_start < 0 or t_start + span > len(transcript):
        return None

    window = transcript[t_start:t_start + span]
    blocks = blocks_from_tmap(tmap, t_start, span)
    if cls != "junction_edited" and len(blocks) > 1:
        return None
    if cls == "junction_edited" and len(blocks) < 2:
        return None

    # Editing, in transcript orientation, before any structural change.
    if cls.endswith("_edited"):
        rate = rng.uniform(0.25, 0.55)
        edited, edit_offs = apply_edits(rng, window, rate)
        if len(edit_offs) < 5:
            return None
    else:
        edited, edit_offs = window, []

    # transcript offset -> genomic coordinate, for the window as drawn
    win_gcoords = tmap[t_start:t_start + span]

    # Structural variation: turn the window into the read, and record where
    # each surviving transcript offset ended up along the read.
    #   t_to_r maps an offset in `edited` to an offset in the read, or None.
    t_to_r = {i: i for i in range(span)}
    read = edited
    cigar_extra = None

    if cls == "del_edited":
        cut = rng.randint(20, READ_LEN - 20)
        read = edited[:cut] + edited[cut + del_len:]
        t_to_r = {}
        for i in range(span):
            if i < cut:
                t_to_r[i] = i
            elif i >= cut + del_len:
                t_to_r[i] = i - del_len
            else:
                t_to_r[i] = None
        cigar_extra = ("D", cut, del_len)

    elif cls == "ins_edited":
        ins_len = rng.randint(2, 5)
        cut = rng.randint(20, READ_LEN - 20)
        ins = "".join(rng.choices("ACGT", k=ins_len))
        read = edited[:cut] + ins + edited[cut:READ_LEN - ins_len]
        t_to_r = {}
        for i in range(span):
            if i < cut:
                t_to_r[i] = i
            elif i < READ_LEN - ins_len:
                t_to_r[i] = i + ins_len
            else:
                t_to_r[i] = None
        cigar_extra = ("I", cut, ins_len)

    elif cls == "clip_edited":
        tail_len = rng.randint(8, 15)
        tail = "".join(rng.choices("ACGT", k=tail_len))
        if rng.random() < 0.5:
            read = tail + edited[:READ_LEN - tail_len]
            t_to_r = {i: (i + tail_len if i < READ_LEN - tail_len else None)
                      for i in range(span)}
            cigar_extra = ("S5", tail_len, 0)
        else:
            read = edited[:READ_LEN - tail_len] + tail
            t_to_r = {i: (i if i < READ_LEN - tail_len else None)
                      for i in range(span)}
            cigar_extra = ("S3", tail_len, 0)

    if cls == "snp_read":
        read = list(read)
        for i in rng.sample(range(READ_LEN), rng.randint(2, 4)):
            read[i] = rng.choice([x for x in "ACGT" if x != read[i]])
        read = "".join(read)

    read = read[:READ_LEN]
    if len(read) != READ_LEN:
        return None

    # Sequencing error, applied last so it is never mistaken for an edit.
    read, err_offs = sequencing_errors(rng, read, err_rate)

    # Emit sense or antisense.
    antisense = rng.random() < 0.5
    final = revcomp(read) if antisense else read

    def to_read_offset(t_off):
        r = t_to_r.get(t_off)
        if r is None or r >= READ_LEN:
            return None
        return (READ_LEN - 1 - r) if antisense else r

    sites = []
    for t_off in edit_offs:
        r = to_read_offset(t_off)
        if r is None:
            continue          # the edit fell inside the deletion or past the end
        pre_rc_offset = (READ_LEN - 1 - r) if antisense else r
        if pre_rc_offset in err_offs:
            continue          # a sequencing error landed on it; not scorable
        sites.append({
            "read_id": rid,
            "genomic_pos": win_gcoords[t_off] + 1,   # 1-based
            "read_offset": r,
            "read_base": final[r],
            "genome_base": genome[win_gcoords[t_off]],
        })

    if cls.endswith("_edited") and len(sites) < 5:
        return None
    if not passes_composition(final):
        return None

    record = {
        "read_id": rid,
        "class": cls,
        "chrom": "chrS",
        "start": blocks[0][0] + 1,
        "blocks": ",".join(f"{s + 1}-{e}" for s, e in blocks),
        "cigar_vs_genome": cigar_from_blocks(blocks),
        "structural": "" if cigar_extra is None else ":".join(map(str, cigar_extra)),
        "antisense": int(antisense),
        "n_edits": len(sites),
        "seq": final,
    }
    return record, sites


def junction_offsets(tmap):
    """Transcript offsets at which the genomic coordinate jumps."""
    out = []
    for t in range(1, len(tmap)):
        if tmap[t] != tmap[t - 1] + 1:
            out.append(t)
    return out


# --- main --------------------------------------------------------------------

MIX = {
    "flat_edited": 0.30,
    "junction_edited": 0.14,
    "del_edited": 0.10,
    "ins_edited": 0.10,
    "clip_edited": 0.12,
    "flat_unedited": 0.16,
    "snp_read": 0.08,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--seed", type=int, default=20260918)
    ap.add_argument("--reads", type=int, default=6000)
    ap.add_argument("--error-rate", type=float, default=0.002)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    genome, exons, elem_span = build_genome(rng)
    transcript, tmap = transcript_of(genome, exons)
    g_to_t = {g: t for t, g in enumerate(tmap)}

    # HISAT2's known-splicesite format: the flanking base on each side of the
    # intron, both 0-based.
    with open(os.path.join(args.outdir, "splicesites.txt"), "w") as f:
        for (_, e_prev), (s_next, _) in zip(exons, exons[1:]):
            f.write(f"chrS\t{e_prev - 1}\t{s_next}\t+\n")

    with open(os.path.join(args.outdir, "genome.fa"), "w") as f:
        f.write(">chrS\n")
        for i in range(0, len(genome), 60):
            f.write(genome[i:i + 60] + "\n")

    records, sites = [], []
    for cls, frac in MIX.items():
        want = int(round(args.reads * frac))
        made = 0
        tries = 0
        while made < want and tries < want * 200:
            tries += 1
            rid = f"{cls}_{made:05d}"
            got = draw_read(rng, genome, transcript, tmap, g_to_t, elem_span,
                            cls, rid, args.error_rate)
            if got is None:
                continue
            rec, st = got
            records.append(rec)
            sites.extend(st)
            made += 1
        if made < want:
            print(f"warning: only made {made}/{want} of {cls}")

    rng.shuffle(records)

    with open(os.path.join(args.outdir, "reads.fastq"), "w") as f:
        for r in records:
            f.write(f"@{r['read_id']}\n{r['seq']}\n+\n{'I' * len(r['seq'])}\n")

    cols = ["read_id", "class", "chrom", "start", "blocks", "cigar_vs_genome",
            "structural", "antisense", "n_edits"]
    with open(os.path.join(args.outdir, "truth_reads.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(records, key=lambda x: x["read_id"]):
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    scols = ["read_id", "genomic_pos", "read_offset", "read_base", "genome_base"]
    with open(os.path.join(args.outdir, "truth_sites.tsv"), "w") as f:
        f.write("\t".join(scols) + "\n")
        for s in sorted(sites, key=lambda x: (x["read_id"], x["read_offset"])):
            f.write("\t".join(str(s[c]) for c in scols) + "\n")

    by_class = {}
    for r in records:
        by_class[r["class"]] = by_class.get(r["class"], 0) + 1
    print(f"genome 300 kb, {len(exons)} exons, elements at {elem_span[0]} and {elem_span[1]}")
    print(f"{len(records)} reads, {len(sites)} edited bases with known coordinates")
    for c in CLASSES:
        n = by_class.get(c, 0)
        ns = sum(1 for s in sites if s["read_id"].startswith(c))
        print(f"  {c:16s} {n:5d} reads  {ns:6d} sites")


if __name__ == "__main__":
    main()

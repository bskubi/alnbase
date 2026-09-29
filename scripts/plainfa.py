#!/usr/bin/env python3
"""Print a region of a FASTA with no index involved: chr19:3034740-3034780.

Streams the file and counts bases, so it cannot be fooled by a stale .fai or
by anything a tool's own index believes. Slow, but it is the arbiter.
"""
import sys

fa, region = sys.argv[1], sys.argv[2]
name, span = region.split(":")
a, b = (int(x.replace(",", "")) for x in span.split("-"))

seen = 0
out = []
inside = False
widths = set()
n_records = 0
with open(fa, "rb") as fh:
    for raw in fh:
        if raw.startswith(b">"):
            rec = raw[1:].split()[0].decode()
            if rec == name:
                n_records += 1
                inside = True
                seen = 0
            elif inside:
                break
            else:
                inside = False
            continue
        if not inside:
            continue
        line = raw.strip()
        widths.add(len(line))
        start, end = seen, seen + len(line)
        if end >= a and start < b:          # 1-based inclusive overlap
            lo = max(a - 1 - start, 0)
            hi = min(b - start, len(line))
            out.append(line[lo:hi].decode())
        seen = end
        if seen >= b:
            break

print(f">{region}")
print("".join(out))
print(f"# records named {name}: {n_records}", file=sys.stderr)
print(f"# line widths seen (last is the final short line): {sorted(widths)[-3:]}", file=sys.stderr)

"""Turn GLORI's pileup into a per-position survivor count for each strand.

GLORI writes one line per (position, strand) holding a comma-separated list of the bases it
kept. All twenty input reads cover all twenty positions, so a position where fewer than ten
reads survive on a strand is a position the trim removed.
"""
import sys
from collections import defaultdict

path, trim = sys.argv[1], int(sys.argv[2])
kept = defaultdict(lambda: {"+": 0, "-": 0})
for line in open(path):
    chrom, pos, strand, ref, bases, acounts = line.rstrip("\n").split("\t")
    kept[int(pos)][strand] = len(bases.split(","))

positions = range(1, 21)
print("  pos       " + "".join(f"{p:>3}" for p in positions))
for strand, label in (("+", "forward"), ("-", "reverse")):
    print(f"  {label:<10}" + "".join(f"{kept[p][strand]:>3}" for p in positions))

# The read is 20 bp and aligns 20M, so read offset == reference offset on the forward
# strand and 20 - offset on the reverse. Both strands should lose `trim` bases at each end.
def lost(strand):
    return sum(1 for p in positions if kept[p][strand] == 0)

print(f"\n  with --trim-head {trim} --trim-tail {trim}, {2 * trim} of 20 positions should be empty on each strand")
print(f"    forward: {lost('+')} empty")
print(f"    reverse: {lost('-')} empty")
if lost("+") != lost("-"):
    print(f"    -> STRAND ASYMMETRY: forward keeps {lost('-') - lost('+')} position(s) the reverse strand drops")

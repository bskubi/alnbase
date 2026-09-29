#!/usr/bin/env bash
# The .aref on-disk layout, decoded by hand, and how damaged or
# wrong-version indexes are refused.
set -uo pipefail
ALNBASE=${ALNBASE:-alnbase}
here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
cp "$here"/ref.fa "$work"; cd "$work"
run() { echo "\$ alnbase $*"; "$ALNBASE" "$@" 2>&1; echo "(exit $?)"; }

"$ALNBASE" index ref.fa ref.aref
echo "== ref.fa"; cat ref.fa
echo "== decode ref.aref"
python3 - <<'PY'
import struct
b = open('ref.aref', 'rb').read()
n = len(b)
n_contigs, version, magic = struct.unpack('<QQQ', b[n-24:])
print(f"file size {n}")
print(f"trailer: n_contigs={n_contigs} version={version} magic={magic.to_bytes(8,'little')!r}")
ENTRY = 48  # four u64s, then the contig's 16 MD5 bytes
table = n - 24 - ENTRY * n_contigs
print(f"data body: bytes [0,{table})  {b[:table]!r}")
for i in range(n_contigs):
    e = b[table+ENTRY*i:table+ENTRY*(i+1)]
    seq_size, seq_off, name_size, name_off = struct.unpack('<QQQQ', e[:32])
    name = b[name_off:name_off+name_size].decode()
    seq = list(b[seq_off:seq_off+seq_size])
    print(f"contig {i}: name={name!r}@{name_off} seq_off={seq_off} len={seq_size} "
          f"bytes={['0b{:08b}'.format(x) for x in seq]} md5={e[32:].hex()}")
PY
echo "(bit 0 A, 1 C, 2 G, 3 T, 4 GAP, 5 PAD, 6 CLIP, 7 SKIP; an index holds bases only: N=0b1111, R=A|G, Y=C|T)"
echo "== the stored digests are the @SQ M5 values samtools dict computes"
samtools dict ref.fa | awk -F'\t' '$1=="@SQ"{print $2, $4}'

echo "== version 4 (newer), versions 2 and 1 (older) and version 0 are all refused"
python3 -c "
import struct
b = bytearray(open('ref.aref','rb').read()); n = len(b)
for v in (4, 2, 1, 0):
    b[n-16:n-8] = struct.pack('<Q', v); open(f'v{v}.aref','wb').write(b)
b = bytearray(open('ref.aref','rb').read()); b[-1] ^= 0xff; open('badmagic.aref','wb').write(b)
open('truncated.aref','wb').write(open('ref.aref','rb').read()[:-30])
open('short.aref','wb').write(b'\0' * 10)
"
run info v4.aref
run info v2.aref
run info v1.aref
run info v0.aref
run info badmagic.aref
run info truncated.aref
run info short.aref
run info nonexistent.aref

echo "== an empty FASTA gives a valid, empty index (a 24-byte trailer)"
: > empty.fa
run index empty.fa empty.aref
run info empty.aref

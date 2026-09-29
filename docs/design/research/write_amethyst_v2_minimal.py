"""Minimal writer for an Amethyst v2.0.0 ("amethyst2.0.0") HDF5 file with h5py.

Layout (see singlecell-formats.md for citations):

    /metadata/version                  scalar string dataset = "amethyst2.0.0"
    /<context>/<barcode>/1             compound: chr S10, pos <i8, c <i8, t <i8   (bp-level obs)
    /<context>/<barcode>/<window_name> compound: chr S10, start <i8, end <i8,
                                                 c <i8, t <i8, c_nz <i8, t_nz <i8

Rules the readers actually rely on:
  * rows of /1 must be grouped by chr (contiguous) -- amethyst::indexChr stores
    (start=min(row index), count=N) per chr and later reads that hyperslab;
    sort by (chr, pos) like every writer does.
  * field names matter, field order does not (facet uses chr,pos,c,t; premethyst/facet
    calls2h5 use chr,pos,t,c; ScaleMethyl/legacy add pct). Readers access by name.
  * chr is fixed-width S10: names longer than 10 bytes are silently truncated by numpy.
  * windows: [start, end) half-open, only windows with >=1 observation written,
    sorted by chr,start,end; c/t = sums, c_nz/t_nz = #positions with c>0 / t>0.
  * compression: gzip level 6 (facet) or 9 (premethyst/ScaleMethyl); no explicit chunks.

Verified: run with h5py 3.6.0 / numpy 1.21.5 / HDF5 1.10.7 (see __main__ self-check).
Not verified against R/rhdf5 (rhdf5 not installed here).
"""
import numpy as np
import h5py

VERSION = "amethyst2.0.0"
OBS_DTYPE = [("chr", "S10"), ("pos", "<i8"), ("c", "<i8"), ("t", "<i8")]
WIN_DTYPE = [("chr", "S10"), ("start", "<i8"), ("end", "<i8"),
             ("c", "<i8"), ("t", "<i8"), ("c_nz", "<i8"), ("t_nz", "<i8")]


def write_obs(f, context, barcode, rows, name="1", level=6):
    arr = np.array(rows, dtype=OBS_DTYPE)
    assert all(len(r[0]) <= 10 for r in rows), "chr name > 10 bytes would be truncated"
    arr = np.sort(arr, order=["chr", "pos"])
    f.create_dataset(f"/{context}/{barcode}/{name}", data=arr,
                     compression="gzip", compression_opts=level)
    return arr


def uniform_windows(obs, size, offset=0):
    """Non-overlapping [start, start+size) windows, start = (pos-offset)//size*size+offset."""
    start = (obs["pos"] - offset) // size * size + offset
    keys = {}
    for chr_, s, c, t in zip(obs["chr"], start, obs["c"], obs["t"]):
        k = (chr_, int(s))
        a = keys.setdefault(k, [0, 0, 0, 0])
        a[0] += c; a[1] += t; a[2] += int(c > 0); a[3] += int(t > 0)
    out = [(k[0], k[1], k[1] + size, *v) for k, v in sorted(keys.items())]
    return np.array(out, dtype=WIN_DTYPE)


def main(path="minimal_amethyst_v2.h5"):
    with h5py.File(path, "w") as f:
        f.create_dataset("/metadata/version", data=VERSION)  # variable-length UTF-8 scalar
        cells = {
            "CATCTG": [("chr1", 16082, 1, 0), ("chr1", 16139, 1, 0), ("chr2", 10, 0, 1),
                       ("chr1", 116141, 0, 1)],
            "GGATCC": [("chr1", 500, 2, 1)],
        }
        for bc, rows in cells.items():
            for ctx in ("CG", "GCH"):  # any context name is just a group name
                obs = write_obs(f, ctx, bc, rows)
                win = uniform_windows(obs, 100000)
                f.create_dataset(f"/{ctx}/{bc}/100000", data=win,
                                 compression="gzip", compression_opts=6)

    # self-check
    with h5py.File(path, "r") as f:
        assert f["/metadata/version"].asstr()[()] == VERSION
        d = f["/CG/CATCTG/1"][:]
        assert d.dtype.names == ("chr", "pos", "c", "t")
        assert list(d["chr"]) == [b"chr1", b"chr1", b"chr1", b"chr2"]
        w = f["/CG/CATCTG/100000"][:]
        print(d)
        print(w)
    print("ok:", path)


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])

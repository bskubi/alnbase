#!/usr/bin/env python3
"""Run one named simulation profile and record what produced the output.

    simulate.py lambda-illumina --out $SCRATCH/sims/lambda-illumina

Profiles live in ../sims/profiles.toml. The point of naming them is
attribution: every figure the validation suite reports depends on the
simulation behind it, so the output directory carries a provenance.json saying
which profile, which parameters and which reference it came from, and the
comparison scripts can carry that name into their own reports.

Needs `bsreadsim` on PATH (the `bsreadsim` micromamba environment).
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 and older
    import tomli as tomllib

PROFILES = Path(__file__).resolve().parent.parent / "sims" / "profiles.toml"


def read_fasta(fasta):
    """The first contig's name, and the digest BSReadSim records for the bases.

    Checked rather than assumed because pointing a profile at the wrong copy of
    a reference would silently invalidate everything measured from it, without
    failing. The digest covers the uppercased bases alone, so it is the same
    whatever the header says or however the lines are wrapped -- which is why
    the contig name is checked separately: a correct sequence under a different
    name still fails to line up with everything already measured.
    """
    contig, digest = None, hashlib.sha256()
    with open(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if contig is None:
                    contig = line[1:].split()[0].decode()
            else:
                digest.update(line.strip().upper())
    return contig, digest.hexdigest()


def simulation_sha256(profile, ref):
    """An identity for the simulation that two machines can agree on.

    BSReadSim's own `configuration_sha256` covers the output directory path, so
    the same simulation run in two places gets two hashes and the hash cannot be
    used to ask "is this the input those numbers came from?". This covers the
    sequence, the technology and the parameters, and nothing else.

    Two runs sharing this hash produce the same reads. They do not produce
    byte-identical BAMs: BSReadSim stamps each run with a fresh UUID read group.
    """
    identity = {
        "technology": profile["technology"],
        "reference_sequence_sha256": ref["sequence_sha256"],
        "args": dict(sorted(profile["args"].items())),
    }
    canonical = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def load_profile(name):
    with open(PROFILES, "rb") as fh:
        conf = tomllib.load(fh)
    if name not in conf["profile"]:
        known = ", ".join(sorted(conf["profile"]))
        sys.exit(f"no profile `{name}` in {PROFILES}\nknown profiles: {known}")
    profile = conf["profile"][name]
    return profile, conf["reference"][profile["reference"]]


def find_reference(ref, name, given, ref_dir):
    if given:
        path = Path(given)
    else:
        path = Path(ref_dir) / f"{name}.fa"
    if not path.exists():
        sys.exit(
            f"reference `{name}` not found at {path}\n"
            f"fetch {ref['accession']} and name its contig `{ref['contig']}`:\n"
            f"  curl -sL '{ref['url']}' \\\n"
            f"    | sed '1s/.*/>{ref['contig']}/' > {path}\n"
            f"(verified: that command reproduces the pinned digest)"
        )
    contig, seen = read_fasta(path)
    if seen != ref["sequence_sha256"]:
        sys.exit(
            f"{path} is not the pinned `{name}` sequence\n"
            f"  expected {ref['sequence_sha256']}\n"
            f"  found    {seen}"
        )
    if contig != ref["contig"]:
        sys.exit(
            f"{path} holds the right sequence under the wrong name\n"
            f"  expected contig `{ref['contig']}`, found `{contig}`\n"
            f"rename it: sed -i '1s/.*/>{ref['contig']}/' {path}"
        )
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("profile")
    ap.add_argument("--out", required=True, help="new output directory")
    ap.add_argument("--reference", help="FASTA to use instead of <ref-dir>/<name>.fa")
    ap.add_argument("--ref-dir", default="ref", help="where references live (default: ./ref)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true", help="print the command and stop")
    args = ap.parse_args()

    profile, ref = load_profile(args.profile)
    fasta = find_reference(ref, profile["reference"], args.reference, args.ref_dir)

    out = Path(args.out)
    if out.exists():
        sys.exit(f"{out} exists; BSReadSim requires a new directory")

    cmd = ["bsreadsim", "run", profile["technology"],
           "--reference", str(fasta), "--output", str(out)]
    for key, value in profile["args"].items():
        cmd += [f"--{key}", str(value)]
    cmd += ["--threads", str(args.threads), "--format", "bam", "--save-truth"]

    if args.dry_run:
        print(" ".join(cmd))
        return
    if not shutil.which("bsreadsim"):
        sys.exit("bsreadsim is not on PATH; activate the `bsreadsim` environment")

    print(" ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)

    # BSReadSim's manifest is the authority on what actually ran; provenance.json
    # only adds what it cannot know, namely which named profile asked for it.
    manifest = json.loads((out / "sim.manifest.json").read_text())
    provenance = {
        "profile": args.profile,
        "summary": profile["summary"],
        "simulation_sha256": simulation_sha256(profile, ref),
        "profiles_toml_sha256": hashlib.sha256(PROFILES.read_bytes()).hexdigest(),
        "reference": {
            "name": profile["reference"],
            "accession": ref["accession"],
            "path": str(fasta.resolve()),
            "sequence_sha256": ref["sequence_sha256"],
        },
        "args": profile["args"],
        "command": cmd,
        "bsreadsim_configuration_sha256": manifest["details"]["configuration_sha256"],
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"{out}/provenance.json written for profile `{args.profile}`", file=sys.stderr)


if __name__ == "__main__":
    main()

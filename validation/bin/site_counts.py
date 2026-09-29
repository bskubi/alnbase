#!/usr/bin/env python3
"""Turn each tool's per-site methylation output into one common table.

    site_counts.py methyldackel PREFIX       OUT.parquet   # PREFIX_CpG.bedGraph, _CHG, _CHH
    site_counts.py bismark-cov  FILE.cov.gz  OUT.parquet   # CG only (bismark2bedGraph default)
    site_counts.py calls        CALLS.parquet OUT.parquet          # from call_observations.py
    site_counts.py truth        OBSERVATIONS.parquet OUT.parquet   # from truth_observations.py

Every output has the same columns:

    context   CG, CHG or CHH
    contig
    pos       0-based position of the cytosine (for a bottom-strand cytosine, the G
              on the top strand)
    m, u      methylated and unmethylated observations

Conventions of the inputs, which this script is the one place to know about:
  - MethylDackel bedGraph: a track line, then chrom, 0-based start, end, percent, m, u.
  - Bismark .cov: chrom, 1-based start, end, percent, m, u; no strand column.
  - Per-read calls from call_observations.py (alnbase hit rows, Bismark extractor
    call files); calls with kept = false are not counted.
  - BSReadSim truth: the observations listed by truth_observations.py, counted by site.
"""

import argparse

import duckdb


def from_methyldackel(prefix: str):
    parts = [
        f"""select '{ctx}' as context, column0 as contig, column1 as pos, column4 as m, column5 as u
            from read_csv('{prefix}_{name}.bedGraph', skip=1, header=false, delim='\t')"""
        for ctx, name in (("CG", "CpG"), ("CHG", "CHG"), ("CHH", "CHH"))
    ]
    return duckdb.sql(" union all ".join(parts))


def from_bismark_cov(path: str):
    return duckdb.sql(f"""
        select 'CG' as context, column0 as contig, column1 - 1 as pos, column4 as m, column5 as u
        from read_csv('{path}', header=false, delim='\t')""")


def from_calls(calls: str):
    return duckdb.sql(f"""
        select context, contig, pos, count(*) filter (where methylated) as m,
               count(*) filter (where not methylated) as u
        from '{calls}' where kept group by all""")


def from_truth(observations: str):
    return duckdb.sql(f"""
        select context, contig, pos, count(*) filter (where methylated) as m,
               count(*) filter (where not methylated) as u
        from '{observations}' where context is not null group by all""")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("kind", choices=["methyldackel", "bismark-cov", "calls", "truth"])
    p.add_argument("input")
    p.add_argument("out")
    a = p.parse_args()
    rel = {
        "methyldackel": lambda: from_methyldackel(a.input),
        "bismark-cov": lambda: from_bismark_cov(a.input),
        "calls": lambda: from_calls(a.input),
        "truth": lambda: from_truth(a.input),
    }[a.kind]()
    rel.order("contig, pos, context").write_parquet(a.out)
    print(duckdb.sql(f"select context, count(*) sites, sum(m) m, sum(u) u from '{a.out}' group by 1 order by 1"))


if __name__ == "__main__":
    main()

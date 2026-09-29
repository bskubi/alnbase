#!/usr/bin/env python3
"""Turn each caller's per-read methylation calls into one common table.

    call_observations.py alnbase 'DIR/*.parquet' OUT.parquet [--keep SQL]
    call_observations.py alnbase-bismark 'DIR/*.parquet' OUT.parquet
    call_observations.py bismark-extractor DIR OUT.parquet --bam ALIGNED.bam

Every output has one row per call, with the same columns:

    qname, mate (1 or 2), contig
    pos         0-based position of the cytosine (for a bottom-strand cytosine, the G
                on the top strand)
    context     CG, CHG or CHH
    methylated
    kept        false for a call the caller's own filters would discard (see --keep)
    aln_pos, aln_reverse, aln_cigar    where the read was aligned

Conventions of the inputs, which this script is the one place to know about:
  - alnbase hit rows from a query file whose query names are `<context>_methylated`
    or `<context>_unmethylated` (presets/methyldackel/queries.toml), extracted with
    `-F qname,is_last_in_template,pos,is_reverse,cigar` plus any columns --keep uses.
    `CpG` in a query name is read as CG. --keep is a DuckDB expression over the hit
    columns, for a caller's filters that are not part of the call itself (MAPQ, flags,
    base quality; presets/methyldackel/recipe.toml, [filter.hit]). Calls failing it are kept in
    the table with kept = false, so that scoring can say which truth they would match.
  - alnbase hit rows from presets/bismark/queries.toml, whose query names are Bismark's
    call letters (Z/z CG, X/x CHG, H/h CHH, U/u unknown context; upper case
    methylated), with the same -F columns. U/u calls are left out, as the Bismark
    extractor leaves them out.
  - Bismark methylation extractor call files (CpG_*.txt, CHG_*.txt, CHH_*.txt): a
    version line, then read name, `+` methylated / `-` unmethylated, contig, 1-based
    position, call letter. They name neither the mate nor the alignment, so both come
    from the BAM the extractor read: the mate is the one whose aligned interval contains
    the call (mates of the pairs compared here do not overlap).
"""

import argparse

import duckdb
import pyarrow as pa
import pysam


def from_alnbase(con, glob: str, keep: str):
    return con.sql(f"""
        select qname, if(is_last_in_template, 2, 1) as mate, ref_name as contig, refr_pos as pos,
               replace(regexp_extract(name, '^([A-Za-z]+)_', 1), 'CpG', 'CG') as context,
               name like '%\\_methylated' escape '\\' as methylated,
               coalesce(({keep}
               ), false) as kept,
               pos as aln_pos, is_reverse as aln_reverse, cigar as aln_cigar
        from read_parquet('{glob}') where name is not null""")


def from_alnbase_bismark(con, glob: str):
    return con.sql(f"""
        select qname, if(is_last_in_template, 2, 1) as mate, ref_name as contig, refr_pos as pos,
               case upper(name) when 'Z' then 'CG' when 'X' then 'CHG' when 'H' then 'CHH' end as context,
               name = upper(name) as methylated, true as kept,
               pos as aln_pos, is_reverse as aln_reverse, cigar as aln_cigar
        from read_parquet('{glob}') where upper(name) in ('Z', 'X', 'H')""")


def from_bismark_extractor(con, directory: str, bam: str):
    records = {k: [] for k in ("qname", "mate", "contig", "aln_pos", "aln_end", "aln_reverse", "aln_cigar")}
    with pysam.AlignmentFile(bam) as f:
        for r in f:
            if r.is_unmapped:
                continue
            row = (r.query_name, 2 if r.is_read2 else 1, r.reference_name, r.reference_start, r.reference_end,
                   r.is_reverse, r.cigarstring)
            for col, value in zip(records.values(), row):
                col.append(value)
    con.register("records", pa.table(records))
    return con.sql(f"""
        with calls as (
            select column0 as qname, column2 as contig, column3 - 1 as pos,
                   case upper(column4) when 'Z' then 'CG' when 'X' then 'CHG' when 'H' then 'CHH' end as context,
                   column1 = '+' as methylated
            from read_csv(['{directory}/CpG_*.txt', '{directory}/CHG_*.txt', '{directory}/CHH_*.txt'],
                          skip=1, header=false, delim='\t',
                          columns={{'column0': 'varchar', 'column1': 'varchar', 'column2': 'varchar',
                                    'column3': 'bigint', 'column4': 'varchar'}}))
        select qname, mate, contig, pos, context, methylated, true as kept, aln_pos, aln_reverse, aln_cigar
        from calls join records using (qname, contig)
        where pos >= aln_pos and pos < aln_end""")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("kind", choices=["alnbase", "alnbase-bismark", "bismark-extractor"])
    p.add_argument("input")
    p.add_argument("out")
    p.add_argument("--keep", default="true", help="alnbase: DuckDB expression a call must satisfy")
    p.add_argument("--bam", help="bismark-extractor: the BAM the extractor read")
    a = p.parse_args()

    con = duckdb.connect()
    if a.kind == "alnbase":
        rel = from_alnbase(con, a.input, a.keep)
    elif a.kind == "alnbase-bismark":
        rel = from_alnbase_bismark(con, a.input)
    else:
        if not a.bam:
            p.error("bismark-extractor needs --bam")
        rel = from_bismark_extractor(con, a.input, a.bam)
    rel.order("contig, pos, qname, mate").write_parquet(a.out)
    print(con.sql(f"""select context, count(*) calls, count(*) filter (where kept) kept,
                             count(distinct (qname, mate, pos)) distinct_observations
                      from '{a.out}' group by 1 order by 1"""))


if __name__ == "__main__":
    main()

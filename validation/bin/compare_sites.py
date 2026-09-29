#!/usr/bin/env python3
"""Compare two per-site count tables (from site_counts.py) site by site.

    compare_sites.py --label NAME A.parquet B.parquet OUT_PREFIX

Joins on context, contig and position, and classifies every site found in either:

    agree        both tables have the site, with the same m and u
    count_diff   both have it, with different counts
    a_only       only A has it
    b_only       only B has it

Writes OUT_PREFIX.sites.parquet (every site, both tables' counts, the category) and
OUT_PREFIX.metrics.parquet (one row per context and metric), and prints a summary.
Methylation-level agreement is also reported, as the mean absolute difference in
m / (m + u) over sites both tables have.
"""

import argparse

import duckdb


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--label", required=True, help="name of this comparison, e.g. methyldackel-vs-alnbase")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("out_prefix")
    args = p.parse_args()

    con = duckdb.connect()
    con.sql(f"""
        create table sites as
        select coalesce(a.context, b.context) as context, coalesce(a.contig, b.contig) as contig,
               coalesce(a.pos, b.pos) as pos,
               a.m as a_m, a.u as a_u, b.m as b_m, b.u as b_u,
               case when a.pos is null then 'b_only'
                    when b.pos is null then 'a_only'
                    when a.m = b.m and a.u = b.u then 'agree'
                    else 'count_diff' end as category
        from '{args.a}' a full outer join '{args.b}' b using (context, contig, pos)""")
    con.sql(f"copy (select '{args.label}' as comparison, * from sites order by contig, pos, context) "
            f"to '{args.out_prefix}.sites.parquet'")

    con.sql("""
        create table per_context as
        select context,
               count(*)::double as sites,
               count(*) filter (where category = 'agree')::double as agree,
               count(*) filter (where category = 'count_diff')::double as count_diff,
               count(*) filter (where category = 'a_only')::double as a_only,
               count(*) filter (where category = 'b_only')::double as b_only,
               coalesce(sum(abs(a_m - b_m) + abs(a_u - b_u)) filter (where category = 'count_diff'), 0)::double
                   as abs_count_difference,
               avg(abs(a_m / (a_m + a_u) - b_m / (b_m + b_u)))
                   filter (where a_m + a_u > 0 and b_m + b_u > 0) as mean_abs_level_difference
        from sites group by context""")
    metrics = con.sql("unpivot per_context on columns(* exclude (context)) into name metric value value")
    con.sql(f"copy (select '{args.label}' as comparison, context, metric, value::double as value "
            f"from metrics) to '{args.out_prefix}.metrics.parquet'")
    print(args.label)
    print(con.sql(f"""pivot (select * from '{args.out_prefix}.metrics.parquet') on metric using first(value)
                      group by context order by context"""))


if __name__ == "__main__":
    main()

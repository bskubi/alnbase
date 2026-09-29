#!/usr/bin/env python3
"""Score a caller's per-read methylation calls against BSReadSim truth, call by call.

    score_against_truth.py TRUTH.parquet CALLS.parquet OUT_PREFIX

TRUTH.parquet comes from truth_observations.py, CALLS.parquet from
call_observations.py. Calls with kept = false (discarded by the caller's own filters)
count as not called.

Calls and truth are joined on read name, mate and position. Each row gets a category:

    agree              same context, same methylation state
    state_differs      same context, methylation state differs
    context_differs    both have a cytosine here, in different contexts
    call_only          the caller has a call where the truth has no measured cytosine
    truth_only         a measured cytosine the caller did not call

and a likely cause, from the truth bits and from where the read was placed:

    filtered           the caller's own filters discarded this call (MAPQ, flags, base quality)
    wrong_locus        the aligner placed the read on another strand or more than 20 bases
                       from where it was simulated
    realigned          the aligner placed the read at its simulated locus, but with its
                       start or CIGAR (indels, ends) differing from the simulation
    sequencing_error   the truth marks a sequencing error at this base
    not_converted      an unmethylated cytosine that bisulfite did not convert
    variant            a simulated variant changed this base
    variant_nearby     a simulated variant within 2 bases changed the context
    deletion_nearby    a deletion in the read within 2 bases changed the context
    insertion_nearby   an insertion in the read within 2 bases changed the context
    none               none of the above

The causes are checked in that order and the first that applies is reported.

Writes OUT_PREFIX.calls.parquet (every row) and OUT_PREFIX.summary.parquet (counts by
context, category and cause), and prints the summary.
"""

import argparse

import duckdb


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("truth")
    p.add_argument("calls")
    p.add_argument("out_prefix")
    a = p.parse_args()

    con = duckdb.connect()
    con.sql(f"create table calls as select * from '{a.calls}'")
    con.sql(f"create table truth as select * from '{a.truth}'")
    # Where each read was simulated and where it was aligned, to recognise misplaced reads.
    con.sql("""
        create table placement as
        select qname, mate,
               case when sim_pos = aln_pos and sim_reverse = aln_reverse and sim_cigar = aln_cigar then 'exact'
                    when sim_reverse = aln_reverse and abs(sim_pos - aln_pos) <= 20 then 'realigned'
                    else 'wrong_locus' end as placement
        from (select distinct qname, mate, sim_pos, sim_reverse, sim_cigar from truth)
        join (select distinct qname, mate, aln_pos, aln_reverse, aln_cigar from calls) using (qname, mate)""")
    con.sql("""
        create table scored as
        select qname, mate, contig, pos, c.context as call_context, t.context as truth_context,
               c.methylated as call_methylated, t.methylated as truth_methylated,
               c.kept, t.converted, t.variant, t.error, t.variant_nearby, t.deletion_nearby, t.insertion_nearby, placement,
               case when c.context is null then 'truth_only'
                    when t.context is null then 'call_only'
                    when c.context <> t.context then 'context_differs'
                    when c.methylated <> t.methylated then 'state_differs'
                    else 'agree' end as category,
               case when not coalesce(c.kept, true) then 'filtered'
                    when placement in ('wrong_locus', 'realigned') then placement
                    when t.error then 'sequencing_error'
                    when not t.methylated and not t.converted then 'not_converted'
                    when t.variant then 'variant'
                    when t.variant_nearby then 'variant_nearby'
                    when t.deletion_nearby then 'deletion_nearby'
                    when t.insertion_nearby then 'insertion_nearby'
                    else 'none' end as cause
        from (select qname, mate, contig, pos, if(kept, context, null) as context, methylated, kept from calls) c
        full outer join truth t using (qname, mate, contig, pos)
        left join placement using (qname, mate)
        where c.context is not null or t.context is not null""")
    con.sql(f"copy scored to '{a.out_prefix}.calls.parquet'")
    con.sql(f"""
        copy (select coalesce(truth_context, call_context) as context, category, cause, count(*) as n
              from scored group by all order by all) to '{a.out_prefix}.summary.parquet'""")
    print(con.sql(f"""
        pivot (select category, cause, n from '{a.out_prefix}.summary.parquet')
        on category using sum(n) group by cause order by cause"""))


if __name__ == "__main__":
    main()

"""Handing a writer one group of rows at a time.

Almost every output format is a `COPY ... TO`, and needs no Python at all. The
exception is a format that is not one file: an Amethyst H5 holds one dataset
per block, and a per-block folder of `.cov` files holds one file per block. Those
need a loop, and a loop that pulls a whole-genome result into memory first is
not a loop anybody can run.

So the division is: the format's `.sql` file does the aggregation *and the
sort*, declaring which column it is cut on; this module walks the sorted stream
and yields one group at a time. Nothing larger than the group currently open is
ever held, and the slices are zero copy.

    con.execute("SELECT ... ORDER BY block, chr, pos")
    for (block,), table in groups(con.to_arrow_reader(8192), ["block"]):
        write(block, table)

The sort is the whole contract, and it is checked twice. `check_sorted` reads
the query plan before a byte is written, so a format whose `ORDER BY` does not
lead with its cut key fails at preflight rather than after half its files are
on disk. `groups` then checks again while streaming, because a plan can be read
correctly and still be about the wrong query, and because a key reappearing
after its group closed would otherwise mean a silently truncated output.
"""

from __future__ import annotations

import json
import re


class NotGrouped(Exception):
    """A group reappeared after it closed: the query did not sort by the key."""


def groups(reader, keys, *, batch_rows=8192):
    """Yield `(key tuple, arrow Table)` from a stream ordered by `keys`.

    The stream must be sorted so that each group's rows are contiguous. Groups
    are found by watching the key columns change, so a group may span any
    number of Arrow batches and one batch may hold the end of one group and the
    start of the next.

    `batch_rows` is accepted for symmetry with `to_arrow_reader` and is unused
    here: the reader decides its own batch size, and this walks whatever it is
    handed.
    """
    if isinstance(keys, str):
        keys = [keys]
    import pyarrow as pa

    seen: set = set()
    held: list = []
    current = None

    def key_at(batch, i, idx):
        return tuple(batch.column(j)[i].as_py() for j in idx)

    def close(key, parts):
        if key in seen:
            raise NotGrouped(
                f"group {key} reappeared after it closed; the query must "
                f"ORDER BY {', '.join(keys)} first"
            )
        seen.add(key)
        return key, pa.Table.from_batches(parts)

    for batch in reader:
        idx = [batch.schema.get_field_index(k) for k in keys]
        if -1 in idx:
            missing = [k for k, j in zip(keys, idx) if j == -1]
            raise KeyError(f"cut key(s) {missing} not in the query's output")
        start = 0
        for i in range(len(batch)):
            key = key_at(batch, i, idx)
            if current is None:
                current = key
            elif key != current:
                held.append(batch.slice(start, i - start))
                yield close(current, held)
                held, current, start = [], key, i
        if start < len(batch):
            held.append(batch.slice(start))
    if held:
        yield close(current, held)


def check_sorted(con, sql: str, keys) -> list:
    """Refuse a streaming format whose `ORDER BY` does not lead with `keys`.

    Reads the sort keys out of the optimized plan. The plan names source
    columns rather than output aliases, so a format must not rename its cut
    key: a query selecting `block AS barcode` and cut on `barcode` is refused
    here, which is a better outcome than discovering it at row one.
    """
    if isinstance(keys, str):
        keys = [keys]
    con.execute("SET explain_output='optimized_only'")
    plan = json.loads(con.execute(f"EXPLAIN (FORMAT json) {sql}").fetchall()[0][1])

    found: list = []

    def bare(key: str) -> str:
        """`db.schema.table.column` -> `column`.

        The qualifier is of any depth: a parquet scan gives `sites.block`, a base
        table `memory.main.hit.block`. An expression matches nothing and is left
        exactly as the plan gave it, so it can never accidentally equal a
        declared cut key.
        """
        key = key.removesuffix(" ASC").removesuffix(" DESC").strip()
        # DuckDB quotes a part that collides with one of its keywords, so a
        # column named `sample` comes back as `s."sample"`. The quotes are
        # spelling, not identity.
        m = re.fullmatch(r'(?:"?\w+"?\.)*"?(\w+)"?', key)
        return m.group(1) if m else key

    def walk(node):
        if node.get("name") == "ORDER_BY":
            found.extend(bare(k) for k in node.get("extra_info", {}).get("Order By", []))
        for child in node.get("children", []):
            walk(child)

    walk(plan if isinstance(plan, dict) else plan[0])
    if not found:
        raise ValueError(f"streaming format has no ORDER BY; it needs {keys}")
    if found[: len(keys)] != list(keys):
        raise ValueError(f"ORDER BY leads with {found[: len(keys)]}, not {list(keys)}")
    return found

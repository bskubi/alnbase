#!/usr/bin/env python3
"""Generate an alnbase query file that reproduces Bismark's XM methylation call.

    python generate.py > queries.toml

Run alnbase with `--insertions emit` (Bismark sees inserted bases; see rule 3).

Bismark (v0.25.1) makes its call at alignment time (`bismark`,
`extract_corresponding_genomic_sequence_single_end` and `methylation_call`;
alnbase-validation/docs/research/methylation-tools.md §1.2). In alnbase's walk every read is
oriented along its conversion strand, so a cytosine call is always "reference C,
read C (methylated) or T (unconverted)", and its context is the columns after it.

Bismark's context rules, and how each is written here:

 1. The context is the next two characters of the read-projected reference: the
    reference bases the read aligns to, in read order, plus two reference bases
    past the end of the alignment. alnbase supplies the bases past the end as pad
    columns, so a context column may be a read base or a pad (alias `o`).
 2. A deleted reference base is skipped, not used as context. alnbase shows it as a
    deletion column (`.` on the read side). So each context position may be
    preceded by 0..MAX_DELETION deletion columns, and every combination is a
    separate pattern in the query. Longer deletions are not covered.
 3. An inserted (or soft-clipped) read base contributes an unknown character `X`.
    With `--insertions emit` alnbase shows it as a column with a gap on the
    reference side, which matches neither G nor H below, so the call becomes U/u.
    (Soft clips produce no alnbase columns; Bismark aligns end-to-end by default.)
 4. Next character G: CG (Z/z). N or X: unknown (U/u). Otherwise the character
    after it decides: G is CHG (X/x), N or X is unknown (U/u), anything else is
    CHH (H/h). Here "anything else" is H (A, C or T); an IUPAC ambiguity code in the
    reference is treated as unknown, where Bismark would use it as H.

All patterns in a query must have the same width, so shorter variants are padded
with trailing `~` columns, which match anything.
"""

MAX_DELETION = 4
WIDTH = 3 + 2 * MAX_DELETION
STATES = {"methylated": "C", "unconverted": "T"}
# Bismark's XM letters: uppercase methylated, lowercase unconverted.
LETTERS = {"CG": "Z", "CHG": "X", "CHH": "H", "unknown": "U"}


def pattern(read_cells: list[str], refr_cells: list[str]) -> tuple[str, str]:
    """Pad two rows of cells to WIDTH with columns that match anything."""
    pad = WIDTH - len(read_cells)
    return "".join(read_cells) + "~" * pad, "".join(refr_cells) + "~" * pad


def cytosine_then_context(base: str, context: list[tuple[int, str]]) -> tuple[str, str]:
    """A cytosine read as `base`, then each context column preceded by deletions.

    `context` is a list of (deletions before this context column, reference code).
    """
    read, refr = [base], ["C"]
    for deletions, code in context:
        read += ["."] * deletions + ["o"]
        refr += ["~"] * deletions + [code]
    return pattern(read, refr)


def section(title: str) -> str:
    return f"\n# {'-' * 70}\n# {title}\n# {'-' * 70}\n"


def main() -> None:
    out = [__doc__.replace("\n", "\n# ").join(["# ", ""]).rstrip("# \n") + "\n"]
    out.append('\n[alias]\no = "{N_}"   # a read base or a pad past the read end: not a deletion\n')
    deletion_counts = range(MAX_DELETION + 1)

    for state, base in STATES.items():
        out.append(section(f"cytosines read as {base} ({state})"))
        shapes = {
            # CG: the first context column is G; the second does not matter.
            "CG": [[(d1, "G")] for d1 in deletion_counts],
            "CHG": [[(d1, "H"), (d2, "G")] for d1 in deletion_counts for d2 in deletion_counts],
            "CHH": [[(d1, "H"), (d2, "H")] for d1 in deletion_counts for d2 in deletion_counts],
            # Any C with two context columns, whatever they hold.
            "any": [[(d1, "~"), (d2, "~")] for d1 in deletion_counts for d2 in deletion_counts],
        }
        names = {}
        for context, variants in shapes.items():
            names[context] = []
            for i, variant in enumerate(variants):
                name = f"{context}_{base}_{i}"
                read, refr = cytosine_then_context(base, variant)
                deletions = ", ".join(str(d) for d, _ in variant)
                out.append(f'[pat.{name}]   # deletions before each context column: {deletions}\n'
                           f'read = "{read}"\nrefr = "{refr}"\n')
                names[context].append(name)

        for context in ("CG", "CHG", "CHH"):
            letter = LETTERS[context] if base == "C" else LETTERS[context].lower()
            out.append(f"[query.{letter}]   # {context}, {state}\n"
                       f'where = "{" or ".join(names[context])}"\n')
        letter = "U" if base == "C" else "u"
        known = " or ".join(names["CG"] + names["CHG"] + names["CHH"])
        out.append(f"[query.{letter}]   # context unknown (N or an inserted base), {state}\n"
                   f'where = "({" or ".join(names["any"])}) and not ({known})"\n')

    out.append(section("the tag: Bismark's letters, '.' where no call is made"))
    out.append('[tag.xm.bases]\nfill = "."\n' +
               "".join(f'{c} = "{c}"\n' for c in "ZzXxHhUu"))
    print("".join(out), end="")


if __name__ == "__main__":
    main()

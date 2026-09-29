# IUPAC ambiguity codes in the reference, beyond N

`n_context` pins that `N` in the reference leaves a context undetermined. `N` is only the
widest of the fifteen IUPAC codes, and the rule that produces that answer is not "N is
special" but subset matching: a reference value matches a pattern column only when its base
set is a **subset** of the column's set (`docs/reference/01-reference-and-codes.md` §1.1,
§4). So some ambiguity codes resolve a context and others do not, and which is which
follows from the sets alone.

`chrP` is five reference triplets followed by nine A, read by one 24M record whose own base
is `A` under every ambiguous position:

| reference | the ambiguous base | ⊆ `H` = {A,C,T}? | call | 0-based |
|---|---|---|---|---|
| `CYG` | `Y` = {C,T} | yes | `CHG_met` | 0 |
| `CRG` | `R` = {A,G} | no | `Cx_met` | 3 |
| `CYT` | `Y` = {C,T} | yes | `CHH_met` | 6 |
| `CSA` | `S` = {C,G} | no | `Cx_met` | 9 |
| `CWG` | `W` = {A,T} | yes | `CHG_met` | 12 |

The two that fail do so twice over: `R` and `S` are no more subsets of `G` than of `H`, so
neither a CG nor a CHG nor a CHH can be claimed and the context falls to `Cx`. `K`, `B`,
`D`, `V` and `N` behave the same way for the same reason, and `M` = {A,C} behaves like `Y`
and `W`.

This is a deliberate choice, not an inevitability: a caller could reasonably refuse every
ambiguous base, as treating `N` specially amounts to. Subset matching instead calls the
context exactly when the reference constrains it enough to, which is the same rule that
governs ambiguity in SEQ, and the ambiguous base is never itself called — `S` includes `C`
but is not a subset of `C`, so no cytosine is claimed at reference 11.

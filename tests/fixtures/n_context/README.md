# An N in the reference next to a cytosine

`chrN` puts `CAG`, `CAA`, `CG`, `CNG` and `CNN` next to one another. The first three are
controls that must come out CHG, CHH and CG. In the last two the reference does not say
what the context is, and `context.toml` is IUPAC-strict -- N is not H -- so both are `Cx`.

The choice is real and worth pinning rather than assuming: MethylDackel treats N as H and
calls the same two cytosines CHG and CHH, and Bismark calls them U, its own unknown code.
Assembly gaps and telomeric N runs put long stretches of N in every real reference, so a
comparison that does not know which rule each side used will see the difference and
mistake it for an error.

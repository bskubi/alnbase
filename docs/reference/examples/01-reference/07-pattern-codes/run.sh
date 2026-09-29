#!/usr/bin/env bash
# How pattern rows in a query file are parsed into code sets, checked with
# `query --trace READ@REFR` (a synthetic one-column-per-character alignment;
# no BAM or index needed). Each case prints the query file, the trace input,
# and whether the query matched or failed to compile.
set -uo pipefail
ALNBASE=${ALNBASE:-alnbase}
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT; cd "$work"
case_() { # case_ "description" "toml" "READ@REFR"
  printf '%s\n' "$2" > q.toml
  echo "=== $1"; sed 's/^/    /' q.toml; echo "    trace: $3"
  "$ALNBASE" query --query-file q.toml --trace "$3" --trace-grid false 2>&1 \
    | grep -E "matches ending at|no match|rror|warning" | sed 's/^ */    -> /' | head -3
}

case_ "alternative spellings: U=T, '-'=GAP on the pattern side; Z=GAP, X=PAD" \
'[query.q]
read = "UZ"
refr = "-X"' 'T.@._'

# A query fires on any window it matches, junction columns alone included, so the
# junction cases below need no aligned base (0.1.3 to 0.1.6 required one).
case_ "junction: ',' alone and 'J' are the same code" \
'[query.q]
read = "J,"
refr = "~~"' ',,@AA'

case_ "clip: ':' alone and 'L' are the same code; it matches a clip observation" \
'[query.q]
read = "L:"
refr = "~~"' '::@AA'

case_ "~ includes a clip" \
'[query.q]
read = "~"
refr = "~"' ':@A'

case_ "% is not a code (it was ERR, removed in 0.1.10)" \
'[query.q]
read = "%"
refr = "~"' 'A@A'

case_ "observed lowercase/U in a trace are folded (the trace decoder is case-insensitive)" \
'[query.q]
read = "N"
refr = "N"' 'U@a'

case_ "a group in a row counts its braces as columns for the row-width check" \
'[query.q]
read = "{C.}~"
refr = "CG"' 'CG@CG'

case_ "... so a group only parses when the other row is as many CHARACTERS wide" \
'[query.q]
read = "{C.}"
refr = "{NJ}"' '.@A'

case_ "lowercase is accepted INSIDE a group ({c.} = {C.})" \
'[query.q]
read = "{c.}"
refr = "{~~}"' 'C@A'

case_ "the practical form: an alias (values are parsed like group bodies, braces optional)" \
'[alias]
f = "c."
[query.q]
read = "f~"
refr = "CG"' '.G@CG'

case_ "TRAP: an alias letter inside braces is NOT the alias; {j} is J = junction" \
'[alias]
j = "A"
[query.q]
read = "{j}"
refr = "{~}"' ',@A'

case_ "lowercase outside a group must be a declared alias" \
'[query.q]
read = "c"
refr = "C"' 'C@C'

case_ "uppercase letters cannot be aliases" \
'[alias]
C = "{C.}"
[query.q]
read = "C"
refr = "~"' 'C@C'

case_ "E is not a code" \
'[query.q]
read = "E"
refr = "C"' 'C@C'

case_ "digits are rejected in rows, so '0' (empty set) and counts are unusable" \
'[query.q]
read = "C2"
refr = "CC"' 'CC@CC'

case_ "a comma inside a group is rejected" \
'[query.q]
read = "{A,C}"
refr = "~~~~~"' 'A@A'

case_ "relational codes go on one side only" \
'[query.q]
read = "="
refr = "="' 'A@A'

case_ "relational + a side that holds no plain base: dead column, compile error" \
'[query.q]
read = "."
refr = "="' 'C@C'

case_ "relational + IUPAC: the other side's set still applies (warning names the unreachable part)" \
'[query.q]
read = "Y"
refr = "="' 'C@C'

case_ "... and Y@= does not fire on a T over C" \
'[query.q]
read = "Y"
refr = "="' 'T@C'

case_ "a relational code cannot be an alias value" \
'[alias]
j = "="
[query.q]
read = "j"
refr = "~"' 'A@A'

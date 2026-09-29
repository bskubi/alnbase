#!/usr/bin/env bash
# Render every set of API documentation in this repository into site/.
#
#   ./make-docs.sh          build both, then print where to open them
#   ./make-docs.sh rust     the alnbase crate only
#   ./make-docs.sh agg      the alnbase-agg Python package only
#
# site/ holds generated pages and is ignored by git. Hand-written documentation
# lives in docs/ and is committed. Keeping the two apart makes it clear which
# files a person may edit: nothing under site/ survives the next build.
#
#   site/index.html   a page that links to both
#   site/rust/        the alnbase crate, from rustdoc
#   site/agg/         the alnbase-agg package, from pdoc
#
# Read the rendered pages rather than the source. A doc comment can look correct
# in the file and still render as an empty section, a wall of text, or a link
# that goes nowhere.

set -euo pipefail

cd "$(dirname "$0")"

REPO=$PWD
SITE=$REPO/site
WHICH=${1:-all}

build_rust() {
    echo "== alnbase (rustdoc)"
    # The alias in .cargo/config.toml supplies --document-private-items. alnbase
    # is a binary crate, so without it almost every page comes out empty.
    cargo docs
    # cargo writes under its own target directory, which is often outside the
    # repository. Link to it rather than copying, so the pages never go stale.
    local target=${CARGO_TARGET_DIR:-$REPO/target}
    rm -rf "$SITE/rust"
    ln -s "$target/doc" "$SITE/rust"
}

build_agg() {
    echo "== alnbase-agg (pdoc)"
    local pdoc=${PDOC:-$REPO/python/.venv/bin/pdoc}
    if [[ ! -x $pdoc ]]; then
        echo "no pdoc at $pdoc" >&2
        echo "build the environment first:" >&2
        echo "    cd python && uv venv --python 3.12 .venv" >&2
        echo "    cd python && VIRTUAL_ENV=.venv uv pip install -e . pdoc" >&2
        return 1
    fi

    # Every submodule is named on the command line. pdoc walks a package's
    # submodules on its own only when the package does not define `__all__`, and
    # __init__.py does define one, because the package exports less than the set
    # of modules in it. Naming the modules keeps both: `__all__` still says what
    # the package exports, and the documentation still covers every module.
    local modules=(alnbase_agg)
    local path name
    for path in $(ls python/alnbase_agg/*.py | sort); do
        name=$(basename "$path" .py)
        [[ $name == __init__ ]] && continue
        modules+=("alnbase_agg.$name")
    done

    # -t points at a template that overrides pdoc's `is_public`, so that the
    # internal helpers are documented too. See the comment in the template.
    (cd python && "$pdoc" "${modules[@]}" -t pdoc-templates -o "$SITE/agg")
}

mkdir -p "$SITE"

case $WHICH in
    all) build_rust; build_agg ;;
    rust) build_rust ;;
    agg) build_agg ;;
    *) echo "usage: $0 [all|rust|agg]" >&2; exit 2 ;;
esac

cat > "$SITE/index.html" <<'HTML'
<!doctype html>
<meta charset="utf-8">
<title>alnbase API documentation</title>
<style>
  body { font: 16px/1.6 system-ui, sans-serif; max-width: 34rem; margin: 4rem auto;
         padding: 0 1rem; color: #222; }
  h1 { font-size: 1.4rem; }
  li { margin: .9rem 0; }
  code { background: #f2f2f2; padding: .1rem .3rem; border-radius: 3px; }
  p.note { color: #555; font-size: .9rem; }
</style>
<h1>alnbase API documentation</h1>
<ul>
  <li><a href="rust/alnbase/index.html">alnbase</a> &mdash; the Rust binary that
      walks alignments and matches queries.</li>
  <li><a href="agg/alnbase_agg.html">alnbase-agg</a> &mdash; the Python package
      that compiles recipes into DuckDB SQL over hit tables.</li>
</ul>
<p class="note">These pages are generated. Rebuild them with
<code>./make-docs.sh</code>. The hand-written guides are in <code>docs/</code>.</p>
HTML

echo
echo "open $SITE/index.html"

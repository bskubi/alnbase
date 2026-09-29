#!/usr/bin/env bash
# Fetch the hyper-editing pipeline at its published commit and make it runnable
# off the authors' cluster.
#
#   setup.sh WORKDIR
#
# Everything here is portability only. No behaviour changes: the patched tree
# produces the same output the original would have produced on the machine it
# was written for. Behaviour changes live in run.sh, which swaps whole stages.
#
# What has to be fixed before the tree runs at all:
#   1. Four Perl scripts carry `#!/private/apps/bin/perl516`, a path that exists
#      on one cluster.
#   2. The five shell scripts have no shebang, and are invoked as programs.
#   3. Nothing in the tree is executable, and the scripts call each other by
#      bare name, so they have to be on PATH.
#   4. filter_sam.pl opens `unique_simple_repeats.txt` by relative path, so it
#      only works with the repository as the working directory.
#   5. The FASTQ extraction step calls `bam2fastx`, from an old TopHat
#      distribution. A four-line wrapper over `samtools fastq` stands in.
set -euo pipefail

COMMIT=5eee27943b399f9703fda2903fb3de25bd0f942a
WORK="${1:?usage: setup.sh WORKDIR}"
HE="$WORK/Hyper-editing"

mkdir -p "$WORK"
if [ ! -d "$HE/.git" ]; then
  git clone -q https://github.com/hagitpt/Hyper-editing.git "$HE"
fi
git -C "$HE" checkout -q "$COMMIT"
git -C "$HE" clean -qfd
git -C "$HE" checkout -q -- .

# 1. shebangs that name a cluster-local perl
sed -i '1s|^#!/private/apps/bin/perl516 -w$|#!/usr/bin/perl -w|' \
  "$HE"/analyse_mm.pl "$HE"/detect_ue.pl "$HE"/Get_orig_read.pl "$HE"/sort_R_read.pl

# 2. shell scripts with no shebang at all
for s in pre_unmapped.sh TransRun.sh analyse_mm.sh run_hyper_editing.sh TransformIndexBWA_genome.sh; do
  if [ "$(head -c 2 "$HE/$s")" != '#!' ]; then
    printf '#!/usr/bin/env bash\n%s' "$(cat "$HE/$s")" > "$HE/$s.tmp"
    mv "$HE/$s.tmp" "$HE/$s"
  fi
done

# 3. executable, and `mkdir` that does not abort a rerun
sed -i 's/^mkdir /mkdir -p /' "$HE/run_hyper_editing.sh"
chmod +x "$HE"/*.pl "$HE"/*.sh

# 4. the simple-repeat list, by absolute path
sed -i "s|^my \$sim_rep_file = \"unique_simple_repeats.txt\";|my \$sim_rep_file = \"$HE/unique_simple_repeats.txt\";|" \
  "$HE/filter_sam.pl"

# 5. bam2fastx stand-in: the pipeline calls it as `bam2fastx --fastq -o OUT IN`
cat > "$HE/bam2fastx" <<'SHIM'
#!/usr/bin/env bash
# Stands in for TopHat's bam2fastx. Only the one call shape the pipeline uses.
set -euo pipefail
out=""; args=()
while [ $# -gt 0 ]; do
  case "$1" in
    --fastq) shift ;;
    -o) out="$2"; shift 2 ;;
    *) args+=("$1"); shift ;;
  esac
done
samtools fastq -n "${args[@]}" > "$out"
SHIM
chmod +x "$HE/bam2fastx"

echo "$COMMIT" > "$WORK/hyper-editing-commit.txt"
echo "patched tree: $HE"

#!/usr/bin/env bash
# Build the design notes to PDF with pandoc, into design/build/.
#
# Usage: design/build-pdf.sh [--proof] [NOTE.md ...]
# With no notes given, every note is built. --proof builds NOTE-proof.pdf, a
# copy for proofreading with wide margins, line numbers and the commit.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="$here/build"

proof=false
if [[ "${1:-}" == "--proof" ]]; then
  proof=true
  shift
fi

if [[ $# -gt 0 ]]; then
  notes=("$@")
else
  notes=(LANGUAGE.md EXAMPLES.md)
fi

# The generated appendices must match the sections before a build.
python3 "$here/appendices.py" --check

# GitHub identifiers keep the anchors the notes already link to.
reader="gfm+tex_math_dollars+pipe_tables"

mkdir -p "$out"

# The authors are listed once, in myst.yml.
authors=()
while IFS= read -r name; do
  authors+=(--metadata "author=$name")
done < <(sed -n 's/^ *- name: //p' "$here/myst.yml")

layout=(
  --variable geometry:margin=2cm
  --variable fontsize=10pt
)
suffix=""
if $proof; then
  # The commit in the footer ties each mark to the text that was printed.
  version="$(git -C "$here" rev-parse --short HEAD), $(date +%Y-%m-%d)"
  if [[ -n "$(git -C "$here" status --porcelain -- "${notes[@]}")" ]]; then
    version="$version, with uncommitted changes"
  fi
  printf '\\newcommand{\\proofversion}{%s}\n' "$version" >"$out/proof-version.tex"
  layout=(
    --variable geometry:a4paper
    --variable geometry:left=2.5cm
    --variable geometry:right=5.5cm
    --variable geometry:top=2.5cm
    --variable geometry:bottom=3cm
    --variable fontsize=11pt
    --include-in-header "$out/proof-version.tex"
    --include-in-header "$here/pandoc/proof-header.tex"
  )
  suffix="-proof"
fi

for note in "${notes[@]}"; do
  name="${note%.md}$suffix"
  pandoc "$here/$note" \
    --from "$reader" \
    --standalone \
    --shift-heading-level-by=-1 \
    --toc --toc-depth=2 \
    --pdf-engine=xelatex \
    "${layout[@]}" \
    --include-in-header "$here/pandoc/pdf-header.tex" \
    --variable mainfont="Helvetica Neue" \
    --variable monofont="Menlo" \
    --variable colorlinks=true \
    "${authors[@]}" \
    --metadata note="$note" \
    --lua-filter "$here/pandoc/notes.lua" \
    --output "$out/$name.pdf"
  echo "built $out/$name.pdf"
done

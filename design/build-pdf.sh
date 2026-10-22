#!/usr/bin/env bash
# Build the design notes to PDF with pandoc, into design/build/.
#
# Usage: design/build-pdf.sh [NOTE.md ...]
# With no notes given, every note is built.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="$here/build"

if [[ $# -gt 0 ]]; then
  notes=("$@")
else
  notes=(DESIGN.md LANGUAGE.md EXAMPLES.md DESIGN_CHOICES.md FINAL_TODOS.md)
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

for note in "${notes[@]}"; do
  name="${note%.md}"
  pandoc "$here/$note" \
    --from "$reader" \
    --standalone \
    --shift-heading-level-by=-1 \
    --toc --toc-depth=2 \
    --pdf-engine=xelatex \
    --include-in-header "$here/pandoc/pdf-header.tex" \
    --variable geometry:margin=2cm \
    --variable mainfont="Helvetica Neue" \
    --variable monofont="Menlo" \
    --variable fontsize=10pt \
    --variable colorlinks=true \
    "${authors[@]}" \
    --metadata note="$note" \
    --lua-filter "$here/pandoc/notes.lua" \
    --output "$out/$name.pdf"
  echo "built $out/$name.pdf"
done

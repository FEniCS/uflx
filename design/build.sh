#!/usr/bin/env bash
# Build the design notes with pandoc.
#
# Usage: design/build.sh [html|pdf|all] [NOTE.md ...]
# Output goes to design/build/. With no notes given, every note is built.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="$here/build"
format="${1:-html}"
shift || true

if [[ $# -gt 0 ]]; then
  notes=("$@")
else
  notes=(LANGUAGE.md DESIGN.md DESIGN_CHOICES.md EXAMPLES.md)
fi

# GitHub identifiers keep the anchors the notes already link to.
reader="gfm+tex_math_dollars+pipe_tables"

build_html() {
  mkdir -p "$out"
  cp "$here/pandoc/style.css" "$out/style.css"
  for note in "${notes[@]}"; do
    local name="${note%.md}"
    pandoc "$here/$note" \
      --from "$reader" \
      --to html5 \
      --standalone \
      --shift-heading-level-by=-1 \
      --toc --toc-depth=3 \
      --mathjax \
      --css style.css \
      --metadata link-ext=html \
      --lua-filter "$here/pandoc/notes.lua" \
      --output "$out/$name.html"
    echo "built $out/$name.html"
  done
}

build_pdf() {
  mkdir -p "$out"
  for note in "${notes[@]}"; do
    local name="${note%.md}"
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
      --metadata link-ext=pdf \
      --metadata note="$note" \
      --lua-filter "$here/pandoc/notes.lua" \
      --output "$out/$name.pdf"
    echo "built $out/$name.pdf"
  done
}

case "$format" in
  html) build_html ;;
  pdf) build_pdf ;;
  all) build_html; build_pdf ;;
  *) echo "usage: $0 [html|pdf|all] [NOTE.md ...]" >&2; exit 2 ;;
esac

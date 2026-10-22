#!/usr/bin/env bash
# Build the design notes with pandoc.
#
# Usage: design/build.sh [html|pdf|all] [NOTE.md ...]
# Output goes to design/build/. With no notes given, every note is built.
# README.md is the landing page, index.html, and has no PDF.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="$here/build"
format="${1:-html}"
shift || true

# Every note, with its label in the navigation bar of the HTML pages.
all_notes=(README.md DESIGN.md LANGUAGE.md EXAMPLES.md DESIGN_CHOICES.md FINAL_TODOS.md)
labels=("UFLx design notes" "Design" "Language" "Examples" "Design choices" "TODOs")

if [[ $# -gt 0 ]]; then
  notes=("$@")
else
  notes=("${all_notes[@]}")
fi

# The generated appendices must match the sections before a build.
python3 "$here/appendices.py" --check

# GitHub identifiers keep the anchors the notes already link to.
reader="gfm+tex_math_dollars+pipe_tables"

html_name() {
  if [[ "$1" == README.md ]]; then echo index; else echo "${1%.md}"; fi
}

# Write the navigation bar of one page, with the page itself marked.
nav_html() {
  local current="$1"
  echo '<nav class="notes">'
  for k in "${!all_notes[@]}"; do
    local note="${all_notes[$k]}"
    local mark=""
    [[ "$note" == "$current" ]] && mark=' aria-current="page"'
    echo "<a href=\"$(html_name "$note").html\"$mark>${labels[$k]}</a>"
  done
  echo '</nav>'
}

build_html() {
  mkdir -p "$out"
  cp "$here/pandoc/style.css" "$out/style.css"
  local nav
  nav="$(mktemp)"
  trap 'rm -f "$nav"' RETURN
  for note in "${notes[@]}"; do
    nav_html "$note" > "$nav"
    local name page=(--toc --toc-depth=3)
    name="$(html_name "$note")"
    [[ "$note" == README.md ]] && page=(--metadata landing=true)
    pandoc "$here/$note" \
      --from "$reader" \
      --to html5 \
      --standalone \
      --shift-heading-level-by=-1 \
      "${page[@]}" \
      --include-before-body "$nav" \
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
    [[ "$note" == README.md ]] && continue
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

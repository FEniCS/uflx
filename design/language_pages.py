"""Locate the language pages and join them for the combined PDF."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAGES = tuple(
    HERE / "language" / f"{name}.md"
    for name in (
        "rules",
        "expressions",
        "topology",
        "geometry",
        "forms",
        "evaluation",
        "rewriting",
        "substitution",
        "differentiation",
        "lowering",
        "optimisation",
        "compiler",
        "extensions",
        "reference",
    )
)
REFERENCE = PAGES[-1]
LINK = re.compile(r"\]\(([^\s)]+)\)")


def headings(text: str) -> list[tuple[str, str]]:
    """Return heading text and page-local identifiers outside code fences."""
    counts: dict[str, int] = {}
    result = []
    fence = ""
    for line in text.splitlines():
        match = re.match(r"^(`{3,}|~{3,})", line)
        if match:
            token = match.group()
            if not fence:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = ""
            continue
        match = re.match(r"^#{1,6} (.+)$", line)
        if fence or not match:
            continue
        title = match.group(1)
        base = re.sub(r"[^\w\- ]", "", re.sub(r"[`*]", "", title).strip().lower())
        base = base.replace(" ", "-")
        count = counts.get(base, 0)
        counts[base] = count + 1
        result.append((title, base if count == 0 else f"{base}-{count}"))
    return result


def pdf_source(note: Path) -> str:
    """Join the language or prepare another note's links to the combined PDF."""
    note = note.resolve()
    landing = HERE / "language.md"
    combined = note.samefile(landing)
    if combined:
        note = landing
    # Keep the original chapter headings, not the extra website page titles.
    bodies = {page: page.read_text().split("\n", 1)[1].lstrip() for page in PAGES}
    intro = landing.read_text().split("## Contents\n", 1)[0]
    counts: dict[str, int] = {}
    targets: dict[tuple[Path, str], str] = {}
    for title, anchor in headings(intro):
        counts[anchor] = counts.get(anchor, 0) + 1
    for page, body in bodies.items():
        local = headings(page.read_text())[1:]
        for title, local_anchor in local:
            # Page-local suffixes must be recomputed across the whole PDF.
            base = re.sub(r"[^\w\- ]", "", re.sub(r"[`*]", "", title).strip().lower())
            base = base.replace(" ", "-")
            count = counts.get(base, 0)
            counts[base] = count + 1
            targets[page, local_anchor] = base if count == 0 else f"{base}-{count}"
        targets[page, ""] = targets[page, local[0][1]]
        targets[page, headings(page.read_text())[0][1]] = targets[page, ""]

    def links(text: str, source: Path) -> str:
        def replace(match: re.Match[str]) -> str:
            target = match.group(1)
            if "://" in target or target.startswith("mailto:"):
                return match.group()
            filename, separator, anchor = target.partition("#")
            destination = (source.parent / filename).resolve() if filename else source
            if destination in PAGES:
                anchor = targets[destination, anchor]
                destination = landing
                separator = "#"
            if combined and destination == landing:
                return f"](#{anchor or headings(intro)[0][1]})"
            if not filename and destination == source:
                return match.group()
            relative = os.path.relpath(destination, note.parent)
            return f"]({relative}{separator}{anchor})"

        return LINK.sub(replace, text)

    if combined:
        return (
            links(intro, landing)
            + "\n"
            + "\n".join(links(body, page) for page, body in bodies.items())
        )
    return links(note.read_text(), note)


def main() -> None:
    """Write Markdown prepared for the existing PDF pipeline."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("note", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(pdf_source(args.note))


if __name__ == "__main__":
    main()

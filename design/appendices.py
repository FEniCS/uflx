"""Generate the appendices of LANGUAGE.md from the sections.

Run ``python3 design/appendices.py`` to rewrite the generated regions, or
``python3 design/appendices.py --check`` to fail when a region is stale or a
section disagrees with an appendix.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
LANGUAGE = HERE / "LANGUAGE.md"
EXAMPLES = HERE / "EXAMPLES.md"

# Classes of a language extension are not part of the core appendices.
EXCLUDED_CHAPTERS = {"Extending the core"}
FREE_FUNCTION_HEADER = ["Free function", "Overload", "Takes", "Returns"]
PYTHON_TYPES = {"int", "float", "str", "tuple", "bool", "None"}
TREE_COLUMN = 42


@dataclass
class ClassEntry:
    """A class described by a signature table."""

    name: str
    signature: str
    chapter: str
    anchor: str
    abstract: bool
    rows: dict[str, str]
    lead: str
    extends: list[str] = field(default_factory=list)


@dataclass
class FreeFunction:
    """A free function, from a class table or a free-function table."""

    signature: str
    chapter: str
    anchor: str
    overload: str
    takes: str
    returns: str


@dataclass
class LegacyNote:
    """A Legacy UFL subsection and the section it belongs to."""

    title: str
    link: str
    text: str


def slugify(heading: str) -> str:
    """Return the GitHub anchor of a heading."""
    text = re.sub(r"[`*]", "", heading).strip().lower()
    text = re.sub(r"[^\w\- ]", "", text)
    return text.replace(" ", "-")


def split_row(line: str) -> list[str]:
    """Split a pipe-table row into cells, keeping escaped pipes."""
    cells = re.split(r"(?<!\\)\|", line.strip().strip("|"))
    return [c.strip() for c in cells]


def code_names(text: str) -> list[str]:
    """Return the class names in backticks, with or without arguments."""
    return re.findall(r"`([A-Z]\w*)(?:\([^`]*\))?`", text)


def property_names(text: str) -> list[str]:
    """Return the property names that start the items of a property row."""
    names = re.findall(r"(?:^|, |and )`([a-z_][a-z_0-9]*)`", text)
    return [n for n in names if n not in PYTHON_TYPES]


SYMBOLS = {r"\otimes": "⊗", r"\psi": "ψ"}


def plain(text: str) -> str:
    """Strip Markdown and math markers for the plain-text tree."""
    for command, symbol in SYMBOLS.items():
        text = text.replace(command, symbol)
    return re.sub(r"[`$\\]", "", text)


class Document:
    """The headings, tables and Legacy UFL notes of one Markdown file."""

    def __init__(self, path: Path):
        """Parse the file at path."""
        self.path = path
        self.lines = path.read_text().split("\n")
        self.classes: list[ClassEntry] = []
        self.functions: list[FreeFunction] = []
        self.legacy: list[LegacyNote] = []
        self.tables: dict[str, list[list[list[str]]]] = {}
        self._parse()

    def _parse(self) -> None:
        chapter = ""
        section = ""
        anchor = ""
        used: dict[str, int] = {}
        in_code = False
        legacy: tuple[str, str] | None = None
        legacy_text: list[str] = []
        prose: list[str] = []
        link_prefix = "" if self.path == LANGUAGE else f"{self.path.name}"

        def flush_legacy() -> None:
            if legacy is not None:
                for paragraph in "\n".join(legacy_text).strip().split("\n\n"):
                    text = " ".join(paragraph.split())
                    # The chapter list of worked examples may follow the last note.
                    if not text.startswith("Worked examples for this chapter"):
                        self.legacy.append(LegacyNote(legacy[0], legacy[1], text))

        i = 0
        while i < len(self.lines):
            line = self.lines[i]
            if line.startswith("```"):
                in_code = not in_code
            heading = None if in_code else re.match(r"^(#{1,4}) (.+)$", line)
            if heading:
                level, title = len(heading.group(1)), heading.group(2)
                slug = slugify(title)
                n = used.get(slug, 0)
                used[slug] = n + 1
                slug = slug if n == 0 else f"{slug}-{n}"
                if title == "Legacy UFL":
                    flush_legacy()
                    legacy = (section, f"{link_prefix}#{slugify(section)}")
                    legacy_text = []
                else:
                    if level <= 3:
                        flush_legacy()
                        legacy = None
                    if level == 2:
                        chapter = title
                    if level in (2, 3):
                        section = title
                        prose = []
                    if level <= 4:
                        anchor = slug
                i += 1
                continue
            if legacy is not None:
                legacy_text.append(line.strip())
            is_table = (
                not in_code
                and line.startswith("|")
                and i + 1 < len(self.lines)
                and re.match(r"^\|\s*-", self.lines[i + 1])
            )
            if is_table:
                header = split_row(line)
                rows = []
                j = i + 2
                while j < len(self.lines) and self.lines[j].startswith("|"):
                    rows.append(split_row(self.lines[j]))
                    j += 1
                self._table(header, rows, chapter, anchor, " ".join(prose))
                self.tables.setdefault(chapter, []).append([header, *rows])
                i = j
                continue
            if not in_code:
                prose.append(line)
            i += 1
        flush_legacy()

    def _table(self, header, rows, chapter, anchor, lead) -> None:
        if header == FREE_FUNCTION_HEADER:
            for signature, overload, takes, returns in rows:
                self.functions.append(
                    FreeFunction(signature.strip("`"), chapter, anchor, overload, takes, returns)
                )
            return
        if len(header) != 2 or header[1] != "" or not code_names(header[0]):
            return
        signatures = re.findall(r"`([A-Z]\w*(?:\([^`]*\))?)`", header[0])
        values = {r[0]: r[1] for r in rows if len(r) == 2}
        for k, signature in enumerate(signatures):
            name = re.match(r"\w+", signature).group(0)
            entry = ClassEntry(
                name=name,
                signature=signature,
                chapter=chapter,
                anchor=anchor,
                abstract="(" not in signature,
                rows=values,
                lead=lead,
                extends=code_names(values.get("Extends", "")),
            )
            self.classes.append(entry)
            functions = re.findall(r"`(\w+\([^`]*\))`", values.get("Free function", ""))
            if k < len(functions):
                overload = re.search(r"overload (`[^`]+`)", values["Free function"])
                self.functions.append(
                    FreeFunction(
                        functions[k],
                        chapter,
                        anchor,
                        overload.group(1) if overload else "",
                        values.get("Takes", ""),
                        values.get("Returns", f"`{name}`"),
                    )
                )


def summary_of(entry: ClassEntry) -> str:
    """Return the summary of one class, from a shared or a split row."""
    summary = entry.rows.get("Summary", "")
    parts = re.findall(r"`(\w+)`: ([^;]+)", summary)
    if parts:
        return dict(parts).get(entry.name, "").strip()
    return summary


def core(classes: list[ClassEntry]) -> list[ClassEntry]:
    """Return the classes of the core, without those of language extensions."""
    return [c for c in classes if c.chapter not in EXCLUDED_CHAPTERS]


def class_hierarchy(classes: list[ClassEntry]) -> str:
    """Render the class hierarchy, one tree per chapter."""
    out = []
    chapters = []
    for c in core(classes):
        if c.extends and c.chapter not in chapters:
            chapters.append(c.chapter)
    for chapter in chapters:
        members = [c for c in core(classes) if c.chapter == chapter and c.extends]
        names = {c.name for c in members}
        out.append(f"# {chapter}")
        roots: list[tuple[str, list[ClassEntry]]] = []
        for c in members:
            if len(c.extends) == 1 and c.extends[0] in names:
                continue
            key = ", ".join(c.extends)
            for root, group in roots:
                if root == key:
                    group.append(c)
                    break
            else:
                roots.append((key, [c]))

        def walk(group: list[ClassEntry], prefix: str) -> None:
            for k, c in enumerate(group):
                last = k == len(group) - 1
                label = f"{prefix}{'└── ' if last else '├── '}{c.name}"
                summary = plain(summary_of(c))
                line = f"{label.ljust(TREE_COLUMN)}{summary}" if summary else label
                out.append(line.rstrip())
                children = [m for m in members if m.extends == [c.name]]
                walk(children, prefix + ("    " if last else "│   "))

        for root, group in roots:
            out.append(root)
            walk(group, "")
        out.append("")
    return "```\n" + "\n".join(out).rstrip() + "\n```\n"


def table(header: list[str], rows: list[list[str]]) -> str:
    """Render a padded pipe table."""
    widths = [max(len(r[k]) for r in [header, *rows]) for k in range(len(header))]
    widths = [max(w, 3) for w in widths]

    def line(cells):
        return "| " + " | ".join(c.ljust(w) for c, w in zip(cells, widths, strict=True)) + " |"

    sep = "| " + " | ".join("-" * w for w in widths) + " |"
    return "\n".join([line(header), sep, *(line(r) for r in rows)]) + "\n"


def free_functions(functions: list[FreeFunction]) -> str:
    """Render the table of free functions, grouped by chapter."""
    rows = []
    previous = None
    for f in functions:
        if f.chapter in EXCLUDED_CHAPTERS:
            continue
        group = f.chapter if f.chapter != previous else ""
        previous = f.chapter
        rows.append([group, f"[`{f.signature}`](#{f.anchor})", f.overload, f.takes, f.returns])
    header = ["Group", "Free function", "`AbstractExpression` overload", "Takes", "Returns"]
    return table(header, rows)


def properties(classes: list[ClassEntry]) -> str:
    """Render the stored, derived and abstract properties of each class."""
    rows = []
    for c in core(classes):
        stored = property_names(c.rows.get("Properties", ""))
        derived = property_names(c.rows.get("Derived", ""))
        provides = property_names(c.rows.get("Provides", ""))
        if not (stored or derived or provides):
            continue
        cells = [", ".join(f"`{n}`" for n in names) for names in (stored, derived, provides)]
        rows.append([f"`{c.name}`", *cells])
    return table(["Class", "Stored", "Derived", "Provides"], rows)


def legacy_ufl(notes: list[LegacyNote]) -> str:
    """Render the Legacy UFL notes as one bullet each."""
    return "\n".join(f"- **[{n.title}]({n.link}).** {n.text}" for n in notes) + "\n"


def opening_parents(entry: ClassEntry) -> list[str] | None:
    """Return the parents that the opening sentence of a section names."""
    m = re.search(
        r"subclass of (the base object|\[?`(\w+)`\]?(?:\([^)]*\))?)(?: and of `(\w+)`)?", entry.lead
    )
    if not m:
        return None
    if m.group(1) == "the base object":
        return ["Base"]
    return [p for p in (m.group(2), m.group(3)) if p]


def checks(lang: Document) -> list[str]:
    """Return the disagreements between the sections and the hand-written appendices."""
    problems = []
    for c in core(lang.classes):
        if c.name != "Base" and not c.extends:
            problems.append(f"{c.name}: the table has no Extends row")
        named = opening_parents(c)
        if named and c.extends and named != c.extends and len(lang_entries(lang, c)) == 1:
            problems.append(f"{c.name}: Extends {c.extends}, the opening sentence says {named}")
    expressions = {
        c.name for c in core(lang.classes) if not c.abstract and is_expression(lang, c.name)
    }
    known = {c.name for c in lang.classes}
    # The first table of each appendix covers every node class.
    for title in ("Appendix: Static type table", "Appendix: Pass table"):
        tables = lang.tables.get(title, [])
        first = tables[0][1:] if tables else []
        listed = {n for r in first for n in code_names(r[0])}
        for name in sorted(expressions - listed):
            problems.append(f"{title}: no row for {name}")
        for rows in tables:
            for name in sorted({n for r in rows[1:] for n in code_names(r[0])} - known):
                problems.append(f"{title}: {name} is not a class of the sections")
    return problems


def lang_entries(lang: Document, entry: ClassEntry) -> list[ClassEntry]:
    """Return the classes that share the table of an entry."""
    return [c for c in lang.classes if c.rows is entry.rows]


def is_expression(lang: Document, name: str) -> bool:
    """Return whether a class extends AbstractExpression."""
    parents = {c.name: c.extends for c in lang.classes}
    seen = set()
    todo = [name]
    while todo:
        n = todo.pop()
        if n == "AbstractExpression":
            return True
        if n not in seen:
            seen.add(n)
            todo.extend(parents.get(n, []))
    return False


def render(lang: Document, examples: Document) -> dict[str, str]:
    """Return the text of each generated region."""
    return {
        "legacy-ufl": legacy_ufl(lang.legacy + examples.legacy),
        "class-hierarchy": class_hierarchy(lang.classes),
        "free-functions": free_functions(lang.functions),
        "properties": properties(lang.classes),
    }


def fill(text: str, regions: dict[str, str]) -> str:
    """Replace the body of each generated region of text."""
    for name, body in regions.items():
        start = rf"<!-- generated: {re.escape(name)}[^>]*-->\n"
        end = rf"<!-- end generated: {re.escape(name)} -->"
        pattern = re.compile(f"({start})(.*?)({end})", re.S)
        if not pattern.search(text):
            raise SystemExit(f"LANGUAGE.md has no generated region {name}")
        text = pattern.sub(lambda m: m.group(1) + body + m.group(3), text)
    return text


def main() -> int:
    """Rewrite or check the generated appendices."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="fail instead of rewriting")
    args = parser.parse_args()
    lang = Document(LANGUAGE)
    examples = Document(EXAMPLES)
    text = LANGUAGE.read_text()
    new = fill(text, render(lang, examples))
    problems = checks(lang)
    for p in problems:
        print(f"warning: {p}", file=sys.stderr)
    if args.check:
        if new != text:
            print("the appendices are stale, run design/appendices.py", file=sys.stderr)
            return 1
        return 1 if problems else 0
    if new != text:
        LANGUAGE.write_text(new)
        print("updated the appendices of LANGUAGE.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())

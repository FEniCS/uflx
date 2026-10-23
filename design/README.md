# UFLx design notes

UFLx is a minimal core language for finite element variational forms,
embedded in Python. These notes specify the language and argue for its
design. The notes are written for users and developers who know a language
for variational forms, such as UFL.

## Contents

Readers may want to start with the **[Worked examples](EXAMPLES.md)**.

- **[Core UFLx language definition](LANGUAGE.md).**
- **[Worked examples](EXAMPLES.md)**.

## AI usage statement

The design notes in this directory have been co-created with Claude Opus 5.5
and Claude Fable 5.1. Both models have presumably been trained on a large
corpus of material on modern differential geometry, algorithms and functional
programming, and on the source code of Firedrake, DOLFINx, UFL, FFCx, TSFC,
GEM and many other packages. That said, Claude Code was instructed not to
read any of those materials or source code directly. In the judgement of the
authors, the quality of the notes and of the resulting design is
substantially higher than the authors could have reached without an LLM.

## Building

The site builds with [Jupyter Book](https://jupyterbook.org) into
`design/_build/html/`, and the PDFs build with pandoc and xelatex into
`design/build/`.

```
cd design && jupyter-book build --html   # the site
cd design && jupyter-book start          # serve the site locally
design/build-pdf.sh                      # every PDF
design/build-pdf.sh LANGUAGE.md
```

# UFLx design notes

UFLx is a minimal core language for finite element variational forms,
embedded in Python. These notes specify the language and argue for its
design. The notes are written for users and developers who know a language
for variational forms, such as UFL.

## Contents

Readers may want to start with the **[Worked examples](EXAMPLES.md)**.

- **[The design of the UFLx core](DESIGN.md).**
- **[Core UFLx language definition](LANGUAGE.md).**
- **[Worked examples](EXAMPLES.md)**.
- **[Design choices](DESIGN_CHOICES.md).**

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

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

The notes build to HTML and PDF with pandoc, into `design/build/`. The PDF
needs xelatex.

```
design/build.sh            # html
design/build.sh pdf
design/build.sh all LANGUAGE.md
```

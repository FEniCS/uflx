# Core UFLx language definition

This note defines UFLx objects and the rules for constructing and evaluating
them as local kernels. The finite element package supplies concrete cell and
interface data, dispatches kernels and assembles their contributions. The note
includes the core's rewrite passes, but leaves numerical algorithms and backend
implementation to the compiler.

The library classes of the examples and further worked examples are in
[examples/index.md](examples/index.md). The note assumes knowledge of differential
geometry and functional analysis, and explains some advanced concepts of
differential geometry in more detail.

This document is written for users and developers who know domain-specific
languages for variational forms, such as UFL and ngsolve.

## Contents

- [Language rules](language/rules.md)
- [Expressions](language/expressions.md)
- [Topological meshes](language/topology.md)
- [Geometry and calculus](language/geometry.md)
- [Functions and forms](language/forms.md)
- [Rewriting and lowering](language/evaluation.md)
  - [Rewrite machinery](language/rewriting.md)
  - [Substitution](language/substitution.md)
  - [Differentiation](language/differentiation.md)
  - [Applications and mappings](language/lowering.md)
- [Extensions](language/extensions.md)
- [Reference](language/reference.md)

The PDF build combines these pages into `language.pdf`.

# Worked examples for the UFLx core

These examples accompany [language.md](../language.md). Each states its
library-supplied inputs locally and uses core names without import lists.
The code constructs symbolic expressions; numerical results require later
passes or a solver.

Start with [Poisson](forms.md#poisson) for variational forms or the
[affinely mapped triangle](integration.md#an-affinely-mapped-triangle)
for geometry. The
[FEniCS layer](forms.md#a-classic-fenics-layer-for-the-poisson-problem)
shows the familiar UFL syntax.

## Topics

- **[Tensor operations](tensors.md).** Construct tensor operations
  from indices and represent symmetry in stored data.
- **[Geometry and integration](integration.md).** Map reference
  cells, integrate over edges and cells, and introduce quadrature and cut-cell rules.
- **[Facets and coupled meshes](coupling.md).** Use traces and side
  maps for facet forms, fluid–solid coupling and embedded networks.
- **[Variational forms](forms.md).** Build constrained, complex,
  mixed and nonlinear forms, and a familiar FEniCS interface.
- **[Surfaces and moving geometry](surfaces.md).** Write tangential
  operators, oriented fluxes and shape derivatives on fixed or deformed geometry.
- **[Error estimation](estimation.md).** Express local residual
  solves and goal-oriented error indicators.

The [reference library](library.md) collects the longer helper
implementations: simplices, static topologies, side tables and the FEniCS layer.

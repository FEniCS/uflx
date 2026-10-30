# Worked examples for the UFLx core

These examples accompany [LANGUAGE.md](LANGUAGE.md). Each states its
library-supplied inputs locally and uses core names without import lists.
The code constructs symbolic expressions; numerical results require later
passes or a solver.

For an introduction to the core constructions, start with the
[affinely mapped triangle](EXAMPLES_INTEGRATION.md#an-affinely-mapped-triangle).
For a familiar entry point from UFL, start with the
[FEniCS layer](EXAMPLES_FORMS.md#a-classic-fenics-layer-for-the-poisson-problem).

## Topics

- **[Tensor operations](EXAMPLES_TENSORS.md).** Construct tensor operations
  from indices and represent symmetry in stored data.
- **[Geometry and integration](EXAMPLES_INTEGRATION.md).** Map reference
  cells, integrate over edges and cells, and introduce quadrature and cut-cell rules.
- **[Facets and coupled meshes](EXAMPLES_COUPLING.md).** Use traces and side
  maps for facet forms, fluid–solid coupling and embedded networks.
- **[Variational forms](EXAMPLES_FORMS.md).** Build constrained, complex,
  mixed and nonlinear forms, and a familiar FEniCS interface.
- **[Surfaces and moving geometry](EXAMPLES_SURFACES.md).** Write tangential
  operators, oriented fluxes and shape derivatives on fixed or deformed geometry.
- **[Error estimation](EXAMPLES_ESTIMATION.md).** Express local residual
  solves and goal-oriented error indicators.

The [reference library](EXAMPLES_LIBRARY.md) collects the longer helper
implementations: simplices, static topologies, side tables and the FEniCS layer.

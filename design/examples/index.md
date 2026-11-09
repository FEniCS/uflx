# Worked examples for the UFLx core

These examples accompany [language.md](../language.md). Each states its
library-supplied inputs locally and uses core names without import lists.
The code constructs symbolic expressions; numerical results require later
passes or a solver. Topologies describe fixed local configurations; the
package supplies concrete cell data, dispatches kernels and assembles results.

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
- **[Physical polygon meshes](physical_meshes.md).** Sketch a language extension
  and compiler for physical polynomial DG spaces.
- **[Facets and coupled meshes](coupling.md).** Use traces and side
  maps for facet forms, fluid–solid coupling and embedded networks.
- **[Variational forms](forms.md).** Build constrained, complex,
  mixed and nonlinear forms, and a familiar FEniCS interface.
- **[Surfaces and moving geometry](surfaces.md).** Write tangential
  operators, oriented fluxes and shape derivatives on fixed or deformed geometry.
- **[Error estimation](estimation.md).** Express local residual
  solves and goal-oriented error indicators.
- **[Patch problems](patches.md).** Build contributions to a DG patch problem
  with a patch-wide multiplier and a zero-mean correction.

The [reference library](library.md) collects the longer helper
implementations: simplices, static topologies, side data and the FEniCS layer.

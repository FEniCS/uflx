# Final TODOs for variational forms in LANGUAGE.md

What the core language cannot yet express when defining variational forms, in
the suggested order of work.

## Blocking

Nothing. Integrals over cells and facets use the induced density in
`integral(f, D)`, and `normal(D, E, s)` gives the outward normal.

## Common forms

1. **Interpolation.** Add the interpolant of an expression in a space, as a
   language extension, and its adjoint on the dual space. Interpolation from
   another mesh needs the cell that contains each point, a runtime gather of
   the cells.
1. **Conditionals and elementary functions.** Add `conditional` and `max` (for
   upwinding, limiters and switching materials) and more elementary functions,
   such as `exp` and `sin`.

## Derivatives

1. **Shape derivatives.** A language extension, with a node for the
   derivative when the coordinates of a mesh move, lowered after the integrals
   and the pull-backs are expanded to a `DirectionalDerivative` with respect to
   the coordinates (see PASSES_DRAFT.md). ALE needs the same extension, with
   the parameterisation of the pieces built from a displacement function.
   COORDINATE_FUNCTIONS_DRAFT.md proposes the coordinates as a function on the
   topology, which makes the variable of the extension a function.
1. **Static real types.** A complex function that appears only in real-valued
   combinations gives a misleading derivative.
1. **Literals in derivatives.** The derivative pass produces many `Zero` and
   `Integer(1)` terms. State whether a later pass or the form compiler folds
   them.

## Geometry

1. **One transformation law for the geometry.** The maps of the table
   "Representation on the reference cell" are instances of the transformation
   law of a tensor density, `TransformationLaw(variances, weight)` (see
   COORDINATE_FUNCTIONS_DRAFT.md). The rows of the core that use the Jacobian
   are instances of the same law: the pull-back of a gradient by a cell adds a
   covariant axis, the normal of a side normalises the covariant push-forward
   of the reference normal, the density of `expand_integrals` transforms the
   reference measure with the weight -1, and the coordinates are the
   push-forward of the coordinate function by the identity. Restate these rows,
   the table and the Piola extension through the law, so that the Jacobian,
   its inverse, its determinant and the density occur only in the lowering of
   the law.
1. **Cancellations of the law.** With the weight explicit, a pass can see
   cancellations that a structural comparison of graphs misses: the
   divergence of a contravariant field of weight 1 is the reference divergence
   over the determinant, the integral of a density of weight 1 is the integral
   on the reference cell, and the flux of a contravariant field through a
   facet is the reference flux, with no geometry. State whether the core, a
   language extension or the form compiler finds them.
1. **The density of a point.** The row $m = 0$ of `density` gives
   `Integer(1)`, which the general row already gives: the determinant of a
   matrix of value shape `(0, 0)` is `Integer(1)` by the Leibniz formula, with
   the one empty permutation, so the row only spares a `square_root` of 1 that
   a compiler folds. Remove the row, or keep it and say why.
1. **Hellan–Herrmann–Johnson.** The table "Representation on the reference
   cell" stops at Regge. Add the double contravariant map of Basix,
   J û Jᵀ / (det J)² for k = n, and with the density squared for k < n.

## Tensors

1. **`Inverse` with an `AbstractExtent`.** The requirement that $k$ is at most $n$ cannot be
   checked when an extent is symbolic.
1. **Ellipsis operations of a library.** `outer`, `dot` and `inner` build
   results with an `EllipsisIndex`, which the reduction cannot split per
   axis. Derive one `Index` per axis, or state that the results do
   not reduce.
1. **Symmetry of a form.** Detecting that a bilinear form is symmetric is left
   to a pass, which is not specified.

## Compilers and JAX

1. **Complex conventions.** The antilinear position 0 and the Wirtinger parts
   differ from the conventions of `jax.grad` and `jax.vjp`. Document the
   conjugate at the boundary.
1. **Non-smooth points.** `abs`, `square_root` and `ImaginaryPart` give NaN
   in derivatives, and the core has no `conditional` to guard them.
1. **Lowering.** Fuse an `IndexSum` of a `Product` into one contraction, lower
   a repeated index as a diagonal gather, and gather with in-bounds indices,
   because tables are never negative.

## Implementation

1. **Pattern matching.** Derive `__match_args__` from the construction
   arguments in `Base.__init_subclass__`, so that passes can use `match`.
   Dataclasses conflict with the equality, the hash and the interning of
   `Base`.

## Notes

1. **A check of reductions.** The fix for `EllipsisIndex` was traced by hand
   for `sigma` and `C - I` only. Add a script that reduces every composition
   of the core.
1. **Generated content.** Decide whether the new sections need the marker
   for LLM generated content.
1. **A duplicate name.** "Symmetric tensors" in EXAMPLES.md uses `sigma`
   twice.
1. **Formatting.** mdformat escapes the math of EXAMPLES.md, as for
   LANGUAGE.md.

## Already expressible

- The Laplacian, as `integral` over `D` of
  `index_sum(gradient(u, D)[i] * gradient(v, D)[i], (i,))`, with the
  Euclidean frame of `RD(n)`.
- The volume and the size $h$ of each cell, and the area of each facet, as an
  `Integral` over the reference cell of the density, with the free cell index.
- Mixed forms, such as Stokes and mixed Poisson with Raviart–Thomas spaces, as
  an argument that is a `tuple` of functions, `form(e, ((v, q), (u, p)))`. The members
  may be on different domains.
- H(curl) and H(div) spaces, as vector-valued spaces with a Piola map, with
  the normal trace through `normal`.
- Real spaces and Lagrange multipliers, as a space with the domain `None` and a
  `Function`.
- Subdomain and marked boundary terms, through topologies with sides of
  codimension 0 and 1.
- Point and vertex terms, through topologies of vertices, with `integral` over
  a decomposed domain of dimension 0.

## Outside the core

- Quadrature degree and kernel options belong to the form compiler.
- Non-local operators on product domains are left to language extensions.
- Time stepping.

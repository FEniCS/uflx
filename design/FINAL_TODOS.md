# Final TODOs for variational forms in LANGUAGE.md

What the core language cannot yet express when defining variational forms, in
the suggested order of work.

## Blocking

Nothing. Integrals over cells and facets use the induced density in
`integral(f, D)`, and `normal(D, E, s)` gives the outward normal.

## Common forms

1. **Linearisation.** `DirectionalDerivative` gives the derivative of an
   expression. Add the Gateaux transformation and its free function, which adds
   an argument to a form, and decide whether it applies to a form, so that a
   tangent can be built from a residual. For shape derivatives and ALE, the
   geometry must depend on a function, such as a parameterisation of the
   pieces from a displacement function on the domain, because a form binds
   only functions and the parameterisation is an expression of runtime data.
1. **Interpolation.** Add the interpolant of an expression in a space, as a
   language extension, and its adjoint on the dual space. Interpolation from
   another mesh needs the cell that contains each point, a runtime gather of
   the cells.
1. **Conditionals and elementary functions.** Add `conditional` and `max` (for
   upwinding, limiters and switching materials) and more elementary functions,
   such as `exp` and `sin`.

## Derivatives

1. **`DirectionalDerivative` through `PullBack`, `Form` and `Action`.** See
   PASSES_DRAFT.md. Through `PullBack` it is needed for
   shape derivatives. Through `Form` and `Action` it is needed to differentiate
   after binding, such as a Hessian-vector product of a bound vector-Jacobian
   product.
1. **Missed occurrences.** `directional_derivative` finds `g` by equality, so
   an expression equal in value to `g` but built differently is not
   differentiated. Decide on a checking pass or a stated warning.
1. **Static real types.** A complex function that appears only in real-valued
   combinations gives a misleading derivative.
1. **Literals in derivatives.** `apply_derivatives` produces many `Zero` and
   `Integer(1)` terms. State whether a later pass or the form compiler folds
   them.

## Tensors

1. **`transpose` in the core.** The examples use `transpose(F) @ F`, but
   `transpose` is a function of a library.
1. **`Inverse` with an `AbstractExtent`.** The requirement that $k$ is at most $n$ cannot be
   checked when an extent is symbolic.
1. **Ellipsis operations of a library.** `outer`, `dot` and `inner` build
   results with an `EllipsisIndex`, which `reduce_component_tensors` cannot
   split per axis. Derive one `Index` per axis, or state that the results do
   not reduce.
1. **Symmetry of a form.** Detecting that a bilinear form is symmetric is left
   to a pass, which is not specified.

## Compilers and JAX

1. **The compiler contract.** What a backend receives decides how a lowering,
   such as to JAX, is written.
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

# Final TODOs for variational forms in LANGUAGE.md

What the core language cannot yet express when defining variational forms, in
the suggested order of work.

## Blocking

Nothing. Integrals over cells and facets use the induced density in
`integral(f, D)`, and `normal(side_map(D, E, s))` gives the outward normal.

## Common forms

1. **Interpolation.** Add the interpolant of an expression in a space, as a
   language extension, and its adjoint on the dual space. Interpolation from
   another mesh needs the cell that contains each point, a runtime gather of
   the cells.
1. **Conditionals and elementary functions.** Add `conditional` and `max` (for
   upwinding, limiters and switching materials) and more elementary functions,
   such as `exp` and `sin`.

## Derivatives

1. **Shape derivatives of spaces on a mesh.** A function in a space on a mesh,
   such as Hermite, is read through `cell_domain(D, j, e)`, which neither the
   lowering nor the substitution of a derivative enters. Its dependence on the
   coordinate function is therefore not differentiated by `ShapeDerivative`.
1. **Static real types.** A complex function that appears only in real-valued
   combinations gives a misleading derivative.
1. **Literals in derivatives.** The derivative pass produces many `Zero` and
   `Integer(1)` terms. State whether a later pass or the form compiler folds
   them.

## Tensors

1. **`Inverse` with an `AbstractExtent`.** The requirement that $k$ is at most $n$ cannot be
   checked when an extent is symbolic.
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

## Passes

1. **Termination of the lowering stage.** Write the argument that repeating
   `apply_gradients` and `expand_pull_backs` ends. The toy takes at most three
   rounds.
1. **Rows applied again.** Decide whether `expand_pull_backs` still applies
   its rows again to each `PullBack` that a rewrite creates, now that the
   lowering stage runs to a fixed point.
1. **Supports in the prose.** Find the prose of LANGUAGE.md that says domain
   for the support of an expression, now that a support may be a topology.

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
1. **Formatting.** mdformat escapes the math of EXAMPLES.md, as for
   LANGUAGE.md.

## Already expressible

- The Laplacian, as `integral` over `D` of
  `index_sum(grad_u[i] * grad_v[i], (i,))`, with
  `grad_u = gradient(push_forward(u, D), D)`, in the Euclidean frame of
  `RD(n)`.
- The volume and the size $h$ of each cell, and the area of each facet, as an
  `Integral` over the reference cell of the density, with the free cell index.
- Mixed forms, such as Stokes and mixed Poisson with Raviart–Thomas spaces, as
  an argument that is a `tuple` of functions, `form(e, ((v, q), (u, p)))`. The members
  may be on different topologies.
- H(curl) and H(div) spaces, as spaces on a topology with a covariant or a
  contravariant `TransformationLaw`, with the normal trace through `normal`.
- Shape derivatives and ALE meshes, through the coordinate function of
  `image(x, Q)` and the shape derivative extension.
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

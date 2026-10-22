# Final TODOs for variational forms in LANGUAGE.md

What the core language cannot yet express when defining variational forms, in
the suggested order of work.

## Blocking

Nothing. Integrals over cells and facets use the induced density in
`integral(f, D)`, and `normal(D, E, s)` gives the outward normal.

## Common forms

1. **Surface gradient and conormal.** `components(du, D)` is refused when $D$ has a lower
   dimension than the ambient coordinate domain. Specify the tangential
   projection (the pseudo-inverse of a non-square Jacobian) for surfaces and
   shells, and the in-surface conormal of the edges of a surface mesh, which
   `normal` refuses.
1. **Linearisation.** Add the derivative with respect to a function or a
   `tuple` of functions, which adds a slot to a form. For shape derivatives and
   ALE, the geometry must depend on a function, such as a parameterisation of
   the pieces from a displacement function on the domain, because a form binds
   only functions and the parameterisation is an expression of runtime data.
1. **Interpolation.** Add the interpolant of an expression in a space, as a
   language extension, and its adjoint on the dual space. Interpolation from
   another mesh needs the cell that contains each point, a runtime gather of
   the cells.
1. **Conditionals and elementary functions.** Add `conditional` and `max` (for
   upwinding, limiters and switching materials) and more elementary functions,
   such as `exp` and `sin`.

## Already expressible

- The Laplacian, as `integral` over `D` of
  `index_sum(components(du, D)[i] * components(dv, D)[i], (i,))`, with the
  Euclidean frame of `RD(n)`.
- The volume and the size $h$ of each cell, and the area of each facet, as an
  `Integral` over the reference cell of the density, with the free cell index.
- Mixed forms, such as Stokes and mixed Poisson with Raviart–Thomas spaces, as
  a slot of a `tuple` of functions, `form(e, ((v, q), (u, p)))`. The members
  may be on different domains.
- H(curl) and H(div) spaces, as spaces of the degree
  `(1,)` and `(n - 1,)`, with the normal trace through `normal`.
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

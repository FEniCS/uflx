# Final TODOs for variational forms in LANGUAGE.md

What the core language cannot yet express when defining variational forms, in
the suggested order of work.

## Blocking

1. **Integral over a decomposed domain.** The section is a placeholder. Define
   the integral of `omega` over `D` as the sum over the groups $j$ of
   `index_sum(integral(s_j * pull_back(omega, cell_domain(D, j)), K_j), (cell_index(T, j),))`,
   with the orientation sign `s_j = det J / |det J|` for $k = n$.
1. **Volume form or density on a decomposed domain.** `volume_form` accepts only
   coordinate domains. Add `volume_form(D)` for $k = n$, so that a 0-form
   integrates as the legacy `f * dx`. For $k < n$ (surfaces, and facets for
   `ds` and `dS`) decide between the induced density `sqrt(det(J^T J))` and an
   orientation from a side. See PR 106's `INTEGRATION.md` §1–2 for the case for
   densities.
1. **Facet normal.** Define the outward normal of a cell of a topology with
   sides, taking a side name, as the outward normal of the cell on that side.
   The two sides give opposite normals. Needed for Neumann data, Nitsche terms
   and DG fluxes.

## Common forms

4. **Surface gradient.** `components(du, D)` is refused when $D$ has a lower
   dimension than the ambient coordinate domain. Specify the tangential
   projection (the pseudo-inverse of a non-square Jacobian) for surfaces and
   shells.
1. **Mixed and block forms.** A slot of a `Form` binds one function. Stokes is
   bilinear in `(v, q)` and `(u, p)`, so decide between a slot that binds a
   `tuple` of functions and a mixed space.
1. **Linearisation.** Add the derivative with respect to a function, which adds
   an argument to a form.
1. **Conditionals and elementary functions.** Add `conditional` and `max` (for
   upwinding, limiters and switching materials) and more elementary functions,
   such as `exp` and `sin`.

## Already expressible

- The Laplacian, once item 2 exists, as
  `index_sum(components(du, D)[i] * components(dv, D)[i], (i,))` times the
  volume form, with the Euclidean frame of `RD(n)`.
- The volume and the size $h$ of each cell, as an `Integral` over $K_j$ of the
  pulled-back volume form with the free cell index. The facet size waits on
  item 2.
- Real spaces and Lagrange multipliers, as a space with the domain `None` and a
  `Function`.
- Subdomain and marked boundary terms, through topologies with sides of
  codimension 0 and 1.
- Point and vertex terms, through topologies of vertices, once item 1 covers
  $k = 0$.

## Outside the core

- Quadrature degree and kernel options belong to the form compiler.
- Non-local operators on product domains are left to language extensions.
- Time stepping.

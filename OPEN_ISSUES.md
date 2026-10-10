# Open issues

Things noticed while working on the domain, parametrization and integral branches
(PRs #100, #102, #103) and while reviewing `uflx/geometry.py`. Local scratch file,
deliberately not committed. Last updated 2026-10-10.

## The geometry model

- [ ] **`g⁻¹` is not named.** `MetricTensor` now exists, but raising and lowering
  indices wants the inverse metric, and `Matrix.compute_inverse` still forms `JᵀJ` on
  its own for the pseudo-inverse. Worth having once something needs it; the sharing
  argument I first gave for naming `g` was weak, since the DAG already dedups
  structurally equal subexpressions.
- [ ] **Orientation is not expressible.** `UnitNormal` fixes the normal's direction but
  not its sign, since nothing knows which side of a facet its cell is on. UFL carries
  `CellOrientation` for this. Needed before an outward normal means anything.
- [ ] **No tangent basis.** The columns of `J` are one, and `MetricTensor` says how to
  measure with them, but nothing names them.
- [ ] **The volume form is not an object.** In differential geometry the measure *is*
  `√det g dX`. Here `Measure` carries `dim`/`codim`/`boundary_only` and knows nothing of
  geometry, while `Integral.pull_back_to_entity` injects `abs(JacobianDeterminant(...))`
  into the integrand by hand.
- [ ] **Nothing second order is expressible**, so no second fundamental form, no
  curvature, no Christoffel symbols, no Kirchhoff-Love shells.
  `AbstractParametrization` exposes only `value` and `jacobian`; this was the explicit
  cost of taking the matrix chain rule, and reversing it means revisiting that interface.
- [ ] **"Push forward" names two different operations.** `geometry.py`'s
  `PushedForwardPoint`/`PulledBackPoint` are `φ` and `φ⁻¹` on *points*; `maps.py`'s
  `PushedForward`/`PulledBack` are element value maps on *fields*, which are the
  pushforward of a vector field and the pullback of a form. Covariant Piola is `J⁻ᵀ` and
  contravariant Piola is `J/det J`, both expressible from `geometry.py`, but nothing
  states the connection.

## `uflx/geometry.py` specifics

- [ ] **`SpatialCoordinate`'s signature changed**, from `SpatialCoordinate(dim)` to
  `SpatialCoordinate(domain, point=None)`. It is re-exported from `uflx/__init__.py`, so
  this is a public break to call out. `codegeneration/main` also pattern-matched
  `SingleSpatialCoordinate` for quadrature substitution; that node still exists but now
  carries a domain and a point too.
- [ ] **`component` re-expands the whole quantity each call**, so reading a 3x2 Jacobian
  componentwise builds its full matrix six times. Waste, not wrongness: value equality
  dedups the graph.

## Exterior facet integrals (`ds`)

The geometry is already there. A facet's map is the reference facet included into the
cell followed by the cell's own map, which is a `ComposedParametrization`, and a facet
domain built that way gives a Jacobian, a metric, a tangential projector, a measure
factor and a normal with no new machinery. Four other things are missing.

- [ ] **Reference geometry, which UFLx cannot supply.** `AbstractEntity` has
  `sub_entities` and `sub_entity_vertices`, which are combinatorics: it knows facet 0 of
  a triangle is vertices `[1, 2]` but not where they are, and `EntityDomain`'s docstring
  says coordinates belong to whoever defines the elements. So the inclusion
  `s -> (1 - s, s)` has to come from outside. Basix has `sub_entity_geometry`; UFL makes
  it the `ReferenceFacetJacobian` terminal and lets FFCx fill it in. **This is the fork
  that shapes the rest**: a new abstract hook on `AbstractEntity`, a consumer-supplied
  `AbstractParametrization` per facet, or unexpandable terminals as in UFL.
- [ ] **The measure means nothing yet.** `dim`, `codim` and `boundary_only` are stored,
  put in `init_args` and `__repr__`, and read by nothing. `ds` and `dS` do not exist;
  `dx = Measure(codim=0)` is the only one. Adding `ds` today would change no behaviour,
  because `pull_back_to_entity` never consults the measure to decide what it is pulling
  back onto.
- [ ] **A second fan-out axis.** A kernel is generated per (cell type, local facet
  index), since the inclusion differs per facet. That generalises the cell-type fan out
  in #103.
- [ ] **`dS` needs the facet seen from both sides**, with a transition between the two
  cells' coordinates, plus the `'+'`/`'-'` restriction the spec promises
  (`LANGUAGE.md:189`). This is the one place an atlas-like structure with genuine
  transition maps is earned.

## Language gaps

- [ ] **`Coefficient.diff` raises `NotImplementedError`** (`functions.py:495`), so
  expanding `grad` of a general Lagrange coefficient fails in
  `EntityGrad.expand_geometry`. Pre-existing; it is why the mixed-mesh tests use a mass
  form.
- [ ] **`IntegralSum` does not canonicalise.** UFL's `Form` merges integrands sharing a
  domain, integral type and subdomain id via `group_form_integrals`; `IntegralSum` keeps
  terms in order and compares order-sensitively, so `a + a` stays two terms. Grouping
  belongs in `simplify`, beside the commutative-operand sorting already there.
- [ ] **`Integral * Integral` is unimplemented**, though `LANGUAGE.md:161` offers
  `(u * dx) * (v * dx)` as valid UFLx.
- [ ] **`FunctionSpace.elements` does two jobs** — several elements on one cell, and one
  element per cell type — and `__init__` requires them all to share an
  `ambient_value_shape`, so the tuple cannot express a Taylor-Hood space.
- [ ] **Multi-mesh forms are indistinguishable.** Domain equality is structural, so two
  different meshes sharing a coordinate element compare equal and
  `Integral`'s domain-agreement assert passes silently. UFL solved this with `ufl_id`.

## Naming and layering

- [ ] **`ParametrizedDomain` and `parametrized_domain` carry the generic names** while
  being specifically the finite element flavour. `FiniteElementParametrizedDomain` would
  need `finite_element_parametrized_domain()` under the class/factory rule in
  DESIGN_CHOICES.md, which is a mouthful.
- [ ] **`uflx/parametrizations.py` is not the plural of a class it contains**, bending
  the plural-module convention.

## Housekeeping

- [ ] **`.gitignore` has `docs/_build/` but the directory is `doc/_build/`**, so a docs
  build leaves untracked output.
- [ ] **Stale build artifacts** in `external/basix_uflx/build/` and
  `external/codegeneration/build/` hold old copies of removed code. Gitignored, so
  harmless, but they turn up in greps.
- [ ] **The branch stack is four deep**: `main` → #100 → #102 → #103 →
  `jhale/differential-geometry-concepts`, the last pushed with no PR of its own. A
  squash merge of any link gives `main` a new SHA for a commit the children already
  contain, which is what already required rebuilding `codegeneration/main`.

## Done

On `jhale/differential-geometry-concepts`, oldest first. The commit messages carry the
reasoning; `git log` is the record, not this list.

- `2253c82` Tidy the geometry quantities before extending them. Immutable nodes,
  one `__repr__`, one pairing helper, `ValueError` for caller mistakes.
- `2027c48` Name the metric a parametrization induces. `MetricTensor` is `g = JᵀJ`, and
  on a manifold `abs(√det g)` is `JacobianDeterminant` structurally.
- `f509a4f` Name the tangential projector. **Bug**: `J @ J⁻¹` claimed `Identity(gdim)` on
  a manifold, where the product has rank tdim.
- `6f734ef` Reduce a composed map's Jacobian before inverting it. **Bug**: no form over a
  composed domain could be expanded, so the parabola had no measure.
- `c5ab0e0` Give a domain of codimension one its unit normal, with no reference geometry
  needed and the sign left to the consumer.
- `43323ec` Ask once which domain an expression's functions are on. `extract_domain` in
  `functions.py`, shared by `Integral.domain` and the gradient pull back.
- `03c3f28` Give the spatial coordinate a domain, so `x` can be expanded at all, and say
  in the point classes why one takes a domain, the other a map, and why pulling back is
  a terminal.

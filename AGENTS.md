# AGENTS.md

Guidance for agents (Claude Code and others) working in this repository.

## What this is

UFLx is an experimental, minimal reimplementation of UFL (Unified Form Language): a
symbolic language embedded in Python for writing finite-element variational forms (e.g.
`inner(grad(v), grad(u)) * dx(mesh)`). It does not own meshes, boundary markers,
coefficient values, or assembly — those belong to a surrounding solver environment (e.g.
DOLFIN/DOLFINx). See [LANGUAGE.md](LANGUAGE.md) for the full language spec and
[DESIGN_CHOICES.md](DESIGN_CHOICES.md) for naming/API conventions — read both before
adding to the object model. [INTEGRATION.md](INTEGRATION.md) is the argument behind the
measure and integral model, written against differential geometry, and lists what is not
built.

## Commands

Install for development (editable):

```
pip install -e ".[ci,lint,docs]"
```

Lint (must pass with zero diff/errors, matches CI):

```
ruff check .
ruff format --check .          # use `ruff format .` to auto-fix
pyrefly check                  # type checker, config in [tool.pyrefly]
```

Run the core test suite:

```
pytest test/
pytest test/test_forms.py::test_name    # single test
pytest -n auto test/                    # parallel, matches CI (needs pytest-xdist)
```

`external/basix_uflx` is a separate installable package excluded from the main `uflx`
package; install and test it independently:

```
pip install external/basix_uflx
pytest external/basix_uflx/test
```

Build docs:

```
cd doc && make html
```

## Architecture

Expression trees (the symbolic core) are not plain Python object graphs — they are
backed by an explicit DAG in `uflx/graphs/` (`Graph` subclasses `networkx.DiGraph`;
nodes are `GraphNode`s exposing `successors` and `init_args`). Whole-tree rewrites
(substituting terminals, extracting real/imaginary parts, applying push-forward/
pull-back maps) go through `uflx/algorithms/replace.py`, which walks the graph and calls
`reconstruct_node` (`uflx/algorithms/reconstruct.py`) to rebuild any node whose
successors changed, by re-invoking `node.__class__(*args)` with replaced `init_args`.
Anything that needs to transform an expression (`complex.py`, `maps.py`, `geometry.py`)
is built on this `replace` primitive rather than ad hoc tree-walking — follow that
pattern for new transformations.

The object model is a layered chain of abstract base classes, each file named after its
plural concept with `Abstract*` base classes (concrete classes are defined by consumers,
e.g. `test/conftest.py`'s `LagrangeElement`, or the `basix_uflx` extension):

- `entities.py` — `AbstractEntity`: topological mesh entities (point/interval/
  triangle/.../hexahedron), defined recursively via their sub-entities.
  `sub_entity_inclusion` and `facet_outward_normal` are the two pieces of reference
  geometry an entity is asked for — the map carrying a sub-entity's coordinates into the
  cell's, and which way is out of the cell across a facet — since the combinatorics say
  which vertices a facet has and not where they are, nor which side the interior is on.
  Their defaults raise; a consumer supplies them.
- `finite_elements.py` — `AbstractFiniteElement` / `AbstractMappedFiniteElement`: what
  basis functions look like on a cell; depends on `entities` and `maps`.
- `maps.py` — `AbstractValueMap` (e.g. `IdentityValueMap`): push-forward/pull-back of
  field values between an entity's coordinates and ambient coordinates, implemented as
  graph rewrites. Which map a field needs follows from the degree of the form it
  represents — identity for a 0-form, `J⁻ᵀ` for a 1-form, `J / det J` for an (n-1)-form
  — and the module docstring pairs each with the `geometry.py` quantity it is. Applying
  a parametrization to a *point* is a different operation, named `ImagePoint` and
  `PreimagePoint` in `geometry.py`.
- `domains.py` — `AbstractDomain`: a set you can integrate over. A coordinate domain's
  points are tuples of numbers (`RD`, `EntityDomain`); a parametrized domain is the
  image of a map out of one. `AbstractParametrization` is that map, evaluated by `value`
  and differentiated by `jacobian`; `IdentityParametrization` is an entity domain's and
  `RD`'s own. `AbstractParametrizedDomain` is what geometry and a measure ask for: a map
  out of a parameter region, a reference cell being one kind of region and a region of
  `RD` another, so being integrable over does not imply being made of cells. The actual
  mesh stays external to UFLx. A parametrization says where it lands (`target`) and not
  only in how many coordinates, which is what lets a domain built from one report its
  `coordinate_domain`: the inclusion of a facet lands in the cell's coordinates where a
  mesh's own map lands in the ambient ones. `UnionDomain` is a domain presented as
  several parametrized pieces, which is what the facets of a mesh are. A domain holds
  its `parametrizations` rather than keying them by parameter region, since several
  parametrizations may share one; needing exactly one is a question `geometry.py` asks
  of the point a quantity is read at. `LabelledDomain` is a domain and an opaque label a
  consumer selects it by — a marked subdomain being a domain and not an id on a measure
  — delegating everything geometric, and `unlabelled` takes the labels off, which an
  integral does before asking whether a function belongs to the domain it is over, and
  which anything asking a structural question of a domain does first: `ComposesADomain`
  is the protocol for a domain presented as another with a map applied, and a label is
  not one.
- `function_spaces.py` — `AbstractFunctionSpace`: standard FE spaces (domain + element),
  constant spaces (shape + scalar type), or non-FE spaces (domain + shape, no element).
  Do not construct `Argument`/`Coefficient` directly from an element — they must come
  from a `FunctionSpace`.
- `variables.py` — `AbstractVariable`: a point of a domain without naming one, which an
  integral's dummy variable is. A variable in a coordinate domain has coordinates, as
  `PointComponent`s, and that is where expanding a parametrization written in closed
  form bottoms out; one on a mesh or on a union has none and says to pull back or split
  first. It sits beside `domains.py` rather than in `functions.py` because a domain's
  parametrizations take a variable and the geometry built on them does too; the two are
  peers and refer to each other, which is why `domains.py` names this module only in its
  type hints. `ReadAtAVariable` is the protocol for a node told later where it is read.
- `functions.py` / `basis_functions.py` — `AbstractFunction` and basis functions
  evaluated at points, in either an entity's coordinates or ambient coordinates.
  `extract_domain` is which domain an expression is read on, read off the parts of it
  that implement `ReadOnADomain`: a function from its space, a geometric quantity from
  the domain it differentiates, a restriction from the domain it restricts to. A part
  already positioned says nothing, which is the same rule an integral's own check uses.
- `parametrizations.py` — `FiniteElementParametrization` and `ParametrizedDomain`, where
  a map is described by a finite element per cell, plus
  `ComposedParametrization`/`ComposedDomain`, which apply two maps in turn via the chain
  rule, and `AffineParametrization`, a closed form map a consumer can state.
  `sub_entity_domain(cell, dim, index)` is one sub-entity of a cell as a domain: the
  sub-entity's coordinate domain carried into the cell's by the inclusion the cell
  supplies, a composition like any other. `sub_entities_of(domain, codim)` builds the
  union whose pieces are a sub-entity's inclusion followed by the cell's own map, which
  at codimension one is the domain an exterior facet integral is over; a label marks
  each piece, so it survives the split onto one. `SubEntityInclusion` is the affine
  inclusion that says which sub-entity it is, and `included_sub_entity` reads that back
  off a piece or a pulled back integral's geometry, which is what identifies a kernel
  along with the cell type. Evaluating an element map is an interpolation sum, so this
  sits above `function_spaces`/`basis_functions` rather than in `domains`, which must
  not know how a parametrization is described.
- `expressions.py` — `AbstractExpression`: the base of every symbolic node
  (`BinaryOperator`, `UnaryOperator`, terminals). What it carries is `value_shape`,
  `successors` and `init_args`, and nothing else: what a node can answer about itself
  alone. The rest is asked of a node where it makes sense to ask, through a protocol
  rather than an attribute every node has to have -- which domain an expression is read
  on is `ReadOnADomain` and `functions.extract_domain`, being real or complex is
  `complex.ComplexValued`, and neither is something a `Sum` can say on its own. There
  are no free indices; index notation is not built.
- `operators.py`, `tensors.py`, `points.py`, `geometry.py`, `complex.py` — concrete
  operators (`inner`, `grad`, ...), tensor/vector construction, point sets, geometric
  quantities (spatial coordinates, Jacobians, ...), and complex-number support
  (`re`/`im` via the `ComplexValued` protocol), all built as `AbstractExpression`
  subclasses / graph rewrites over them.
- `restrictions.py` — `Restriction`/`restrict`: a field read at the point another
  domain's parametrization gives, `ι*u`, which is how a field on a mesh appears in an
  integral over its facets. It lives in the integrand and not in the measure, which is
  why a gradient under one keeps the cell's inverse Jacobian rather than the facet's
  tangential one. `evaluate_at` is the shared pass that tells an expression where it is
  read, over the `ReadAtAVariable` nodes (`variables.py`) -- a function, a geometric
  quantity or a restriction, each built before there is a point and told later. A field
  is read through one map and never two: restricting along the composition is what a
  second restriction would be, so a sub-entity of a sub-entity is reached with
  `sub_entities_of` at the cell's codimension. On an interface the map is named by a
  side, `restrict(u, gamma.side("-"))`, and `jump`/`avg` are the two sides' difference
  and mean.
- `interfaces.py` — `Interface`: a domain two others meet along, with a *presentation*
  of it per side, which is the domain parametrized by that side's own map, a facet piece
  or a face of a box. The atlas `UnionDomain` says it is not. `Side` is the interface as
  one side presents it, a domain in its own right that a restriction or an
  `OutwardNormal` takes; its `chart_point` is where a point of the interface sits in the
  side's own region, the integral's variable on the primary side and the chart
  transition `ψ₋⁻¹(ψ₊(X))` on any other, written as a `PreimagePoint` of an `ImagePoint`
  with the side's map named by the side (`SideParametrization`). An interface of unions
  is `MadeOfPieces`, one pair of pieces at a time. There is no second spelling: UFL's
  `dS` is a mesh's facets presented twice, and a mesh meeting a submesh of its entities
  is the entities presented by each, a mixed-domain form being each field restricted to
  its side. The protocols the geometry and the integrals ask of a side are in
  `domains.py`, since this module sits above `geometry.py`.
- `legacy.py` — the legacy UFL spellings (`dx`, `ds`, `dS`, their `(3)` forms, `jump`,
  `avg` and `facet_normal`), which infer the domain and insert the restrictions to build
  the native form. A side is named by `plus(u)`/`minus(u)`, a placeholder the measure
  resolves once it has built the interface.
  `add_plus_minus_restriction(AbstractExpression)` puts a `__call__` on the class named
  so UFL's `u("+")` works, the one place anything reaches into the core's classes, done
  on request and naming the class rather than on import so that import order decides
  nothing and the change is written where it is made; `remove_plus_minus_restriction`
  undoes it, and a type checker accepts only `plus(u)`. Opt in by importing it; nothing
  in `uflx` imports it, since inference is what INTEGRATION.md note 2 argues against.
- `integrals.py` — `Measure`, a domain paired with the density on it, built by
  `dx(domain)`; `Integral` (`expr * measure`), which takes its domain from its measure;
  and `IntegralSum`, which adding integrals gives. An integral over a domain of several
  cell types is pulled back as one integral per cell type, since each has its own
  coordinate domain; `Integral.on_cell_type` and `split_by_cell_type` do that, carrying
  the measure onto the restricted domain. `split_by_piece` is the same fan out over the
  pieces of a union, run first since a piece names one parametrization and so one cell
  type. `IntegralSum` splits its own terms both ways, since the pull back walks a sum
  leaves first and so reaches each integral with no chance to split it on the way. The
  pull back is onto the parameter region of the domain's one parametrization, which is
  what lets an interface pull back like a facet. There is no `ds` or `dS`: an exterior
  facet integral is this measure over `sub_entities_of(mesh, 1)` with the fields
  restricted, and an interior one is it over the interface of those facets with
  themselves, `interface(facets, (("+", facets), ("-", facets)))`, with the fields
  restricted to its sides, which [INTEGRATION.md](INTEGRATION.md) §3 and §5 argue for.

Symbolic (Gateaux) differentiation and expression/integral transformation between domain
configurations are the two core transformation procedures UFLx provides on top of this
object model (see LANGUAGE.md §4) — both are lazy, expressed as graph operators rather
than eagerly evaluated.

## Python style

Follows FEniCS project conventions (see [dolfinx](https://github.com/FEniCS/dolfinx)'s
`AGENTS.md`), adapted for this pure-Python, MIT-licensed package:

- **Formatting/linting**: `ruff check` and `ruff format --check`, configured in
  `pyproject.toml`. Line length 100, 4-space indent. Rule set includes pydocstyle (`D`,
  Google convention), pycodestyle, pyflakes, isort, pyupgrade,
  flake8-import-conventions, NumPy-specific rules, and a few more (`RUF`, `FLY`, `LOG`,
  `ISC`). Run the formatter/linter locally before calling a change done — don't rely on
  CI to catch formatting.

- **Import order** (ruff's isort, default grouping): future → standard-library →
  third-party (`networkx`, `pytest`) → first-party (`uflx`) → local-folder.

- **Docstrings**: Google style (`Args:`, `Returns:`, etc.), required on modules, classes
  and methods (`D` rules are enforced, `__init__.py` is exempted from unused imports
  only — not docstrings).

- **Type hints**: used throughout and checked with `pyrefly` (`[tool.pyrefly]` in
  `pyproject.toml`); not currently enforced by a ruff annotation rule, but new public
  functions/methods should still be annotated to match the rest of the codebase.

- **File header**: most files under `uflx/` (the core object-model layer —
  `entities.py`, `finite_elements.py`, `domains.py`, `function_spaces.py`,
  `functions.py`, `variables.py`, `expressions.py`, `operators.py`, `integrals.py`,
  `basis_functions.py`, `__init__.py`) start with:

  ```python
  # Copyright (C) <year> <author(s)>
  #
  # This file is part of UFLx (https://www.fenicsproject.org)
  #
  # SPDX-License-Identifier:    MIT
  ```

  followed by a module docstring. Some newer/infrastructure files (`uflx/graphs/` and
  its `algorithms/` submodule, `maps.py`, `tensors.py`, `points.py`, `geometry.py`,
  `complex.py`, `utils.py`) currently omit it — match the header convention for new
  files in the core object-model layer; add your name/year rather than replacing
  existing authors when making a substantive change to a file that has one.

## Conventions (see DESIGN_CHOICES.md)

- Module files are named as the plural of what they contain (e.g. `domains.py`)
  precisely so the module (`uflx.domains`) and a same-named factory function
  (`uflx.domain`) can coexist.
- A class's factory/constructor function has the same name as the class, lowercased
  (mid-word capitals become `_x`, e.g. a hypothetical `FooX` class pairs with
  `foo_x()`).
- `__init__` methods do not choose a value for the caller; add a factory or named
  instance instead (as `dx(domain)` is a factory for a measure). A `None` default for a
  field not yet known is fine.

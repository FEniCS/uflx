# The design of the UFLx core

Authors: Jack S. Hale, Matthew Scroggs, Michal Habera

This note is the argument behind [LANGUAGE.md](LANGUAGE.md). LANGUAGE.md is the
specification and is normative. This note explains the six ideas the
specification is built on, the decisions that shape every section of it, and
what is deliberately left out. Read this note first, then use LANGUAGE.md as a
reference. Every heading below names the sections of LANGUAGE.md that carry
the rules.

## In one sentence

> UFLx is a typed symbolic language for writing variational forms.

A UFLx program is Python code that builds a graph of immutable objects. The
graph is the form. Nothing is evaluated or simplified while the graph is
built. A small number of passes, each doing one thing, rewrite the graph into
the shape a form compiler consumes, and everything a form compiler needs is
in the graph.

## Reading guide

| If you want to                               | Read in LANGUAGE.md                                   |
| -------------------------------------------- | ----------------------------------------------------- |
| See a whole problem written in the language  | A tour: the Poisson problem                           |
| Know what an object may and may not do       | Rules of the language                                 |
| Look up one node class                       | The section named after the class, and the appendices |
| Follow a form from the user to a kernel      | Core transformations, Concept: kernels                |
| Add something without changing the core      | Extending the core                                    |

## Idea 1: the graph is the form, and only passes change it

Every object extends one base class and obeys ten rules (The base object
`Base`). The ones that matter for a reader:

- An object is built from a tuple of construction arguments, is immutable, and
  is equal to another object exactly when the classes and the construction
  arguments are equal. Equality is structural, so `a + b` and `b + a` differ.
- The edges of the graph are the base objects among the construction
  arguments. There is no information outside the graph.
- Validation happens in constructors and nowhere else. A free function such as
  `transpose` only reads static properties of its arguments and builds a fixed
  composition; it never inspects structure and never simplifies.
- A pass is a pure function on a graph. Substitution renames bound indices and
  bound functions to avoid capture, as in any language with binders.

The consequence that organises the whole document is in Node class or free
function. A concept gets a node class only if no fixed composition of
existing classes builds it, or if a numerical algorithm the user chooses
evaluates it, or if a core transformation must meet it as written. Everything
else is a free function. So `transpose` is a composition, `Division` is a
node because `f * g ** -1` rounds differently, and `Integral` over a mesh is
a node because the lowering should happen in a pass, not at construction.

## Idea 2: an expression is an array-valued differential form with a static type

Every expression carries six static properties, derived at construction and
never changed (Expressions and static type):

| Property        | What it records                                                  |
| --------------- | ---------------------------------------------------------------- |
| value shape     | the array axes, such as `(3,)` for a vector                      |
| free indices    | the open indices, in order of first occurrence, each with extent |
| scalar type     | real or complex                                                  |
| support         | the domain the expression lives on, or none for a constant       |
| degree          | the differential-form degree of each slot, such as `(1,)`        |
| argument spaces | the space of each bound argument, for a variational form         |

Two choices here remove most of the operators of legacy UFL.

First, tensor algebra is index notation. `a[i] * b[i]` is the pointwise
product with a free index, `index_sum(a[i] * b[i], (i,))` is the
contraction, and `component_tensor(A[k, l], (l, k))` is the transpose. An
`EllipsisIndex` covers any number of axes, so `inner`, `dot` and `outer` are
one composition each for every rank. Summation is never implicit. Scope is
lexical: an index refers to the innermost `ComponentTensor` or `IndexSum`
that binds it.

Second, the calculus is exterior calculus. An expression of degree `(k,)` is
a k-form. `ExteriorDerivative` is the only derivative, `Product` is the wedge
product, `TensorProduct` adds a slot for double forms such as a metric or a
Regge field, and `Components` turns the first slot into an array axis in the
coordinate frame. The gradient of a function `u` on a domain `D` is then
`components(exterior_derivative(u, D), D)`, and the same line is the surface
gradient on a surface in three dimensions, because `Components` projects the
ambient frame onto the tangent space. There are no dimension-specific
operators, which is the point of Basic principles.

## Idea 3: domains are leaves or decompositions, and maps are how things move

```
AbstractDomain                      the support of an expression
├── AbstractCoordinateDomain        a leaf: RD(n), or a reference cell from a library
└── DecomposedDomain                the union of the images of families of maps
          │ pieces
          ▼
      DomainMap(P, psi, Q)          psi: P -> Q, with the image a cell, a facet, a point
          ▲
          │ second operand
      PullBack(f, map)              moves an expression on Q to P
```

The core has no way to write a set except as the image of a map, and a chain
of maps must start somewhere. It starts at a coordinate domain, where the
coordinates are canonical: the coordinate function is the identity, the
volume form is `dX_0 ^ ... ^ dX_{n-1}`, and `Integral` is the Lebesgue
integral (Domains). The core supplies `RD(n)`; a library supplies reference
cells.

Everything else is a `DecomposedDomain`: a union of families of `DomainMap`s
from a domain into a domain. A mesh is one piece per cell type, each a family
of maps from the reference cell indexed by the cell index, and the decomposed
domain then has a topology. A set of sensor points is a family of maps from
`RD(0)`. A facet domain is a family of maps from the reference facet. Pieces
may themselves have a decomposed parameter domain, so macro elements are
nothing special.

Expressions never move between domains except by `PullBack`. A function on a
mesh pulled back by a cell map is the function on that cell. The pull-back is
the only place where a Jacobian can appear, and it appears there as a
composition: the Jacobian is `components(exterior_derivative(psi, P), P)`,
its inverse or pseudo-inverse is the one node `Inverse`, and the density of a
cell of any dimension is `square_root` of the determinant of `J^T J`. One
rule serves every dimension and every codimension (Expanding pull-backs).

Orientation is a consequence rather than a datum. A coordinate domain is
oriented by its coordinate order, a map transports that orientation, and
a 0-form integrates against a density that needs none. Where a k-form over a
surface needs consistent orientations, the mesh library supplies a sign per
cell as an ordinary expression (Concept: orientation of a domain).

## Idea 4: a facet is a cell seen from a side

The core does not number facets, edges or vertices of a cell. A topology of
sub-entities names, for each of its cells, the cell it is seen from on each
side and the pairing of vertices, as two tables from the mesh library (Sides
of a topology). That one construction gives boundary facets, interior facets
with two sides, marked subdomains, submesh interfaces, periodic boundaries,
vertices for point loads and edges where three surface triangles meet, each
as a topology whose groups list their sides.

Restriction is then not a node. `restriction(u, E, "+")` is the pull-back of
`u` by a lazy `SideMap`, and a pass rewrites it on each piece as two nested
pull-backs: first by the cell on that side, selected through the table by
`gathered_cell_domain`, then by the affine inclusion of the reference facet
into that cell through the paired vertices (Restriction). Two restrictions of
one field land on the same reference facet with the same parameter, so a jump
is a subtraction. The outward normal of the cell on a side is the inverse
Jacobian applied to the reference facet normal, which the paired vertices
determine (Normal). The interior facet example in LANGUAGE.md follows one
`dS` form through all of this to a kernel.

## Idea 5: a function is a terminal, and a form binds it

A `Function` is a member of a space and nothing more. It is not a test
function, a trial function or a coefficient by itself. A `Form` binds some of
the functions in an expression as its arguments, by position, and the same
function may be bound in one form and fixed in another, which is what the
Gateaux derivative and the action of a form need (Function, Form).

The form also declares the linearity of each argument. Binding `v` says the
expression is linear in `v`; binding `conjugate(v)` says it is linear in the
conjugate, so antilinear in `v`. The argument spaces record this through the
complex conjugate space `conjugate_space(V)`, which is `V` itself when `V` is
real. So a load, a function in `dual(V)`, adds to a residual only when their
linearities agree, and there is no convention about which argument is
conjugated (Argument spaces, ConjugateSpace, DualSpace). For real problems
none of this is visible.

Mixed problems need no mixed space. An argument may be a tuple of functions
in different spaces, on different domains, with different degrees, and a
library assembles the tuple as one system or in blocks (Form).

## Idea 6: six lazy transformations do all the work

```
form(integral(f, D), (v, u))                  what the user wrote
  │ expand_integrals          one Integral per piece, an index sum over cells, the density
  │ expand_pull_backs         pull-backs to the terminals; Components becomes K^T times reference components
  │ apply_derivatives         d to the terminals; d of the coordinate is dX
  │ expand_components         components of dX are the identity; interior products
  │ apply_actions             substitute into bound arguments
  │ reduce_component_tensors  component_tensor(e, (i,))[j] becomes e with i := j
  ▼
form(index_sum(integral(g, K), (c,)), (v, u)) g on the reference cell K, reading cell data by gathers
```

Each transformation evaluates one node class, has a rule per operand class,
and leaves the node where no rule applies (Core transformations, Appendix:
Pass table). They terminate and may run together. What remains after them is
exactly what a library owns: a function pulled back by a cell, the exterior
derivative of such a pull-back, and the `Components` of that derivative,
which are the basis functions and their reference gradients that the library
tabulates. Pull-backs nest rather than compose, so a function pulled back by a
cell and then by a facet inclusion is still "this function, on this cell, at
these facet points".

The result is the shape of a kernel without a kernel interface being
designed. The index sum over the cell index is the assembler's loop. The
gathered indices through the cell index, such as the vertex positions and the
side tables, are the data the assembler packs. The rest is the kernel, and
the inputs of the kernel are the leaves that remain (Concept: kernels).

## The tour

The Poisson problem on a mesh of two triangles, from LANGUAGE.md:

```python
K = ReferenceSimplex(2)
T = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
c = cell_index(T, 0)

X = coordinate(K)
k, w = index("k"), index("w")
p = component_tensor(y[gathered_index(cell_vertices, (c, w)), k], (c, w, k))
x = component_tensor(
    p[c, 0, k] + (p[c, 1, k] - p[c, 0, k]) * X[0] + (p[c, 2, k] - p[c, 0, k]) * X[1],
    (k,),
)
D = decomposed_domain((domain_map(K, x, RD(2)),))

u, v, f = function(V, "u"), function(V, "v"), function(V, "f")

i = index("i")
grad_u = components(exterior_derivative(u, D), D)
grad_v = components(exterior_derivative(v, D), D)
a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), D), (v, u))
L = form(integral(f * v, D), (v,))
```

Every idea above is visible. The mesh is a family of affine maps from the
reference triangle, built by gathers through the cell index (ideas 1 and 3).
The gradient is the exterior derivative turned into an array (idea 2). The
functions are plain members of a space, and the forms bind them (idea 5).
`integral(f, D)` is one lazy node over the mesh, which the passes turn into a
sum over cells of integrals over the reference triangle (idea 6). On a mesh
of triangles and quadrilaterals, only `D` changes.

## The decisions that shape everything

Each of these is argued in the named section of LANGUAGE.md. They are listed
here so that a reader knows which choices are load-bearing.

- **No simplification at construction**, so the graph keeps what the user
  wrote and a form compiler can use it directly (No simplification during
  construction).
- **A fixed composition, never a node class**, unless no composition exists,
  an algorithm is the user's choice, or a transformation must meet the
  concept as written (Node class or free function).
- **Index notation and ellipsis indices** instead of tensor operators, so
  every tensor operation is one composition for every rank (Index,
  EllipsisIndex).
- **Exterior calculus with one derivative**, and `Components` as the bridge
  to arrays, so there are no dimension-specific operators (Degree,
  ExteriorDerivative, Components).
- **One kind of domain with two forms**, leaf and decomposed, and `DomainMap`
  as the only way between them, so every node that takes a domain takes any
  domain (AbstractDomain, DomainMap).
- **`Integral` is lazy over a decomposed domain**, so the user's integral
  survives to the compiler and the lowering lives in one pass (Integral,
  Expanding integrals).
- **The core numbers nothing**: no facets, no sub-entities, no argument
  numbers. Side tables and paired vertices replace facet numbers, and
  positions in `arguments` replace argument numbers (Sides of a topology,
  Form).
- **Pull-backs nest**, so the cell a function is read on is never lost and a
  library sees only pull-backs by cells (Expanding pull-backs).
- **Linearity is declared, not inferred**, and recorded as a conjugate space,
  so complex forms are typed and real forms are untouched (ConjugateSpace).
- **Constants have no support**, and a constant with free indices is written
  with `Zero`, so a literal is the same object on every domain (Support,
  Zero).
- **Extension is by passes and leaves**, never by properties on core nodes,
  so a new question such as units touches no core class (Extending the core).

## What is outside the core, on purpose

- Mesh data: coordinates, connectivity, side tables, marked regions, point
  location and the inverse of a cell map are tables and arrays that a mesh
  library supplies as leaves (Sides of a topology, Points and point
  evaluation).
- Finite elements: how a function is represented on a cell, including the
  Piola maps, is the library's choice behind the pull-back by a cell
  (AbstractFunctionSpace).
- Quadrature: a pass with the rule as an argument, not a node (Example:
  quadrature).
- Elementary functions beyond `SquareRoot`, conditionals, and `min` and
  `max`: language extensions (Elementary functions).
- The Gateaux derivative: planned as a core transformation, with only its
  shape stated so far (Gateaux derivative).
- Ragged indexing, such as a different number of quadrature points per cut
  cell, and interfaces that are not whole sub-entities: not yet covered
  (Extents, Sides of a topology).

## Correspondence with legacy UFL

| Legacy UFL                               | UFLx                                                    |
| ---------------------------------------- | ------------------------------------------------------- |
| `Mesh`                                   | a `DecomposedDomain` with a topology                    |
| `SpatialCoordinate(mesh)`                | `coordinate(D)`                                         |
| `grad(u)`                                | `components(exterior_derivative(u, D), D)`              |
| `dx`, `ds`, `dS`, `dx(1)`                | `integral(f, D)` over the cells, a facet topology, an interior facet topology, a marked topology |
| `u('+')`                                 | `restriction(u, E, "+")`                                |
| `FacetNormal(mesh)`                      | `normal(D, E, "+")`                                     |
| `TestFunction`, `TrialFunction`          | functions bound by `form(e, (v, u))`                    |
| `Coefficient`                            | a function not bound by the form                        |
| `Cofunction`, `Coargument`               | a fixed and a bound function in `dual(V)`               |
| `MixedElement`, `split`                  | an argument that is a tuple of functions                |
| `action`, `adjoint`                      | `action(f, position, g)`, rebinding with `conjugate`     |
| `derivative`                             | the Gateaux derivative, not yet covered                 |
| `Real` space, `Constant`                 | a space with no domain, an expression with no support   |

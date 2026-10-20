# Core UFLx language definition

Authors: Jack S. Hale, Matthew Scroggs, Michal Habera

This note describes the smallest core language behind UFLx. It intentionally
avoids implementation details, compiler algorithms, and broad redesign ideas.
It assumes knowledge of differential geometry, and functional analysis. Some
advanced differential geometry concepts are explained in more detail.

This document is aimed at users and developers with extensive knowledge of
domain specific languages for variational forms (e.g. UFL, ngsolve).

## AI usage statement

This document has been co-created with Claude Opus 5.5 which has presumably
been trained on a large corpus of material related to modern differential
geometry, algorithms, functional programming and the source code of Firedrake,
DOLFINx, UFL, FFCx, TSFC, GEM and many other packages. That said, Claude Code
was specifically instructed not to look at any materials or source code
directly. The quality of the document and resulting design is substantially
higher than what could be achieved through 'old-fashioned' means. Sections with
fully LLM generated content are explicitly marked.

## Contents

- [Scope](#scope)
  - [What UFLx is and isn't](#what-uflx-is-and-isnt)
  - [Basic principles](#basic-principles)
  - [A tour: the Poisson problem](#a-tour-the-poisson-problem)
- [Rules of the language](#rules-of-the-language)
  - [Basic implementation principles](#basic-implementation-principles)
  - [The base object `Base`](#the-base-object-base)
  - [Leaves and subclassing](#leaves-and-subclassing)
  - [Properties and free functions](#properties-and-free-functions)
  - [No simplification during construction](#no-simplification-during-construction)
  - [Passes](#passes)
  - [Kinds of code that produce objects](#kinds-of-code-that-produce-objects)
  - [Enforcement](#enforcement)
- [Coordinate domains](#coordinate-domains)
  - [AbstractCoordinateDomain](#abstractcoordinatedomain)
  - [RD](#rd)
  - [Example: the unit simplex](#example-the-unit-simplex)
- [Expressions and static type](#expressions-and-static-type)
  - [`AbstractExpression`](#abstractexpression)
  - [Index](#index)
  - [EllipsisIndex](#ellipsisindex)
  - [Extents](#extents)
  - [AbstractIndexValue](#abstractindexvalue)
  - [GatheredIndex](#gatheredindex)
  - [Free indices](#free-indices)
  - [Scalar type](#scalar-type)
  - [Support](#support)
  - [Degree](#degree)
  - [Expressions as differential forms](#expressions-as-differential-forms)
  - [Node class or free function](#node-class-or-free-function)
- [Algebra](#algebra)
  - [Literals](#literals)
  - [Zero](#zero)
  - [Arithmetic](#arithmetic)
  - [Sum](#sum)
  - [Product](#product)
  - [Power](#power)
  - [Negation and subtraction](#negation-and-subtraction)
  - [Division](#division)
  - [Absolute](#absolute)
  - [Conjugate](#conjugate)
  - [TensorProduct](#tensorproduct)
  - [Elementary functions (very incomplete)](#elementary-functions-very-incomplete)
  - [Tensor](#tensor)
  - [Indexed](#indexed)
  - [ComponentTensor and IndexSum](#componenttensor-and-indexsum)
  - [Tensor algebra with EllipsisIndex](#tensor-algebra-with-ellipsisindex)
  - [Transpose and matrix product](#transpose-and-matrix-product)
  - [Inverse](#inverse)
  - [Operator overloads](#operator-overloads)
  - [Numbers in expressions](#numbers-in-expressions)
- [Calculus on a coordinate domain](#calculus-on-a-coordinate-domain)
  - [Coordinate](#coordinate)
  - [Valid pair](#valid-pair)
  - [ExteriorDerivative](#exteriorderivative)
  - [Components](#components)
  - [Volume form](#volume-form)
  - [Integral](#integral)
- [Maps between domains](#maps-between-domains)
  - [Domains](#domains)
  - [ParameterisedDomain](#parameteriseddomain)
  - [Example: a paraboloid](#example-a-paraboloid)
  - [Example: an affinely mapped triangle](#example-an-affinely-mapped-triangle)
  - [PullBack](#pullback)
  - [Concept: orientation of a domain](#concept-orientation-of-a-domain)
  - [Concept: integral over a domain](#concept-integral-over-a-domain)
  - [Example: orientation of an edge and of a triangle](#example-orientation-of-an-edge-and-of-a-triangle)
  - [Example: integrating over an edge](#example-integrating-over-an-edge)
  - [Example: integrating over the boundary of a triangle](#example-integrating-over-the-boundary-of-a-triangle)
  - [Points and point evaluation](#points-and-point-evaluation)
  - [Example: quadrature](#example-quadrature)
  - [Example: integrating over a cut cell](#example-integrating-over-a-cut-cell)
  - [Example: the centroid of a triangle](#example-the-centroid-of-a-triangle)
- [Topology](#topology)
  - [AbstractCellType](#abstractcelltype)
  - [Example: the simplex as a cell type](#example-the-simplex-as-a-cell-type)
  - [AbstractTopology](#abstracttopology)
  - [Sides of a topology](#sides-of-a-topology)
  - [Example: a static topology](#example-a-static-topology)
- [Reference cells](#reference-cells)
  - [Reference cell](#reference-cell)
  - [Example: the unit simplex as a reference cell](#example-the-unit-simplex-as-a-reference-cell)
- [Decomposed domains](#decomposed-domains)
  - [DecomposedDomain](#decomposeddomain)
  - [Restriction](#restriction)
  - [Integral over a decomposed domain](#integral-over-a-decomposed-domain)
  - [Normal](#normal)
  - [Concept: kernels](#concept-kernels)
- [Function spaces](#function-spaces)
  - [AbstractFunctionSpace](#abstractfunctionspace)
  - [Function](#function)
- [Variational forms](#variational-forms)
  - [Argument spaces](#argument-spaces)
  - [DualSpace](#dualspace)
  - [Form](#form)
  - [Action](#action)
  - [Example: the Laplace-Beltrami operator](#example-the-laplace-beltrami-operator)
  - [Example: the flux through a surface](#example-the-flux-through-a-surface)
- [Core transformations](#core-transformations)
  - [Expanding pull-backs](#expanding-pull-backs)
  - [Applying derivatives](#applying-derivatives)
  - [Expanding components](#expanding-components)
  - [Applying actions](#applying-actions)
- [Extending the core](#extending-the-core)
  - [Language extensions](#language-extensions)
  - [Example: interpolation](#example-interpolation)
  - [New questions](#new-questions)
  - [New operations](#new-operations)
  - [New subclasses](#new-subclasses)
  - [New transformations](#new-transformations)
- [Appendix: Class hierarchy](#appendix-class-hierarchy)
- [Appendix: Free functions](#appendix-free-functions)
- [Appendix: Properties](#appendix-properties)
- [Appendix: Static type table](#appendix-static-type-table)
- [Appendix: Pass table](#appendix-pass-table)

## Scope

### What UFLx is and isn't

In one sentence:

> UFLx is a typed symbolic language for writing variational forms.

The UFLx core language (henceforth the core) should not use finite element or
numerical analysis terminology - differential geometry, and functional
analysis, and their intersection, without finite element (numerical analysis)
concepts, should be all that is required. This draft of the scheme requires
knowledge of calculus on intrinsic manifolds, i.e. a coordinate-free
representation. This provides a very compact and dimension-free representation
of many concepts already in UFL.

Finite element and numerical analysis concepts, when necessary, will be added
on top of the core, or provided by external libraries (Basix being the
'in-house' example). Basix and DOLFINx and their assumptions should not be
visible in the core.

UFLx is *not* a symbolic language for reasoning about problems in functional
analysis or differential geometry. Hence, the core should only contain the
minimal concepts and relations necessary for writing a wide range of
variational forms as code for subsequent algorithmic transformation.

A UFLx program is ordinary Python code that constructs symbolic objects and
connects them through a graph. The primary results are usually variational
forms, such as bilinear forms, linear forms and functionals. Python code will
be used to motivate examples but a full concrete implementation is out of
scope.

UFLx does not handle concrete data or numerical algorithms; e.g. actual meshes,
boundary markers, numerical values of coefficients, or assembly procedures.
Those belong to the surrounding problem-solving environment, such as DOLFINx,
Firedrake, FIAT and Basix. UFLx instead specifies a set of leaf nodes, and the
external library provides the concrete data behind the leaf nodes.

### Basic principles

- All core object definitions should be closely related to their strict
  mathematical definition.
- Extra definitions, properties or terminology outside of the established
  mathematical fields noted must be avoided.
- Strictness and expressivity, not brevity, is key in the core - brevity can be
  added later with higher-level APIs. The cost of writing code against a
  precise lower-level API has massively decreased with LLMs, and so the end
  library could potentially be more concise than it was in the past.
- The mathematical fields provide a rich foundation for writing variational
  problems in a uniform way across geometric dimensions, i.e. 'special cases'
  for particular dimensions do not exist, and so should not exist in the core.
- The core does not contain algorithms for the convenience of a form
  compiler. A form compiler states the forms that the form compiler accepts.
  A pass outside the core must lower a form to those forms and must produce
  only node classes of the core and of the
  [language extensions](#language-extensions) that the form compiler accepts.
  For example, a very basic form compiler might only accept integrals over a
  reference cell, such as a unit simplex, as in each cell term of a legacy UFL
  `dx` integral.

> Legacy note: TSFC lowers a UFL form to GEM, which repeats much of the
> algebra of UFL, partly because UFL simplifies during construction. The UFLx
> core does not simplify during construction, so my impression is that a form
> compiler could use a lowered UFLx form directly, adding only leaves for
> tables and kernel inputs. I have not checked this with the TSFC or FFCx
> developers.

### A tour: the Poisson problem

The bilinear form $\int_\Omega \nabla u \cdot \nabla v \,\mathrm{d}x$ and the
linear form $\int_\Omega f v \,\mathrm{d}x$, on a square divided into two
triangles, touch most chapters of the core. Each step links to the section
that defines it.

```python
# Topology: two triangles, and the mesh as a decomposed domain.
K = ReferenceSimplex(2)
T = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())

# The index of the cells of type 0 of T, the only cell type, K.
c = cell_index(T, 0)

# The affine map x of each triangle, gathered from the vertex positions y
# through the table cell_vertices. The index c of x gives D the topology T.
X = coordinate(K)
k, w = index("k"), index("w")
p = component_tensor(y[gathered_index(cell_vertices, (c, w)), k], (c, w, k))
x = component_tensor(
    p[c, 0, k] + (p[c, 1, k] - p[c, 0, k]) * X[0] + (p[c, 2, k] - p[c, 0, k]) * X[1],
    (k,),
)
D = decomposed_domain((parameterised_domain(K, x, RD(2)),))

# Functions in a space of a library, such as Lagrange of polynomial degree 1.
u, v, f = function(V, "u"), function(V, "v"), function(V, "f")

# The integrands on D, with the gradient in the frame of RD(2).
i = index("i")
grad_u = components(exterior_derivative(u, D), D)
grad_v = components(exterior_derivative(v, D), D)
a_D = index_sum(grad_u[i] * grad_v[i], (i,))
L_D = f * v

# The forms bind the arguments.
a = form(integral(a_D, D), (v, u))
L = form(integral(L_D, D), (v,))
```

- **Topology and geometry.** `T` is a [static
  topology](#example-a-static-topology) with the single cell type `K`, and `c`
  the index of the cells of type 0, the cells of type `K`. `y` and
  `cell_vertices` are runtime data of the mesh, an array of the vertex
  positions and an [`AbstractIndexValue`](#abstractindexvalue) table, read
  through a [gathered index](#gatheredindex). `D` is the mesh, a [decomposed
  domain](#decomposeddomain) of one family of triangles, and `D.topology` is
  `T`, read from the index `c` of `x`.
- **Functions.** `V` is an [`AbstractFunctionSpace`](#abstractfunctionspace) of
  a library on `D`, with the degree `()` and the value shape `()`. `u`, `v` and
  `f` are [functions](#function). None is a test or trial function by itself.
- **Calculus.** The gradient is the [exterior
  derivative](#exteriorderivative), a 1-form, turned into an array by
  [`Components`](#components). `i` is an [index](#index), and
  [`IndexSum`](#componenttensor-and-indexsum) contracts.
- **Integration.** [`integral(f, D)`](#integral-over-a-decomposed-domain) of
  a 0-form is the sum over the cell index `c` of [`Integral`](#integral)s over
  the reference triangle `K`, each of the [pull-back](#pullback) to `K` times
  the density $\lvert \det J \rvert$ of the cell.
- **Forms.** [`form`](#form) binds `v` in slot 0 and `u` in slot 1. `f` stays
  fixed. `a` and `L` are expressions with the [argument
  spaces](#argument-spaces) `(V, V)` and `(V,)`.

Nothing is evaluated during construction. The [core
transformations](#core-transformations) then rewrite the graph:
`apply_derivatives` and `expand_pull_backs` move the derivatives and the
pull-backs to the leaves, and turn the components on `D` into the inverse
Jacobian times the components on `K`. `pull_back(u, cell_domain(D, 0))`
remains, because the library of `V` gives `u` on each cell. A form compiler
then splits the graph at the cell index (see [Concept: kernels](#concept-kernels)).
The `IndexSum` over `c` becomes the loop of the assembler, the gathers of the
vertex positions become the packing of the kernel inputs, and the rest is the
kernel.

On a mesh of triangles and quadrilaterals, only the domain changes. The
topology has two cell types, and the domain one piece for each:

```python
Kt, Kq = ReferenceSimplex(2), ReferenceQuadrilateral()
M = ...  # a topology with the cell types (Kt, Kq)

# The affine maps of the triangles, with the free index cell_index(M, 0), and
# the bilinear maps of the quadrilaterals, with cell_index(M, 1).
x_tri, x_quad = ..., ...

D = decomposed_domain(
    (parameterised_domain(Kt, x_tri, RD(2)), parameterised_domain(Kq, x_quad, RD(2)))
)
```

`V`, `u`, `v`, `f`, the integrands and the forms are unchanged. `u` is one
function on `D`, and the library of `V` gives `u` on each cell of each type.
Each integral over `D` is now a sum of two terms, one `IndexSum` over the
triangles on `Kt` and one over the quadrilaterals on `Kq`, and a form compiler
makes a kernel for each. `ReferenceQuadrilateral` stands for a reference cell of a
library.

Point values use the same `integral`. A set of points in the mesh is a
decomposed domain of dimension 0, one family of points with an ordinary index:

```python
# m sensor points in D, with the coordinates q, the measured values d and the
# loads g, runtime arrays of the language extension of the cut cell example.
m = runtime_extent("m")
q = runtime_array("q", (m, 2))
d = runtime_array("d", (m,))
g = runtime_array("g", (m,))
j = index("j")

# The position of point j in RD(2), with no support and the free index j.
q_j = component_tensor(q[j, k], (k,))

# The family of points, each a copy of RD(0) placed in D at q_j.
points = parameterised_domain(RD(0), q_j, D)

# The points as a domain, one piece with the ordinary index j.
S = decomposed_domain((points,))

# A load g[j] at each point, the legacy UFL PointSource, as a linear form.
L_point = form(integral(g[j] * v, S), (v,))

# The misfit to the measurements, a functional for data assimilation.
misfit = integral((u - d[j]) ** Integer(2), S)
```

`runtime_extent` and `runtime_array` are not in the core (see [Example:
integrating over a cut cell](#example-integrating-over-a-cut-cell)). The
integral of a 0-form over a point is the value at the point, so `misfit`
is $\sum_j (u(q_j) - d_j)^2$. The index `j` of `d[j]` and `g[j]` is the index
of the family, so each datum pairs with its point. `S` has no topology, because
`j` is not a cell index, and the library of `V` locates the cell of each point
to evaluate `u` and `v` there.

## Rules of the language

### Basic implementation principles

- Explicit is better than implicit.

- Boilerplate and error-prone patterns should be minimised.

- UFLx should be deterministic on the same machine using the same byte-code.
  Serialisation/pickling should be considered early on.

- Performance is *important* and so any data structures and algorithms in the
  core should be carefully chosen and tested for scalability.

- The minimal language should only need a relatively small number of object
  kinds. When designing the language, consider whether objects can be composed,
  rather than extended. Consider what could be layered on top of the core.

### The base object `Base`

The base object in UFLx must be extended by all objects that are included in
the graph. The graph is a directed acyclic graph (DAG). An edge means "is
defined from".

Everything needed to define the variational problem must be recoverable from
the graph. No information is held outside of the graph.

Objects must satisfy the following rules:

1. A construction argument must be one of a small curated set of types: a base
   object, or exactly an `int`, a finite `float`, a `str`, or a `tuple` whose
   items are of these types. A subclass of these Python types is not accepted,
   so a `bool` is not an `int` and a named tuple is not a `tuple`. Construction
   arguments are therefore deeply immutable, have a hash that is deterministic
   across processes on the same machine, and are compared by type and value.
   `1`, `1.0` and `True` are distinct, `0.0` and `-0.0` are distinct, and every
   value is equal to itself. Construction raises `TypeError` for an argument of
   any other type, and `ValueError` for a `float` that is not finite.
2. Construction must depend only on the construction arguments, not on global
   state.
3. The state of an object is immutable after construction. Every result that
   the object produces after construction depends only on the class and the
   construction arguments.
4. An object must be reconstructible from its class and construction arguments,
   and the reconstructed object must be equal to the original. The base object
   records the construction arguments, and the property
   `construction_arguments` returns the construction arguments, in order, as a
   `tuple`. The constructor of a subclass receives the construction arguments,
   validates the construction arguments, and may store data derived from the
   construction arguments.
5. The edges of an object must be exactly the base objects among its
   construction arguments, including the base objects nested in `tuple`s, in
   order.
6. Any order that the core or a pass relies on must be deterministic across
   processes. A pass, and any other walk of the graph such as hashing,
   serialisation or printing, follows the edges of each object in order (rule
   5). Where unrelated objects must be sorted, the sort uses a total order
   defined by the class and the construction arguments alone.
7. Equality (`==`) must be structural and return a `bool`. Whether `is` is
   structural is implementation dependent. An implementation in which `is` is
   structural must intern objects, so that equal objects are the same object.
   Equality does not decide mathematical equivalence, so $a + b$ and $b + a$
   are not equal. Mathematical equivalencies are applied in later graph
   algorithms.
8. All validation of an object happens in the constructor of the class of the
   object. This rule is a hard rule. A free function never validates and never
   converts arguments. A free function only reads derived properties of
   arguments and constructs objects. An invalid argument is rejected by the
   constructor of the first object that the argument makes invalid. Unless the
   section of a class states otherwise, the constructor raises `TypeError` for
   an argument of a type or a class that the constructor does not accept, and
   `ValueError` for an argument of an accepted type or class with a value that
   the constructor does not accept.
9. A constructor has no default values, so every construction argument is
   written at construction. A free function may supply a default value.
10. A concrete class states an `identifier`, a `str` with the public name of
    the class qualified by the package name, such as `"uflx.Sum"` for the core
    or `"library.UnitSimplex"` for a library. The hash of rule 1, the total
    order of rule 6 and serialisation use the identifier and the construction
    arguments, and never the location of the class in the code. Construction
    raises `TypeError` for a class with no identifier. The package name keeps
    the identifiers of different libraries apart, and a library keeps the
    identifiers of the classes of the library distinct.

> Implementation note: Python randomises the hash of a `str` for each process,
> so rule 1 needs a custom hash for `str` in an implementation in Python.

Hence two objects are equal if and only if the two objects have the same class
and equal construction arguments.

Because objects are immutable, an object need only be validated once, at
construction. Every object that exists is therefore valid.

An object cannot be changed. To obtain a modified object, construct a new
object and reconstruct every object that depends on the original object, up to
the root of the graph.

The above design implies that construction arguments should be small and cheap
to hash and compare, and the graph implementation must be performant.

### Leaves and subclassing

Requires: [The base object `Base`](#the-base-object-base).

See also: [RD](#rd), [ParameterisedDomain](#parameteriseddomain), [Language
extensions](#language-extensions).

Only leaves are subclassed. A leaf is an object that the core cannot build
from other objects, such as the coordinate domain `RD(n)` or a function whose
values an external library supplies. An object built from other objects is a
single class of the core or a free function that returns objects of the core,
never a subclass. For example, every domain is a `ParameterisedDomain`, and an
affinely mapped triangle is a free function that returns a
`ParameterisedDomain`. A language extension is the one exception (see
[Language extensions](#language-extensions)).

A concrete class of the core is final. A library subclasses only an abstract
class of the core, such as `AbstractCoordinateDomain`, so a pass that
dispatches on a concrete class of the core never meets a subclass of the
class.

A parent class never constructs a subclass. A property of a parent class may
return an object whose class is a subclass, because the property only
returns an object that already exists. A subclass may narrow the return
type of a property.

### Properties and free functions

Requires: [The base object `Base`](#the-base-object-base).

See also: [`AbstractExpression`](#abstractexpression),
[ParameterisedDomain](#parameteriseddomain), [Coordinate](#coordinate),
[Transpose and matrix product](#transpose-and-matrix-product).

A property reads and a free function builds. Data that an object stores, or
derives from the class and the construction arguments of the object, is a
property of the object, such as `value_shape` of an expression or
`parameter_domain` of a domain. A property never constructs an object. Any
operation that constructs objects is a pure free function, such as
`coordinate(P)` or `transpose(A)`. A class method, where
needed, immediately calls a free function.

Free functions are not objects in the graph. A free function builds or returns
objects in the graph. Operator overloads call the free functions.

A constructor takes the construction arguments by position only, in order, and
raises `TypeError` for an argument by name. So an object has one `tuple` of
construction arguments. A free function may take arguments by name, such as
`parameterised_domain(P=T, psi=psi, Q=RD(2))`.

A stored property returns a construction argument. A derived property returns
data derived from the class and the construction arguments. A concrete class
has the derived properties of the abstract classes of the concrete class.

### No simplification during construction

Requires: [The base object `Base`](#the-base-object-base).

Simplifications and rewriting should be completely avoided during user graph
construction, maintaining the user's intent. Only later selectable algorithmic
passes should re-write. A good example is `(a + b) + c` in floating point
arithmetic where `a`, `b` and `c` are unknown terminals (i.e. to be provided
later). The user may have very good reason to have `(a + b)` performed before
adding the answer to `c`.

### Passes

Requires: [The base object `Base`](#the-base-object-base).

See also: [Index](#index), [ComponentTensor and
IndexSum](#componenttensor-and-indexsum).

A pass is a pure free function that walks a graph after construction of the
graph. The result of a pass depends only on the arguments of the pass. A pass
returns a value, such as a unit, or a new graph. A pass builds a new graph by
constructing an object for each node that the pass replaces, and by
reconstructing every object that depends on a replaced node, up to the root
(see [The base object](#the-base-object-base)). A pass may cache results for
each object, because objects are immutable.

Substitution avoids capture. A `ComponentTensor` or an `IndexSum` binds
indices (see [ComponentTensor and
IndexSum](#componenttensor-and-indexsum)), and a [`Form`](#form) binds
functions. If a pass substitutes an expression
into the expression of a binding node, and a bound index of the binding node is
free in the substituted expression, the pass renames the bound index. The new
name is the derived index `Index(("uflx.rename", node, m))`, or the derived
`EllipsisIndex` with the same label, or the function
`Function(V, ("uflx.rename", node, m))` in the space $V$ of a bound function,
with `node` the binding node before the substitution and `m` the position of
the bound index or function. For example, replacing
`c` by `b[i]` in `index_sum(a[i] * c, (i,))` gives
`index_sum(a[j] * b[i], (j,))`, with `j` the renamed index, and the result has
the free index `i`. Without the renaming, the result would be
`index_sum(a[i] * b[i], (i,))`, with no free index.

A pass does one thing, such as differentiating or lowering. Passes compose,
and the user of the core chooses which passes run and in what order.

### Kinds of code that produce objects

Requires: [The base object `Base`](#the-base-object-base), [Properties and free
functions](#properties-and-free-functions), [Passes](#passes).

See also: [Transpose and matrix product](#transpose-and-matrix-product),
[Indexed](#indexed), [New questions](#new-questions).

Three kinds of code produce objects:

| Kind                                               | When the code runs              | Inspects arguments to choose what to build                                | Validates                     |
| -------------------------------------------------- | ------------------------------- | ------------------------------------------------------------------------- | ----------------------------- |
| Constructor                                        | Whenever an object is made      | No, a constructor inspects arguments only to validate                     | Yes, and only constructors do |
| Construction free function, such as `transpose`    | While the user builds the graph | Only the static type, to build one fixed composition for each static type | No                            |
| Pass, such as differentiation, lowering or `units` | After construction              | Yes, and a pass may rewrite the graph                                     | No                            |

A construction free function fails when a constructor rejects the
arguments, for example `transpose` of a vector raises `IndexError` from the
constructor of `Indexed`. A pass fails with `TypeError` for a node class that
the pass has no rule for, or with an error of the pass, such as inconsistent
units.

### Enforcement

Requires: [The base object `Base`](#the-base-object-base).

Enforcement must be implemented to help library developers - look at tools such
as `dataclass`, deterministic hashing, hash-consing, `WeakValueDictionary`,
`ruff` and `semgrep`.

https://github.com/semgrep/semgrep

This rules should be put in place at the start of the project.

It is not clear at the moment how best to do this - it may be that custom
static analysis e.g. (semgrep) is more robust than hacking Python to be 'well
behaved' in some places, for example checking the purity of functions is not
possible in Python.

## Coordinate domains

We begin with differential geometry because the function spaces and operators
of UFLx are defined on a domain. The function spaces and operators require the
topological structure of the domain, which is the dimension, boundary,
orientation and decomposition into cells of the domain. The core needs no
metric beyond the Euclidean metric of the ambient coordinate domain, which
gives the components of a differential form on a domain of lower dimension
(see [Components](#components)) and the unit normal (see [Normal](#normal)).

### AbstractCoordinateDomain

Requires: [The base object `Base`](#the-base-object-base).

See also: [Support](#support), [Coordinate](#coordinate),
[ParameterisedDomain](#parameteriseddomain).

An abstract base class for coordinate domains `AbstractCoordinateDomain`.

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$ whose
points are named by the coordinates of the points themselves. A coordinate
domain is a parameter space. A coordinate domain is the parameter domain of a
domain, the ambient coordinate domain of a domain, and the support of an
expression. A coordinate domain is not a
`ParameterisedDomain`. The domain whose points are the points of a coordinate domain
$P$ is `parameterised_domain(P, coordinate(P), P)` (see
[ParameterisedDomain](#parameteriseddomain)).

A new coordinate domain, such as a reference cell, is a subclass of
`AbstractCoordinateDomain`, because the core cannot build the set from smaller
objects.

| `AbstractCoordinateDomain` |                                              |
| -------------------------- | -------------------------------------------- |
| Provides                   | `dimension`, an `int` $n \ge 0$              |
| Defined by                 | A subclass, such as `RD` or a reference cell |

Examples:

- The unit simplex in $R^n$ and $R^n$ itself are coordinate domains.
- A face of codimension 1 of the unit simplex in $R^n$ is *not* a coordinate
  domain, because the face has empty interior in $R^n$. The face is the image
  of the unit simplex in $R^{n-1}$ under a map, so the face is a parameterised
  domain (defined below).

### RD

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain).

See also: [Example: a paraboloid](#example-a-paraboloid).

A concrete subclass of `AbstractCoordinateDomain` for the whole of $R^n$.

| `RD(n)`    |                                                                 |
| ---------- | --------------------------------------------------------------- |
| Properties | `dimension`, an `int` $n \ge 0$                                 |
| Raises     | `TypeError` if $n$ is not an `int`, and `ValueError` if $n < 0$ |

Example:

- `RD(2)` is the parameter domain of the paraboloid example below.

### Example: the unit simplex

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain).

`UnitSimplex` is not part of the core. The examples use `UnitSimplex` as a
concrete `AbstractCoordinateDomain`.

Define the unit simplex in $R^n$,

$\{p \in R^n : p_i \ge 0 \text{ for every } i, \ \sum_i p_i \le 1\}$,

as a concrete `AbstractCoordinateDomain`.

```python
class UnitSimplex(AbstractCoordinateDomain):
    identifier = "library.UnitSimplex"

    def __init__(self, n: int):
        if type(n) is not int:
            raise TypeError("Dimension must be an int.")
        if n < 0:
            raise ValueError("Dimension must be non-negative.")
        self._n = n

    @property
    def dimension(self) -> int:
        return self._n
```

`UnitSimplex(1)` is the unit interval $[0, 1]$. `UnitSimplex(2)` is the unit
triangle with vertices $(0, 0)$, $(1, 0)$ and $(0, 1)$.

## Expressions and static type

### `AbstractExpression`

Requires: [The base object `Base`](#the-base-object-base).

Every `AbstractExpression` derives the value shape, the free indices, the
scalar type, the support, the degree and the [argument
spaces](#argument-spaces) of the expression. The value shape is
a `tuple` of [extents](#extents), one for each axis. An expression with value
shape `()` is a scalar.

### Index

Requires: [The base object `Base`](#the-base-object-base), [Properties and free
functions](#properties-and-free-functions).

See also: [Indexed](#indexed), [Transpose and matrix
product](#transpose-and-matrix-product).

> LLM generated content.

Index notation writes tensor algebra as fixed compositions, so that operators
such as `transpose` and `matrix_product` are free functions rather than node
classes.

A concrete subclass of the base object for a free index. An `Index` is not an
expression.

| `Index(label)` |                                                                                                 |
| -------------- | ----------------------------------------------------------------------------------------------- |
| Free function  | `index(label)`                                                                                  |
| Properties     | `label`, a `str` for a user index, or a `tuple` that contains graph objects for a derived index |
| Raises         | `TypeError` if the label is neither a `str` nor a `tuple`                                       |

An `Index` has no extent of its own, and takes the extent of the axis that the
`Index` indexes.

The user names user indices. Two occurrences of one `Index` are the same
index, so `a[i] * b[i]` is the component-wise product (see [Free
indices](#free-indices)), and `A[i, i]` is the diagonal entry (see
[Indexed](#indexed)). An index that the user means to be new needs a new name.

A free function that needs fresh indices derives the indices from the
operands. The label of a derived index is the identifier of the free function,
followed by the arguments of the free function and a number. The identifier is
the public name of the free function, qualified by the package name, such as
`"uflx.transpose"` for the core or `"library.trace"` for a library. For example,
`transpose(A)` uses `Index(("uflx.transpose", A, 0))` and
`Index(("uflx.transpose", A, 1))`. An index whose label contains `A` has an
edge to `A`. Every node that uses an `Index` has an edge to the `Index`. So an index
whose label contains `A` cannot occur inside `A`, because `A` would then
contain `A`, and the graph is acyclic. Hence
a derived index is fresh with respect to the operands, and the free function
does not inspect the operands to find a fresh index. A global counter is not
allowed, because a counter makes construction depend on global state.

A pass that needs fresh indices derives the indices in the same way, with the
identifier of the pass, the arguments of the rule and a number. For example,
`expand_components` writes the new axis of `Components(g, P)` with the index
`Index(("uflx.expand_components", g, P, 0))`. The renamed index of
[Passes](#passes) follows the same rule.

The edges of an index to the objects in the label of the index identify the
index and are not operands. The label keeps the operands as the operands were
when the index was derived, and a later pass may rewrite the operands. So a
pass that computes or collects from the value of an expression does not walk
into the label of an index.

### EllipsisIndex

Requires: [Index](#index).

See also: [Indexed](#indexed), [Tensor algebra with
EllipsisIndex](#tensor-algebra-with-ellipsisindex).

> LLM generated content.

> Legacy UFL: This numpy-like concept does not appear to be in legacy UFL.

A concrete subclass of the base object for an index that covers any number of
axes. A `EllipsisIndex` is not an expression.

| `EllipsisIndex(label)` |                                                           |
| ---------------------- | --------------------------------------------------------- |
| Free function          | `ellipsis_index(label)`                                   |
| Properties             | `label`, with the rules of the label of an `Index`        |
| Raises                 | `TypeError` if the label is neither a `str` nor a `tuple` |

A `EllipsisIndex` covers every axis of an `Indexed` that the other entries do not
cover, so a `EllipsisIndex` may cover no axis. The extent of a `EllipsisIndex` is the
`tuple` of the extents of the covered axes. In every other rule of index
notation, a `EllipsisIndex` behaves as an `Index`.

### Extents

Requires: [The base object `Base`](#the-base-object-base).

See also: [Free indices](#free-indices), [Indexed](#indexed),
[Example: quadrature](#example-quadrature), [Example: integrating over a cut
cell](#example-integrating-over-a-cut-cell).

An extent is the length of an axis of a value shape, or the extent of a free
index. An extent is an `int` $n \ge 0$ or an `AbstractExtent`.

An abstract subclass of the base object `AbstractExtent` for an extent whose
value is not known during construction. An `AbstractExtent` is not an
expression. A concrete `AbstractExtent` is a subclass outside the core, such as
`RuntimeExtent` in the
[cut cell example](#example-integrating-over-a-cut-cell).

An `AbstractExtent` has no arithmetic and no ordering, so the rules of the
core, which compare extents by equality, accept an `AbstractExtent` unchanged.
Two `AbstractExtent`s are the same extent only if the two `AbstractExtent`s are
equal as base objects, so a concrete `AbstractExtent` stores the arguments that
identify the extent, such as a label.

Example:

- With an `AbstractExtent` `q`, an expression `points` of value shape `(q, 2)`
  and an expression `weights` of value shape `(q,)` share the extent `q`, so
  `weights[i] * points[i, 0]` is accepted. With a second `AbstractExtent` `r`,
  an expression `weights` of value shape `(r,)` raises `ValueError` in
  `weights[i] * points[i, 0]`, because the index `i` has two extents.

> Open: ragged indexing. An axis whose length depends on another index, such
> as the quadrature points of each cut cell, cannot be written. A dependent
> extent, an `AbstractExtent` whose label contains the other index, would cover
> it. An index value read from a table is a [GatheredIndex](#gatheredindex).

### AbstractIndexValue

Requires: [The base object `Base`](#the-base-object-base), [Extents](#extents).

See also: [Indexed](#indexed), [Sides of a topology](#sides-of-a-topology).

An abstract subclass of the base object `AbstractIndexValue` for an index whose
value is fixed but not known during construction, such as the local number of
the facet of a cell in a facet integral, which the solver gives for each facet.
An `AbstractIndexValue` is not an expression and not a free index.

| `AbstractIndexValue` |                                                        |
| -------------------- | ------------------------------------------------------ |
| Provides             | `extent`, an extent, and `shape`, a `tuple` of extents |
| Defined by           | A subclass outside the core                            |

An `AbstractIndexValue` with the shape `()` is one value, an `int` $v$ with
$0 \le v <$ `extent`. An `AbstractIndexValue` with a non-empty shape is a table
of such values, such as the cell on one side of each facet, which the mesh
gives, and is read through a [GatheredIndex](#gatheredindex).

An `AbstractIndexValue` with the shape `()` is an entry of [Indexed](#indexed).
The entry covers one axis, as an `int` entry does, and the extent of the entry
must equal the extent of the axis.

> Legacy note: GEM, the intermediate language of TSFC, has `VariableIndex` for
> the same purpose -- not fact checked.

Example:

- With an `AbstractIndexValue` `l` with the extent 3, and the table `starts` of
  value shape `(3, 2)` of the
  [boundary example](#example-integrating-over-the-boundary-of-a-triangle),
  `starts[l, k]` is the start of edge `l`, with the free index `k`. With an
  `AbstractIndexValue` with the extent 4, the same expression raises
  `IndexError`.

### GatheredIndex

Requires: [AbstractIndexValue](#abstractindexvalue), [Index](#index),
[Extents](#extents).

See also: [Indexed](#indexed), [Restriction](#restriction).

A concrete subclass of the base object for an index value read from a table.
A `GatheredIndex` is not an expression and not a free index.

| `GatheredIndex(table, indices)` |                                                                                                                                                                                   |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Free function                   | `gathered_index(table, indices)`                                                                                                                                                  |
| Properties                      | `table`, an `AbstractIndexValue`, and `indices`, an `int` or an `Index` for each axis of the shape of the table                                                                   |
| Static type                     | The extent of the table, and the free indices of the `Index` entries, each with the extent of the axis                                                                            |
| Raises                          | `TypeError` for an entry that is neither exactly `int` nor `Index`, and `IndexError` if the number of entries differs from the number of axes, or for an `int` entry out of range |

The value is the entry of the table at the indices.

A `GatheredIndex` is an entry of [Indexed](#indexed), and covers one axis
whose extent equals the extent of the `GatheredIndex`. The free indices of the
`GatheredIndex` are free in the `Indexed`.

Example, with a table `cells` of the shape `(4,)`, the extent 2 and the values
`(0, 1, 1, 0)`, such as the cell of each of the four boundary edges of the
square divided into two triangles, and an array `a` of value shape `(2,)` with
one value for each triangle:

```python
e = index("e")
g = gathered_index(cells, (e,))

# The value of a on the triangle of each edge, with the free index e of extent 4.
a_e = a[g]

# The array (a[0], a[1], a[1], a[0]) of value shape (4,).
b = component_tensor(a_e, (e,))

# The sum over the edges, 2 * (a[0] + a[1]) after evaluation.
total = index_sum(a_e, (e,))
```

- `gathered_index(cells, (2,))` has no free index and the value 1, so
  `a[gathered_index(cells, (2,))]` is `a[1]`.
- For an array `c` of value shape `(3,)`, `c[g]` raises `IndexError`, because
  the extent 2 of `g` differs from the extent 3 of the axis.

### Free indices

Requires: [`AbstractExpression`](#abstractexpression), [Index](#index),
[EllipsisIndex](#ellipsisindex), [Extents](#extents).

See also: [Sum](#sum), [Tensor](#tensor), [Indexed](#indexed), [ComponentTensor
and IndexSum](#componenttensor-and-indexsum).

Every `AbstractExpression` derives the free indices, a `tuple` of pairs of an
`Index` and an extent, or of a `EllipsisIndex` and a `tuple` of extents. The
order is the order of first occurrence, never the order of a `set`.

Unless a node class states otherwise, the free indices of a node are the free
indices of the edges of the node that are expressions, in edge order, with each
`Index` taken once at the first occurrence. An `Index` that occurs in two edges
must have the same extent in both edges, or construction raises `ValueError`.

The following node classes state otherwise:

- `Indexed`: see [Indexed](#indexed).
- `ComponentTensor` and `IndexSum`: see
  [ComponentTensor and IndexSum](#componenttensor-and-indexsum).
- `PullBack`: see [PullBack](#pullback).

Summation is never implicit. `a[i] * b[i]` has the free index `i` and denotes
the component-wise product. The sum over `i` is `index_sum(a[i] * b[i], (i,))`.
The operands of a node may have different free indices, so `a[i] + b[j]` has
the free indices `i` and `j` and denotes the outer sum.

### Scalar type

Requires: [`AbstractExpression`](#abstractexpression).

See also: [Literals](#literals), [Absolute](#absolute).

The scalar type is real or complex. Unless the [static type
table](#appendix-static-type-table) states otherwise, an expression is complex
if and only if at least one edge is a complex expression.

The scalar type is the arithmetic in which the expression is evaluated, not a
statement about the value. A real `SquareRoot` or a real `Power` whose value is
not a real number, such as `square_root(Integer(-1))`, is undefined.

### Support

Requires: [`AbstractExpression`](#abstractexpression),
[AbstractCoordinateDomain](#abstractcoordinatedomain).

The support is an `AbstractCoordinateDomain` or a
[`DecomposedDomain`](#decomposeddomain). The support is absent for an expression
with no support, such as a literal. Unless the [static type
table](#appendix-static-type-table) states otherwise, the support
of a node is the support of the edges of the node that are expressions with a support.
All such edges must have the same support, or construction raises `ValueError`. A
node with no such edge has no support.

### Degree

Requires: [`AbstractExpression`](#abstractexpression).

See also: [Appendix: Static type table](#appendix-static-type-table), [New
subclasses](#new-subclasses), [ExteriorDerivative](#exteriorderivative).

The degree is a `tuple` with the form degree of each slot of the expression,
and every entry is at least 1. An expression with the degree `()` has no
slots and is a 0-form, such as a scalar or an array. The axes of the value
shape are not slots. The degree has no default. Every concrete class states
the degree of the class, and a leaf class of a downstream library states the
degree as part of the contract of the class (see
[New subclasses](#new-subclasses)).

The degree of every node class of the core is in the
[static type table](#appendix-static-type-table). If the
degree of an edge breaks the rule of the node, construction raises
`ValueError`. A rule that requires the degree `()` of an edge is a deliberate
statement, such as for `SquareRoot`, because the square root of a differential
form is not defined.

Examples of the degree, on $R^3$ with the coordinate function $X$, the
coordinates $x, y, z$ and a scalar function $f$:

| Object                                                                     | Value shape | Degree   |
| -------------------------------------------------------------------------- | ----------- | -------- |
| The scalar $f$, or the literal `Integer(1)`                                | `()`        | `()`     |
| The coordinate function `X`                                                | `(3,)`      | `()`     |
| The 1-form $\mathrm{d}f$                                                   | `()`        | `(1,)`   |
| The 1-form $\mathrm{d}x$                                                   | `()`        | `(1,)`   |
| The 2-form $\mathrm{d}x \wedge \mathrm{d}y$                                | `()`        | `(2,)`   |
| The volume form $\mathrm{d}x \wedge \mathrm{d}y \wedge \mathrm{d}z$        | `()`        | `(3,)`   |
| The vector-valued 1-form $\mathrm{d}X$, with one 1-form for each component | `(3,)`      | `(1,)`   |
| The double form $\mathrm{d}x \otimes \mathrm{d}y$, a covariant 2-tensor    | `()`        | `(1, 1)` |

### Expressions as differential forms

Requires: [Support](#support), [Degree](#degree).

See also: [Zero](#zero).

An expression denotes an array-valued differential form, a multi-form, on the
support of the expression, with the degree of the expression. A 0-form is an
array-valued function. An expression with no support and the degree `()` denotes a
single array. An expression with argument spaces other than `()` denotes a
multilinear map from the argument spaces to such differential forms (see
[Variational forms](#variational-forms)). A leaf with a non-empty degree has a
support, except `Zero`. So an
expression with no support and a non-empty degree is built from `Zero`s, such as
`Zero((), (), (1,), ()) + Zero((), (), (1,), ())`, and denotes the zero form of that
degree on every support.

### Node class or free function

Requires: [Properties and free functions](#properties-and-free-functions),
[`AbstractExpression`](#abstractexpression), [Free indices](#free-indices),
[Scalar type](#scalar-type), [Support](#support), [Degree](#degree),
[AbstractCoordinateDomain](#abstractcoordinatedomain).

See also: [Transpose and matrix product](#transpose-and-matrix-product), [Volume
form](#volume-form), [ExteriorDerivative](#exteriorderivative),
[Division](#division), [Elementary functions (very
incomplete)](#elementary-functions-very-incomplete),
[Tensor algebra with
EllipsisIndex](#tensor-algebra-with-ellipsisindex),
[Language extensions](#language-extensions).

A concept must be a free function, not a node class, if the free function can
build the concept as one fixed composition of existing node classes for each
static type of the arguments. The static type of an expression is the value
shape, the free indices, the scalar type, the support, the degree and the
argument spaces, and the
static type of a domain is the stored properties of the domain, such as the
dimension of a coordinate domain or the pieces of a decomposed domain. A free function may read
the static type to choose the composition, and never reads the structure of
an argument, such as the class, the edges or the values. Otherwise the concept
gets its own node class. A concept that is evaluated by a numerical algorithm
that the user supplies also gets its own node class, even if a fixed
composition exists, because the composition would fix the algorithm. For
example, the transpose is always `component_tensor(A[k, l], (l, k))`, and
`volume_form(P)` is the product of the $n$ coordinate 1-forms of a coordinate
domain of dimension $n$, so both are free functions. No composition of
existing nodes differentiates, so `ExteriorDerivative` is a node class.
`Division` and `SquareRoot` are evaluated by a numerical algorithm that the
user supplies, so each is a node class. A free
function that builds a fixed composition is a definition, not a
simplification. The free function does not read the structure of the
operands, so the free function never removes structure that the user built.
A `EllipsisIndex` makes most tensor algebra one fixed composition for every
value shape (see [EllipsisIndex](#tensor-algebra-with-ellipsisindex)), so a free function rarely
needs to read a value shape.

A node class is part of the core only if a definition or a rule of the core
needs the node class, or if the node class is an operation of the scalar
field, such as conjugation. A rule for the node class itself does not count.
Otherwise the node class belongs to a
[language extension](#language-extensions). For example, the core writes a
rational number as the `Division` of two `Integer`s, so `Division` is part of
the core.

## Algebra

### Literals

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices), [Scalar type](#scalar-type), [Support](#support),
[Degree](#degree).

See also: [Division](#division).

A literal is a number written down in an expression, e.g. `Integer(2) * u`. A
literal's value is part of the expression itself. Every literal has value shape
`()`, no free indices, no support and the degree `()`.

The literals are concrete subclasses of `AbstractExpression`:

```
Integer(2)
Float(3.5)
ImaginaryUnit()
```

Each literal constructor accepts only exact types. Construction raises
`TypeError` for a value of another type, and `ValueError` for a `float` that is
not finite, as in rule 1 of [the base object](#the-base-object-base):

- `Integer(v)` takes an `int` $v$. A `bool` is rejected, although `bool` is a
  subclass of `int`.
- `Float(v)` takes a finite `float` $v$. `Float(-0.0) != Float(0.0)`.
- `ImaginaryUnit()` takes no arguments and denotes $i$.

A rational number is the `Division` of two `Integer`s, such as
`Integer(1) / Integer(3)`. In floating-point arithmetic, the division of two
integers rounds once, the same as the conversion of a rational number, so the
core has no rational literal. A complex number is a composition such as
`a + ImaginaryUnit() * b`, so the core has no complex literal.

### Zero

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices), [Support](#support), [Degree](#degree).

See also: [Indexed](#indexed), [Applying derivatives](#applying-derivatives).

| `Zero(value_shape, free_indices, degree, argument_spaces)` |                                                                                                                                                                                    |
| ---------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Properties                                                 | `value_shape`, `free_indices`, `degree` and `argument_spaces`, as in [Extents](#extents), [Free indices](#free-indices), [Degree](#degree) and [Argument spaces](#argument-spaces) |
| Static type                                                | The stored properties, with no support                                                                                                                                             |
| Raises                                                     | `TypeError` for an argument of another form, and `ValueError` for a negative extent, an index that occurs twice, or a slot degree less than 1                                      |

`Zero` is the zero of any static type. `apply_derivatives` produces `Zero`
(see [Applying derivatives](#applying-derivatives)). No fixed composition gives
a zero with free indices or with a non-empty degree, because
`Integer(0) * a[i]` has edges to `a`, and $0 \, \mathrm{d}x$ has edges to
$\mathrm{d}x$. So `Zero` is a node class. For example,
`Zero((), (), (1,), ())` is the zero 1-form, `Zero((3,), (), (), ())` is the
zero array of value shape `(3,)`, and `Zero((), (), (), (V,))` is the zero
linear form on a space `V`.

### Arithmetic

Requires: [`AbstractExpression`](#abstractexpression), [Free
indices](#free-indices).

See also: [ComponentTensor and IndexSum](#componenttensor-and-indexsum), [Tensor
algebra with EllipsisIndex](#tensor-algebra-with-ellipsisindex).

`Sum`, `Product`, `Division` and `Power` are concrete subclasses of
`AbstractExpression` that each store two expressions. Arithmetic is scalar.
Both operands must have value shape `()`, or construction raises
`ValueError`, and the result has value shape `()`. Free indices pass through
arithmetic (see [Free indices](#free-indices)), so tensor arithmetic is
written in index notation. For example, a scalar `a` times a tensor `T` of any
rank is `component_tensor(a * T[r], (r,))`, with `r` a
[EllipsisIndex](#ellipsisindex).

| Class      | Free function | Overload |
| ---------- | ------------- | -------- |
| `Sum`      | `add`         | `+`      |
| `Product`  | `multiply`    | `*`      |
| `Division` | `divide`      | `/`      |
| `Power`    | `power`       | `**`     |

### Sum

Requires: [Arithmetic](#arithmetic), [Degree](#degree).

See also: [ExteriorDerivative](#exteriorderivative).

| `Sum(first, second)` |                                                                 |
| -------------------- | --------------------------------------------------------------- |
| Free function        | `add(first, second)`, with the overload `+`                     |
| Requires             | The operands of [Arithmetic](#arithmetic), with the same degree |
| Raises               | `ValueError` if a requirement fails                             |

`Sum` adds two differential forms of the same degree. For example,
$\mathrm{d}x + \mathrm{d}y$ has the degree `(1,)`, and $f + \mathrm{d}x$
raises `ValueError`.

### Product

Requires: [Arithmetic](#arithmetic), [Free indices](#free-indices),
[Support](#support), [Degree](#degree),
[AbstractCoordinateDomain](#abstractcoordinatedomain).

See also: [TensorProduct](#tensorproduct), [ComponentTensor and
IndexSum](#componenttensor-and-indexsum),
[ExteriorDerivative](#exteriorderivative).

| `Product(first, second)` |                                                                                                                |
| ------------------------ | -------------------------------------------------------------------------------------------------------------- |
| Free function            | `multiply(first, second)`, with the overload `*`                                                               |
| Requires                 | The operands of [Arithmetic](#arithmetic). If both operands have a non-empty degree, each has exactly one slot |
| Static type              | The degrees added slot by slot, where the degree `()` adds 0 to every slot                                     |
| Raises                   | `ValueError` if a requirement fails                                                                            |

`Product` is the wedge product $\alpha \wedge \beta$. For two 0-forms, the
wedge product is the ordinary product, and $f \wedge \omega = f \omega$. With
two slots in an operand, the interior product of [Expanding
components](#expanding-components) removes a slot from one operand only, and
the slots of the operands no longer align. A slot of the result may have a
degree greater than the dimension of the support, and then the result is the
zero form of that degree.

The wedge product is graded commutative,
$\alpha \wedge \beta = (-1)^{kl} \beta \wedge \alpha$ for a $k$-form $\alpha$
and an $l$-form $\beta$. Equality is structural, so `Product` never reorders
the operands, and a later pass may apply the identity. `Product` never adds a
slot, so the product of two 1-forms is a 2-form with the degree `(2,)`, not a
double form with the degree `(1, 1)`.

No fixed composition of the other node classes multiplies, so `Product` is a
node class. Every node class combines the free indices of operands that
differ, by the general rule of [Free indices](#free-indices), unless the node
class states otherwise. So a contraction is an `IndexSum` of any node class,
such as `index_sum(a[i] * b[i], (i,))` or `index_sum(a[i] / b[i], (i,))`.

Examples, with a scalar $f$, the 1-forms $\mathrm{d}x$ and $\mathrm{d}y$ and a
double form $S$ with the degree `(1, 1)`:

| Product                                         | Degree              |
| ----------------------------------------------- | ------------------- |
| $f \cdot f$                                     | `()`                |
| $f \, \mathrm{d}x$                              | `(1,)`              |
| $\mathrm{d}x \wedge \mathrm{d}y$                | `(2,)`              |
| $\mathrm{d}x \wedge \mathrm{d}x$, the zero form | `(2,)`              |
| $f \, S$                                        | `(1, 1)`            |
| $\mathrm{d}x \wedge S$                          | raises `ValueError` |

### Power

Requires: [Arithmetic](#arithmetic), [Support](#support).

See also: [Applying derivatives](#applying-derivatives).

| `Power(first, second)` |                                                                             |
| ---------------------- | --------------------------------------------------------------------------- |
| Free function          | `power(first, second)`, with the overload `**`                              |
| Requires               | The operands of [Arithmetic](#arithmetic), with an exponent with no support |
| Raises                 | `ValueError` if a requirement fails                                         |

The derivative of a power whose exponent has a support needs a logarithm,
which the core does not have.

### Negation and subtraction

Requires: [Sum](#sum), [Product](#product), [Literals](#literals).

Negation and subtraction are free functions, not node classes:

- `negate(f)`, with the overload `-f`, is `multiply(Integer(-1), f)`.
- `subtract(f, g)`, with the overload `f - g`, is `add(f, negate(g))`.

Both compositions are exact in floating-point arithmetic, so neither
composition changes the arithmetic that the user wrote.

### Division

Requires: [Arithmetic](#arithmetic), [Product](#product), [Power](#power),
[Literals](#literals).

Division is a node class, because `multiply(f, power(g, Integer(-1)))`
rounds twice and differs from `f / g` in floating-point arithmetic.

| `Division(first, second)` |                                                                                        |
| ------------------------- | -------------------------------------------------------------------------------------- |
| Free function             | `divide(first, second)`, with the overload `/`                                         |
| Requires                  | The operands of [Arithmetic](#arithmetic), with the degree `()` for the second operand |
| Raises                    | `ValueError` if a requirement fails                                                    |

The quotient has the degree of the first operand, so `dx / Integer(2)` is a
1-form.

### Absolute

Requires: [`AbstractExpression`](#abstractexpression), [Scalar
type](#scalar-type).

| `Absolute(f)` |                                           |
| ------------- | ----------------------------------------- |
| Free function | `absolute(f)`, with the overload `abs(f)` |
| Requires      | A scalar $f$                              |
| Static type   | Value shape `()`, real                    |
| Raises        | `ValueError` if a requirement fails       |

For a complex argument, `Absolute` is the modulus.

### Conjugate

Requires: [`AbstractExpression`](#abstractexpression), [Scalar
type](#scalar-type), [Degree](#degree).

See also: [Variational forms](#variational-forms).

| `Conjugate(f)` |                                              |
| -------------- | -------------------------------------------- |
| Free function  | `conjugate(f)`                               |
| Requires       | A scalar differential form $f$ of any degree |
| Raises         | `ValueError` if a requirement fails          |

`Conjugate` is the complex conjugate, acting on the coefficients of the form.
No fixed composition of the other node classes conjugates, so `Conjugate` is a
node class. A sesquilinear variational form, such as the $L^2$ inner product
$\int u \bar{v}$, needs `Conjugate`.

### TensorProduct

Requires: [Arithmetic](#arithmetic), [Product](#product), [Free
indices](#free-indices), [Degree](#degree).

See also: [ExteriorDerivative](#exteriorderivative),
[ParameterisedDomain](#parameteriseddomain).

| `TensorProduct(first, second)` |                                                                                 |
| ------------------------------ | ------------------------------------------------------------------------------- |
| Free function                  | `tensor_product(first, second)`, with no overload                               |
| Requires                       | Two operands with value shape `()`, as in [Arithmetic](#arithmetic)             |
| Static type                    | The degree `first.degree + second.degree`, and the free indices as in `Product` |
| Raises                         | `ValueError` if a requirement fails                                             |

`TensorProduct` is the tensor product $\alpha \otimes \beta$ of two
differential forms. An operand with the degree `()` adds no slot, so
$f \otimes \omega = f \omega$.

The tensor product is not commutative, and `TensorProduct` never reorders the
operands. No fixed composition of the other node classes adds a slot, so
`TensorProduct` is a node class. `Product` adds degrees within a slot, and
`TensorProduct` adds slots.

Examples, with the setup of [ExteriorDerivative](#exteriorderivative) on
`RD(3)` and `i = index("i")`:

| Expression                                        | Degree              | Meaning                                                    |
| ------------------------------------------------- | ------------------- | ---------------------------------------------------------- |
| `tensor_product(dx, dy)`                          | `(1, 1)`            | $\mathrm{d}x \otimes \mathrm{d}y$                          |
| `index_sum(tensor_product(dX[i], dX[i]), (i,))`   | `(1, 1)`            | The Euclidean metric                                       |
| `tensor_product(dx * dy, dz)`                     | `(2, 1)`            | $(\mathrm{d}x \wedge \mathrm{d}y) \otimes \mathrm{d}z$     |
| `tensor_product(dx, dy) * tensor_product(dy, dz)` | raises `ValueError` | Each operand of `Product` with a slot has exactly one slot |
| `tensor_product(dX, dX)`                          | raises `ValueError` | `dX` has value shape `(3,)`                                |

For a parameterised domain $D$ with the parameter domain $P$ and
`dpsi = exterior_derivative(D.parameterisation, P)`,
`index_sum(tensor_product(dpsi[i], dpsi[i]), (i,))` is the metric that $\psi$
induces on $P$, the pull-back of the Euclidean metric.

### Elementary functions (very incomplete)

Requires: [`AbstractExpression`](#abstractexpression), [Power](#power),
[Division](#division), [Literals](#literals).

> UFLx core is likely to contain very few elementary functions - a language
> extension e.g. `uflx.elementary` will implement them later.

`SquareRoot` is the only elementary function so far.

| `SquareRoot(f)` |                                     |
| --------------- | ----------------------------------- |
| Free function   | `square_root(f)`                    |
| Requires        | A scalar $f$ with the degree `()`   |
| Static type     | Value shape `()`                    |
| Raises          | `ValueError` if a requirement fails |

`SquareRoot` is a node class, not the composition `power(f, Integer(1) /
Integer(2))`. IEEE 754 requires a correctly rounded square root but does not
require a correctly rounded power. So the two differ in floating-point
arithmetic, and the result of `power` depends on the math library.

### Tensor

Requires: [`AbstractExpression`](#abstractexpression), [Degree](#degree).

See also: [Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative).

A concrete subclass of `AbstractExpression` that stacks expressions along a
new first axis.

| `Tensor(entries)` |                                                                                                           |
| ----------------- | --------------------------------------------------------------------------------------------------------- |
| Free function     | `tensor(entries)`                                                                                         |
| Properties        | `entries`, a non-empty `tuple` of `AbstractExpression`s with the same value shape `s` and the same degree |
| Static type       | Value shape `(n,) + s`, with $n$ the number of entries                                                    |
| Raises            | `ValueError` if a requirement fails                                                                       |

A matrix is a `Tensor` of `Tensor`s. A `Tensor` of differential forms is an
array-valued differential form, so the new axis is an axis of the value shape
and never a slot.

A value shape may contain the extent 0, for example the value shape `(0,)` of
the coordinate function of `RD(0)`. `Tensor` cannot build an expression with
the extent 0, because an empty `tuple` has no entry to give `s`. An array with
an extent 0 has no entries, so the array is the `Zero` of that value shape.

Examples, with `zero, one, two = Integer(0), Integer(1), Integer(2)`:

- `tensor((one, two))` has value shape `(2,)` and no support.
- With the 1-forms $\mathrm{d}x$ and $\mathrm{d}y$ on `RD(2)`,
  `tensor((dx, dy))` has value shape `(2,)` and the degree `(1,)`, and is the
  vector-valued 1-form $\mathrm{d}X$.
- `tensor((tensor((one, zero)), tensor((zero, one))))` has value shape
  `(2, 2)` and no support.
- With `r0 = Float(0.0)` and `r1 = Float(1.0)`,
  `tensor((tensor((r1, r0)), tensor((r0, r1))))` has `Float` entries, so the
  matrix is not equal to `tensor((tensor((one, zero)), tensor((zero, one))))`.
- `tensor((tensor((one, two)),))` has value shape `(1, 2)`.
- `tensor((one, tensor((one, two))))` raises `ValueError`, because the value
  shapes of the entries differ.
- `tensor((one, dx))` raises `ValueError`, because the degrees of the entries
  differ.

### Indexed

Requires: [`AbstractExpression`](#abstractexpression), [Index](#index),
[EllipsisIndex](#ellipsisindex), [Extents](#extents),
[AbstractIndexValue](#abstractindexvalue), [Free indices](#free-indices),
[Support](#support), [Degree](#degree).

See also: [Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative),
[Example: a paraboloid](#example-a-paraboloid).

> LLM generated content.

A concrete subclass of `AbstractExpression` that selects one component of a
tensor-valued expression.

| `Indexed(expression, indices)` |                                                                                                                                                                                                                                                 |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Free function                  | `indexed(expression, indices)`, with the overload `expression[indices]`                                                                                                                                                                         |
| Properties                     | `expression`, of rank $r$, and `indices`, a `tuple` of `int`, [`Index`](#index), [`EllipsisIndex`](#ellipsisindex), [`AbstractIndexValue`](#abstractindexvalue) and [`GatheredIndex`](#gatheredindex) entries, with at most one `EllipsisIndex` |
| Static type                    | Value shape `()`, and the free indices of the expression followed by the `Index` and `EllipsisIndex` entries and the free indices of the `GatheredIndex` entries, with each `Index` taken once                                                  |
| Raises                         | See the table below                                                                                                                                                                                                                             |

Indexing is full indexing only. Each `int`, each `Index`, each
`AbstractIndexValue` and each `GatheredIndex` entry covers one axis, and a `EllipsisIndex` entry covers every axis that the other entries do not
cover. Without a `EllipsisIndex`, the number of entries must equal $r$. With a
`EllipsisIndex`, the number of other entries must be at most $r$. Each `int` entry
$i_j$ must index an axis with an `int` extent $n_j$ and must satisfy
$0 \le i_j < n_j$, where `(n_0, ..., n_{r-1})` is the value
shape of the expression. Each `AbstractIndexValue` and each `GatheredIndex`
entry must have the extent $n_j$. Each `Index` entry $i_j$ takes the extent $n_j$. A
`EllipsisIndex` takes the `tuple` of the extents of the axes that the `EllipsisIndex`
covers.

An `Index` may occur in two or more entries, or may already be free in the
expression. All the occurrences are the same index, so the indexed axes take
the same value. `A[i, i]` is the diagonal entry of a matrix `A`, and
`index_sum(A[i, i], (i,))` is the trace. All the occurrences must have the same
extent. A `EllipsisIndex` may also already be free in the expression, with the
same rule and the `tuple` of extents of the `EllipsisIndex`. So with a
`EllipsisIndex` `r` that is free with the extents `(2, 3)` in an expression `B`
of value shape `(2, 3)`, `B[r]` is the block diagonal of `B`.

The rank $r$ counts the axes of the value shape only. The indices select axes of the value
shape and never a slot, so the slots of a differential form cannot be indexed.

`X[0]` and `A[0, 1]` are `indexed(X, (0,))` and `indexed(A, (0, 1))`.
Construction of `Indexed` raises:

| Error        | Cause                                                                                                       |
| ------------ | ----------------------------------------------------------------------------------------------------------- |
| `TypeError`  | An entry that is neither exactly `int`, `Index`, `EllipsisIndex`, `AbstractIndexValue` nor `GatheredIndex`. |
| `ValueError` | No entries.                                                                                                 |
| `IndexError` | A number of entries that does not cover exactly the rank. An expression with value shape `()` has rank 0.   |
| `ValueError` | More than one `EllipsisIndex` entry.                                                                        |
| `IndexError` | A negative `int` entry, or an `int` entry $i_j \ge n_j$.                                                    |
| `IndexError` | An `int` entry for an axis whose extent is not an `int`.                                                    |
| `IndexError` | An `AbstractIndexValue` or `GatheredIndex` entry whose extent differs from the extent of the axis.          |
| `ValueError` | An index with two extents, in two entries or in an entry and the expression.                                |

A negative index is rejected rather than converted, so the stored indices are
always the indices that the user wrote.

Examples:

- For the coordinate function $X$ of `RD(3)`, `X[2]` has value shape `()`.
- For a matrix `A` of value shape `(2, 2)`, `A[1, 0]` has value shape `()`,
  and `A[1]` raises `IndexError`.
- `X[-1]` and `X[True]` raise errors.
- `tensor((Integer(1) - X[0], X[0]))` has value shape `(2,)` and the support
  `RD(2)`.
- `tensor((tensor((X[0], X[1])),))` has value shape `(1, 2)`.
- For the vector-valued 1-form `dX = tensor((dx, dy))` on `RD(2)`, `dX[0]` is
  the 1-form $\mathrm{d}x$, with value shape `()` and the degree `(1,)`.
- In the paraboloid example, the height `X[0]**Integer(2) + X[1]**Integer(2)`
  is an expression on `RD(2)` with value shape `()`.

### ComponentTensor and IndexSum

Requires: [Index](#index), [EllipsisIndex](#ellipsisindex), [Free
indices](#free-indices), [Degree](#degree), [Indexed](#indexed).

See also: [ExteriorDerivative](#exteriorderivative).

`ComponentTensor` and `IndexSum` are concrete subclasses of
`AbstractExpression` that bind free indices.
`ComponentTensor(expression, indices)` turns the indices into tensor axes.
`IndexSum(expression, indices)` sums over the indices.

| `ComponentTensor(expression, indices)`, `IndexSum(expression, indices)` |                                                                                                                                                                                                                   |
| ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Free function                                                           | `component_tensor(expression, indices)`, `index_sum(expression, indices)`                                                                                                                                         |
| Properties                                                              | `expression`, a scalar, and `indices`, a non-empty `tuple` of distinct `Index` and `EllipsisIndex` entries free in the expression                                                                                 |
| Static type                                                             | The free indices of the expression without the indices. `ComponentTensor` has the value shape of the extents of the indices in order, with every extent of a `EllipsisIndex`, and `IndexSum` the value shape `()` |
| Raises                                                                  | `ValueError` if a requirement fails                                                                                                                                                                               |

The indices bind axes of the value shape, so a `ComponentTensor` of a
differential form is an array-valued differential form, and an `IndexSum` of
differential forms is a sum of differential forms of the same degree.

One `IndexSum` holds a whole contraction, such as the sum over both indices of
`A[i, j] * B[i, j]`. A pass that chooses the order of summation then finds the
contraction in one node, not in a chain of nested nodes.

Examples, with `i, j = index("i"), index("j")`, two vectors `a` and `b` of value
shape `(3,)`, and a matrix `A` of value shape `(2, 2)`:

- `A[i, j]` has value shape `()` and the free indices `((i, 2), (j, 2))`.
- `a[i] * b[i]` has value shape `()` and the free indices `((i, 3),)`.
- `index_sum(a[i] * b[i], (i,))` is the inner product of `a` and `b`, with no
  free indices.
- `index_sum(A[i, j] * A[i, j], (i, j))` is the sum of the squares of the
  entries of `A`, with no free indices.
- `component_tensor(A[i, j], (j, i))` has value shape `(2, 2)` and is the
  transpose of `A`, written with user indices.
- With the vector-valued 1-form `dX` of value shape `(3,)` and the degree
  `(1,)`, `index_sum(a[i] * dX[i], (i,))` is the 1-form
  $a_0 \mathrm{d}x + a_1 \mathrm{d}y + a_2 \mathrm{d}z$, with the degree
  `(1,)`.
- `a[i] + b[j]` has the free indices `((i, 3), (j, 3))`, the outer sum of `a`
  and `b`.
- `A[i, i]` is the diagonal entry, with the free indices `((i, 2),)`.
  `index_sum(A[i, i], (i,))` is the trace of `A`, and
  `component_tensor(A[i, i], (i,))` is the diagonal of `A`, for every extent of
  `A`, including an `AbstractExtent`.

> Note: `opt_einsum` can directly plan a contract path for a given contraction,
> but assumes dense operators.

### Tensor algebra with EllipsisIndex

Requires: [EllipsisIndex](#ellipsisindex), [Indexed](#indexed), [ComponentTensor
and IndexSum](#componenttensor-and-indexsum), [Sum](#sum), [Product](#product),
[Conjugate](#conjugate).

A `EllipsisIndex` writes an operation on tensors of every rank as one fixed
composition. The operations in the table below are not part of the core. A
library may provide the operations as free functions (see
[New operations](#new-operations)). With derived ellipsis indices `r` and `s`
and an ordinary index `k`:

| Operation          | Fixed composition                                              |
| ------------------ | -------------------------------------------------------------- |
| `scale(a, T)`      | `component_tensor(a * T[r], (r,))`                             |
| `tensor_sum(A, B)` | `component_tensor(A[r] + B[r], (r,))`                          |
| `outer(u, v)`      | `component_tensor(u[r] * v[s], (r, s))`                        |
| `dot(u, v)`        | `component_tensor(index_sum(u[r, k] * v[k, s], (k,)), (r, s))` |
| `inner(A, B)`      | `index_sum(A[r] * conjugate(B[r]), (r,))`                      |

`A[r] + B[r]` raises `ValueError` if the value shapes of `A` and `B` differ,
because the extents of `r` differ. The overload `A + B` calls `add`, so
`A + B` accepts only scalars (see [Arithmetic](#arithmetic)).

### Transpose and matrix product

Requires: [Index](#index), [Indexed](#indexed), [ComponentTensor and
IndexSum](#componenttensor-and-indexsum), [Product](#product).

Transpose and matrix product are free functions built with derived indices:

- `transpose(A)` is `component_tensor(A[k, l], (l, k))`, with
  $k$ = `Index(("uflx.transpose", A, 0))` and
  $l$ = `Index(("uflx.transpose", A, 1))`.
- `matrix_product(A, B)`, with the overload `A @ B`, is
  `component_tensor(index_sum(A[i, k] * B[k, j], (k,)), (i, j))`, with $i$, $k$
  and $j$ the derived indices `Index(("uflx.matrix_product", A, B, m))` for
  $m = 0, 1, 2$.

An operand that does not have rank 2, including a scalar operand, is rejected
with `IndexError` by the constructor of `Indexed`. Unequal inner extents in
`matrix_product` are rejected with `ValueError` by the constructor of `Product`,
because the index $k$ has two extents. So the overload `@` accepts only rank 2
operands. Contractions of other ranks, such as a matrix-vector product, are
written in index notation.

Examples, with a matrix `A` of value shape `(2, 2)`:

- `transpose(A)` is not equal to `component_tensor(A[i, j], (j, i))` for user
  indices `i` and `j`, because the index labels differ.

### Inverse

Requires: [`AbstractExpression`](#abstractexpression), [Degree](#degree),
[Division](#division), [Elementary
functions](#elementary-functions-very-incomplete), [Transpose and matrix
product](#transpose-and-matrix-product).

See also: [Expanding pull-backs](#expanding-pull-backs), [Applying
derivatives](#applying-derivatives).

| `Inverse(A)`  |                                                                                         |
| ------------- | --------------------------------------------------------------------------------------- |
| Free function | `inverse(A)`                                                                            |
| Requires      | $A$ with value shape `(n, k)`, $k \le n$, and the degree `()`, and $A$ real for $k < n$ |
| Static type   | Value shape `(k, n)`                                                                    |
| Raises        | `ValueError` if a requirement fails                                                     |

`Inverse(A)` is the pseudo-inverse $A^+ = (A^T A)^{-1} A^T$ of a matrix $A$
of rank $k$, the left inverse of $A$ with $A^+ A$ the identity. For $k = n$,
$A^+$ is the inverse matrix $A^{-1}$. If $A$ has rank less than $k$,
`Inverse(A)` is undefined. Construction does not check the rank of $A$.

For $k = n$, the adjugate of $A$ divided by the determinant of $A$ is a
fixed composition, and for $k < n$, `inverse(transpose(A) @ A) @ transpose(A)`
is a fixed composition. `Inverse` is a node class, not those compositions,
because the compositions differ in floating-point arithmetic from a
factorisation of $A$, as for `Division` and `SquareRoot`. A pass chooses how
to evaluate `Inverse`, such as in closed form for a small matrix, from
$\det(A^T A)$, which the volume density of a domain of lower dimension also
needs. One node class for every $k \le n$ gives one rule for every dimension,
such as the rule for `Components` in
[Expanding pull-backs](#expanding-pull-backs).

Examples, with `one, two, three = Integer(1), Integer(2), Integer(3)`:

- With the Jacobian matrix
  `J = tensor((tensor((two, one)), tensor((one, three))))` of the [affinely
  mapped triangle](#example-an-affinely-mapped-triangle), `inverse(J)` has
  value shape `(2, 2)`, and a later pass evaluates `inverse(J)` to
  $\frac{1}{5} \begin{pmatrix} 3 & -1 \\ -1 & 2 \end{pmatrix}$.
- With the Jacobian matrix `Je = tensor((tensor((one,)), tensor((one,))))` of
  value shape `(2, 1)`, of the edge from $(0, 0)$ to $(1, 1)$, `inverse(Je)`
  has value shape `(1, 2)`, and a later pass evaluates `inverse(Je)` to
  $\begin{pmatrix} 1/2 & 1/2 \end{pmatrix}$.
- `inverse(tensor((one, two)))` raises `ValueError`, because the value shape
  `(2,)` is not the value shape of a matrix.
- `inverse(transpose(Je))` with value shape `(1, 2)` raises `ValueError`,
  because $k > n$.
- With the coordinate coframe `dX` of `RD(2)`, `inverse(tensor((dX, dX)))`
  raises `ValueError`, because the degree is `(1,)`.

### Operator overloads

Requires: [`AbstractExpression`](#abstractexpression), [Sum](#sum),
[Product](#product), [Division](#division), [Power](#power), [Negation and
subtraction](#negation-and-subtraction), [Absolute](#absolute),
[Indexed](#indexed), [Transpose and matrix
product](#transpose-and-matrix-product).

See also: [Numbers in expressions](#numbers-in-expressions), [Appendix: Free
functions](#appendix-free-functions).

Every overload in the [table of free functions](#appendix-free-functions) is a method of `AbstractExpression` (see
[`AbstractExpression`](#abstractexpression)). No other class defines
arithmetic, matrix or indexing overloads. Overloading should be considered
carefully outside of standard mathematical meanings.

An operator overload calls the free function with the operands unchanged, and
the constructor validates (see [Numbers in expressions](#numbers-in-expressions)).
`AbstractExpression` defines no reflected overloads, such as `__radd__`.

`__getitem__` is the one exception. Python passes `f[i]` to `__getitem__` as the
single entry `i`, and passes `f[i, j]` as the `tuple` `(i, j)`. `__getitem__`
converts a single entry into a `tuple` of one entry and calls `indexed`. The
conversion exists because of Python indexing syntax.

This base class is sketched in pseudo-code according to the principles outlined
above. Later classes will not be sketched in detail - `AbstractExpression` is
the only class with overloads.

```
class AbstractExpression(Base, ABC):
    """Abstract base class for expressions."""
    @property
    @abstractmethod
    def value_shape(self) -> tuple[int | AbstractExtent, ...]:
        """The value shape of the expression."""

    @property
    @abstractmethod
    def free_indices(
        self,
    ) -> tuple[
        tuple[
            Index | EllipsisIndex,
            int | AbstractExtent | tuple[int | AbstractExtent, ...],
        ],
        ...,
    ]:
        """The free indices and extents, in order of first occurrence."""

    @property
    @abstractmethod
    def scalar_type(self) -> Literal["real", "complex"]:
        """The scalar type of the expression."""

    @property
    @abstractmethod
    def support(self) -> AbstractCoordinateDomain | DecomposedDomain | None:
        """The support of the expression, or None for an expression with no support."""

    @property
    @abstractmethod
    def degree(self) -> tuple[int, ...]:
        """The differential form degree of each slot, with every entry at least 1."""

    @property
    @abstractmethod
    def argument_spaces(
        self,
    ) -> tuple[AbstractFunctionSpace | tuple[AbstractFunctionSpace, ...], ...]:
        """The space of each slot of a form, or () for an expression that is not a form."""

    # Operator overloading must be in the class.
    def __add__(self, other: AbstractExpression) -> AbstractExpression:
        """Add."""
        return add(self, other)

    def __sub__(self, other: AbstractExpression) -> AbstractExpression:
        """Subtract."""
        # Calls the free function subtract.
        ...

    def __mul__(self, other: AbstractExpression) -> AbstractExpression:
        """Multiply."""
        # Calls the free function multiply.
        ...

    def __matmul__(self, other: AbstractExpression) -> AbstractExpression:
        """Matrix multiply."""
        # Calls the free function matrix_product.
        ...

    def __truediv__(self, other: AbstractExpression) -> AbstractExpression:
        """Divide."""
        # Calls the free function divide.
        ...

    def __neg__(self) -> AbstractExpression:
        """Negate."""
        # Calls the free function negate.
        ...

    def __abs__(self) -> AbstractExpression:
        """Absolute value."""
        # Calls the free function absolute.
        ...

    def __pow__(self, exponent: AbstractExpression) -> AbstractExpression:
        """Raise to a power."""
        # Calls the free function power.
        ...

    def __getitem__(self, indices: object) -> AbstractExpression:
        """Index."""
        # Converts a single entry into a tuple of one entry, then calls indexed.
        ...

    def __iter__(self) -> NoReturn:
        """Raise TypeError, because expressions are not iterable."""
        # Without __iter__, Python iterates through __getitem__.
        ...

# A free function
def add(first: AbstractExpression, second: AbstractExpression) -> Sum:
    """Add two scalar differential forms of the same degree.

    Args:
        first: The first summand.
        second: The second summand.

    Returns:
        The sum of the two expressions, with the summands in the given order.

    Raises:
        TypeError: From the constructor of Sum, if a summand is not an
            expression.
        ValueError: From the constructor of Sum, if a summand does not have
            value shape (), if an index has two extents, if the summands
            have different supports, or if the degrees differ.
    """
    # The constructor of Sum validates.
    return Sum(first, second)
```

### Numbers in expressions

Requires: [Literals](#literals), [Operator overloads](#operator-overloads),
[Tensor](#tensor), [Indexed](#indexed), [RD](#rd).

See also: [Coordinate](#coordinate).

Every number in an expression is written with literals, such as `Integer(2)`,
`Float(0.5)`, `Integer(1) / Integer(3)` or `ImaginaryUnit()`. A Python `int`, `float`,
`complex`, `Fraction` or NumPy scalar is never converted and is rejected with
`TypeError`. The type of every literal is a choice that the user makes
explicitly.

Examples, with `X = coordinate(RD(2))`:

- `X[0]**Integer(2)`, `Integer(2) * X[0]` and `Integer(1) - X[0]` are
  accepted.
- `X[0]**2`, `2 * X[0]`, `1 - X[0]` and `0.5 * X[0]` raise `TypeError`.
- `tensor((1, 2))` raises `TypeError`. `tensor((Integer(1), Integer(2)))` is
  accepted.
- `2 + X[0]` raises `TypeError`, because `AbstractExpression` defines no
  reflected overloads.

## Calculus on a coordinate domain

### Coordinate

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain),
[`AbstractExpression`](#abstractexpression), [Support](#support),
[Degree](#degree).

See also: [ExteriorDerivative](#exteriorderivative),
[ParameterisedDomain](#parameteriseddomain).

A concrete subclass of `AbstractExpression` for the coordinate function of a
coordinate domain.

Every concrete coordinate domain $P$ provides a coordinate function
`Coordinate(P)`. The coordinate function $X$ of $P$ is the terminal expression
that maps each point of $P$ to the $n$-tuple of that point:

$X : P \to R^n$, $X(p) = p$.

| `Coordinate(P)` |                                                          |
| --------------- | -------------------------------------------------------- |
| Free function   | `coordinate(P)`                                          |
| Properties      | `domain`, a coordinate domain or a decomposed domain $P$ |
| Static type     | Value shape `(n,)`, real                                 |

The coordinate function is a 0-form.
The 1-forms $\mathrm{d}X_i$ come from the exterior derivative, not from
`Coordinate`.

Only coordinate domains and decomposed domains provide a coordinate function.
The coordinate function of a decomposed domain $D$ maps each point of $D$ to
the coordinates of the point in the ambient coordinate domain of $D$, with the
value shape `(n,)` for the dimension $n$ of the ambient coordinate domain (see
[DecomposedDomain](#decomposeddomain)). The coordinates of the points of a
domain are given by the parameterisation of the domain (see
[ParameterisedDomain](#parameteriseddomain)), which is built from the
coordinate function of the parameter domain.

Examples:

- `coordinate(RD(3))` is `Coordinate(RD(3))`, with value shape `(3,)` and the
  support `RD(3)`.
- `coordinate(RD(3)) == coordinate(RD(3))` is `True`, and
  `coordinate(RD(2)) == coordinate(RD(3))` is `False`.
- In the paraboloid example, `X = coordinate(RD(2))` has value shape `(2,)`.

### Valid pair

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain),
[Support](#support).

An expression and a coordinate domain $(f, P)$ are a valid pair if $P$ is an
`AbstractCoordinateDomain` or a [`DecomposedDomain`](#decomposeddomain) and $f$
has the support $P$ or no support. Otherwise
construction raises `TypeError` for $P$ or `ValueError` for $f$. An $f$ with
no support is a constant on $P$.

### ExteriorDerivative

Requires: [Valid pair](#valid-pair), [Coordinate](#coordinate),
[Degree](#degree), [Product](#product), [Indexed](#indexed).

See also: [Applying derivatives](#applying-derivatives).

A concrete subclass of `AbstractExpression` for the exterior derivative
$\mathrm{d}\omega$ of a differential form $\omega$ on a coordinate domain $P$.

| `ExteriorDerivative(f, P)` |                                                                     |
| -------------------------- | ------------------------------------------------------------------- |
| Free function              | `exterior_derivative(f, P)`                                         |
| Properties                 | `expression` $f$ and `domain` $P$                                   |
| Requires                   | A [valid pair](#valid-pair) $(f, P)$, and $f$ with at most one slot |
| Static type                | The value shape of $f$                                              |
| Raises                     | `TypeError` for $P$ and `ValueError` for $f$ if a requirement fails |
| Evaluated by               | [`apply_derivatives`](#applying-derivatives)                        |

An $f$ with no support is a constant on $P$.

The exterior derivative acts on differential forms with at most one slot. If
$f$ has the degree `()`, the exterior derivative has the degree `(1,)`. If $f$
has the degree `(k,)`, the exterior derivative has the degree `(k + 1,)`. If
$f$ has two or more slots, construction raises `ValueError`, because the
exterior derivative of a double form needs a connection. A degree greater than
the dimension of $P$ gives the zero form, as for [`Product`](#product).

The exterior derivative acts on each component of an array-valued differential
form, so the axes of the value shape pass through unchanged.

No fixed composition of the other node classes differentiates, so
`ExteriorDerivative` is a node class. `ExteriorDerivative` is lazy, and
[`apply_derivatives`](#applying-derivatives) evaluates the derivative. The
exterior derivative of `Coordinate(P)` is the vector-valued 1-form
$\mathrm{d}X$ and remains.

Examples on $R^3$:

```python
P = RD(3)
X = coordinate(P)
x, y, z = X[0], X[1], X[2]

def d(f):
    return exterior_derivative(f, P)

# The coordinate 1-forms, each with value shape () and the degree (1,).
dx, dy, dz = d(x), d(y), d(z)

# The coordinate coframe, with value shape (3,) and the degree (1,).
dX = d(X)
```

A pass that applies the rules of `apply_derivatives` gives the familiar
results of vector calculus:

| Expression                                   | Degree | After a pass                                                               |
| -------------------------------------------- | ------ | -------------------------------------------------------------------------- |
| `d(x * y * z)`                               | `(1,)` | $yz \, \mathrm{d}x + xz \, \mathrm{d}y + xy \, \mathrm{d}z$, the gradient  |
| `d(-y * dx + x * dy)`                        | `(2,)` | $2 \, \mathrm{d}x \wedge \mathrm{d}y$, the curl $(0, 0, 2)$                |
| `d(x * dy * dz + y * dz * dx + z * dx * dy)` | `(3,)` | $3 \, \mathrm{d}x \wedge \mathrm{d}y \wedge \mathrm{d}z$, the divergence 3 |
| `d(d(x * y * z))`                            | `(2,)` | The zero 2-form, because the curl of a gradient is zero                    |
| `d(d(-y * dx + x * dy))`                     | `(3,)` | The zero 3-form, because the divergence of a curl is zero                  |
| `d(Integer(1))`                              | `(1,)` | The zero 1-form, because a constant has no gradient                        |

The second row is the reason that half the integral of
$-y \, \mathrm{d}x + x \, \mathrm{d}y$ around a closed curve in the plane is
the enclosed area, by Stokes' theorem.

- `d(dx * dy * dz)` has the degree `(4,)` and is the zero form, because
  `RD(3)` has dimension 3.
- With a double form $S$ with the degree `(1, 1)`, such as a metric, `d(S)`
  raises `ValueError`.
- `exterior_derivative(X, RD(2))` raises `ValueError`, because `X` has the
  support `RD(3)`.

### Components

Requires: [Valid pair](#valid-pair), [ExteriorDerivative](#exteriorderivative),
[Degree](#degree), [ComponentTensor and
IndexSum](#componenttensor-and-indexsum).

See also: [PullBack](#pullback), [Expanding components](#expanding-components),
[ParameterisedDomain](#parameteriseddomain).

A concrete subclass of `AbstractExpression` that turns the first slot of a
differential form into an axis of the value shape, with the coordinate frame of
a coordinate domain $P$ of dimension $n$, or of the ambient coordinate domain
of a decomposed domain $P$.

| `Components(f, P)` |                                                                                                                           |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------- |
| Free function      | `components(f, P)`                                                                                                        |
| Properties         | `expression` $f$ and `domain` $P$, with $n$ the dimension of $P$, or of the ambient coordinate domain of a decomposed $P$ |
| Requires           | A [valid pair](#valid-pair) $(f, P)$, and $f$ with at least one slot                                                      |
| Static type        | Value shape `f.value_shape + (n,)`                                                                                        |
| Raises             | `ValueError` if a requirement fails                                                                                       |
| Evaluated by       | [`expand_components`](#expanding-components)                                                                              |

For a decomposed domain $P$ of dimension $k$, the frame is the coordinate
frame of the ambient coordinate domain of $P$, of dimension $n$. For $k = n$,
that frame is a frame of $P$. For $k < n$, such as a surface in $R^3$, the
coordinate vector fields $\partial_j$ of the ambient coordinate domain are not
tangent to $P$, and $\partial_j$ in the definition below is the orthogonal
projection of $\partial_j$ onto the tangent space of $P$, with the Euclidean
metric of the ambient coordinate domain. The array is then the array of the
extension of $f$ that is zero on the normal space of $P$, such as the surface
gradient of a function on $P$.

With $\partial_j$ the coordinate vector fields of $P$, dual to the coordinate
1-forms $\mathrm{d}X_j$, the entry $j$ of the new axis is the interior product
of $\partial_j$ with the first slot of $f$,

$$
\mathrm{Components}(f, P)_j(v_2, \dots, v_k) = f(\partial_j, v_2, \dots, v_k).
$$

If the first slot of $f$ has the degree $k > 1$, the first slot of the result has the
degree $k - 1$. If the first slot of $f$ has the degree 1, the result does not
have that slot. For example, the degree `(1,)` gives `()`, the degree `(2,)`
gives `(1,)`, and the degree `(1, 1)` gives `(1,)`.

No fixed composition of the other node classes turns a slot into an axis, so
`Components` is a node class. The reverse direction is a fixed composition,
`index_sum(v[i] * dX[i], (i,))` for a vector `v` of value shape `(n,)`.

The components depend on the coordinates of $P$. So a pull-back does not
commute with `Components`. `PullBack(Components(f, Q), D)` is the array of the
components of $f$ on $Q$, evaluated at $\psi(p)$. If $P$ has the dimension of
$Q$, `expand_pull_backs` moves the pull-back inside `Components` with the
inverse of the Jacobian matrix (see
[Expanding pull-backs](#expanding-pull-backs)). Otherwise
`PullBack(Components(f, Q), D)` remains until `expand_components` evaluates the
`Components` (see [Expanding components](#expanding-components)).

With `Components`, the arrays of standard differential calculus are fixed
compositions. With the examples of [ExteriorDerivative](#exteriorderivative)
on `RD(3)`, and `f = x * y * z`:

| Expression                              | Value shape | Degree | Meaning                                                      |
| --------------------------------------- | ----------- | ------ | ------------------------------------------------------------ |
| `components(d(f), P)`                   | `(3,)`      | `()`   | The gradient $(yz, xz, xy)$, after evaluation                |
| `components(d(components(d(f), P)), P)` | `(3, 3)`    | `()`   | The Hessian of $f$                                           |
| `components(dX, P)`                     | `(3, 3)`    | `()`   | The identity matrix                                          |
| `components(dx * dy, P)`                | `(3,)`      | `(1,)` | $(\mathrm{d}y, -\mathrm{d}x, 0)$                             |
| `components(components(dx * dy, P), P)` | `(3, 3)`    | `()`   | The antisymmetric matrix of $\mathrm{d}x \wedge \mathrm{d}y$ |

- `components(x, P)` raises `ValueError`, because `x` has no slot.
- For a parameterised domain $D$ with the parameter domain $P$,
  `components(exterior_derivative(D.parameterisation, P), P)` is the Jacobian
  matrix of value shape `(m, n)`, with $m$ the dimension of the ambient
  coordinate domain of $D$.
- For a decomposed domain `M` of dimension 2 in `RD(3)` and a function `u` on
  `M`, `components(exterior_derivative(u, M), M)` is the surface gradient of
  `u`, with value shape `(3,)`, and
  `components(exterior_derivative(coordinate(M), M), M)` is the orthogonal
  projection onto the tangent plane of `M`, with value shape `(3, 3)`.

### Volume form

Requires: [Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative),
[Product](#product), [Literals](#literals).

See also: [Integral](#integral), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

The free function `volume_form(P)` is the coordinate volume form of an
`AbstractCoordinateDomain` $P$ of dimension $n$,

$$
\mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}.
$$

With `dX = exterior_derivative(coordinate(P), P)`, `volume_form(P)` is
`dX[0] * ... * dX[n - 1]`, nested to the left as Python evaluates `*`, so
`volume_form(RD(3))` is `(dX[0] * dX[1]) * dX[2]`. The degree is `(n,)`. For $n = 0$,
`volume_form(P)` is `Integer(1)`. The free function reads the dimension of
$P$, so the free function is one fixed composition for each dimension. The
volume form fixes the orientation of $P$ (see [Integral](#integral)). The
volume form needs no metric. A Riemannian volume form is a density times
`volume_form(P)`.

Examples:

- `volume_form(RD(3))` has the degree `(3,)` and is
  $\mathrm{d}x \wedge \mathrm{d}y \wedge \mathrm{d}z$.
- With `T = UnitSimplex(2)`, `integral(volume_form(T), T)` is the area $1/2$
  of `T`.
- `volume_form(RD(0))` is `Integer(1)`.

### Integral

Requires: [Valid pair](#valid-pair), [Volume form](#volume-form),
[Support](#support), [Degree](#degree).

See also: [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

At this point we have everything we need to integrate an expression over a
single domain. Note we do not need to explicitly introduce any
topological notions of the domain to build a basic integral.

`Integral(f, P)` is the integral of a top form $f$ over a coordinate domain
$P$ of dimension $n$, with the orientation
$\mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}$ of the coordinates $X$
of $P$. With $f = f_{0 \dots n-1} \, \mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1}$,

$$
\int_P f = \int_P f_{0 \dots n-1} \,\mathrm{d}X,
$$

with the Lebesgue measure $\mathrm{d}X$ on $P$. For $n = 0$, the top form is a
0-form, and the integral is the value of $f$ at the single point of $R^0$.

| `Integral(f, P)` |                                                                                          |
| ---------------- | ---------------------------------------------------------------------------------------- |
| Free function    | `integral(f, P)`                                                                         |
| Properties       | `expression` $f$ and `domain` $P$, an `AbstractCoordinateDomain` of dimension $n$        |
| Requires         | A [valid pair](#valid-pair) $(f, P)$, with $f$ of the degree `(n,)`, or `()` for $n = 0$ |
| Static type      | The value shape and the free indices of $f$                                              |
| Raises           | `TypeError` for $P$ and `ValueError` for $f$ if a requirement fails                      |

An array-valued differential form integrates component by component. The
free function `integral(f, D)` also accepts a decomposed domain, as a sum of
integrals over the parameter domains of the pieces (see [Integral over a
decomposed domain](#integral-over-a-decomposed-domain)). `Integral` is lazy, and a later
pass might evaluate the integral, for example by quadrature.

Note that an `Integral` is an expression, so integrals can be added,
multiplied and used wherever an expression is accepted.

#### Examples

With `T = UnitSimplex(2)` and `X = coordinate(T)`:

- `integral(X[0] * volume_form(T), T)` is $\int_T X_0 \,\mathrm{d}X$. A later
  pass evaluates the integral to $1/6$.
- `integral(X[0], T)` raises `ValueError`, because `X[0]` is a 0-form, not a
  2-form.
- `integral(X[0] * volume_form(T), T) * integral(X[1] * volume_form(T), T)` is a
  product of two integrals, with value shape `()` and no support.

## Maps between domains

### Domains

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain), [Node class or
free function](#node-class-or-free-function).

See also: [ParameterisedDomain](#parameteriseddomain).

A domain $D$ is the image $\psi(P)$ of a coordinate domain $P$ under a
parameterisation $\psi$ into a coordinate domain $Q$. The core has two classes
of domain, [`ParameterisedDomain`](#parameteriseddomain), which stores $P$,
$\psi$ and $Q$, and [`DecomposedDomain`](#decomposeddomain), a union of
families of parameterised domains, such as the cells of a mesh.

A particular kind of domain, such as an affinely mapped triangle or a point, is
a free function that returns a `ParameterisedDomain`, not a subclass. Every
kind of domain is a fixed composition of $P$, $\psi$ and $Q$, so every kind of
domain is a free function (see
[Node class or free function](#node-class-or-free-function)).

### ParameterisedDomain

Requires: [Domains](#domains),
[AbstractCoordinateDomain](#abstractcoordinatedomain),
[Coordinate](#coordinate), [Support](#support), [Valid pair](#valid-pair), [Free indices](#free-indices),
[Scalar type](#scalar-type), [Degree](#degree).

See also: [PullBack](#pullback), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

A concrete subclass of `Base` for a domain $D$, the image of a parameter
domain $P$ under a parameterisation $\psi$ into an ambient coordinate domain
$Q$.

| `ParameterisedDomain(P, psi, Q)` |                                                                                                                                                                                                                                 |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Free function                    | `parameterised_domain(P, psi, Q)`                                                                                                                                                                                               |
| Properties                       | `parameter_domain` $P$, an `AbstractCoordinateDomain` of dimension $k$, `parameterisation` $\psi$, and `ambient_domain` $Q$, an `AbstractCoordinateDomain` or a [`DecomposedDomain`](#decomposeddomain)                         |
| Requires                         | A [valid pair](#valid-pair) $(\psi, P)$, $\psi$ real with value shape `(n,)` for the dimension $n$ of $Q$, or of the ambient coordinate domain of a decomposed $Q$, the degree `()` and the argument spaces `()`, and $k \le n$ |
| Static type                      | `free_indices`, the free indices of $\psi$                                                                                                                                                                                      |
| Raises                           | `TypeError` for $Q$, and `ValueError` if another requirement fails                                                                                                                                                              |

$\psi$ maps each point of $P$ to a point of $Q$, and the image $\psi(P)$ is
$D$. The coordinates of $Q$ name the points of $D$. The 'topological dimension'
of $D$ is $k$, `D.parameter_domain.dimension`, and the 'geometric dimension' is
$n$, `D.ambient_domain.dimension`. Construction does not check that $\psi(P)$
lies in $Q$.

A $\psi$ with free indices makes $D$ a family of domains, with one domain for
each value of the free indices. For example, with an array `x` of value shape
`(m, n)` and indices `i` and `k`,
`parameterised_domain(RD(0), component_tensor(x[i, k], (k,)), RD(n))` is a
family of $m$ points in `RD(n)`. A family of domains writes a set of points
once, such as the points of a quadrature rule.

$Q$ is usually `RD(n)`. A domain in a smaller $Q$, such as an edge of the unit
triangle in the unit triangle, or a cell of a decomposed domain in the
decomposed domain, pulls back an expression on $Q$ (see
[PullBack](#pullback)).

Examples:

- A sphere in $R^3$ has geometric dimension 3, topological dimension 2 and the
  ambient coordinate domain `RD(3)`.
- With `T = UnitSimplex(2)`, `parameterised_domain(T, coordinate(T), T)` is the
  unit triangle as a domain, with $k = n = 2$.
- The hypotenuse of the unit triangle is a parameterised domain with $k = 1$
  and $n = 2$. $P$ is the unit interval in $R^1$, and
  $\psi(X) = (1 - X_0, X_0)$. $Q$ is `RD(2)` for the hypotenuse in the plane,
  and `UnitSimplex(2)` for the hypotenuse in the unit triangle.
- With `X = coordinate(RD(2))` and `h = X[0]**Integer(2) + X[1]**Integer(2)`,
  `parameterised_domain(RD(2), tensor((X[0], X[1], h)), RD(3))` is the
  paraboloid, with $k = 2$ and $n = 3$.
- A physical triangle in $R^2$ is a parameterised domain with $k = n = 2$. $P$
  is the unit triangle, $\psi$ is the affine map from $P$ onto the physical
  triangle, and $Q$ is `RD(2)`. See the affinely mapped triangle example.
- A point $p$ in $R^2$ is a parameterised domain with $k = 0$ and $n = 2$. $P$
  is `RD(0)`, $\psi$ is the constant `tensor((p0, p1))`, with no support, and $Q$
  is `RD(2)`.

### Example: a paraboloid

Requires: [RD](#rd), [Coordinate](#coordinate), [Tensor](#tensor),
[Indexed](#indexed), [ParameterisedDomain](#parameteriseddomain).

The paraboloid over `RD(2)` collects the examples of the sections above.

```python
# AbstractCoordinateDomain: RD(2), no edges.
P = RD(2)

# AbstractExpression: Coordinate(RD(2)), value shape (2,).
X = coordinate(P)

# AbstractExpression: Tensor((Indexed(X, (0,)), Indexed(X, (1,)), Sum(...))),
# value shape (3,).
psi = tensor((X[0], X[1], X[0]**Integer(2) + X[1]**Integer(2)))

# ParameterisedDomain(RD(2), psi, RD(3)), k = 2, n = 3.
D = parameterised_domain(P, psi, RD(3))
```

### Example: an affinely mapped triangle

Requires: [Example: the unit simplex](#example-the-unit-simplex),
[Coordinate](#coordinate), [Tensor](#tensor), [Components](#components),
[ParameterisedDomain](#parameteriseddomain).

The physical triangle with vertices $v_0 = (1, 0)$, $v_1 = (3, 1)$ and
$v_2 = (2, 3)$ is the image of the unit triangle under an affine map.

Important point is that there is no topology here, just geometry.

The affine map from `UnitSimplex(2)` onto the physical triangle is

$$\psi(X) = v_0 + X_0 (v_1 - v_0) + X_1 (v_2 - v_0)
= (1 + 2 X_0 + X_1, X_0 + 3 X_1).$$

```python
# AbstractCoordinateDomain: UnitSimplex(2), no edges.
T = UnitSimplex(2)

# AbstractExpression: Coordinate(UnitSimplex(2)), value shape (2,).
X = coordinate(T)

# AbstractExpression: Tensor(..., ...), value shape (2,).
one, two, three = Integer(1), Integer(2), Integer(3)
psi = tensor((one + two * X[0] + X[1], X[0] + three * X[1]))

# ParameterisedDomain(UnitSimplex(2), psi, RD(2)), k = n = 2.
K = parameterised_domain(T, psi, RD(2))
```

The Jacobian matrix of $\psi$ is the constant matrix

$$J = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix},$$

so $\det J = 5$, and the area of the physical triangle is $5$ times the area
$1/2$ of the unit triangle.

The physical triangle is a parameterised domain with $k = n$, not a coordinate
domain, although the physical triangle has non-empty interior in $R^2$. The
class records how the domain is built.

### PullBack

Requires: [ParameterisedDomain](#parameteriseddomain), [Valid pair](#valid-pair),
[Support](#support),
[Free indices](#free-indices), [Degree](#degree),
[ExteriorDerivative](#exteriorderivative), [Example: an
affinely mapped triangle](#example-an-affinely-mapped-triangle).

See also: [Expanding pull-backs](#expanding-pull-backs).

A concrete subclass of `AbstractExpression` for the pull-back $\psi^* f$ of a
differential form $f$ on $Q$ by the parameterisation $\psi$ of a domain $D$,
with $P$ the parameter domain and $Q$ the ambient coordinate domain of $D$.

| `PullBack(f, D)` |                                                                                                            |
| ---------------- | ---------------------------------------------------------------------------------------------------------- |
| Free function    | `pull_back(f, D)`                                                                                          |
| Properties       | `expression` $f$, and `domain` $D$, a `ParameterisedDomain` or a [`SideMap`](#restriction)                 |
| Requires         | A [valid pair](#valid-pair) $(f, Q)$, and an equal extent for an `Index` free in both $f$ and $D$          |
| Static type      | The value shape of $f$, and the free indices of $f$ followed by those of $D$, with each `Index` taken once |
| Raises           | `TypeError` for $D$, and `ValueError` if another requirement fails                                         |
| Evaluated by     | [`expand_pull_backs`](#expanding-pull-backs)                                                               |

The second operand is a domain, not an expression, because a domain states
$Q$. The pull-back by a `SideMap` is a [restriction](#restriction).

For a 0-form $f$, the pull-back is the composition $f \circ \psi$, the
expression on $P$ with the value $f(\psi(p))$ at each point $p$ of $P$. For a
differential form with slots, the pull-back acts on each slot,

$$
(\psi^* f)_p(v_1, \dots, v_k) = f_{\psi(p)}(\mathrm{d}\psi \, v_1, \dots, \mathrm{d}\psi \, v_k),
$$

so `PullBack` accepts any degree, including a double form such as a metric.
The axes of the value shape are not transformed. An axis of an expression on
$Q$ is an array in the coordinate frame of $Q$, so the pull-back keeps the
frame of $Q$. A slot with a degree greater
than the dimension of $P$ gives the zero form, so a 2-form pulled back to a
curve, and a 1-form pulled back to a point, are zero.

`PullBack` is lazy, and [`expand_pull_backs`](#expanding-pull-backs)
evaluates the pull-back. A `PullBack` of a downstream terminal or of a
[`Function`](#function) remains, and the library of the terminal evaluates the
terminal at $\psi(p)$ and transforms the slots.

A domain $E$ in the parameter domain of $D$ has the image $\psi(E)$, a subset
of $D$. The image is the parameterised domain with the parameterisation
`pull_back(D.parameterisation, E)` on the parameter domain of $E$, in the
ambient coordinate domain of $D$, such as `E` in the example below.

Example, with `T`, `X`, `psi`, `one` and `K` of the affinely mapped triangle:

```python
# The unit interval, the parameter domain of the edge, with coordinate t.
S = UnitSimplex(1)
t = coordinate(S)[0]

# The edge of T from (1, 0) to (0, 1), as a domain in T.
e = parameterised_domain(S, tensor((one - t, t)), T)

# PullBack(f, e), with the support UnitSimplex(1), value shape () and the degree ().
f = X[0] * X[1]
fe = pull_back(f, e)

# PullBack(dx, e), with the support UnitSimplex(1) and the degree (1,).
dx = exterior_derivative(X[0], T)
dxe = pull_back(dx, e)

# The edge of K from (3, 1) to (2, 3), as a domain in RD(2).
E = parameterised_domain(S, pull_back(psi, e), RD(2))
```

- A later pass evaluates `fe` to $(1 - t) t$ and the parameterisation of `E`
  to $(3 - t, 1 + 2 t)$.
- A later pass evaluates `dxe` to $\mathrm{d}(1 - t) = -\mathrm{d}t$, because
  $X_0 = 1 - t$ on the edge.
- `pull_back(dx * exterior_derivative(X[1], T), e)` has the degree `(2,)` and
  is the zero form, because the edge has dimension 1.
- The Jacobian matrix of the parameterisation of `E` is
  $\begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix} \begin{pmatrix} -1 \\ 1
  \end{pmatrix} = \begin{pmatrix} -1 \\ 2 \end{pmatrix}$, so the edge of `K`
  has the length $\sqrt{5}$.
- `pull_back(f, parameterised_domain(S, tensor((one - t, t)), RD(2)))`
  raises `ValueError`, because that domain lies in `RD(2)`, not in `T`.

### Concept: orientation of a domain

Requires: [ParameterisedDomain](#parameteriseddomain),
[Volume form](#volume-form), [PullBack](#pullback).

An orientation is a choice of direction, i.e. which way an edge runs, or which
side of a surface is up. A coordinate domain $P$ is oriented by the order of
the coordinates of $P$, as fixed by the volume form of $P$ (see
[Volume form](#volume-form)). A domain $D$ = `ParameterisedDomain(P, psi, Q)`
takes the orientation of $P$ through $\psi$. For example, an edge runs in the
direction in which $\psi$ moves as the parameter increases. For a map $r$ of
$P$ onto $P$, $\psi$ and $\psi \circ r$ parameterise the same set. The two
orientations are equal if the Jacobian determinant of $r$ is positive, and
opposite if the Jacobian determinant of $r$ is negative. For example,
$r(t) = 1 - t$ reverses an edge, and $r(X) = (X_1, X_0)$ reverses a triangle.
The orientation needs no metric.

If $k = n$, $D$ also has the orientation of $Q$. With $J$ the Jacobian matrix
of $\psi$,

$$
\psi^* (\mathrm{d}Y_0 \wedge \dots \wedge \mathrm{d}Y_{n-1})
= \det J \, \mathrm{d}X_0 \wedge \dots \wedge \mathrm{d}X_{n-1},
$$

with $Y$ the coordinates of $Q$ and $X$ the coordinates of $P$. The two
orientations agree where $\det J > 0$, and differ where $\det J < 0$. The core
has no node class for the determinant, because $\det J$ is the coefficient of
`pull_back(volume_form(Q), D)`.

If $k < n$, such as for a surface in $R^3$, $Q$ gives $D$ no orientation. Each
piece of a decomposed domain of dimension $k < n$ takes the orientation of the
parameter domain of the piece, and the core does not relate the orientations
of two cells. An integral of a $k$-form over such a domain, such as in Stokes'
theorem on a surface, needs the cells to be oriented consistently, as does a
term that turns a $k$-form into a 0-form or a $(k - 1)$-form into a vector
field, such as the scalar curl of a 1-form on a surface. The mesh library
provides a
consistent orientation where one is needed, such as an expression with the
sign $+1$ or $-1$ for each cell, and the core needs nothing more (see
[Example: the flux through a surface](#example-the-flux-through-a-surface)).
An integral of a 0-form uses a density and needs no orientation (see
[Integral over a decomposed domain](#integral-over-a-decomposed-domain)).

### Concept: integral over a domain

Requires: [PullBack](#pullback), [Integral](#integral),
[Absolute](#absolute), [Division](#division),
[Concept: orientation of a domain](#concept-orientation-of-a-domain).

The integral of a $k$-form $\omega$ on $Q$ over a domain $D$ =
`ParameterisedDomain(P, psi, Q)` of dimension $k$, with the orientation of $D$,
is

$$
\int_D \omega = \int_P \psi^* \omega,
$$

which is `integral(pull_back(omega, D), P)`. Reversing the orientation of $D$
changes the sign of the integral.

An integral that must not depend on the orientation, such as the legacy UFL
integral `f * dx` over a cell whose orientation is arbitrary, multiplies by
the sign $s = \det J / |\det J|$ for $k = n$. Then
$\int_P s \, \psi^* (f \, \mathrm{d}Y_0 \wedge \dots \wedge \mathrm{d}Y_{n-1})
= \int_P (f \circ \psi) \, |\det J| \, \mathrm{d}X$.

### Example: orientation of an edge and of a triangle

Requires: [Concept: orientation of a domain](#concept-orientation-of-a-domain),
[Concept: integral over a domain](#concept-integral-over-a-domain),
[Components](#components), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

The examples use `T`, `X`, `one`, `S`, `t` and `e` of the
[PullBack](#pullback) example, and `K` of the
[affinely mapped triangle example](#example-an-affinely-mapped-triangle).

The edge `e`, with $\psi(t) = (1 - t, t)$, runs from vertex 1 to vertex 2 of
the unit triangle. The edge with $\psi(t) = (t, 1 - t)$ is the same set with
the opposite orientation, so the integral of $x \, \mathrm{d}y$ changes sign:

```python
x, y = X[0], X[1]
dy = exterior_derivative(y, T)

# The edge e reversed, the same set with the opposite orientation.
e_reversed = parameterised_domain(S, tensor((t, one - t)), T)

# A later pass evaluates the two integrals to 1/2 and -1/2.
integral(pull_back(x * dy, e), S)
integral(pull_back(x * dy, e_reversed), S)
```

On `e_reversed`, $x = t$ and $\mathrm{d}y$ pulls back to $-\mathrm{d}t$, so
the integral is $\int_0^1 -t \, \mathrm{d}t = -1/2$.

The triangle `K` has $\det J = 5$, so the orientation of `K` agrees with the
orientation of `RD(2)`:

```python
# The pull-back of the volume form of RD(2), det J dX_0 ^ dX_1, a 2-form on T.
w = pull_back(volume_form(RD(2)), K)

# The 2-form as the antisymmetric matrix [[0, det J], [-det J, 0]] of value shape (2, 2).
W = components(components(w, T), T)

# Entry [0, 1] follows the order of dX_0 ^ dX_1, so the entry is det J and not -det J.
# A later pass evaluates det J to 5.
det_J = W[0, 1]

# The sign of det J, 1 here.
s = det_J / absolute(det_J)

# The area of K. A later pass evaluates the integral to 5/2, whatever the sign of det J.
integral(s * w, T)
```

Swapping $v_1$ and $v_2$ gives $\det J = -5$ and $s = -1$. The integral of `w`
alone is then $-5/2$, and the integral of `s * w` is still $5/2$.

### Example: integrating over an edge

Requires: [PullBack](#pullback), [Integral](#integral),
[Components](#components), [Elementary functions (very
incomplete)](#elementary-functions-very-incomplete).

The integral over the edge `e` of the [PullBack](#pullback) example is the
integral over `S` of the pull-back by `e`. The examples use `T`, `X`, `one`,
`S`, `t`, `e` and `dx` of that example.

```python
x, y = X[0], X[1]
dy = exterior_derivative(y, T)
```

A 1-form integrates over the edge without a metric. On the edge,
$x = 1 - t$ and $\mathrm{d}y$ pulls back to $\mathrm{d}t$, so

```python
# A later pass evaluates the integral to 1/2.
integral(pull_back(x * dy, e), S)
```

is $\int_0^1 (1 - t) \, \mathrm{d}t = 1/2$.

The flux of a vector field $u$ through the edge is the integral of the 1-form
$u_0 \, \mathrm{d}y - u_1 \, \mathrm{d}x$, which also needs no metric. For
$u = X$, the 1-form pulls back to $(1 - t) \, \mathrm{d}t + t \, \mathrm{d}t$:

```python
# A later pass evaluates the integral to 1.
integral(pull_back(x * dy - y * dx, e), S)
```

The flux is 1. The divergence of $X$ is 2 and the area of `T` is $1/2$, so
the flux of $X$ out of `T` is 1. The 1-form $x \, \mathrm{d}y - y \,
\mathrm{d}x$ is zero on the other two edges of `T`, so the whole flux passes
through `e`. This is the divergence theorem, without a normal vector and
without a metric.

The integral of a scalar function with respect to arc length needs the length
element $|\psi'| \, \mathrm{d}t$, which comes from the Euclidean metric of
`T`. The metric enters as an expression that the user writes:

> This example writes the length element by hand, for a curve in the
> Euclidean plane, because `Integral` over a coordinate domain takes the
> density from the user. The general concept is the volume form of the metric
> that $\psi$ induces, $\sqrt{\det(J^T J)} \, \mathrm{d}t$. The free function
> `integral(f, D)` over a decomposed domain builds that density for each piece
> (see [Integral over a decomposed
> domain](#integral-over-a-decomposed-domain)).

```python
# The Jacobian matrix of the edge, with value shape (2, 1).
J = components(exterior_derivative(e.parameterisation, S), S)

# The length element of the edge.
ds = square_root(J[0, 0] * J[0, 0] + J[1, 0] * J[1, 0]) * volume_form(S)

# A later pass evaluates the integral to the length sqrt(2) of e.
integral(ds, S)

# A later pass evaluates the integral to sqrt(2) / 6.
integral(pull_back(x * y, e) * ds, S)
```

### Example: integrating over the boundary of a triangle

Requires: [ParameterisedDomain](#parameteriseddomain),
[ComponentTensor and IndexSum](#componenttensor-and-indexsum),
[Example: integrating over an edge](#example-integrating-over-an-edge).

The three edges of `T` are one family of domains, with a free index. The
example uses the names of
[Example: integrating over an edge](#example-integrating-over-an-edge).

```python
zero = Integer(0)
i, k = index("i"), index("k")

# The start and the end vertex of each edge, counterclockwise.
starts = tensor((tensor((zero, zero)), tensor((one, zero)), tensor((zero, one))))
ends = tensor((tensor((one, zero)), tensor((zero, one)), tensor((zero, zero))))

# The family of the three edges of T, with the free index i.
psi_edges = component_tensor(starts[i, k] + t * (ends[i, k] - starts[i, k]), (k,))
edges = parameterised_domain(S, psi_edges, T)

# A later pass evaluates the integral to 1/2.
index_sum(integral(pull_back(x * dy, edges), S), (i,))
```

The result is the integral of $x \, \mathrm{d}y$ over the boundary of `T`,
which is $\int_T \mathrm{d}x \wedge \mathrm{d}y = 1/2$ by Stokes' theorem.
Without the `IndexSum`, the integral has the free index `i`, and the value
for each `i` is the integral over one edge.

The orientation of each edge comes from the order of `starts` and `ends`. Every
domain of a family has the same parameter domain, here `S`.

### Points and point evaluation

Requires: [RD](#rd), [ParameterisedDomain](#parameteriseddomain),
[PullBack](#pullback), [Integral](#integral), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

A point is a domain! A point in a coordinate domain $Q$ of dimension $n$ is
`parameterised_domain(RD(0), x, Q)`, with $x$ a constant expression of value
shape `(n,)`. The core needs no node class and no free function for a point.

The value of a 0-form $f$ of any value shape at a point $p$ is
`integral(pull_back(f, p), RD(0))`, because the integral of a 0-form over
`RD(0)` is the value at the single point (see [Integral](#integral)).
`pull_back(f, p)` is an expression on `RD(0)`, and the integral has no support, so
the integral combines with any expression.

A point in `T` is given in the coordinates of `T`. For a point $x$ of `K`, the
surrounding environment supplies $\psi^{-1}(x)$. For a non-affine $\psi$, the
surrounding environment supplies a numerical approximation of $\psi^{-1}(x)$,
for example by Newton's method.

Examples, with `T`, `X`, `psi` and `K` of the affinely mapped triangle, and
`f = X[0] * X[1]`:

```python
third = Integer(1) / Integer(3)

# ParameterisedDomain: the centroid of T, as a point in T.
c = parameterised_domain(RD(0), tensor((third, third)), T)

# Integral: the value of f at the centroid, with no support and value shape ().
# A later pass evaluates the integral to 1/9.
fc = integral(pull_back(f, c), RD(0))

# ParameterisedDomain: the centroid of K, as a point in RD(2). A later pass
# evaluates pull_back(psi, c) to (2, 4/3).
C = parameterised_domain(RD(0), pull_back(psi, c), RD(2))
```

- With a second point
  `d = parameterised_domain(RD(0), tensor((third, Integer(2) / Integer(3))), T)`,
  `integral(pull_back(f, c), RD(0)) + integral(pull_back(f, d), RD(0))` is the action of
  a sum of two Dirac measures on $f$.
- `parameterised_domain(RD(0), tensor((third, third)), RD(2))` is the point
  $(1/3, 1/3)$ in the plane. With that point as `q`, `pull_back(f, q)` raises
  `ValueError`, because `q` lies in `RD(2)`, not in `T`.
- `parameterised_domain(RD(0), X, T)` raises `ValueError`, because `X` has the
  support `T`, not `RD(0)`.

### Example: quadrature

> This is very much inspired by what GEM/TSFC does but expresses it in core
> primitives. Clearly by combining these ideas UFLx has quite some additional
> expressivity. Perhaps having this in UFLx will make TSFC simpler, if it
> eventually targets UFLx.

Requires: [Points and point evaluation](#points-and-point-evaluation),
[ParameterisedDomain](#parameteriseddomain), [PullBack](#pullback),
[Integral](#integral), [Components](#components), [Tensor](#tensor),
[ComponentTensor and IndexSum](#componenttensor-and-indexsum),
[Extents](#extents).

This example is a 'look ahead' to expressing the concept of quadrature
integration in the core language.

A quadrature rule replaces the integral by the inner product of the weights
with the values at the points. With `T`, `X` and `f = X[0] * X[1]` of
[Points and point evaluation](#points-and-point-evaluation), the edge midpoint
rule is:

```python
zero, half = Integer(0), Integer(1) / Integer(2)
sixth = Integer(1) / Integer(6)
i, k = index("i"), index("k")

points = tensor((tensor((half, zero)), tensor((half, half)), tensor((zero, half))))
weights = tensor((sixth, sixth, sixth))

# A family of three points in T, with the free index i.
xs = parameterised_domain(RD(0), component_tensor(points[i, k], (k,)), T)

q = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

The points are one domain with the free index `i`, so `f` occurs once in `q`,
and a form compiler can evaluate `q` as a loop over the points.

A later pass evaluates `q` to $1/24$, which is
`integral(f * volume_form(T), T)`, because the rule is exact for $f$. The two
expressions are not equal, so quadrature is a pass, and the rule is an
argument of the pass. The pass finds the coefficient $c$ of a top form
$c \, \mathrm{d}X_0 \wedge \mathrm{d}X_1$ on `T` as
`components(components(c * volume_form(T), T), T)[0, 1]`. The result is an
expression of the core.

After the passes, `q` is

```python
index_sum(weights[i] * (points[i, 0] * points[i, 1]), (i,))
```

`expand_pull_backs` replaces `X` by the parameterisation of `xs`. A pass
replaces the integral over `RD(0)` of a 0-form by the 0-form, and a pass
replaces `component_tensor(points[i, k], (k,))[0]` by `points[i, 0]`. A form
compiler could generate:

```python
weights = (1 / 6, 1 / 6, 1 / 6)
points = ((1 / 2, 0.0), (1 / 2, 1 / 2), (0.0, 1 / 2))


def kernel() -> float:
    q = 0.0
    for i in range(3):
        q += weights[i] * (points[i][0] * points[i][1])
    return q
```

The `IndexSum` over `i` is the loop over the points, and the literals are
tables.

With `points` and `weights` given by leaves of a
[language extension](#language-extensions), of value shapes `(q, 2)` and
`(q,)` for an `AbstractExtent` `q`, the same expression is runtime quadrature
(see [Example: integrating over a cut
cell](#example-integrating-over-a-cut-cell)).
The number of points is known only at assembly, and the loop bound `q` is an
input of the kernel.

### Example: integrating over a cut cell

Requires: [Example: quadrature](#example-quadrature), [Extents](#extents).

In CutFEM, a level set $\phi$ cuts a cell, and an integral covers only the
part of the cell where $\phi < 0$. The solver computes a quadrature rule for
that part of each cut cell, so the number of points differs from cell to cell.

> The points and the weights are leaves of a small
> [language extension](#language-extensions) for kernel inputs:
> `RuntimeExtent(label)`, an `AbstractExtent`,
> `RuntimeArray(label, value_shape)`, an array such as quadrature data, the
> coefficients of a function on a cell, or a constant, and
> `RuntimeIndex(label, extent)`, an `AbstractIndexValue`, such as the local
> number of a facet.

With `T`, `f`, `i` and `k` of the quadrature example:

```python
# These are language extensions, not the core.
q = runtime_extent("q")

points = runtime_array("points", (q, 2))
weights = runtime_array("weights", (q,))

# The family of q points in the part of T where phi < 0.
# At this point we have no Function (phi) or if().
xs = parameterised_domain(RD(0), component_tensor(points[i, k], (k,)), T)

a = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

For all the cells of a topology, see the open item on ragged indexing in
[Extents](#extents).

### Example: the centroid of a triangle

> LLM came up with this example.

Requires: [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle), [PullBack](#pullback),
[Integral](#integral), [Division](#division),
[ComponentTensor and IndexSum](#componenttensor-and-indexsum),
[Points and point evaluation](#points-and-point-evaluation).

An integral is an expression, so the centroid of `K` is a ratio of integrals.
The example uses `T` and `K` of the affinely mapped triangle.

```python
Y = coordinate(RD(2))
vol = volume_form(RD(2))
k = index("k")

# A later pass evaluates the area of K to 5/2.
area = integral(pull_back(vol, K), T)

# A later pass evaluates the centroid of K to (2, 4/3).
centroid = component_tensor(integral(pull_back(Y[k] * vol, K), T) / area, (k,))
```

The pull-back of `vol` by `K` gives the factor $\det J = 5$, without a
determinant. The centroid is the parameterisation of the point `C` of
[Points and point evaluation](#points-and-point-evaluation), which is the
image of the centroid of `T`. The two expressions are not equal, and a later
pass evaluates both to $(2, 4/3)$.

> At this point, we can do no more on geometry before we do some work in
> topology (cells, facets, families, orientations) and functions. Then, we can
> return to geometry later.

## Topology

### AbstractCellType

Requires: [The base object `Base`](#the-base-object-base).

See also: [Reference cell](#reference-cell),
[AbstractTopology](#abstracttopology).

An abstract subclass of the base object `AbstractCellType` for the type of a
cell, with numbered vertices.

| `AbstractCellType` |                                                                         |
| ------------------ | ----------------------------------------------------------------------- |
| Provides           | `dimension`, an `int` $n \ge 0$, and `num_vertices`, an `int` $m \ge 1$ |
| Defined by         | A subclass outside the core, such as `Simplex`                          |

The core does not number the facets or the other sub-entities of a cell type.
A sub-entity of a cell is named by the vertex numbers of the cell that the
sub-entity contains, such as the vertex numbers of a facet, and the mesh or a
library such as Basix decides which sets of vertices are sub-entities (see
[Sides of a topology](#sides-of-a-topology)).

A topology of polygons, as for the virtual element method, has a cell type for
each number of vertices. A cell type whose number of vertices is known only at
run time is not yet covered.

### Example: the simplex as a cell type

> LLM generated content.

Requires: [AbstractCellType](#abstractcelltype).

`Simplex(n)` is the simplex of dimension $n \ge 0$ as a cell type, with
$n + 1$ vertices. `Simplex(0)` is a point.

```python
class Simplex(AbstractCellType):
    identifier = "library.Simplex"

    def __init__(self, n: int):
        if type(n) is not int:
            raise TypeError("Dimension must be an int.")
        if n < 0:
            raise ValueError("Dimension must be non-negative.")
        self._n = n

    @property
    def dimension(self) -> int:
        return self._n

    @property
    def num_vertices(self) -> int:
        return self._n + 1
```

### AbstractTopology

Requires: [The base object `Base`](#the-base-object-base), [Index](#index),
[Extents](#extents), [Free indices](#free-indices),
[AbstractCellType](#abstractcelltype).

See also: [Sides of a topology](#sides-of-a-topology),
[Restriction](#restriction),
[Integral over a decomposed domain](#integral-over-a-decomposed-domain),
[DecomposedDomain](#decomposeddomain).

An abstract subclass of the base object `AbstractTopology` for a set of cells
of given cell types, without coordinates.

| `AbstractTopology` |                                                                                             |
| ------------------ | ------------------------------------------------------------------------------------------- |
| Provides           | `cell_types`, `num_cells` and `sides`, and `side_tables(i, side)` for a topology with sides |
| Derived            | `dimension`, the dimension of the cell types                                                |
| Defined by         | A subclass outside the core, such as the topology of the mesh of a solver                   |

- `cell_types` is a `tuple` of `AbstractCellType`s $C_0, C_1, \dots$ of the
  same dimension, one for each group of cells. A cell of type $j$, a cell of
  group $j$, has the vertices of $C_j$, numbered as in $C_j$.
  The cell types are distinct, except that two groups with different sides may
  have the same cell type (see [Sides of a topology](#sides-of-a-topology)).
- `num_cells` is a `tuple` with the number of cells of each type, an `int`,
  such as for a macro element, or an `AbstractExtent`, such as for the mesh of
  a solver. For the mesh of a solver, the solver loops over the cells, and a
  kernel computes the terms of one cell.
- `sides` is `()`, or the sides of [Sides of a topology](#sides-of-a-topology).

Two topologies are the same topology only if the two topologies are equal as
base objects, so a concrete `AbstractTopology` stores the arguments that
identify the topology, such as a label.

The free function `cell_index(T, j)` returns the derived index
`Index(("uflx.cell_index", T, j, 0))` for the cells of type $j$, numbered from
0. The index has the extent `T.num_cells[j]`. The index differs between types
and between topologies, even if the numbers of cells are equal. Restriction
and the integral over a topology use the index, so user code rarely needs the
index.

By the rules of [free indices](#free-indices), an expression with the free
index `cell_index(T, j)` has one value for each cell of type $j$.

> Legacy note: In legacy UFL, a `Mesh` holds both the topology and the
> coordinate element, and the measure `dx` hides the sum over the cells. In
> UFLx, the topology, the parameterisations of the pieces of a decomposed
> domain and the `IndexSum` over `cell_index(T, j)` are separate.

Examples:

- A square divided into two triangles has the cell types `(K,)`, for the
  triangle cell type `K`, and the numbers of cells `(2,)`.
- The Alfeld split of a triangle, a macro element, has the numbers of cells
  `(3,)`.
- A mesh of a solver with triangles and quadrilaterals has two cell types, and
  two `AbstractExtent`s as the numbers of cells.

### Sides of a topology

Requires: [AbstractCellType](#abstractcelltype),
[AbstractTopology](#abstracttopology).

See also: [Restriction](#restriction).

A topology $S$ whose cells are sub-entities of cells of other topologies has
sides. `S.sides` has one entry for each cell type of $S$. Entry $i$ is a
`tuple` of one or more sides with distinct names, each a triple `(name, T, j)`
of a side name, a `str`, a topology $T$ and a cell type number $j$ of $T$. Two
entries of the same cell type of $S$ differ, so a group of cells of $S$ is
fixed by the cell type and the sides together. Each
cell of type $i$ of $S$ is seen from one cell of type $j$ of $T$ on each side,
and is a whole sub-entity of that cell, the image of cell type $i$ of $S$ under
the affine map through the paired vertices (see [Restriction](#restriction)).
The sub-entity may be the cell itself, with the vertices in the order of the
cell, a facet, or a sub-entity of lower dimension, such as an edge of a
tetrahedron or a vertex. The concrete mesh pairs each vertex of a
cell of $S$ with a vertex of the cell on each side, the same vertex for a
matching facet, and a vertex one period away for a periodic boundary. A
topology whose cells are not seen from other cells has the sides `()`.

A concrete topology with sides also provides `side_tables(i, side)`, a pair of
[`AbstractIndexValue`](#abstractindexvalue) tables `(cells, vertices)` for the
side named `side` of group $i$, with the side `(side, T, j)`:

- `cells`, of the shape `(S.num_cells[i],)` and the extent `T.num_cells[j]`.
  Entry $e$ is the cell of type $j$ of $T$ on the side of the cell $e$ of $S$.
- `vertices`, of the shape `(S.num_cells[i], m)`, with $m$ the number of
  vertices of cell type $i$ of $S$, and the extent the number of vertices of
  cell type $j$ of $T$. Entry $(e, v)$ is the local vertex of the cell on the
  side that is paired with vertex $v$ of the cell $e$.

[Restriction](#restriction) reads the tables, and picks the value of an
expression on a cell of $S$ by the side name. The core does not check the
sides or the tables.

The concrete mesh decides which sub-entities are the cells of $S$, and from
which cells on each side each cell of $S$ is seen. For example, the cells of
$S$ may be the facets of a marked boundary or the cells of a marked subdomain,
and a facet of two cells may be seen from one cell only. The sides may have
different topologies, such as a mesh and a submesh, or two meshes that meet.
The topology of a side may itself have sides, so the cells of $S$ may be
facets of facets.

For a sub-entity of lower dimension than a facet, such as a vertex of a
triangle mesh, several cells contain the sub-entity, and the concrete mesh
chooses the one cell on a side. A function that is continuous across the cells
has the same restriction from each of the cells, and the restriction of any
other function depends on the choice.

The concrete mesh also names the sides, and decides which cell is on which
side. For example, the facets of an interface between fluid and solid
triangles of the one cell type 0 of `T` have the sides
`(("fluid", T, 0), ("solid", T, 0))`, and the cell on the side `"fluid"` is the
fluid triangle.

> Legacy note: The legacy UFL measures `ds` and `dS` integrate over a topology
> of facets with one and two sides, and legacy UFL names the two sides `'+'`
> and `'-'`. The legacy UFL measure `dx(1)` integrates over a topology of the
> cells of a subdomain, and the legacy UFL measure `dP` integrates over a
> topology of vertices.

A cell of $S$ that is not a whole sub-entity of the cell on a side, as with
hanging nodes, a mortar interface or a line through the cells of a 3D mesh, is
not yet covered. For such a side, the inclusion of
[Restriction](#restriction) is expected to be $\phi_c^{-1} \circ \psi_e$, with
$\psi_e$ the map of the cell $e$ of $S$ into the ambient coordinate domain and
$\phi_c^{-1}$ the inverse of the map of the cell $c$ on the side, a terminal
that a downstream library provides. Which parts of a cell of $S$ meet which
cells of a side stays outside the language.

Examples, with a topology `T` of triangles with the cell types `(K,)`,
`I = Simplex(1)` and `P = Simplex(0)`:

- A topology `B` of the boundary facets of `T` has the cell types `(I,)` and
  the sides `((("inside", T, 0),),)`. `cell_index(B, 0)` is the index of the
  facets.
- A topology of the facets shared by two triangles has the cell types `(I,)`
  and the sides `((("+", T, 0), ("-", T, 0)),)`.
- A topology of the corners of the boundary has the cell types `(P,)` and the
  sides `((("+", B, 0), ("-", B, 0)),)`, since each corner is seen from two
  facets of `B`.
- A topology of the triangles of a marked subdomain of `T` has the cell types
  `(K,)` and the sides `((("inside", T, 0),),)`.
- A topology of the vertices of `T`, such as for point loads, has the cell
  types `(P,)` and the sides `((("inside", T, 0),),)`, with one triangle chosen
  for each vertex.
- With a submesh `M` of the boundary facets of `T`, with the cell types
  `(I,)`, a topology of the facets has the sides
  `((("mesh", T, 0), ("submesh", M, 0)),)`, and each facet is a cell of `M`.
- A topology of the left boundary facets of a periodic `T` has the sides
  `((("left", T, 0), ("right", T, 0)),)`, with each vertex paired with the
  vertex one period to the right.
- For a surface mesh `U` of triangles in $R^3$, a topology of the edges where
  three triangles meet has the sides `((("a", U, 0), ("b", U, 0), ("c", U, 0)),)`.
- For a topology `M` of triangles and quadrilaterals with the cell types
  `(Kt, Kq)`, a topology of all the facets shared by two cells has three
  groups of the cell type `I`, with the cell types `(I, I, I)` and the sides
  `((("+", M, 0), ("-", M, 0)), (("+", M, 0), ("-", M, 1)), (("+", M, 1), ("-", M, 1)))`,
  for the facets between two triangles, between a triangle and a
  quadrilateral, and between two quadrilaterals. An integral over the topology
  has one term for each group, and a form compiler may generate one kernel for
  each group or one kernel for all of them.
- For a topology `U3` of prisms with the cell types `(W,)`, a topology of the
  boundary facets has the cell types `(K, Q)`, with the quadrilateral cell type
  `Q`, and the sides `((("inside", U3, 0),), (("inside", U3, 0),))`.

### Example: a static topology

> LLM generated content.

Requires: [AbstractTopology](#abstracttopology),
[Example: the simplex as a cell type](#example-the-simplex-as-a-cell-type),
[Sides of a topology](#sides-of-a-topology).

`StaticTopology(cell_types, cells, sides)` is a topology whose cells are known
during construction. Entry $j$ of `cells` is a `tuple` with the vertices of
each cell of type $j$, as vertex numbers of the topology. Local vertex $i$ of
cell $c$ of type $j$ is the vertex `cells[j][c][i]` of the topology. A
`StaticTopology` with sides has one topology on every side, uses the vertex
numbers of that topology, pairs each vertex with the same vertex, and computes
the side tables from the vertex numbers. For a
cell seen from two cells of the same group, the cell with the lower number is
on the first side, and otherwise the group of each side decides.
For a sub-entity of lower dimension than a facet, the cell with the lowest
number that contains the sub-entity is chosen.

> This is primarily useful for testing UFLx ideas later, but also covers
> Alfeld splits.

```python
class StaticTopology(AbstractTopology):
    identifier = "library.StaticTopology"

    def __init__(
        self,
        cell_types: tuple[AbstractCellType, ...],
        cells: tuple[tuple[tuple[int, ...], ...], ...],
        sides: tuple[tuple[tuple[str, AbstractTopology, int], ...], ...],
    ):
        if len(cells) != len(cell_types):
            raise ValueError("Need one tuple of cells for each cell type.")
        self._cell_types = cell_types
        self._cells = cells
        self._sides = sides

    @property
    def cell_types(self) -> tuple[AbstractCellType, ...]:
        return self._cell_types

    @property
    def num_cells(self) -> tuple[int, ...]:
        return tuple(len(c) for c in self._cells)

    @property
    def sides(self) -> tuple[tuple[tuple[str, AbstractTopology, int], ...], ...]:
        return self._sides
```

With `K = Simplex(2)` and `I = Simplex(1)`:

```python
# A square divided into two triangles, with the numbers of cells (2,).
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())

# The Alfeld split of a triangle at vertex 3, with the numbers of cells (3,).
alfeld = StaticTopology((K,), (((0, 1, 3), (1, 2, 3), (2, 0, 3)),), ())

# The facets of the boundary of the square, each seen from one triangle.
boundary = StaticTopology(
    (I,), (((0, 1), (1, 3), (3, 2), (2, 0)),), ((("inside", square, 0),),)
)

# The facet of the diagonal of the square, seen from two triangles.
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
```

A facet of a cell is a set of vertex numbers of the topology, found through the
facets of `K`. In `square`, the facet $\{1, 2\}$ is a facet of two cells and
the other four facets are facets of one cell. In `alfeld`, the three facets
that contain vertex 3 are facets of two cells, and the three other facets are
facets of one cell. `boundary` and `diagonal` have the numbers of cells `(4,)`
and `(1,)`.

`StaticTopology` does not check the sides. A `StaticTopology` has no
coordinates. The positions of the vertices come from the
[coordinates of a decomposed domain](#decomposeddomain).

## Reference cells

### Reference cell

Requires: [AbstractCoordinateDomain](#abstractcoordinatedomain), [RD](#rd),
[Tensor](#tensor), [ParameterisedDomain](#parameteriseddomain),
[Inverse](#inverse), [AbstractCellType](#abstractcelltype).

See also: [AbstractTopology](#abstracttopology).

An abstract subclass of `AbstractCoordinateDomain` and of `AbstractCellType`
for a reference cell `AbstractReferenceCell`, a cell type whose vertices have
coordinates.

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Provides                | `vertices`, a constant expression with no support and value shape `(m, n)`, with row $v$ the coordinates of vertex $v$, for dimension $n$ and $m$ vertices |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](#example-the-unit-simplex-as-a-reference-cell)                                                   |

The free function `affine_map(F, points)` returns the affine map on a reference
cell $F$ of dimension $k$ that sends vertex $v$ of $F$ to row $v$ of `points`
for $v \le k$, with `points` of value shape `(k + 1, m)`. With the vertices
$V_v$ of $F$, the matrix $E$ with the columns $V_{d+1} - V_0$, and the rows
$p_v$ of `points`, the map is

$$
X \mapsto p_0 + \sum_{d < k} (p_{d+1} - p_0) \, \big(E^{-1} (X - V_0)\big)_d,
$$

an expression on $F$ with value shape `(m,)`, a fixed composition for each $k$
with [Inverse](#inverse). For $k = 0$ the map is the constant $p_0$. The first
$k + 1$ vertices of $F$ must not lie in a hyperplane of $F$.

A sub-entity of $K$ with the vertex numbers $(w_0, w_1, \dots)$ and the
reference cell $F$ is `parameterised_domain(F, affine_map(F, points), K)`, with
row $v$ of `points` the vertex $w_v$ of $K$. The core does not check that the
image is a sub-entity of $K$.

> A downstream library such as Basix provides the vertex numbers of the facets,
> so a facet is an ordinary
> `ParameterisedDomain` in the form, and an integral over a facet is the
> integral of a pull-back, as in the [boundary
> example](#example-integrating-over-the-boundary-of-a-triangle).

Examples:

- The unit triangle has the vertices $(0, 0)$, $(1, 0)$ and $(0, 1)$. The
  sub-entity with the vertex numbers $(1, 2)$ and the unit interval as $F$ is
  the hypotenuse, with $\psi(t) = (1 - t, t)$, as in the
  [ParameterisedDomain](#parameteriseddomain) examples.
- The unit interval has the vertices $(0)$ and $(1)$, and the sub-entity with
  the vertex number $(1)$ and the reference point as $F$ is the end point.
- The reference point has the single vertex $()$, the point of $R^0$.

### Example: the unit simplex as a reference cell

> LLM generated content.

Requires: [Example: the unit simplex](#example-the-unit-simplex),
[Example: the simplex as a cell type](#example-the-simplex-as-a-cell-type),
[Literals](#literals), [Tensor](#tensor), [Reference cell](#reference-cell).

`ReferenceSimplex(n)` is the cell type `Simplex(n)` as a reference cell, the
unit simplex in $R^n$. Vertex 0 is the origin, and vertex $v \ge 1$ is the unit
vector $e_{v-1}$. `ReferenceSimplex(0)` is the reference point.

```python
class ReferenceSimplex(Simplex, AbstractReferenceCell):
    identifier = "library.ReferenceSimplex"

    @property
    def vertices(self) -> AbstractExpression:
        n = self.dimension
        if n == 0:
            return Zero((1, 0), (), (), ())
        return tensor(
            tuple(
                tensor(tuple(Integer(1 if k == v - 1 else 0) for k in range(n)))
                for v in range(n + 1)
            )
        )
```

`ReferenceSimplex(2)` and `ReferenceSimplex(1)` have the vertices of the unit
triangle and the unit interval in the [Reference cell](#reference-cell)
examples.

`ReferenceSimplex(2)` and `Simplex(2)` are different cell types, so a topology
of `Simplex(2)` cells has no reference cells. `ReferenceSimplex(2)` and
`UnitSimplex(2)` are the same set but different coordinate domains, so an
expression on one is not an expression on the other. Likewise,
`ReferenceSimplex(0)` and `RD(0)` are different coordinate domains.

## Decomposed domains

### DecomposedDomain

Requires: [AbstractTopology](#abstracttopology),
[Reference cell](#reference-cell), [ParameterisedDomain](#parameteriseddomain),
[PullBack](#pullback), [Free indices](#free-indices).

See also: [Sides of a topology](#sides-of-a-topology),
[AbstractFunctionSpace](#abstractfunctionspace).

A concrete subclass of `Base` for a domain that is a union of pieces, each a
family of parameterised domains.

| `DecomposedDomain(pieces)` |                                                                                                                                                           |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Free function              | `decomposed_domain(pieces)`                                                                                                                               |
| Properties                 | `pieces`, a non-empty `tuple` of [`ParameterisedDomain`](#parameteriseddomain)s with the same ambient coordinate domain $Q$ and the same dimension $k$    |
| Derived                    | `ambient_domain` $Q$, `dimension` $k$, and `topology`, $T$ if the free indices of piece $j$ are `(cell_index(T, j),)` for every $j$, and otherwise `None` |
| Raises                     | `TypeError` for a piece that is not a `ParameterisedDomain`, and `ValueError` if another requirement fails                                                |

A piece with free indices is a family, one domain for each value of the free
indices, and $D$ is the union of all the domains of all the pieces. A mesh has
one piece for each cell type $j$, `parameterised_domain(K_j, x_j, Q)`, with the
reference cell $K_j$ and the coordinates $x_j$ of the cells, with the free
index `cell_index(T, j)`. Then the topology of $D$ is $T$. A piece with an
ordinary `Index`, such as a cloud of points, gives no topology, and a single
`ParameterisedDomain` with no free indices is a decomposed domain of one
piece. The core does not check that the domains meet only in sets of measure
zero.

A decomposed domain $D$ is a domain, not a coordinate domain. The points of $D$
are named by the coordinates of $Q$, and $D$ may have a lower dimension than
$Q$, such as a surface in $R^3$. Two decomposed domains are the same only if
the two decomposed domains are equal as base objects, so expressions on two
decomposed domains have different supports and do not combine (see
[Support](#support)).

The free function `cell_domain(D, j)` returns piece $j$ with the ambient domain
$D$, `parameterised_domain(P_j, psi_j, D)`, with the free indices of piece $j$.
An expression $f$ with the support $D$ is defined through the pieces.
`pull_back(f, cell_domain(D, j))` is an expression on $P_j$, the value of $f$
on each domain of piece $j$. The pull-back of `coordinate(D)` is $\psi_j$, and
the pull-back of a function in a space on $D$ is evaluated by the library of
the space (see [PullBack](#pullback)).

[Restriction](#restriction) moves an expression on $D$ to the cells of a
topology with sides.

> Legacy note: In legacy UFL, a `Mesh` is a decomposed domain, and
> `SpatialCoordinate(mesh)` is `coordinate(D)`.

Example, with `K = ReferenceSimplex(2)`, the topology
`square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())` of the
[static topology example](#example-a-static-topology) on reference cells, an
array `y` of value shape `(4, 2)` with the positions $(0, 0)$, $(1, 0)$,
$(0, 1)$ and $(1, 1)$ of the vertices, and a table `cell_vertices` of the shape
`(2, 3)`, the extent 4 and the values `((0, 1, 2), (1, 3, 2))`, the vertices of
each triangle:

```python
X = coordinate(K)
c = cell_index(square, 0)
k, v = index("k"), index("v")

# The position of each local vertex of each triangle, of value shape (2, 3, 2),
# gathered from the vertex positions.
p = component_tensor(y[gathered_index(cell_vertices, (c, v)), k], (c, v, k))

# The affine map of each triangle.
x = component_tensor(
    p[c, 0, k] + (p[c, 1, k] - p[c, 0, k]) * X[0] + (p[c, 2, k] - p[c, 0, k]) * X[1],
    (k,),
)
D = decomposed_domain((parameterised_domain(K, x, RD(2)),))
```

- `D.topology` is `square`, because `x` has the free index `c`.
  `cell_domain(D, 0)` is the family of the two triangles, with the free index
  `c`.
- `pull_back(coordinate(D), cell_domain(D, 0))` is `x`, which a later pass
  evaluates on triangle 1 to $(1 - X_1, X_0 + X_1)$.
- With `x` of value shape `(3,)` and `RD(3)`, the decomposed domain is a
  surface of dimension 2 in `RD(3)`, and `components(du, D)` raises
  `ValueError`, because the frame of `RD(3)` is not a frame of the surface.

For a topology `M` of triangles and quadrilaterals, with the cell types
`(Kt, Kq)`, $D$ has one piece for each cell type, and a function `u` on `D` is
one expression, pulled back through each cell type:

```python
D = decomposed_domain(
    (parameterised_domain(Kt, x_tri, RD(2)), parameterised_domain(Kq, x_quad, RD(2)))
)

# On Kt with the free index cell_index(M, 0), and on Kq with cell_index(M, 1).
u_tri = pull_back(u, cell_domain(D, 0))
u_quad = pull_back(u, cell_domain(D, 1))
```

### Restriction

Requires: [DecomposedDomain](#decomposeddomain),
[Sides of a topology](#sides-of-a-topology), [GatheredIndex](#gatheredindex),
[PullBack](#pullback), [Reference cell](#reference-cell).

See also: [Integral over a decomposed domain](#integral-over-a-decomposed-domain).

The restriction of an expression on a decomposed domain $D$ to the cells of a
topology $S$ with sides, from the side $s$, is the pull-back along the
inclusion of each cell of $S$ into the cell on the side $s$. The side tables of
$S$ give the cell and the vertices of the inclusion, through
[GatheredIndex](#gatheredindex) entries, so the restriction is a pull-back and
needs no node class of its own.

For the group $i$ of $S$ with the side `(s, T, j)`, with $T$ the topology of $D$,
the side tables `(cells, vertices)`, the cell type $F$ of group $i$ and
`e = cell_index(S, i)`, two families of domains with the free index `e` make
the inclusion:

- The cell on the side, $B_i$, is `cell_domain(D, j)` with every `Indexed`
  entry `cell_index(T, j)` of the parameterisation replaced by
  `gathered_index(cells, (e,))`. The library of a function in a space on $D$
  evaluates the pull-back by $B_i$ on the cell `cells[e]`.
- The inclusion into the cell, $A_i$, is
  `parameterised_domain(F, affine_map(F, points), K_j)`, with row $v$ of
  `points` the vertex `gathered_index(vertices, (e, v))` of $K_j$. $A_i$ sends
  each vertex of the cell $e$ to the paired vertex of the cell on the side. For
  a vertex, $A_i$ is the constant vertex, and for the cell itself with the
  vertices in order, $A_i$ is the identity map.

The restriction of $f$ to group $i$ from the side $s$ is
`pull_back(pull_back(f, B_i), A_i)`, an expression on $F$ with the free index
`e`. Each side uses the vertex numbering of the cell $e$, so the restrictions
from two sides are expressions on the same $F$, with the same coordinates, and
the difference is the jump. A cell of $S$ with more than two sides, such as an
edge where three triangles meet, has no jump, and a form uses the restriction
from each side.

| `SideMap(D, E, s)` |                                                                                                         |
| ------------------ | ------------------------------------------------------------------------------------------------------- |
| Free function      | `side_map(D, E, s)`                                                                                     |
| Properties         | `ambient_domain` $D$, `parameter_domain` $E$, a decomposed domain with the topology $S$, and `side` $s$ |
| Requires           | $D$ has a topology, and every group of $S$ has a side named $s$ with the topology of $D$                |
| Raises             | `ValueError` if a requirement fails                                                                     |

A `SideMap` maps the cells of $E$ into $D$ through the side $s$.
`pull_back(f, SideMap(D, E, s))` is an expression with the support $E$, the
restriction of $f$ to all the groups at once. The pull-back of it by
`cell_domain(E, i)` is `pull_back(pull_back(f, B_i), A_i)` (see [Expanding
pull-backs](#expanding-pull-backs)).

The free function `side_domain(D, S)` returns the decomposed domain of $S$ with
the geometry of the first side of each group, with piece $i$
`parameterised_domain(F_i, pull_back(pull_back(coordinate(D), B_i), A_i), Q)`
for the first side, with $F_i$ the cell type of group $i$ and $Q$ the ambient
coordinate domain of $D$. The topology of the result is $S$. The core does not check that the other sides have the same
geometry. The free function `restriction(f, E, s)` returns
`pull_back(f, side_map(D, E, s))` with $D$ the support of $f$, or $f$ if $f$
has no support.

The pull-back gives the trace of each degree. The restriction of a 0-form is
the value on the cell of $S$, and the restriction of the 1-form $\mathrm{d}u$
keeps only the derivatives along the cell of $S$. To keep the whole array of
$\mathrm{d}u$, as for a normal derivative, restrict `components(du, D)`, whose
axes are in the frame of the ambient coordinate domain from every side. The
outward normal of the cell on a side is [`normal(D, E, s)`](#normal).

> Legacy note: The legacy UFL restriction `f('+')` is `restriction(f, E, "+")`
> on the decomposed domain $E$ of the interior facets, whose sides are named
> `'+'` and `'-'`.

Example, with `D` and `square` of the [DecomposedDomain](#decomposeddomain)
example, `I = ReferenceSimplex(1)`, a function `u` on `D`, and the diagonal of
`square`:

```python
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
E = side_domain(D, diagonal)
jump = restriction(u, E, "+") - restriction(u, E, "-")
```

- The side tables from `"+"` are the cells `(0,)` and the vertices `((1, 2),)`,
  and from `"-"` the cells `(1,)` and the vertices `((0, 2),)`, so
  $A_+(t) = (1 - t, t)$ and $A_-(t) = (0, t)$. The coordinates of `E` are
  $(1 - t, t)$.
- If `pull_back(u, cell_domain(D, 0))` is
  `a[c, 0] + a[c, 1] * X[0] + a[c, 2] * X[1]`, the two restrictions pulled
  back to the diagonal are $a_{0,0} + a_{0,1}(1 - t) + a_{0,2} t$ and
  $a_{1,0} + a_{1,2} t$, equal for all $t$ if `u` is continuous at the
  vertices 1 and 2.
- With the diagonal listed as `(2, 1)`, the vertices from `"+"` are
  `((2, 1),)`, $A_+(t) = (t, 1 - t)$, and the restriction of
  $\mathrm{d}x_0$ from `"+"` changes sign, from $-\mathrm{d}t$ to
  $\mathrm{d}t$.
- On the topology of the vertices of `square`, with the cell type
  `ReferenceSimplex(0)`, vertex 3 is seen from cell 1 as local vertex 1, so
  `u` restricts to $a_{1,0} + a_{1,1}$ there.
- For the facets `S` of a triangle and quadrilateral mesh, with the three
  groups of [Sides of a topology](#sides-of-a-topology), the same two lines
  with `E = side_domain(D, S)` give one jump for all three groups.

### Integral over a decomposed domain

Requires: [Integral](#integral),
[Concept: orientation of a domain](#concept-orientation-of-a-domain),
[Concept: integral over a domain](#concept-integral-over-a-domain),
[DecomposedDomain](#decomposeddomain).

The free function `integral(f, D)` for a decomposed domain $D$ of dimension $k$
integrates a 0-form against the density that the pieces induce, and a $k$-form
with the orientation of the pieces. $f$ has the support $D$ or no support, and
the integral is the sum over the pieces $j$ of

```python
# A 0-form f, with the density rho_j of piece j.
index_sum(integral(pull_back(f, cell_domain(D, j)) * rho_j * volume_form(P_j), P_j), indices_j)

# A k-form f, with the orientation of piece j.
index_sum(integral(pull_back(f, cell_domain(D, j)), P_j), indices_j)
```

with $P_j$ the parameter domain and `indices_j` the free indices of piece $j$,
and no `IndexSum` for a piece with no free indices. With the Jacobian matrix
$J_j$ = `components(exterior_derivative(psi_j, P_j), P_j)` of value shape
`(n, k)`, the density is $\rho_j = \sqrt{\det(J_j^T J_j)}$, the length, the
area or the volume of the image of a small region of $P_j$. The determinant of
a $k \times k$ matrix is a fixed composition for each $k$. For $k = n$,
$\rho_j = \lvert \det J_j \rvert$, and for $k = 0$, the empty determinant is 1,
so the integral over a set of points is the sum of the values.

The density needs no orientation, so the cells of a mesh may be oriented in any
way, and an integral of a positive 0-form is positive. The integral of a
$k$-form, such as for Stokes' theorem, uses the orientation of the pieces (see
[Concept: orientation of a domain](#concept-orientation-of-a-domain)).

An $f$ with the support of the ambient domain $Q$ of $D$ is pulled back by
piece $j$ itself in place of `cell_domain(D, j)`, so $D$ may be a set of points
or of curves in a mesh. A free index of $f$ that is also a free index of piece
$j$ takes the value of the member, so `d[i]` pairs a datum with each point of a
family with the index `i`. The free function reads only the pieces of $D$,
which are static, so the integral is one fixed composition for each $D$, and
`Integral` stays a node class over a coordinate domain.

> Legacy note: `f * dx`, `f * ds` and `f * dS` are `integral(f, D)` over the
> decomposed domain of the cells, of the boundary facets and of the interior
> facets, with the fields restricted to the facets in `f`. The three are one
> integral over three domains.

Examples, with `D` of the [DecomposedDomain](#decomposeddomain) example, a
function `u` on `D`, and the decomposed domain `B` of the boundary facets of
`D`, with the side `"inside"`:

- `integral(Integer(1), D)` is the area 1 of the square, the sum of two
  integrals over `K`.
- `integral(u, D)` is the legacy UFL `u * dx`.
- `integral(Integer(1), B)` is the perimeter 4, and
  `integral(restriction(u, B, "inside"), B)` is the legacy UFL `u * ds`.
- `integral(restriction(du, B, "inside"), B)` integrates the 1-form
  $\mathrm{d}u$ along the boundary facets, with the orientation of each facet.

### Normal

Requires: [Restriction](#restriction), [Reference cell](#reference-cell),
[Inverse](#inverse).

A concrete subclass of `AbstractExpression` for the unit outward normal of the
cell on a side, on the cells of a topology with sides.

| `Normal(D, E, s)` |                                                                                                       |
| ----------------- | ----------------------------------------------------------------------------------------------------- |
| Free function     | `normal(D, E, s)`                                                                                     |
| Properties        | `ambient_domain` $D$, `parameter_domain` $E$ and `side` $s$, as for [`SideMap`](#restriction)         |
| Requires          | The requirements of `SideMap(D, E, s)`, and $E$ of the dimension $k - 1$ for the dimension $k$ of $D$ |
| Static type       | Value shape `(n,)` for the dimension $n$ of the ambient coordinate domain, the support $E$, real      |
| Raises            | `ValueError` if a requirement fails                                                                   |
| Evaluated by      | [`expand_pull_backs`](#expanding-pull-backs)                                                          |

`Normal(D, E, s)` is a terminal of the geometry, as `Coordinate` is. The
pull-back by `cell_domain(E, i)`, with $B_i$ and $A_i$ of
[Restriction](#restriction) for the side $s$ and the reference cell $K_j$ of
the cell on the side, is

$$
\nu = \frac{K^T \hat m}{\lvert K^T \hat m \rvert},
$$

with $K$ the [inverse](#inverse) of the Jacobian matrix of $B_i$ at the
points of $A_i$, the pseudo-inverse for $k < n$,
and $\hat m$ a normal of the cell of $E$ in $K_j$, orthogonal to the columns of
the Jacobian matrix of $A_i$ and pointing away from the centroid of the
vertices of $K_j$. Both are fixed compositions for each $k$ and $n$, from the vertices
of $K_j$ and the side tables, so the core needs no numbering of the facets.
The centroid lies inside a convex reference cell, so $\hat m$ points out of
$K_j$. For $k = n$, the normals from two sides of one cell of $E$ are
opposite. For $k < n$, such as on the edges of a surface in $R^3$, $\nu$ is
the unit conormal, tangent to the cell on the side and orthogonal to the cell
of $E$. The conormals from two sides are opposite only where the cells on the
two sides have the same tangent space, so on a surface a jump weights each
side by the normal of that side, as in `u_p * n_p + u_m * n_m`.

`Normal` is the normal of a cell on a side, not the unit normal field of a
domain of dimension $n - 1$, such as of a surface in $R^3$. That field needs
an orientation of the domain (see [Concept: orientation of a
domain](#concept-orientation-of-a-domain)) and belongs to a downstream library.

Examples, with `B` of the [Integral over a decomposed
domain](#integral-over-a-decomposed-domain) examples, the decomposed domain `E`
of the interior facets of `D` with the sides `"+"` and `"-"`,
`grad_u = components(exterior_derivative(u, D), D)` and `i = index("i")`:

```python
n_B = normal(D, B, "inside")

# The flux of grad u out of the square, the legacy UFL dot(grad(u), n) * ds.
flux = integral(index_sum(restriction(grad_u, B, "inside")[i] * n_B[i], (i,)), B)

# The jump of the normal derivative across the diagonal, from the side "+".
n_E = normal(D, E, "+")
jump = restriction(grad_u, E, "+") - restriction(grad_u, E, "-")
penalty = integral(index_sum(jump[i] * n_E[i], (i,)) ** Integer(2), E)
```

### Concept: kernels

Requires: [DecomposedDomain](#decomposeddomain), [GatheredIndex](#gatheredindex),
[Integral over a decomposed domain](#integral-over-a-decomposed-domain).

> This section describes how a downstream library is expected to use the core,
> not a rule of the core.

An integral over a decomposed domain is an `IndexSum` over the cell index of
each group, of an expression on a reference cell that reads arrays through
gathered indices with the cell index. A form compiler splits the integral
there:

- The `IndexSum` over the cell index is the loop of an assembler over the
  cells, or over the cells of a topology with sides.
- A gathered index with the cell index, such as of a map from the local to the
  global degrees of freedom, of the coordinates or of the side tables, is the
  data of one cell, which the assembler packs and gives to a kernel.
- The rest, an expression on a reference cell without the cell index, is the
  kernel.

The inputs of a kernel are the leaves that remain after the split, each with
the gathered index that reads the leaf, such as the coefficients of a function
through a map of degrees of freedom on one side, or a runtime extent. So the
interface of a kernel follows from the expression, and a form compiler may
describe the interface beyond a fixed convention, such as for two sides of
different cell types.

The assembler adds the result of a kernel into a global array through the same
tables, which is outside the language. A whole-array compiler may instead keep
the integral whole, as a batched gather, one kernel for all the cells and a
scatter-add.

A kernel depends on no particular mesh. Before compiling a kernel, a pass
replaces the topologies, the tables, the arrays and the extents of a mesh by
placeholders numbered in the order of appearance, and renames the arguments of
each form by the slot, so one compiled kernel serves every mesh with the same
cell types and spaces, and every naming of the arguments.

## Function spaces

A function is a member of a function space, such as a coefficient or an
argument of a variational form.

### AbstractFunctionSpace

Requires: [DecomposedDomain](#decomposeddomain), [Degree](#degree),
[Scalar type](#scalar-type).

See also: [PullBack](#pullback).

An abstract subclass of the base object `AbstractFunctionSpace` for a space of
functions on a decomposed domain.

| `AbstractFunctionSpace` |                                                                                                                                                             |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Provides                | `domain`, a [`DecomposedDomain`](#decomposeddomain) $D$ or `None`, `degree`, `value_shape`, `scalar_type` and `argument_spaces`, `()` for a space of fields |
| Defined by              | A subclass outside the core, such as a subclass for the finite element spaces of Basix, and [`DualSpace`](#dualspace) in the core                           |

A space with the domain `None` has functions that do not depend on a domain.

The dual space of $V$ is [`dual(V)`](#dualspace), a space of the core.

A function in $V$ is an expression with the support $D$ and with the degree,
the value shape and the scalar type of $V$, defined through the cells of $D$
(see [DecomposedDomain](#decomposeddomain)). The concrete subclass defines the
pull-back of a function in $V$ by `cell_domain(D, j)` for every cell type $j$
of $D$, such as with one finite element for each cell type. The core does not
represent the polynomials, the degrees of freedom or the basis of $V$, which belong to the
concrete subclass. Two spaces are the same space only if the two spaces are
equal as base objects, so a concrete `AbstractFunctionSpace` stores the
arguments that identify the space, such as the finite element.

A space with the domain `None` and the argument spaces `()` is a space of
constants. A function in the space
is an expression with no support and the degree `()`, a single array with the
value shape and the scalar type of $V$. The degree is `()`, because an
expression with no support and a non-empty degree is a `Zero` (see
[Expressions as differential forms](#expressions-as-differential-forms)), and
a constant $k$-form is a constant array times the coordinate $k$-forms of a
domain. The integrals of a form decide the cells on which a constant acts. A
space of constants is needed for a constant unknown, such as a Lagrange
multiplier, because an argument comes from a space. Constant data needs no
space, because an expression with no support, such as a literal or an array
that a kernel receives, is already a constant.

> Legacy note: The legacy UFL `Real` space is a space of constants, and the
> legacy UFL `Constant` is constant data.

The degree and the value shape fix how a function in $V$ transforms between a
cell and the reference cell. The pull-back by a cell transforms each slot and
keeps the frame of each axis of the value shape (see [PullBack](#pullback)). So a space needs no
separate map to the reference cell:

| Space                             | Degree     | Value shape | Legacy UFL map      |
| --------------------------------- | ---------- | ----------- | ------------------- |
| Lagrange                          | `()`       | `()`        | identity            |
| Vector Lagrange                   | `()`       | `(n,)`      | identity            |
| Nédélec, $H(\mathrm{curl})$       | `(1,)`     | `()`        | covariant Piola     |
| Raviart–Thomas, $H(\mathrm{div})$ | `(n - 1,)` | `()`        | contravariant Piola |
| Discontinuous densities           | `(n,)`     | `()`        | $L^2$ Piola         |
| Regge                             | `(1, 1)`   | `()`        | double covariant    |

> Legacy note: In legacy UFL, a `FunctionSpace` holds a mesh and an element,
> and the element states the value shape, the map to the reference cell and
> the Sobolev space. In UFLx, the degree and the value shape replace the map,
> and the space holds a decomposed domain, the legacy mesh.

A mixed problem needs no mixed space, because a slot of a [form](#form) may
bind a `tuple` of functions in several spaces. The Sobolev space of $V$ is not
yet covered. The Sobolev space is expected to be stated through the
[restrictions](#restriction) of the functions of $V$ to the facets of the
topology of $D$ from the two sides of each facet.

Examples, with a decomposed domain `D` of `ReferenceSimplex(2)` cells in
`RD(2)`, so $n = 2$:

- The Lagrange space of polynomial degree 1 on `D` has the degree `()` and the
  value shape `()`.
- The Nédélec and the Raviart–Thomas spaces on `D` both have the degree `(1,)`,
  because $n - 1 = 1$. A Raviart–Thomas function is the vector field of legacy
  UFL contracted with the volume form.
- The discontinuous Lagrange space of polynomial degree 0 on `D`, as densities,
  has the degree `(2,)`.

### Function

Requires: [AbstractFunctionSpace](#abstractfunctionspace),
[`AbstractExpression`](#abstractexpression), [Support](#support),
[Degree](#degree).

See also: [PullBack](#pullback), [Restriction](#restriction),
[Variational forms](#variational-forms).

A concrete subclass of `AbstractExpression` for a function, a member of a
function space.

| `Function(V, label)` |                                                                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| Free function        | `function(V, label)`                                                                                                         |
| Properties           | `space` $V$, an `AbstractFunctionSpace`, and `label`, a `str` or a `tuple` with the rules of the label of an [Index](#index) |
| Static type          | The value shape of $V$, and no free indices                                                                                  |
| Raises               | `TypeError` if a requirement fails                                                                                           |

The label identifies the function, so two functions are the same only if the
spaces and the labels are equal. A function in a [dual space](#dualspace) is an
antilinear form, such as a known load.

A `Function` is a terminal. The core does not represent the values of a
function. The pull-back of a function by a cell of the domain is evaluated by
the library of the space (see [PullBack](#pullback)), so a library chooses how
a function is represented on a cell. The exterior derivative of a function
remains after [applying derivatives](#applying-derivatives).

A function is not an argument or a coefficient by itself. Whether a function is
an argument of a form is stated by the form (see [Form](#form)).

> Legacy note: The legacy UFL `Coefficient`, `TestFunction` and `TrialFunction`
> are all functions. A test function and a trial function are the functions
> that a form states as the arguments in the slots 0 and 1.

Examples, with a decomposed domain `D` and a space `V` on `D` of the degree
`()` and the value shape `()`:

- `u = function(V, "u")` has the support `D`, the degree `()` and the value
  shape `()`, and `function(V, "u") == u` is `True`.
- `exterior_derivative(u, D)` is the 1-form $\mathrm{d}u$ on `D`.
- With the decomposed domain `E` of the interior facets of `D`,
  `restriction(u, E, "+")` is `u` from the side `"+"`.

## Variational forms

A variational form is an expression that is abstracted over some functions, the
arguments, as a multilinear map is over the arguments of the map. A form binds
the arguments, as a [`ComponentTensor`](#componenttensor-and-indexsum) binds
indices, and an action gives a function to a slot of a form, as an
[`Indexed`](#indexed) gives an index to an axis.

### Argument spaces

Requires: [`AbstractExpression`](#abstractexpression),
[AbstractFunctionSpace](#abstractfunctionspace).

Every `AbstractExpression` derives the argument spaces, a `tuple` with the
space of each slot of a form, and `()` for an expression that is not a form.
The space of a slot is an `AbstractFunctionSpace` for a slot of one function,
or a `tuple` of at least two spaces for a slot of a `tuple` of functions, the
product of the spaces. So each slot has one representation. The arity is the number of argument spaces. A functional has the arity
0, a linear form the arity 1 and a bilinear form the arity 2, whatever the
number of functions in each slot.

A node passes on the argument spaces of an operand in which the node is linear
or antilinear, and every other operand must have the argument spaces `()`, or
construction raises `ValueError`. The [static type
table](#appendix-static-type-table) gives the rule of each node class. A leaf
class states the argument spaces as part of the contract of the class, as the
degree.

The argument spaces are static, as the value shape and the degree are, so a
form with a wrong arity is rejected during construction. The arity and the
degree are independent. The degree counts the slots of a differential form,
and the arity counts the functions that a form binds.

### DualSpace

Requires: [Argument spaces](#argument-spaces),
[AbstractFunctionSpace](#abstractfunctionspace), [Function](#function).

See also: [Action](#action).

A concrete subclass of `AbstractFunctionSpace` for the space of antilinear
forms on a space or on the product of spaces.

| `DualSpace(S)` |                                                                                                                                                               |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Free function  | `dual(S)`                                                                                                                                                     |
| Properties     | `space` $S$, an `AbstractFunctionSpace`, or a `tuple` of at least two spaces                                                                                  |
| Provides       | The domain `None`, the degree `()`, the value shape `()`, the argument spaces `(S,)` and the scalar type of $S$, complex if any space of a `tuple` is complex |
| Raises         | `TypeError` for an $S$ that is not a space or a `tuple` of spaces, and `ValueError` for a `tuple` of fewer than two spaces                                    |

A function in `dual(S)` is an antilinear form on $S$, as a form is in the slot
of a function under a `Conjugate`, such as the slot of $v$ in
$\int u \bar v$. So a residual and an assembled load have the same
linearity, and the difference of the two is a form. For a real scalar type,
antilinear and linear coincide. A fixed function in `dual(S)` is a known
antilinear form, such as a load, and a function bound by a [form](#form) is a
variable antilinear form. The pairing of a function `l` in `dual(S)` with a
function `w` in $S$ is `action(l, 0, w)`, antilinear in `w`, and the library
that provides the values of `l` evaluates the pairing.

For a `tuple` $S$, a function in `dual(S)` is one antilinear form on the product
of the spaces, such as the load of a mixed problem, although no single
function is in the product.

`dual(dual(V))` is the space of antilinear forms on `dual(V)`, with the
argument spaces `(dual(V),)`. `dual(dual(V))` is not `V`, because `V` has the
argument spaces `()`. The canonical map takes a function `w` in `V` to
`form(conjugate(action(l, 0, w)), (l,))`, with a function `l` in `dual(V)`,
linear in `w` and antilinear in `l`, and a library may store both with the
same values.

> Legacy note: The legacy UFL `Cofunction` and `Coargument` are a fixed and a
> bound function in `dual(V)`. In legacy UFL, the dual of the dual space is
> $V$. In UFLx, the canonical map is a form, so the spaces keep distinct
> argument spaces.

Examples, with the spaces `V` and `Q` and the functions `v` in `V` and `q` in
`Q` of the [Form](#form) examples:

- `l = function(dual(V), "l")` has the argument spaces `(V,)`, and
  `action(l, 0, v)` is the pairing $\langle l, v \rangle$.
- `m = function(dual((V, Q)), "m")` has the argument spaces `((V, Q),)`, and
  `action(m, 0, (v, q))` is the pairing of `m` with `(v, q)`.
- With a function `w` in `V` and a function `l` in `dual(V)`,
  `form(conjugate(action(l, 0, w)), (l,))` has the argument spaces
  `(dual(V),)`, the static type of a function in `dual(dual(V))`, and
  evaluates every antilinear form on `V` at `w`.
- With the [interpolation extension](#example-interpolation), a space `W` and
  a function `l` in `dual(W)`, `form(action(l, 0, interpolate(v, W)), (v,))`
  is the adjoint of the interpolation from `V` to `W` applied to `l`, with the
  argument spaces `(V,)`.

### Form

Requires: [Argument spaces](#argument-spaces), [Function](#function).

| `Form(expression, arguments)` |                                                                                                                                                                                                    |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Free function                 | `form(expression, arguments)`                                                                                                                                                                      |
| Properties                    | `expression`, an expression with the argument spaces `()`, and `arguments`, a `tuple` of entries, each a [function](#function) or a `tuple` of at least two functions                              |
| Requires                      | Every function in `arguments` distinct                                                                                                                                                             |
| Static type                   | That of the expression, with the argument spaces the space of each function entry, or the `tuple` of the spaces of the members of a `tuple` entry, in order                                        |
| Raises                        | `TypeError` for an entry that is not a `Function` or a `tuple` of functions, and `ValueError` for a repeated function, a `tuple` of fewer than two functions or an expression with argument spaces |

A form binds the arguments, as a `ComponentTensor` binds indices. The slot of an
entry is the position in `arguments`. The expression is linear or antilinear in
each slot, as the expression makes it, such as antilinear under a `Conjugate`.
For a `tuple` entry, the expression is linear in the members jointly, so a term
may depend on one member only. Construction does not check the linearity. Every
other function in the expression is fixed, so one function may be an argument
of one form and fixed in another.

A slot of a `tuple` of functions gives the forms of mixed problems with no mixed
space. Each member is a function in an ordinary space, so the members may have
different degrees, such as a Raviart–Thomas and a discontinuous function, and
different domains, such as a Lagrange multiplier on the boundary. A library may
assemble the slots as one system or in blocks. The block of a member of each
slot is the form of the expression with the other members replaced by `Zero`.

The expression may have any value shape and support. A form whose expression
is a field is a linear operator, such as the interpolation operator
`form(interpolate(u, W), (u,))` of the [interpolation
extension](#example-interpolation). The support is that of the expression, so
a form of integrals has no support, although the arguments have one. A substitution renames a bound argument to avoid capture, as for indices (see
[Passes](#passes)).

> Legacy note: In legacy UFL, the arguments of a form are the test functions
> and the trial functions in the integrand, with the slot numbers of the
> functions, and a form is not an expression. In UFLx, a function has no slot
> number, the form states the arguments, and a functional is any expression
> with the argument spaces `()`. The legacy UFL `MixedElement` and `split` are
> a slot of a `tuple` of functions.

Examples, with `u`, `v` and a fixed `w` in a space `V` on `D`, and an
expression `a` linear in `u` and `v`:

- `form(a, (v, u))` is a bilinear form with the argument spaces `(V, V)`.
- `form(a, (v,))` is a linear form in `v`, such as a residual.
- `form(a, (v, u)) + form(a, (v,))` raises `ValueError`.
- With a velocity space `V` and a pressure space `Q`, functions `u`, `v` in `V`
  and `p`, `q` in `Q`, and the Stokes expression `e` of the terms
  $\nabla u : \nabla v$, $-p \, \mathrm{div}\, v$ and $-q \, \mathrm{div}\, u$,
  `form(e, ((v, q), (u, p)))` is a bilinear form with the argument spaces
  `((V, Q), (V, Q))`. `form(e, (v, q, u, p))` states a wrong form of arity
  4, because no term of `e` is linear in all four functions.

### Action

Requires: [Form](#form).

| `Action(f, slot, g)` |                                                                                                                                                                                       |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Free function        | `action(f, slot, g)`                                                                                                                                                                  |
| Properties           | `form`, an expression with more than `slot` argument spaces, `slot`, and `argument`, the expression $g$, or a `tuple` of expressions                                                  |
| Requires             | $g$ has the static type of a function in the space of the slot, with the support of that function or no support, or for a `tuple` of spaces, a `tuple` with one such member per space |
| Static type          | That of $f$, with the argument spaces of $f$ without the slot.                                                                                                                        |
| Raises               | `ValueError` if a requirement fails                                                                                                                                                   |
| Evaluated by         | [`apply_actions`](#applying-actions)                                                                                                                                                  |

`Action` is lazy, and an action on a leaf, such as a function in a dual
space, remains for the library of the leaf.

Legacy UFL operations are actions and rebindings:

- The action on the last argument of a bilinear form is `action(f, 1, w)`.
- The adjoint of `form(e, (v, u))` is `form(conjugate(e), (u, v))`.
- A derivative with respect to a function is expected to add an argument.

Examples, with `a`, `u`, `v` and `w` of the [Form](#form) examples and
`l = function(dual(V), "l")`:

- `action(form(a, (v, u)), 1, w)` is a linear form in `v`, which a pass
  evaluates to `form(b, (v,))`, with `b` the expression `a` with `u` replaced
  by `w`.
- `action(action(form(a, (v, u)), 1, w), 0, w)` is a functional.
- `action(l, 0, w)` is the pairing of `l` with `w`, with no support.
- `action(form(a, (v, u)) + 2 * form(a, (v, u)), 1, w)` moves into both
  operands of the sum, and into the form of the product.
- With the Stokes form of the [Form](#form) examples and fixed functions
  `w` in `V` and `r` in `Q`, `action(form(e, ((v, q), (u, p))), 1, (w, r))`
  is a linear form in `(v, q)`, the residual at `(w, r)`.

### Example: the Laplace-Beltrami operator

Requires: [A tour: the Poisson problem](#a-tour-the-poisson-problem),
[Components](#components), [Integral over a decomposed
domain](#integral-over-a-decomposed-domain), [Form](#form).

See also: [Expanding pull-backs](#expanding-pull-backs), [Inverse](#inverse).

The Poisson forms of the [tour](#a-tour-the-poisson-problem) on a surface in
$R^3$ are the forms of the Laplace–Beltrami operator. Only the geometry
changes. The vertex positions `y` have value shape `(m, 3)`, and the ambient
coordinate domain is `RD(3)`:

```python
# A triangulated surface in RD(3), with K, T, c and p of the tour.
x = component_tensor(
    p[c, 0, k] + (p[c, 1, k] - p[c, 0, k]) * X[0] + (p[c, 2, k] - p[c, 0, k]) * X[1],
    (k,),
)
M = decomposed_domain((parameterised_domain(K, x, RD(3)),))

# Functions in a space of a library on M.
u, v, f = function(V, "u"), function(V, "v"), function(V, "f")

# The surface gradients, of value shape (3,) in the frame of RD(3).
i = index("i")
grad_u = components(exterior_derivative(u, M), M)
grad_v = components(exterior_derivative(v, M), M)

a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), M), (v, u))
L = form(integral(f * v, M), (v,))
```

- `M` has the dimension 2 and the ambient coordinate domain `RD(3)`, so
  [`Components`](#components) projects the frame of `RD(3)` onto the tangent
  plane of each cell, and `grad_u` is the surface gradient of `u`.
- `integral(..., M)` multiplies the integrand on `K` by the density
  $\rho = \sqrt{\det(J^T J)}$, with $J$ the Jacobian matrix of `x`, of value
  shape `(3, 2)`. The integrands are 0-forms, so the forms need no orientation
  of `M`.
- `expand_pull_backs` rewrites the pull-back of `grad_u` by the cells as
  $K^T \hat g$, with $K$ = `inverse(J)` the pseudo-inverse, of value shape
  `(2, 3)`, and $\hat g$ the components of the exterior derivative of the
  pull-back of `u` on `K`. The integrand of `a` on `K` is then
  $\hat g_u^T K K^T \hat g_v \, \rho$, and $K K^T = (J^T J)^{-1}$ is the
  inverse of the metric that `x` induces.
- With a function `w` in a space on `M` of the value shape `(3,)`,
  `index_sum(components(exterior_derivative(w, M), M)[i, i], (i,))` is the
  surface divergence of `w`.
- With `P = components(exterior_derivative(coordinate(M), M), M)`, the
  projection onto the tangent plane,
  `components(exterior_derivative(grad_u, M), M)` is the tangential derivative $H$ of the surface gradient, which is
  not symmetric in general. The covariant Hessian of `u` is $P H P$.

> Legacy note: Legacy UFL writes `inner(grad(u), grad(v)) * dx` on a mesh of
> topological dimension 2 in $R^3$, and `grad` is the surface gradient because
> of the dimensions of the mesh. In UFLx, the same rule of `Components` gives
> the gradient in every dimension.

### Example: the flux through a surface

Requires: [Example: the Laplace-Beltrami
operator](#example-the-laplace-beltrami-operator), [Concept: orientation of a
domain](#concept-orientation-of-a-domain), [Integral over a decomposed
domain](#integral-over-a-decomposed-domain).

The flux of a vector field through a surface is the integral of a 2-form, so
the flux needs an orientation of the surface. With `M`, `T` and `c` of the
[Laplace-Beltrami example](#example-the-laplace-beltrami-operator), where `M` is
a closed surface such as a triangulated sphere:

```python
# The 2-form of the flux of the position field through a surface in RD(3).
Y = coordinate(RD(3))
dY = exterior_derivative(Y, RD(3))
omega = Y[0] * dY[1] * dY[2] + Y[1] * dY[2] * dY[0] + Y[2] * dY[0] * dY[1]

# The sign of each triangle from the mesh library, +1 where the order of the
# vertices of the triangle agrees with the outward orientation of M.
N = T.num_cells[0]
o = runtime_array("o", (N,))

# The outward flux, 3 times the volume inside M.
flux = integral(o[c] * omega, M)

# The area of M, with no orientation.
area = integral(Integer(1), M)
```

- `omega` has the support `RD(3)`, so `integral(omega, M)` pulls `omega` back
  by the triangles of `M` (see [Integral over a decomposed
  domain](#integral-over-a-decomposed-domain)). Each triangle contributes the
  flux with the orientation of the parameterisation of the triangle, which the
  order of the vertices of the triangle in `T` fixes. If the mesh library
  numbers some triangles clockwise and others anticlockwise, seen from
  outside, the contributions of the two groups have opposite signs, and the
  sum is not the flux.
- `o[c]` has the free index `c` of the triangles, so each triangle takes the
  sign of that triangle. With the signs, every contribution uses the outward
  orientation. A mesh library may instead number the vertices of every
  triangle consistently, and then `integral(omega, M)` is the flux without
  `o`.
- `area` integrates a 0-form against the density $\sqrt{\det(J^T J)}$, which
  does not depend on the order of the vertices, so `area` needs no signs.
- A surface with no consistent orientation, such as a Möbius strip, has no
  signs `o`. The area of the surface is defined, and the flux is not.

`runtime_extent` and `runtime_array` are not in the core (see [Example:
integrating over a cut cell](#example-integrating-over-a-cut-cell)). The same
signs enter a term that turns a 2-form on `M` into a 0-form, such as the
scalar curl of a 1-form on `M`.

## Core transformations

Requires: [Passes](#passes), [PullBack](#pullback),
[ExteriorDerivative](#exteriorderivative), [Components](#components),
[Action](#action).

The core has four transformations, each a pass (see [Passes](#passes)) that
evaluates one lazy node class. Each has a default for the node classes of the
operand, and a table of the operand classes with a rule of their own. The
[pass table](#appendix-pass-table) gives every node class for every core
transformation.

| Core transformation                             | Evaluates            | Default                                                     |
| ----------------------------------------------- | -------------------- | ----------------------------------------------------------- |
| [`expand_pull_backs(f)`](#expanding-pull-backs) | `PullBack`           | Commutes with every node class                              |
| [`apply_derivatives(f)`](#applying-derivatives) | `ExteriorDerivative` | Linear for the linear node classes, the chain rule for `()` |
| [`expand_components(f)`](#expanding-components) | `Components`         | Linear for the linear node classes                          |
| [`apply_actions(f)`](#applying-actions)         | `Action`             | Moves into the operand that carries the slot                |

The linear node classes are `Sum`, `Tensor`, `Indexed`, `ComponentTensor` and
`IndexSum`, linear in every operand, and `Division`, linear in the first
operand. An evaluated node that is linear moves into the operands of a linear
node class in which the node is linear.

A rule rewrites a node of the evaluated class to a mathematically equal
expression of node classes of the core, and never simplifies. The rewritten
expression has the free indices of the node, with the same extents, possibly
in another order, and the argument spaces of the node. The rule dispatches on
the class of the operand and moves the evaluated node towards the leaves. The
pass rewrites inner nodes first, until no rule applies. A node without a rule
is rebuilt on the rewritten edges, and the evaluated node remains on it.

In a `Form`, an evaluated node moves only into the expression, and in an
`Action` only into the form $f$, because the arguments and $g$ are in a slot,
not on the support. For example,
$\psi^* \mathrm{action}(f, s, g) = \mathrm{action}(\psi^* f, s, g)$,
because the action replaces $\psi^* u$ in $\psi^* f$ by $\psi^* g$.

Each core transformation rewrites only the evaluated class, so the four run
separately or together. Evaluation runs the four together until no rule
applies. A `PullBack` moves inside an `ExteriorDerivative`, and never the
reverse, and an action removes an `Action`, so evaluation terminates. A
downstream lazy node class may join evaluation if the rules of the class
rewrite only that class, move the node towards the leaves, and never reverse a
core rule, as moving an `ExteriorDerivative` inside a `PullBack` would.

A terminal is a `Coordinate`, a [`Function`](#function), or an expression with a
support whose class is a leaf class of a downstream library. After evaluation,
every `PullBack`, `ExteriorDerivative`, `Components` and `Action` applies to an
expression that contains a downstream terminal, except the coordinate frame
$\mathrm{d}X$, the components of $\mathrm{d}X$ on a decomposed domain of lower
dimension than the ambient coordinate domain, and a `PullBack` of a constant
by a family of domains. An
expression without a downstream terminal evaluates to a sum of products and
tensor products of 0-forms and the coordinate 1-forms $\mathrm{d}X_i$, inside
any `Form`.

### Expanding pull-backs

Requires: [Core transformations](#core-transformations), [PullBack](#pullback),
[Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative),
[Components](#components), [Inverse](#inverse).

For a domain $D$ with the parameterisation $\psi$, the parameter domain $P$
and the ambient coordinate domain $Q$, $\psi^*$ commutes with every node class
of the core, with these rules:

| Operand                                                                   | Rewrite of $\psi^*$ applied to the operand                                           |
| ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| An expression $c$ with no support                                         | $c$, if every free index of $D$ is free in $c$, and otherwise remains                |
| `Coordinate(Q)`                                                           | $\psi$                                                                               |
| `ExteriorDerivative(f, Q)`                                                | $\mathrm{d}(\psi^* f)$, with the exterior derivative on $P$                          |
| `PullBack(g, D2)`, for a `ParameterisedDomain` `D2`                       | The pull-back of $g$ by the image of $D$ under `D2`                                  |
| `PullBack(f, SideMap(D2, E, s))`, for `cell_domain(E, i)`                 | `pull_back(pull_back(f, B_i), A_i)` (see [Restriction](#restriction))                |
| `Normal(D2, E, s)`, for `cell_domain(E, i)`                               | $K^T \hat m / \lvert K^T \hat m \rvert$ (see [Normal](#normal))                      |
| `PullBack(f, SideMap(D2, E, s))` or `Normal(D2, E, s)`, for any other $D$ | Remains                                                                              |
| `Components(f, Q)`, $P$ of the dimension of $Q$                           | $\sum_i K_{ij} \, \mathrm{Components}(\psi^* f, P)_i$ for entry $j$ of the last axis |
| `Components(f, Q)`, $P$ of lower dimension than $Q$                       | Remains                                                                              |
| `Function` or a downstream terminal                                       | Remains, for the library of the terminal                                             |

The rule for an expression with no support comes before every other rule and
before commuting, whatever the class of the expression, so a
pull-back never moves inside an `Integral`. The integral has no support, and
the integrand stays on the coordinate domain of the integral. A pull-back of a
constant by a family remains, because `PullBack(c, D)` has the free indices of
the family $D$ and $c$ does not. For example, with the family of quadrature
points `xs` of [Example: quadrature](#example-quadrature),
`index_sum(integral(pull_back(Integer(1), xs), RD(0)), (i,))` counts the
points.

The image of $D$ under `D2` is
`parameterised_domain(P, pull_back(D2.parameterisation, D), D2.ambient_domain)`
(see [PullBack](#pullback)). So every `PullBack` that remains applies to `g`
directly, such as a function on a physical cell pulled back to a facet of the
reference cell.

In the rule for `Components`, $J$ is
`components(exterior_derivative(D.parameterisation, P), P)`, of value shape
`(n, k)` with $k$ the dimension of $P$ and $n$ the dimension of the frame of
`Components(f, Q)`, and $K$ is `inverse(J)`, the pseudo-inverse for $k < n$
(see [Inverse](#inverse)). The vector fields $\sum_i K_{ij} \, \partial_i$ on $P$
push forward to the orthogonal projections of the $\partial_j$ of the frame
onto the tangent space of $D$, which are the $\partial_j$ themselves for
$k = n$. For a decomposed $Q$ of dimension $k < n$, `Components(f, Q)` is
defined with those projections (see [Components](#components)), so the rule
holds for $f$ on $Q$. If $P$ has lower dimension than $Q$, the same
composition gives the components of the part of $f$ tangent to $D$, not the
components of $f$, so the pull-back remains.

### Applying derivatives

Requires: [Core transformations](#core-transformations),
[ExteriorDerivative](#exteriorderivative), [Product](#product),
[TensorProduct](#tensorproduct), [Zero](#zero), [Conjugate](#conjugate).

The exterior derivative is linear for the linear node classes and follows the
chain rule for the node classes of the degree `()`, with these rules:

| Operand                                                                     | Rewrite of $\mathrm{d}$ applied to the operand                                                                     |
| --------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `Product` $\alpha \wedge \beta$, with $\alpha$ a $k$-form                   | $\mathrm{d}\alpha \wedge \beta + (-1)^k \alpha \wedge \mathrm{d}\beta$                                             |
| `TensorProduct` $\alpha \otimes \beta$                                      | The rule of `Product`, written with `Product`                                                                      |
| `Division` $a / b$                                                          | $\mathrm{d}a / b - (\mathrm{d}b \wedge a) / (b b)$                                                                 |
| `Inverse` $A^+$                                                             | $-A^+ \, \mathrm{d}A \, A^+ + A^+ (A^+)^T \, \mathrm{d}A^T \, (I - A A^+)$                                         |
| `Absolute` $\lvert a \rvert$                                                | $(\bar{a} \, \mathrm{d}a + a \, \overline{\mathrm{d}a}) / (2 \lvert a \rvert)$                                     |
| `Conjugate` $\bar{\omega}$                                                  | $\overline{\mathrm{d}\omega}$                                                                                      |
| `ExteriorDerivative`, an expression with no support, or a `PullBack` of one | `Zero`, with the value shape, the free indices and the argument spaces of the operand and the degree raised by one |
| `Coordinate(P)`, `Normal`, `Function` or a downstream terminal              | Remains                                                                                                            |
| `Components`, or a `PullBack` of an expression with a support               | Remains                                                                                                            |

`Components` remains because the exterior derivative and the interior product
do not commute. A `PullBack` with a support remains because moving the
exterior derivative inside the pull-back would reverse a rule of
[`expand_pull_backs`](#expanding-pull-backs).

A `TensorProduct` under an exterior derivative has at most one slot, so
$\alpha \otimes \beta = \alpha \wedge \beta$. Written with `TensorProduct`,
the two terms would have different numbers of slots. In the rule for
`Inverse`, the matrix products are in the value shape, with the wedge product
of the 0-forms $A^+$ and $(A^+)^T$ with the 1-forms $\mathrm{d}A$ and
$\mathrm{d}A^T$, and $I$ is the identity matrix of value shape `(n, n)`, the
`Tensor` of `Integer(1)` and `Integer(0)` entries. For $k = n$,
$I - A A^+$ is zero and the rule is $-A^{-1} \, \mathrm{d}A \, A^{-1}$, but
the rule does not drop the term, because the core never simplifies. `Absolute` is not
complex differentiable, and for a real $a$ the rule is the chain rule. The
rule for `Conjugate` holds because the exterior derivative is real. The
exterior derivative of `Coordinate(P)` is the coordinate frame $\mathrm{d}X$.

### Expanding components

Requires: [Core transformations](#core-transformations),
[Components](#components), [Product](#product), [TensorProduct](#tensorproduct),
[Zero](#zero), [Conjugate](#conjugate).

The interior product $\iota_j$ with $\partial_j$ is linear for the linear node
classes, and the new axis stays the last axis of the value shape, with these
rules:

| Operand                                                                                                                              | Rewrite of $\iota_j$ applied to the operand                                                                |
| ------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------- |
| `Product` $\alpha \wedge \beta$, with the first slot of $\alpha$ of the degree $k$                                                   | $\iota_j\alpha \wedge \beta + (-1)^k \alpha \wedge \iota_j\beta$, and $f \, \iota_j\beta$ for a 0-form $f$ |
| `TensorProduct` $\alpha \otimes \beta$                                                                                               | $(\iota_j\alpha) \otimes \beta$ if $\alpha$ has a slot, and $f \, \iota_j\beta$ for a 0-form $f$           |
| $\mathrm{d}X_i$, on a support of the dimension of the frame                                                                          | $\delta_{ij}$                                                                                              |
| `Zero`                                                                                                                               | `Zero` with the new value shape and degree                                                                 |
| `Conjugate` $\bar{\alpha}$                                                                                                           | $\overline{\iota_j \alpha}$                                                                                |
| `Function` or a downstream terminal                                                                                                  | Remains                                                                                                    |
| $\mathrm{d}X_i$ on a decomposed domain of lower dimension than the frame, any other `ExteriorDerivative`, `Components` or `PullBack` | Remains                                                                                                    |

The rule for $\mathrm{d}X_i$ applies to `ExteriorDerivative(Coordinate(P), P)`,
reached through `Indexed`. On a decomposed domain of lower dimension than the
frame, the components of $\mathrm{d}X$ are the orthogonal projection onto the
tangent space, which `expand_pull_backs` evaluates on each piece as $J K$.
The components of the exterior derivative of any
other expression remain, such as the gradient
`components(exterior_derivative(u, P), P)` of a function `u`.

The components of `dX` are the identity matrix, the `Tensor` of `Integer(1)`
and `Integer(0)` entries, or `Zero((0, 0), (), (), ())` for $n = 0$. The rule
for `Conjugate` holds because $\partial_j$ is real.

### Applying actions

Requires: [Core transformations](#core-transformations), [Action](#action),
[Form](#form), [Passes](#passes).

An action moves into the operand that carries the slot, with these rules:

| Operand                                                  | Rewrite of `Action(f, slot, g)` with $f$ the operand                                                                           |
| -------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `Form(e, arguments)`                                     | $e$ with the entry of the slot replaced by $g$, member by member for a `tuple`, in a `Form` of the other entries if any remain |
| `Sum`, `Tensor`                                          | The action of every operand                                                                                                    |
| `Product`, `TensorProduct`                               | The action of the operand that carries the slot, with the slot counted within that operand                                     |
| `Zero`                                                   | `Zero` with the argument spaces without the slot                                                                               |
| `Function` or a downstream leaf, such as in a dual space | Remains, for the library of the leaf                                                                                           |
| `Action`                                                 | Remains                                                                                                                        |

The substitution into a `Form` avoids capture (see [Passes](#passes)). The
pass rewrites inner actions first, so an outer action meets an `Action` only on
a leaf, and the outer action remains.

## Extending the core

Requires: [Passes](#passes), [Scalar type](#scalar-type), [Absolute](#absolute).

Extensibility is a requirement of the core. Downstream libraries must be able
to extend UFLx without changing the core.

> Legacy note: In legacy UFL, properties are stored on the node classes, and
> the set of properties is fixed. Every new property touches every node class,
> so many language features were hacked in, for example complex numbers and
> dual spaces.

The core has a small, fixed set of node classes and an open set of passes. A
downstream library extends UFLx mainly by adding passes, not by adding
properties to the node classes of the core. A new pass needs one rule for each
node class of the core and does not change the core. Only a language extension
adds node classes.

The core supports five kinds of extension:

> LLM generated table, could be made more concise.

| Kind of extension                    | Examples                                                    | Section                                     |
| ------------------------------------ | ----------------------------------------------------------- | ------------------------------------------- |
| New questions asked of an expression | Physical units, polynomial degree, linearity checks         | [New questions](#new-questions)             |
| New subclasses                       | Elements and functions from Basix, new coordinate domains   | [New subclasses](#new-subclasses)           |
| New operations                       | Trace, divergence                                           | [New operations](#new-operations)           |
| New transformations                  | Gateaux derivative                                          | [New transformations](#new-transformations) |
| Language extensions                  | External operators, [interpolation](#example-interpolation) | [Language extensions](#language-extensions) |

A property that the rules of the core depend on is part of the core, not an
extension. For example, the scalar type is part of the core, because
`Absolute` is real for a complex operand. Adding such a property changes the
core.

### Language extensions

Requires: [Leaves and subclassing](#leaves-and-subclassing), [Passes](#passes),
[Core transformations](#core-transformations), [Transpose and matrix
product](#transpose-and-matrix-product).

A language extension adds node classes that are not leaves, in a package
outside the core. A node class of a language extension subclasses
`AbstractExpression`, never a concrete class of the core (see [Leaves and
subclassing](#leaves-and-subclassing)). Each node class of a language
extension provides, as part of the class:

- a rule for each [core transformation](#core-transformations) that keeps
  the node, which produces only node classes of the core or of the language
  extension, and
- a lowering to an equal expression of node classes of the core, if such an
  expression exists.

A pass without a rule for a node class of a language extension lowers the node
and continues, so a pass such as `units` needs no rule for the node class. For
a node class without a lowering, every pass without a rule raises
`TypeError`. A rule must agree with the lowering followed by the rules for the
lowered expression, up to mathematical equality.

UFLx core is likely to include its own language extensions, such as UFL
compatibility layers and a simpler layer based on classical notions of
differential geometry compatible with most engineering and physics syllabus.

For example, a language extension for local external operators adds a node
class for an operation known only by name, such as `exp`, a material law or a
neural network. The rule for `apply_derivatives` differentiates the node class.
The rule for `expand_pull_backs` pulls back the operands, because the
operation is pointwise. The node class has no lowering.

A named composition that must survive a pass, or print as written, is a node
class of a language extension whose lowering is the composition. For example,
a `Transpose` node class lowers to `component_tensor(A[k, l], (l, k))`. The
rule for `expand_pull_backs` moves the pull-back inside `Transpose`, because
the transpose is pointwise, so `Transpose` survives the pull-back.

### Example: interpolation

Requires: [Language extensions](#language-extensions),
[AbstractFunctionSpace](#abstractfunctionspace), [DualSpace](#dualspace),
[Form](#form), [Action](#action).

The interpolation extension ships with UFLx and adds the node class
`Interpolate` for the interpolant of an expression in a space.

| `Interpolate(f, V)` |                                                                                                                                      |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Free function       | `interpolate(f, V)`                                                                                                                  |
| Properties          | `expression` $f$ and `space` $V$, an `AbstractFunctionSpace` with the argument spaces `()`                                           |
| Requires            | $f$ with the value shape and the degree of $V$, no free indices, the support `V.domain` or no support, and a real $f$ for a real $V$ |
| Static type         | That of a function in $V$, with the argument spaces of $f$                                                                           |
| Raises              | `ValueError` if a requirement fails                                                                                                  |
| Lowering            | None                                                                                                                                 |

`interpolate(f, V)` is the function in $V$ with the degrees of freedom of $f$.
The core does not represent the degrees of freedom of $V$, so the library of
$V$ defines the interpolant, as the library defines the pull-back of a
function in $V$. The interpolant is a composition of the core only after a
pull-back by a cell, so the node class has no lowering. The rules for the core
transformations are:

| Core transformation | Rule for `Interpolate(f, V)`                                                                                                                                                           |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `expand_pull_backs` | The pull-back by `cell_domain(V.domain, j)` is the composition that the library of $V$ gives from the pull-back of $f$ by the cell, such as values at points or moments as `Integral`s |
| `apply_derivatives` | The node remains under the exterior derivative                                                                                                                                         |
| `expand_components` | The node remains under `Components`                                                                                                                                                    |
| `apply_actions`     | `interpolate` of the action of $f$                                                                                                                                                     |

The interpolant of a derivative differs from the derivative of the
interpolant, except for some pairs of spaces, so the exterior derivative
remains on the node, as on a function. `Interpolate` is linear in $f$, so the
node passes on the argument spaces of $f$, and an action moves into $f$.

The interpolation operator from a space $U$ to $V$ is the form
`form(interpolate(u, V), (u,))`, with a function `u` in $U$, a form whose
value is a field. The adjoint of the operator takes a function `l` in
`dual(V)` to the form `form(action(l, 0, interpolate(u, V)), (u,))`,
antilinear in `u`, with the static type of a function in `dual(U)`. Binding
`l` as well gives the form of the interpolation matrix. The action on `l`
remains, and the library of `l` evaluates the pairing, the transpose of the
interpolation. So the extension needs no class for the adjoint, because the
core provides [dual spaces](#dualspace), forms and actions.

Interpolation between the spaces of two domains needs the cell of one domain
that contains a point of the other, and is left to a later version of the
extension.

> Legacy note: In legacy UFL, `Interpolate(f, V)` is a `BaseFormOperator`, as
> an external operator is, and the adjoint interpolation takes a `Coargument`.
> In DOLFINx, interpolation happens outside UFL.

Examples, with the decomposed domain `D`, the Lagrange spaces `V` and `V2` of
polynomial degrees 1 and 2 on `D`, a Nédélec space `W` on `D` of the degree
`(1,)`, a function `u` in `V` and a function `l` in `dual(V2)`:

- `interpolate(u * u, V)` has the static type of a function in `V`.
- `interpolate(exterior_derivative(u, D), W)` is the interpolant of
  $\mathrm{d}u$ in `W`.
- `interpolate(u, W)` raises `ValueError`, because the degree of `u` is `()`
  and the degree of `W` is `(1,)`.
- `form(interpolate(u, V2), (u,))` is the prolongation from `V` to `V2`, and
  `form(action(l, 0, interpolate(u, V2)), (u,))` is the restriction of `l` to
  a linear form on `V`, the transpose of the prolongation, as in multigrid.

### New questions

Requires: [Passes](#passes), [Language extensions](#language-extensions),
[Coordinate](#coordinate), [ExteriorDerivative](#exteriorderivative),
[Components](#components).

A downstream library asks a new question of an expression with a pass. The pass
has one rule for each node class of the core and one rule for each leaf class
of the library. The pass needs no rule for a node class of a language extension
that has a lowering (see [Language extensions](#language-extensions)).

For example, a library that checks physical units, like the work of Habera et
al.:

1. Assigns a unit to leaves in a table of the library. For example, with
   `X = coordinate(RD(2))`, the table `{X: "m"}` states that the coordinates
   are in metres. A leaf that is not in the table, such as a literal, is
   dimensionless.
2. Writes passes, such as `units(expression, table)`, which returns the unit of
   any expression. `units` derives the unit of an expression from the units of
   the leaves, in the same way as the core derives the value shape.
   `X[0] * X[1]` is in square metres. The gradient
   `components(exterior_derivative(X[0], RD(2)), RD(2))` is in metres divided
   by metres, so the gradient is dimensionless.
   `X[0] + X[0] * X[1]` is an error of the units library.
3. Calls `units` on an expression whenever the library needs the check.

The core never knows about units. Objects of the core are immutable and are
compared by content, so the leaves are reliable keys of the table.

An error reported by `units` does not make an object invalid. An object is
valid if the constructors of the core accept the object.

The pass never guesses a rule. A `TypeError` from `units` means that the units
library does not know a node class. An error of the units library means that the
expression has inconsistent units.

### New operations

Requires: [Node class or free function](#node-class-or-free-function), [Tensor
algebra with EllipsisIndex](#tensor-algebra-with-ellipsisindex), [Language
extensions](#language-extensions).

A downstream library builds a new operation from existing node classes, as a
free function. As in the core, the free function may read the static type of
the operands, such as a value shape, to choose the composition (see
[Node class or free function](#node-class-or-free-function)). For
example, `determinant(A)` reads the extent $n$ of `A` and returns the Leibniz
formula for an $n \times n$ matrix. So `determinant(A)` requires an `int`
extent $n$ (see [Extents](#extents)). `trace(A)` reads no extent, because
`index_sum(A[k, k], (k,))`, with a derived index `k`, is the trace for every
extent. An operation that a
`EllipsisIndex` expresses, such as `scale(a, T)`, needs no value shape (see
[EllipsisIndex](#tensor-algebra-with-ellipsisindex)).

A library node class that a pass replaces by node classes of the core gives
nothing that the free function does not give, unless a core transformation
must keep the structure of the operation. Such an operation, and an operation
that is not a composition of node classes of the core, needs a language
extension.

### New subclasses

Requires: [Leaves and subclassing](#leaves-and-subclassing), [New
operations](#new-operations), [Language extensions](#language-extensions),
[Coordinate](#coordinate), [Example: an affinely mapped
triangle](#example-an-affinely-mapped-triangle).

Every abstract class is a contract. The contract states what a subclass
provides, so that every core transformation handles the subclass through the
abstract class alone. A downstream library may subclass an
abstract class by fulfilling the contract, without changing the core. Only
leaves are subclassed (see
[Leaves and subclassing](#leaves-and-subclassing)). An
operation built from other expressions is not a subclass but a free function
(see [New operations](#new-operations)) or a node class of a language extension
(see [Language extensions](#language-extensions)). A contract is
part of the class, never a global registry. For example, a
downstream subclass of `AbstractCoordinateDomain`, such as
[`UnitSimplex`](#example-the-unit-simplex), provides the
dimension, and every core transformation handles the coordinate function of
the subclass through `Coordinate`.

Information that a downstream library attaches to the graph goes in one of two
places:

- Information that defines the variational problem goes in the construction
  arguments of a leaf. For example, the space that a function belongs to
  changes the meaning of an expression, so the space is part of the leaf.
- Information that only checks or analyses the variational problem may stay in
  a table of the library, keyed by leaves. For example, units check an
  expression but do not change the meaning of the expression.

### New transformations

Requires: [Language extensions](#language-extensions).

A new transformation is a pass with one rule for each node class of the core
and one rule for each leaf class that the transformation rewrites. For example,
a Gateaux derivative rewrites the functions that the derivative is taken with
respect to, so the Gateaux derivative needs a rule for the leaf classes of
those functions. A new transformation must produce only node classes of the
core. A transformation that needs a new node class needs a language extension.

There should be helpers to check whether a transformation covers all classes
of the core.

## Appendix: Class hierarchy

Every object in the graph extends `Base`. The classes defined in this document
form the following hierarchy.

> LLM generated graphs.

```
# Coordinate domains
Base
└── AbstractCoordinateDomain              parameter space, stores n
    └── RD                                R^n

# Expressions and static type
Base
├── AbstractExpression
├── Index                                 named by a label, not an expression
├── EllipsisIndex                         covers any number of axes
├── AbstractExtent                        extent not known during construction
├── AbstractIndexValue                    index value or table not known during construction
└── GatheredIndex                         index value read from a table

# Algebra
AbstractExpression
├── Integer                               literal
├── Float                                 literal
├── ImaginaryUnit                         literal
├── Zero                                  zero of any shape and degree, no support
├── Sum                                   +
├── Product                               *, the wedge product
├── Power                                 **
├── Division                              /
├── Absolute                              abs
├── Conjugate
├── TensorProduct                         ⊗, concatenates slots
├── SquareRoot
├── Inverse                               inverse, or pseudo-inverse of a tall matrix
├── Tensor                                stacks entries along a new axis
├── Indexed                               full indexing, int, Index or index value entries
├── ComponentTensor                       binds indices to tensor axes
└── IndexSum                              sums over bound indices

# Calculus on a coordinate domain
AbstractExpression
├── Coordinate                            terminal, X on P
├── ExteriorDerivative                    d with respect to P
├── Components                            first slot to an axis, with the frame of P
└── Integral                              integral of a top form over P, no support

# Maps between domains
Base
└── ParameterisedDomain                   stores P, ψ and Q
AbstractExpression
└── PullBack                              pull-back of a differential form by a domain, support P of the domain

# Topology
Base
├── AbstractCellType                      numbered vertices
└── AbstractTopology                      cells without coordinates

# Reference cells
AbstractCoordinateDomain, AbstractCellType
└── AbstractReferenceCell                 coordinates of the vertices

# Decomposed domains
Base
├── DecomposedDomain                      union of families of parameterised domains
└── SideMap                               the cells of a topology with sides into the cells on one side
AbstractExpression
└── Normal                                unit outward normal of the cell on a side

# Function spaces
Base
└── AbstractFunctionSpace                 degree and value shape of functions
    └── DualSpace                         antilinear forms on a space or a product
AbstractExpression
└── Function                              a member of a space, a terminal

# Variational forms
AbstractExpression
├── Form                                  binds functions as arguments
└── Action                                gives an expression to a slot of a form
```

## Appendix: Free functions

> LLM generated table.

| Group                           | Free function                                                   | `AbstractExpression` overload | Takes                                                                                                       | Returns                                                                     |
| ------------------------------- | --------------------------------------------------------------- | ----------------------------- | ----------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| Expressions and static type     | [`index(label)`](#index)                                        |                               | `label: str \| tuple`                                                                                       | `Index`                                                                     |
|                                 | [`ellipsis_index(label)`](#ellipsisindex)                       |                               | `label: str \| tuple`                                                                                       | `EllipsisIndex`                                                             |
|                                 | [`gathered_index(table, indices)`](#gatheredindex)              |                               | `table: AbstractIndexValue`, `indices: tuple`                                                               | `GatheredIndex`                                                             |
| Algebra                         | [`add(f, g)`](#sum)                                             | `+`                           | `f, g: AbstractExpression`                                                                                  | `Sum`                                                                       |
|                                 | [`multiply(f, g)`](#product)                                    | `*`                           | `f, g: AbstractExpression`                                                                                  | `Product`                                                                   |
|                                 | [`power(f, g)`](#power)                                         | `**`                          | `f, g: AbstractExpression`                                                                                  | `Power`                                                                     |
|                                 | [`negate(f)`](#negation-and-subtraction)                        | unary `-`                     | `f: AbstractExpression`                                                                                     | `Product` of `Integer(-1)` and `f`                                          |
|                                 | [`subtract(f, g)`](#negation-and-subtraction)                   | `-`                           | `f, g: AbstractExpression`                                                                                  | `Sum` of `f` and `negate(g)`                                                |
|                                 | [`divide(f, g)`](#division)                                     | `/`                           | `f, g: AbstractExpression`                                                                                  | `Division`                                                                  |
|                                 | [`absolute(f)`](#absolute)                                      | `abs`                         | `f: AbstractExpression`                                                                                     | `Absolute`                                                                  |
|                                 | [`conjugate(f)`](#conjugate)                                    |                               | `f: AbstractExpression`                                                                                     | `Conjugate`                                                                 |
|                                 | [`tensor_product(f, g)`](#tensorproduct)                        |                               | `f, g: AbstractExpression`                                                                                  | `TensorProduct`                                                             |
|                                 | [`square_root(f)`](#elementary-functions-very-incomplete)       |                               | `f: AbstractExpression`                                                                                     | `SquareRoot`                                                                |
|                                 | [`tensor(entries)`](#tensor)                                    |                               | `entries: tuple[AbstractExpression, ...]`                                                                   | `Tensor`                                                                    |
|                                 | [`indexed(f, indices)`](#indexed)                               | `f[i]`                        | `f: AbstractExpression`, `indices: tuple`                                                                   | `Indexed`                                                                   |
|                                 | [`component_tensor(f, indices)`](#componenttensor-and-indexsum) |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`                                                       | `ComponentTensor`                                                           |
|                                 | [`index_sum(f, indices)`](#componenttensor-and-indexsum)        |                               | `f: AbstractExpression`, `indices: tuple[Index, ...]`                                                       | `IndexSum`                                                                  |
|                                 | [`transpose(f)`](#transpose-and-matrix-product)                 |                               | `f: AbstractExpression`                                                                                     | `ComponentTensor` of an `Indexed`                                           |
|                                 | [`matrix_product(f, g)`](#transpose-and-matrix-product)         | `@`                           | `f, g: AbstractExpression`                                                                                  | `ComponentTensor` of an `IndexSum`                                          |
|                                 | [`inverse(f)`](#inverse)                                        |                               | `f: AbstractExpression`                                                                                     | `Inverse`                                                                   |
| Calculus on a coordinate domain | [`coordinate(P)`](#coordinate)                                  |                               | `P: AbstractCoordinateDomain \| DecomposedDomain`                                                           | `Coordinate`                                                                |
|                                 | [`exterior_derivative(f, P)`](#exteriorderivative)              |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain \| DecomposedDomain`                                  | `ExteriorDerivative`                                                        |
|                                 | [`components(f, P)`](#components)                               |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain \| DecomposedDomain`                                  | `Components`                                                                |
|                                 | [`volume_form(P)`](#volume-form)                                |                               | `P: AbstractCoordinateDomain`                                                                               | `Product` of the coordinate 1-forms                                         |
|                                 | [`integral(f, P)`](#integral)                                   |                               | `f: AbstractExpression`, `P: AbstractCoordinateDomain \| DecomposedDomain`                                  | `Integral`, or a sum of `IndexSum`s of `Integral`s for a `DecomposedDomain` |
| Maps between domains            | [`parameterised_domain(P, psi, Q)`](#parameteriseddomain)       |                               | `P: AbstractCoordinateDomain`, `psi: AbstractExpression`, `Q: AbstractCoordinateDomain \| DecomposedDomain` | `ParameterisedDomain`                                                       |
|                                 | [`pull_back(f, D)`](#pullback)                                  |                               | `f: AbstractExpression`, `D: ParameterisedDomain \| SideMap`                                                | `PullBack`                                                                  |
| Topology                        | [`cell_index(T, j)`](#abstracttopology)                         |                               | `T: AbstractTopology`, `j: int`                                                                             | `Index`                                                                     |
| Reference cells                 | [`affine_map(F, points)`](#reference-cell)                      |                               | `F: AbstractReferenceCell`, `points: AbstractExpression`                                                    | The affine map of the vertices                                              |
| Decomposed domains              | [`decomposed_domain(pieces)`](#decomposeddomain)                |                               | `pieces: tuple[ParameterisedDomain, ...]`                                                                   | `DecomposedDomain`                                                          |
|                                 | [`cell_domain(D, j)`](#decomposeddomain)                        |                               | `D: DecomposedDomain`, `j: int`                                                                             | `ParameterisedDomain` of the cells of type `j`                              |
|                                 | [`side_domain(D, S)`](#restriction)                             |                               | `D: DecomposedDomain`, `S: AbstractTopology`                                                                | `DecomposedDomain` of the cells of `S`                                      |
|                                 | [`side_map(D, E, side)`](#restriction)                          |                               | `D, E: DecomposedDomain`, `side: str`                                                                       | `SideMap`                                                                   |
|                                 | [`restriction(f, E, side)`](#restriction)                       |                               | `f: AbstractExpression`, `E: DecomposedDomain`, `side: str`                                                 | `PullBack` of `f` by `side_map(D, E, side)`, or `f` with no support         |
|                                 | [`normal(D, E, side)`](#normal)                                 |                               | `D, E: DecomposedDomain`, `side: str`                                                                       | `Normal`                                                                    |
| Function spaces                 | [`function(V, label)`](#function)                               |                               | `V: AbstractFunctionSpace`, `label: str \| tuple`                                                           | `Function`                                                                  |
| Variational forms               | [`form(expression, arguments)`](#form)                          |                               | `expression: AbstractExpression`, `arguments: tuple[Function \| tuple[Function, ...], ...]`                 | `Form`                                                                      |
|                                 | [`dual(S)`](#dualspace)                                         |                               | `S: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]`                                             | `DualSpace`                                                                 |
|                                 | [`action(f, slot, g)`](#action)                                 |                               | `f: AbstractExpression`, `slot: int`, `g: AbstractExpression \| tuple[AbstractExpression, ...]`             | `Action`                                                                    |

The value-shape and other constraints on the arguments are stated in the section
of each free function. A free function and an operator overload accept only the
types in the "Takes" column.

A free function whose "Returns" column names only a node class `C` returns
`C(...)` of the arguments. Any other free function returns the composition in
the "Returns" column. The sections of the node classes do not
repeat this rule.

## Appendix: Properties

A derived property listed for an abstract class is abstract, unless the table
gives the value. Each concrete class implements an abstract property by the
rules in the section of the concrete class.

> LLM generated table.

| Class                                          | Stored properties                                                                               | Derived properties                                                                   |
| ---------------------------------------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| `Base`                                         | `construction_arguments`                                                                        |                                                                                      |
| `AbstractCoordinateDomain`                     | `dimension`, in every concrete subclass such as `RD`                                            |                                                                                      |
| `AbstractExpression`                           |                                                                                                 | `value_shape`, `free_indices`, `scalar_type`, `support`, `degree`, `argument_spaces` |
| `Index`, `EllipsisIndex`                       | `label`                                                                                         |                                                                                      |
| `AbstractIndexValue`                           | `extent`, `shape`, in every concrete subclass                                                   |                                                                                      |
| `GatheredIndex`                                | `table`, `indices`                                                                              | `extent`, `free_indices`                                                             |
| `Integer`, `Float`                             | `value`                                                                                         |                                                                                      |
| `Zero`                                         | `value_shape`, `free_indices`, `degree`, `argument_spaces`                                      |                                                                                      |
| `Sum`, `Product`, `Power`, `Division`          | `first`, `second`                                                                               |                                                                                      |
| `Absolute`, `Conjugate`                        | `operand`                                                                                       |                                                                                      |
| `TensorProduct`                                | `first`, `second`                                                                               |                                                                                      |
| `SquareRoot`, `Inverse`                        | `operand`                                                                                       |                                                                                      |
| `Tensor`                                       | `entries`                                                                                       |                                                                                      |
| `Indexed`, `ComponentTensor`, `IndexSum`       | `expression`, `indices`                                                                         |                                                                                      |
| `Coordinate`                                   | `domain`                                                                                        |                                                                                      |
| `ExteriorDerivative`, `Components`, `Integral` | `expression`, `domain`                                                                          |                                                                                      |
| `ParameterisedDomain`                          | `parameter_domain`, `parameterisation`, `ambient_domain`                                        | `free_indices`                                                                       |
| `PullBack`                                     | `expression`, `domain`                                                                          |                                                                                      |
| `AbstractCellType`                             | `dimension`, `num_vertices`, in every concrete subclass                                         |                                                                                      |
| `AbstractReferenceCell`                        | `vertices`, in every concrete subclass                                                          |                                                                                      |
| `AbstractTopology`                             | `cell_types`, `num_cells`, `sides`, in every concrete subclass                                  | `dimension`                                                                          |
| `SideMap`                                      | `parameter_domain`, `ambient_domain`, `side`                                                    | `free_indices`, always `()`                                                          |
| `Normal`                                       | `ambient_domain`, `parameter_domain`, `side`                                                    |                                                                                      |
| `DecomposedDomain`                             | `pieces`                                                                                        | `ambient_domain`, `dimension`, `topology`                                            |
| `AbstractFunctionSpace`                        | `domain`, `degree`, `value_shape`, `scalar_type`, `argument_spaces`, in every concrete subclass |                                                                                      |
| `DualSpace`                                    | `space`                                                                                         | `domain`, `degree`, `value_shape`, `scalar_type`, `argument_spaces`                  |
| `Function`                                     | `space`, `label`                                                                                |                                                                                      |
| `Form`                                         | `expression`, `arguments`                                                                       | `argument_spaces`                                                                    |
| `Action`                                       | `form`, `slot`, `argument`                                                                      | `argument_spaces`                                                                    |

## Appendix: Static type table

The degree and the argument spaces of each node class of the core. An operand
not named in the argument spaces column must have the argument spaces `()`.

| Class                                    | Edges that are expressions must have                      | Degree of the node                            | Argument spaces of the node                |
| ---------------------------------------- | --------------------------------------------------------- | --------------------------------------------- | ------------------------------------------ |
| `Integer`, `Float`, `ImaginaryUnit`      | No such edges                                             | `()`                                          | `()`                                       |
| `Zero`                                   | No such edges                                             | The stored degree                             | The stored argument spaces                 |
| `Sum`                                    | Equal degrees                                             | The degree of the operands                    | Those of the operands, which must be equal |
| `Product`                                | The degree `()`, or exactly one slot each                 | The slotwise sum (see [Product](#product))    | Those of the operands, concatenated        |
| `Power`, `Absolute`                      | The degree `()`                                           | `()`                                          | `()`                                       |
| `Division`                               | The degree `()` for the second operand                    | The degree of the first operand               | Those of the first operand                 |
| `Conjugate`                              | Any degree                                                | The degree of the operand                     | Those of the operand                       |
| `TensorProduct`                          | Any degree                                                | The concatenated degrees                      | Those of the operands, concatenated        |
| `SquareRoot`, `Inverse`                  | The degree `()`                                           | `()`                                          | `()`                                       |
| `Tensor`                                 | Equal degrees                                             | The degree of the entries                     | Those of the entries, which must be equal  |
| `Indexed`, `ComponentTensor`, `IndexSum` | Any degree                                                | The degree of the expression                  | Those of the expression                    |
| `Coordinate`                             | No such edges                                             | `()`                                          | `()`                                       |
| `Normal`                                 | No such edges                                             | `()`                                          | `()`                                       |
| `Function`                               | No such edges                                             | The degree of the space                       | Those of the space                         |
| `Form`                                   | Any degree                                                | The degree of the expression                  | The space of each entry of the arguments   |
| `Action`                                 | Any degree                                                | The degree of the form                        | Those of the form without the slot         |
| `ExteriorDerivative`                     | The degree `()` or `(k,)`                                 | `(1,)` for `()`, and `(k + 1,)` for `(k,)`    | Those of the operand                       |
| `Components`                             | At least one slot                                         | The first slot lowered by 1, and removed at 0 | Those of the operand                       |
| `Integral`                               | `(n,)` with $n$ the dimension of $P$, or `()` for $n = 0$ | `()`                                          | Those of the integrand                     |
| `PullBack`                               | Any degree                                                | The degree of the expression                  | Those of the expression                    |

The value shape and the free indices of each node class are in the section
of the class. The support and the scalar type of the node classes that differ from the
default rules of [Support](#support) and [Scalar type](#scalar-type):

| Class                                            | Support of the node                | Scalar type of the node           |
| ------------------------------------------------ | ---------------------------------- | --------------------------------- |
| `ImaginaryUnit`                                  | None                               | Complex                           |
| `Zero`                                           | None                               | Real                              |
| `Absolute`                                       | The default                        | Real                              |
| `Coordinate`, `ExteriorDerivative`, `Components` | The coordinate domain $P$          | The default                       |
| `Integral`                                       | None                               | The default                       |
| `PullBack`                                       | The parameter domain of the domain | The default                       |
| `Normal`                                         | $E$                                | Real                              |
| `Function`                                       | The domain of the space, or none   | The scalar type of the space      |
| `Form`                                           | The support of the expression      | The scalar type of the expression |
| `Action`                                         | The support of the form            | The scalar type of the form       |

## Appendix: Pass table

What each core transformation does with each node class as the operand of the
evaluated node. "In" moves the evaluated node into the operands, by linearity
or commuting. "Chain" is the chain rule. "Rule" is a rule of the transformation. "Slot" moves the
evaluated node into the expression or the form only. "Remains" keeps the
evaluated node on the operand. "Rule, Remains" is a rule for some operands of
the class and "Remains" for the others, as the table of the transformation
states. "—" cannot occur, because the static type rules out the pair.

For `expand_pull_backs`, the entry is for an operand with a support. An
operand with no support takes the rule for an expression with no support,
whatever the class of the operand.

| Class                                    | `expand_pull_backs` | `apply_derivatives` | `expand_components` | `apply_actions` |
| ---------------------------------------- | ------------------- | ------------------- | ------------------- | --------------- |
| `Integer`, `Float`, `ImaginaryUnit`      | Rule                | Rule                | —                   | —               |
| `Zero`                                   | Rule                | Rule                | Rule                | Rule            |
| `Sum`, `Tensor`                          | In                  | In                  | In                  | Rule            |
| `Product`, `TensorProduct`               | In                  | Rule                | Rule                | Rule            |
| `Division`                               | In                  | Rule                | In                  | In              |
| `Power`, `SquareRoot`                    | In                  | Chain               | —                   | —               |
| `Absolute`, `Inverse`                    | In                  | Rule                | —                   | —               |
| `Conjugate`                              | In                  | Rule                | Rule                | In              |
| `Indexed`, `ComponentTensor`, `IndexSum` | In                  | In                  | In                  | In              |
| `Coordinate`                             | Rule                | Remains             | —                   | —               |
| `Normal`                                 | Rule, Remains       | Remains             | —                   | —               |
| `Function`                               | Remains             | Remains             | Remains             | Remains         |
| `Form`                                   | Slot                | Slot                | Slot                | Rule            |
| `Action`                                 | Slot                | Slot                | Slot                | Remains         |
| `ExteriorDerivative`                     | Rule                | Rule                | Rule, Remains       | In              |
| `Components`                             | Rule, Remains       | Remains             | Remains             | In              |
| `Integral`                               | Rule                | Rule                | —                   | In              |
| `PullBack`                               | Rule, Remains       | Rule, Remains       | Remains             | In              |

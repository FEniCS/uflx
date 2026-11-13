# Functions and forms

## Function spaces

A function is a member of a function space, such as a coefficient or an
argument of a variational form.

### `AbstractFunctionSpace`

`AbstractFunctionSpace` describes a space of fields or multilinear maps,
with or without a spatial argument.

| `AbstractFunctionSpace` |                                                                                                                                                       |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `Base` |
| Provides                | `domain`, an [`AbstractDomain`](topology.md#abstractdomain), an [`AbstractTopology`](topology.md#abstracttopology) or `None`, `value_shape`, `scalar_type`, and `argument_spaces`, `()` for a space of fields |
| Defined by              | A subclass outside the core, such as a subclass for the finite element spaces of Basix, and [`DualSpace`](#dualspace) in the core                           |

| Property | Meaning |
| -------- | ------- |
| `domain` | Where the space's fields are defined: an `AbstractTopology`, an `AbstractDomain`, or `None` for values with no spatial argument. |
| `value_shape` | Array shape of a function's output, with `()` for a scalar. |
| `scalar_type` | Real or complex arithmetic for the function's expression representation. |
| `argument_spaces` | Ordered input spaces when a function represents a multilinear map; `()` for a field. |

A space with `domain is None` has functions with no spatial argument.

Ordinary finite-element spaces supply independent cell-local basis coefficients.
The package establishes assembled continuity through DOF maps and transformations:
it gathers inputs and scatters contributions without merging local tensor axes. Explicitly constrained or macro-element spaces may define other semantics.
The package may bind functions in different local configurations to restrictions
of the same field; those associations are external to the core.

A library may supply `transformation_law(V)` for spaces with a
[`TransformationLaw`](geometry.md#transformationlaw). This is not a required space property
or a core overload; the core constructor remains `transformation_law(variances, weight)`.

#### Functions in a space

`Function(space, label)` takes its domain, value shape, scalar type and
argument spaces from `space`. These describe its expression
representation; a form may represent the same space element differently
(see [The double dual](#the-double-dual)). For fields on core domains and
topologies, the concrete subclass supplies evaluation as follows:

- A space whose values on the reference cell do not depend on the geometry is
  on a topology $T$. The concrete subclass defines
  `pull_back(u, cell_map(T, j))` for every local cell slot $j$, with local
  data supplied by the package and an appropriate element for each cell
  type. The field on a mesh is a [`PushForward`](geometry.md#pushforward).
- A space whose values on the reference cell depend on the geometry is on a
  mesh $D$, such as Hermite, Argyris, Morley and the other elements that are
  not Piola equivalent. The concrete subclass defines
  `pull_back(u, cell_map(D, j))`, and reads $x$ from the cell map.
- For a coordinate domain, such as a reference cell, a function is a terminal
  that the concrete subclass evaluates.

A geometry-dependent cell read must expose its coordinate dependencies or
provide a derivative rule before geometry differentiation can proceed.
Independent data and its transport are specified by the library
[reconstruction](extensions.md#geometry-dependent-reconstruction); a mesh-bound
space alone does not determine that transport.

The library evaluates these cell reads and their reference derivatives,
including at supplied reference positions. It also evaluates applications of
functional terminals and spatial transformations of their outputs. Core
lowering leaves these reads for the library; it does not compute basis or
coefficient values.

Extension-defined domains use their own evaluation rules; reference-cell
reads are not required. Geometry dependencies must still be exposed.

For evaluation on a core `Mesh`, the solver locates physical points and
supplies the selected cells and reference positions before an evaluation is
constructed. At interfaces or overlaps, the cell selection is part of defining
the value. Evaluation uses a [cell map](topology.md#cellmap) followed by
[reference-cell evaluation](geometry.md#reference-cell-evaluation).
A coordinate-only `DomainMap` into a mesh is invalid.

The concrete subclass owns the polynomials, degrees of freedom and basis of
the space, and stores the arguments identifying the space, such as the finite
element. Space identity is structural.

A library may provide [cell-local interpolation](extensions.md#cell-local-interpolation)
using the sampling rule and basis of the target space. These are library
data, not additional properties required by `AbstractFunctionSpace`.

#### Spaces of constants

A space with the domain `None` and the argument spaces `()` is a space of
constants. A function in the space is an expression with no domain, a
single array with the value shape and scalar type supplied by the space.

The integrals of a form decide the cells on which a constant acts.

A constant unknown, such as a Lagrange multiplier, needs a space because form
arguments come from spaces. Constant data needs none: a literal or kernel-input
array with no domain is already a constant.

#### Mixed problems

A mixed problem needs no mixed space, because an argument of a [form](#form)
may be a `tuple` of functions in several spaces.

#### Examples

The examples use a topology `T` of `ReferenceSimplex(2)` cells and the mesh
`D = mesh(x, EuclideanSpace(2))`, so $k = n = 2$.

- The Lagrange space of polynomial degree 1 on `T` has the value shape `()`.
- The Nédélec and the Raviart–Thomas spaces on `T` both have the value shape
  `(2,)`, and differ in the transformation law.
- The Hermite space is on `D`.

```{dropdown} Legacy UFL
The legacy UFL `Real` space is a space of constants, and the legacy UFL `Constant` is constant data.

In legacy UFL, a `FunctionSpace` holds a mesh and an element, and the element states the value shape, the map to the reference cell and the Sobolev space. In UFLx, a space may be defined on a topology, with an explicit law for its push-forward to a mesh, or on the mesh itself.
```

### `Function`

`Function` is an externally supplied member of a function space. It is a
field when its argument spaces are `()`, and otherwise a multilinear map
with those argument spaces.

| `Function(space, label)` |                                                                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| Extends              | `AbstractExpression` |
| Free function        | `function(space, label)`                                                                                                         |
| Takes                | `space: AbstractFunctionSpace`, `label: str \| tuple` |
| Properties | `space`, the function space supplying its type and evaluation rules, and `label`, the identifier distinguishing functions in that space; labels follow the rules of an [`Index`](expressions.md#index) |
| Static type          | The value shape, scalar type and argument spaces of `space`, the domain of `space`, and no free indices |

The label identifies the function, so two functions are the same only if the
spaces and the labels are equal. `Function(DualSpace(S), label)` is an
antilinear form, such as a known load (see [DualSpace](#dualspace)).

Multilinearity applies to the argument slots, not to spatial dependence.
Evaluation follows [Functions in a space](#functions-in-a-space).

The [form](#form), rather than the function, decides whether a function is an
argument or fixed data.

For a scalar space `V` on a topology `T`, `u = function(V, "u")` has domain
`T` and value shape `()`.

For a matrix, a library may provide `matrix_space(V, U)` with argument
spaces `(conjugate_space(V), U)`, no domain, scalar value shape, and the
common scalar type of `V` and `U`. For fields `u` in `U` and `v` in `V`:

```python
A = function(matrix_space(V, U), "A")
load = form(apply(A, (conjugate(v), u)), (conjugate(v),))
value = apply(load, (conjugate(v),))
```

`load` has argument spaces `(conjugate_space(V),)`; `value` has none and
represents $v^H A u$ in coefficient notation. No matrix entries are stored
in the expression.

```{dropdown} Legacy UFL
The legacy UFL `Coefficient`, `TestFunction` and `TrialFunction` are all functions. A test function and a trial function are the functions that a form states as the arguments at the positions 0 and 1.
```

## Variational forms

A variational form abstracts an expression over functions to make a multilinear
map. It binds functions as [`ComponentTensor`](expressions.md#componenttensor-and-indexsum)
binds indices; application supplies their values. Partial application selects explicit
binders by identity, not by slot number.

### Argument spaces

Every expression derives a `tuple` of argument spaces in position order, or
`()` for a non-form. Each entry is an `AbstractFunctionSpace` for one function,
or a `tuple` of at least two spaces with one scalar type for a product-space
argument. Different slots may have different scalar types.

Arity counts argument spaces, not functions within each argument: a functional
has arity 0, a linear form 1 and a bilinear form 2.

A node passes on the argument spaces of an operand in which the node is linear
or antilinear, and every other operand must have the argument spaces `()`, or
construction raises `ValueError`, unless its contract accepts a functional
value with a specified signature. The [static type
table](reference.md#appendix-static-type-table) gives the rule of each node class. A terminal
class states the argument spaces as part of the contract of the class, as the
value shape.

Binding an argument declares linearity in $V$ or antilinearity in
`conjugate_space(V)` (see [`Form`](#form)).

Static argument spaces let construction reject incompatible arities and sums
of linear and antilinear forms.

#### Bindable functions

A `Function` is bindable if its argument spaces are `()` or its space is a
`DualSpace`. `Form` and `DirectionalDerivative` bind these functions. A dual of a product is one functional argument; a `tuple`
of binders is one jointly linear product-space argument.

The argument spaces of a functional describe the values it accepts. Binding
that functional adds its space as one argument of the enclosing form, without
adding its own argument spaces again. Arbitrary multilinear operator-valued
arguments remain an open question.

### `ConjugateSpace`

`ConjugateSpace(space)` records conjugation of a complex function space.
It lets a form declare antilinearity in an argument. Writing $V$ for `space`,
the resulting space is $\bar V$.

| `ConjugateSpace(space)` |                                                                                                   |
| ------------------- | ------------------------------------------------------------------------------------------------- |
| Extends             | `AbstractFunctionSpace` |
| Free function       | `conjugate_space(space)` |
| Takes               | `space: AbstractFunctionSpace` |
| Returns             | `space` for a real `space`, `space.space` for a `ConjugateSpace`, and `ConjugateSpace(space)` otherwise |
| Properties | `space`, the original complex function space; it must not itself be a `ConjugateSpace` |
| Derived             | `domain`, `value_shape` and `scalar_type`, those of `space`; `argument_spaces`, those of `space` with each space conjugated, member by member for a `tuple` |
| Raises | The node constructor raises `ValueError` for a real space or a `ConjugateSpace`; the factory handles these cases as stated above |

A linear form on $\bar V$ is an antilinear form on $V$.
The core uses $\bar V$ only as an argument space, to record that a form
is antilinear in an argument. The conjugate
of an argument space is `conjugate_space(V)` for a space $V$ and $V$ for
`ConjugateSpace(V)`, so conjugating twice gives the space back (see
[`Conjugate`](expressions.md#conjugate)).

For `l` with space `DualSpace(V)`, `l` has argument spaces
`(conjugate_space(V),)` and `conjugate(l)` has `(V,)`.
`conjugate_space(dual_space(V))` and `dual_space(conjugate_space(V))` are
not identified structurally.

Space substitution rebuilds conjugate spaces through `conjugate_space`,
preserving these canonical forms.

### `DualSpace`

`DualSpace(space)` describes bounded antilinear functionals: maps that take a
member of `space` and return a scalar. `space` may also be a tuple of spaces,
representing one product-space argument.

Every space $V$ is assumed reflexive. The core identifies $V$ with its
double dual by the canonical map defined
[below](#the-double-dual).

| `DualSpace(space)` | |
| -------------- | - |
| Extends | `AbstractFunctionSpace` |
| Takes | `space: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]` |
| Properties | `space`, the space on whose members the functionals act, or a `tuple` of at least two spaces with the same scalar type; a single space must not be a `DualSpace` |
| Derived | `domain`, `value_shape`, `scalar_type` and `argument_spaces`, below |
| Raises | `ValueError` for a `DualSpace` or an invalid tuple |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `dual_space(space)` | | `space: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]` | `space.space` for a `DualSpace`, otherwise `DualSpace(space)` |

The derived properties describe a function in the dual space:

| Property | Meaning |
| -------- | ------- |
| `domain` | `None`; a functional has no spatial argument. |
| `value_shape` | `()`; applying a functional produces a scalar. |
| `scalar_type` | The scalar type of `space`, shared by all members when it is a tuple. |
| `argument_spaces` | One slot: the conjugate of `space`, or a tuple of the conjugated member spaces when `space` is a tuple. |

The descriptions below write `S` for `space`.

Space substitution rebuilds a dual through `dual_space`, preserving its
canonical form. For a `tuple` $S$, taking the dual twice returns that tuple;
product values remain tuples of members, not a single `Function`.

#### Antilinear forms

A function whose space is `DualSpace(S)` is antilinear on $S$ and linear on
its conjugate space. A residual bound as `(conjugate(v),)` and an assembled
load therefore have the same argument spaces and can be subtracted.

A linear form on a complex $S$ has argument spaces `(S,)` and represents an
element of `dual_space(conjugate_space(S))`. The core therefore fixes no
convention on the argument in which a form is antilinear.

For `l` with space `DualSpace(S)` and `w` in $S$, the pairing is
`apply(l, (conjugate(w),))`: the argument space of `l` is conjugated, and
[applications](#apply) are linear in their arguments.
The library supplying `l` evaluates the pairing.

#### The double dual

`dual_space(dual_space(V))` returns `V`. The canonical identification is

$$
\iota_V(u)(\ell) = \overline{\ell(u)},
$$

with the conjugation absent over real scalars. Thus
`form(conjugate(apply(l, (conjugate(u),))), (conjugate(l),))`
represents $u$ in the double dual.

Space membership and expression representation are separate. A
`Function` in the double dual is a `Function` in `V`. A scalar form with
one argument in `conjugate_space(dual_space(V))` represents an element of
`V` but retains that argument and its scalar output. Assembly may produce
the corresponding field; space normalisation does not change a form's
argument spaces or make it a pointwise field expression.

#### Examples

Let `v` and `q` be fields in spaces `V` and `Q` with the same scalar type.
A function in a dual space pairs with a conjugated argument:

```python
l = function(dual_space(V), "l")
pairing = apply(l, (conjugate(v),))
```

`l` has argument spaces `(conjugate_space(V),)`, and `pairing` is
$\langle l, v \rangle$.

For scalar fields `f` and `v` on `D`, with `v` in `V`, a residual form
and an assembled load have compatible argument spaces:

```python
R = form(integral(f * conjugate(v), D), (conjugate(v),))
residual = R - l
```

Both have argument spaces `(conjugate_space(V),)`. For complex `V`, the
linearly bound form `form(integral(f * v, D), (v,))` has `(V,)`; subtracting
`l` from it raises `ValueError`.

For a product of spaces, the argument is a tuple:

```python
m = function(dual_space((V, Q)), "m")
pairing_mixed = apply(m, ((conjugate(v), conjugate(q)),))
```

`m` has argument spaces `((conjugate_space(V), conjugate_space(Q)),)`.
A form may bind both participants of the pairing:

```python
P = form(pairing, (l, conjugate(v)))
```

`P` has argument spaces `(dual_space(V), conjugate_space(V))`: the pairing
is linear in `l` and antilinear in `v`.

```{dropdown} Legacy UFL
The legacy UFL `Cofunction` and `Coargument` are a fixed and a bound function in `dual_space(V)`. Both languages identify the double dual with the primal space.
```

### `Form`

`Form(expression, arguments)` turns selected functions in `expression` into
inputs of a multilinear map. `expression` is the body to evaluate;
`arguments` declares which functions are inputs and in what order. Other
functions remain fixed data for that form.

| `Form(expression, arguments)` |                                                                                                                                                                                                    |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                       | `AbstractExpression` |
| Free function                 | `form(expression, arguments)`                                                                                                                                                                      |
| Takes                         | `expression: AbstractExpression`, `arguments: tuple[Function \| Conjugate \| tuple[Function \| Conjugate, ...], ...]` |
| Properties | `expression`, the body of the form, and `arguments`, the ordered binding declarations described below |
| Requires                      | Every function in `arguments` distinct and [bindable](#bindable-functions); each tuple entry has one scalar type; no bound function occurs free in any resulting argument space, including inherited slots |
| Static type                   | That of the expression, with the argument spaces of the expression followed by the space `V` of each entry `v`, `conjugate_space(V)` of each entry `conjugate(v)`, or the `tuple` of the spaces of the members of a `tuple` entry, in order |
| Raises                        | `ValueError` for a conjugated expression that is not a function, a repeated or non-bindable function, a tuple with fewer than two members or mixed scalar types, or an argument space depending on a bound function |

Each entry of `arguments` declares one new argument slot:

| Entry | Meaning |
| ----- | ------- |
| A bindable `Function` | The form is linear in this function. |
| `conjugate(function)` for a bindable `Function` | The form is antilinear in this function. |
| A tuple of at least two functions or conjugated functions | One product-space argument; the form is jointly linear in the declared values. |

An entry's binder index is its index in `arguments`; its argument position
is the arity of `expression` plus that index. An entry `v` declares
linearity in `v`; `conjugate(v)` declares linearity in $\bar v$, hence
antilinearity in `v`. They give the same argument space for a real space.
A `tuple` entry declares joint linearity over its members' common scalar
field, so a term may depend on only one member. Output scalar type does not
determine argument linearity: a real argument may have a complex-valued output. The user is responsible for these declarations. Construction checks the
structural requirements above; an optional [linearity check](#checking-linearity)
may help diagnose mistakes.

Binding is lexical: an occurrence refers to the innermost enclosing `Form`
that binds the function, and is fixed otherwise. Thus a function may be an
argument in one form and fixed in another. Substitution renames bound functions
to avoid capture (see [Substitution and binding](substitution.md#substitution-and-binding)).

The expression may itself be a form; new arguments follow its existing
arguments, as in a [derivative form](#derivative-forms). No core pass
flattens `Form(Form(e, a), b)` into `Form(e, a + b)`, because products of forms
would require argument reordering. Assembly reads argument order from the
static type.

A form may depend nonlinearly on free parameters while remaining multilinear
in its declared arguments: mathematically, $a(x;v,u)$. Geometry $x$ can be such
a parameter without occupying an argument slot. Argument-space independence
is checked through construction dependencies, including domains and maps;
labels are opaque. It applies also to slots inherited from `expression`.
Spaces may depend on free geometry parameters, but differentiation in such a
parameter requires the fixed-space representation described under
[geometry dependence](differentiation.md#geometry-preparation).

A dual argument accepts a terminal load or a symbolic form with the same
functional signature. Substituting a symbolic form into a pairing creates
an application of that form (see [Expanding applications](lowering.md#expanding-applications)).

Tuple arguments express mixed problems without mixed spaces. Their members may
have different value shapes and domains, such as a Raviart–Thomas field and a
discontinuous field, or a boundary Lagrange multiplier. A library may assemble
one system or separate blocks. A block retains one member of each argument
and replaces the other members by `Zero`.

The expression may have any value shape and domain. A field-valued form
with one linear argument is a linear operator; a form of integrals has no
domain, even if its arguments do.

Let `u` and `v` be scalar functions in a space `V` on `D`:

```python
a = integral(u * v, D)
L = form(a, (v,))
B = form(a, (v, u))
```

`L` binds `v` and keeps `u` fixed. `B` binds both, giving argument spaces
`(V, V)`. The sum `B + L` raises `ValueError` because their arities differ.

For conjugated arguments and assembled loads, see [DualSpace](#dualspace).

For a mixed problem, let `u`, `v` belong to a velocity space `V`, and `p`, `q`
to a pressure space `Q`. Let `e` be the Stokes expression with terms
$\nabla u : \nabla v$, $-p \, \mathrm{div}\, v$ and $-q \, \mathrm{div}\, u$:

```python
stokes = form(e, ((v, q), (u, p)))
wrong_arity = form(e, (v, q, u, p))
```

`stokes` has argument spaces `((V, Q), (V, Q))`. The declaration in
`wrong_arity` is invalid: no term is linear in all four functions. An optional
linearity check may diagnose this mistake.

#### Checking linearity

Linearity is a mathematical requirement of a form declaration, not a property
that construction must prove. For an entry `v`, the expression must be linear
in `v`; for `conjugate(v)`, it must be antilinear in `v`. Zero satisfies either
convention. A tuple entry declares joint dependence in its members.

A library may provide a structural check of these declarations (see [New
questions](extensions.md#new-questions)). The core does not require a checking pass or a
classification API. An inconclusive check does not establish that a declaration
is invalid. Checks may assume the stated contracts of opaque operations;
libraries need not expose their implementations to prove linearity.

An implementation may also recognise linearity to simplify an expression. For
example, a complex-linear differential has a zero antilinear part. If the
analysis is inconclusive, evaluation uses the general rules. This changes
performance, not the meaning or admissibility of the expression. Users may
construct known derivative forms directly with `Form`.

```{dropdown} Legacy UFL
In legacy UFL, the arguments of a form are the test functions and the trial functions in the integrand, with the argument numbers of the functions, and a form is not an expression. In UFLx, a function has no argument number, the form states the arguments, and a functional is any expression with the argument spaces `()`. The legacy UFL `MixedElement` and `split` are an argument that is a `tuple` of functions.
```

### `Apply`

`Apply(form, arguments)` represents evaluation of a multilinear expression by supplying a
value for every argument slot. `form` may be a `Form`, a compatible sum or
product, or an externally supplied multilinear function. `arguments` contains
the input values in the order given by `form.argument_spaces`.

The name `arguments` has different roles in `Form` and `Apply`: `Form.arguments`
declares bound functions; `Apply.arguments` stores the values supplied to them.

| `Apply(form, arguments)` | |
| -------------------- | - |
| Extends | `AbstractExpression` |
| Free function | `apply(form, arguments)` |
| Takes | `form: AbstractExpression`, `arguments: tuple[AbstractExpression \| tuple[AbstractExpression, ...], ...]` |
| Properties | `form`, the multilinear expression to evaluate, and `arguments`, the ordered tuple of supplied values |
| Requires | `len(arguments) == len(form.argument_spaces) > 0`; each supplied value admissible for the corresponding space below |
| Returns | `form` for `apply(form, ())` when `form.argument_spaces == ()`, otherwise `Apply(form, arguments)` |
| Static type | The value shape, scalar type and domain of `form`; no argument spaces; the free indices of `form` followed by those of all supplied values, each `Index` taken once |
| Raises | `ValueError` for an incorrect argument count, inadmissible value or conflicting free-index extents |

#### Admissible argument values

| Argument space | Admissible value |
| -------------- | ------------------- |
| A space `space` | Same value shape and argument spaces as a function in `space`; scalar type `space.scalar_type` or real; domain `space.domain` or no domain |
| `ConjugateSpace(space)` | The conjugate type of an admissible value for `space`, including argument spaces; supplies the conjugated function value |
| Tuple of spaces | A tuple with one admissible member per space |

A product-space argument occupies one slot. Free indices give families of
applications. Admissibility checks static signatures, not discrete-space
membership. Real values enter complex arguments by complex embedding before
substitution, member by member.

For `Form(e, entries)`, the first `len(e.argument_spaces)` values supply
inherited slots. Substitute the remaining values simultaneously for `entries`,
then apply the resulting body to the inherited values. For a conjugated entry,
substitute the conjugate value for its function. Internal binders freshen to
avoid capture; supplied values are never substituted into one another.

Applications distribute over compatible sums and split across products by
operand arity. Libraries evaluate terminal applications. To leave slots open
on an arbitrary map, bind fresh functions around its full application:

```python
load = form(apply(A, (conjugate(v), u)), (conjugate(v),))
```

Spatial differentiation and pull-back transform the output, preserving the
supplied arguments:

| Expression | Value |
| ---------- | ----- |
| `gradient(apply(A, G), P)` | `apply(gradient(A, P), G)` |
| `pull_back(apply(A, G), psi)` | `apply(pull_back(A, psi), G)` |

Here `G` is the complete argument tuple. For $A(g)=a g$,
$(\nabla A)(g)=\nabla(a g)$, including the spatial derivative of $g$.

### Partial application

`partial(form, mapping)` supplies values for selected arguments of an explicit
outer `Form`, leaving the others open. `form` is the binder to partially
apply; `mapping` associates each selected binding declaration with its supplied
value. Substitution is simultaneous. This is a transformation pass, not a lazy node.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `partial(form, mapping)` | | `form: Form`, `mapping: dict` | The form with selected outer entries supplied, or its body when none remain |

Keys must be complete entries of `form.arguments`, including conjugation and
tuple grouping; values satisfy [argument admissibility](#admissible-argument-values).
Inherited slots and inner binders are not selectors. Other expression types,
including sums, are rejected; use `apply` instead.

Retained entries keep their order and identities. If a retained function
occurs free in any supplied value under the construction walk, raise
`ValueError`; explicitly `rebind` it first. Internal binders may freshen.
Conjugation and promotion follow `Apply`. An empty mapping returns `form`;
supplying every outer entry returns the body, including any inherited slots.

The result keeps `form`'s value shape, scalar type and domain. Its slots are the
inherited slots followed by retained entries. Combine the free indices of `form`
and all supplied values, requiring matching extents. Complete the body's type
before rebuilding the binder, retaining indices even from unused values.

```python
F = form(integral(u * v, D), (v, u))
partial(F, {u: w})  # form(integral(w * v, D), (v,))
partial(F, {u: v})  # ValueError: retained-argument capture
partial(rebind(F, {v: p}), {u: v})  # form(integral(v * p, D), (p,))
apply(F, (w, p))  # integral(p * w, D)
```

For `stokes = form(e, ((v, q), (u, p)))`,
`partial(stokes, {(u, p): (w, r)})` supplies one complete product-space entry.

### Rebinding

`rebind(binder, mapping)` renames selected declarations and their bound
occurrences. `binder` is a form or index-binding node; `mapping` associates
existing declarations with their replacement declarations. Renaming is
simultaneous and preserves the meaning of the bound map. This is a
transformation pass.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `rebind(binder, mapping)` | | `binder: Form \| ComponentTensor \| IndexSum`, `mapping: dict` | The binder with selected declarations renamed |

For `Form`, keys and values are complete argument entries; preserve order,
grouping, member spaces and conjugation. For index binders, keys are declared
indices and values are indices; preserve extents. Renaming respects lexical
scope. Swaps are allowed; internal binders may freshen to prevent capture.

Raise `ValueError` for unknown declarations, changed signatures, duplicate
resulting names, or targets free in the original binder under the construction
walk. Unsupported binders are Python type errors. An empty mapping returns `binder`.

Rebinding preserves alpha-equivalence: corresponding bound occurrences refer
to the same lexical declarations, with free names and parameter signatures
unchanged. It does not change structural equality, interning or structural keys:
differently named binders need not compare equal. Dependent parameter types
are not identified by alpha-renaming.

`replace` substitutes free data; `rebind` renames declarations; `partial`
supplies declared arguments. None overloads Python `()`.

```{dropdown} Legacy UFL
The action on the trial argument of an explicit bilinear form becomes
`partial(F, {u: w})`. Full evaluation uses `apply(F, (v_value, u_value))`.
The adjoint of `form(e, (conjugate(v), u))` is
`form(conjugate(e), (conjugate(u), v))`. Derivative forms bind direction
arguments, as described below.
```

### Derivative forms

`directional_derivative` gives the derivative in a supplied direction. When
the differential is linear in that direction, binding a fresh real direction
function constructs the corresponding derivative form:

```python
J = form(directional_derivative(F, x, a, w), (w,))
```

The form represents the map from a direction `w` to the derivative of `F`
with respect to `x` at `a` in that direction. Declaring it as a form assumes
linearity in `w`; existence of directional derivatives alone does not
establish this.

The new slot follows `F`'s existing slots. The direction must be admissible
for `x` and fresh in `F` and `a` under the construction walk. Complex
directions require the projections or real direction blocks below unless the
differential already has the declared complex linearity.

#### Complex and Wirtinger derivatives

For a complex variable, write $D[h]$ for its real differential at $a$. Its
complex-linear and antilinear contributions are

$$
D_+[h]=\tfrac12(D[h]-iD[ih]),\qquad
D_-[h]=\tfrac12(D[h]+iD[ih]).
$$

They sum to $D[h]$. For $F:\mathbb C\to\mathbb C$, $D_+[h]=F_z h$ and
$D_-[h]=F_{\bar z}\bar h$; setting $h=1$ gives the Wirtinger derivatives.
For fields, the parts are operators. For holomorphic $F$, $D_-=0$ and
$D_+[h]=F'(z)h$; construction does not check holomorphicity.

Bind a fresh complex direction `h`:

```python
i = Complex(Integer(0), Integer(1))
d = directional_derivative(F, z, a, h)
di = directional_derivative(F, z, a, i * h)
J_plus = form((d - i * di) / Integer(2), (h,))
J_minus = form((d + i * di) / Integer(2), (conjugate(h),))
```

Recover the full differential by supplying `k` to the new `J_plus` slot and
`conjugate(k)` to `J_minus`. Add their applications; the forms have different
argument spaces. Each contribution may be complex even when their sum is real.
For complex tuples, rotate every direction member by $i$. Library helpers may
provide these compositions.

#### Mixed real and complex variables

Derivative variables may mix real and complex spaces; each point and direction
member must be admissible for its variable. Form slots instead require one
scalar type per tuple. Use separate blocks, without rotating real variables.
For complex `z`, choose a fresh real direction `q` in a compatible library space:

```python
J_real = form(directional_derivative(F, z, a, q), (q,))
J_imag = form(directional_derivative(F, z, a, i * q), (q,))
```

The slots are real; outputs may be complex. Real-variable derivatives add
further blocks. For real matrices, apply `real_part` or `imaginary_part` to
derivative bodies before binding real test and direction functions.

See [Complex residuals and Hessians](../examples/forms.md#complex-residuals-and-hessians).

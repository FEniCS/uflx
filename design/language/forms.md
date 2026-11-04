# Functions and forms

## Function spaces

A function is a member of a function space, such as a coefficient or an
argument of a variational form.

### `AbstractFunctionSpace`

`AbstractFunctionSpace` describes a space of fields or multilinear maps,
with or without a domain or topology.

| `AbstractFunctionSpace` |                                                                                                                                                       |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `Base` |
| Provides                | `domain`, an [`AbstractDomain`](geometry.md#abstractdomain), an [`AbstractTopology`](geometry.md#abstracttopology) or `None`, `value_shape`, `scalar_type`, and `argument_spaces`, `()` for a space of fields |
| Defined by              | A subclass outside the core, such as a subclass for the finite element spaces of Basix, and [`DualSpace`](#dualspace) in the core                           |

A space with the domain `None` has functions that do not depend on a domain.

A library may supply `transformation_law(V)` for spaces with a
[`TransformationLaw`](geometry.md#transformationlaw). This is not a required space property
or a core overload; the core constructor remains `transformation_law(variances, weight)`.

#### Functions in a space

`Function(V, label)` uses the domain of $V$ and the value
shape, scalar type and argument spaces of $V$. These describe its expression
representation; a form may represent the same space element differently
(see [The double dual](#the-double-dual)). For a field,
the concrete subclass defines the function on the cells:

- A space whose values on the reference cell do not depend on the geometry is
  on a topology $T$. The concrete subclass defines
  `pull_back(u, cell_map(T, j, e))` for every cell group $j$ and every entry
  $e$, such as with one finite element for each cell type. The field on a
  mesh is a [`PushForward`](geometry.md#pushforward).
- A space whose values on the reference cell depend on the geometry is on a
  mesh $D$, such as Hermite, Argyris, Morley and the other elements that are
  not Piola equivalent. The concrete subclass defines
  `pull_back(u, cell_map(D, j, e))`, and reads $x$ from the cell map.
- For a coordinate domain, such as a reference cell, a function is a terminal
  that the concrete subclass evaluates.

The library evaluates these cell reads and their reference derivatives,
including at supplied reference positions. It also evaluates actions on
functional terminals and spatial transformations of their outputs. Core
lowering leaves these reads for the library; it does not compute basis or
coefficient values.

The solver locates physical points and supplies the selected cells and
reference positions before an evaluation is constructed. At interfaces or
overlaps, the cell selection is part of defining the value. A coordinate-only
`DomainMap` into a mesh is invalid (see [A family of points in a
mesh](geometry.md#a-family-of-points-in-a-mesh)).

The concrete subclass owns the polynomials, degrees of freedom and basis of
$V$, and stores the arguments identifying the space, such as the finite
element. Space identity is structural.

A library may provide [cell-local interpolation](extensions.md#cell-local-interpolation)
using the sampling rule and basis of the target space. These are library
data, not additional properties required by `AbstractFunctionSpace`.

#### Spaces of constants

A space with the domain `None` and the argument spaces `()` is a space of
constants. A function in the space is an expression with no domain, a
single array with the value shape and the scalar type of $V$.

The integrals of a form decide the cells on which a constant acts.

A constant unknown, such as a Lagrange multiplier, needs a space because form
arguments come from spaces. Constant data needs none: a literal or kernel-input
array with no domain is already a constant.

#### Mixed problems

A mixed problem needs no mixed space, because an argument of a [form](#form)
may be a `tuple` of functions in several spaces.

#### Examples

The examples use a topology `T` of `ReferenceSimplex(2)` cells and the mesh
`D = image(x, RD(2))`, so $k = n = 2$.

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

| `Function(V, label)` |                                                                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| Extends              | `AbstractExpression` |
| Free function        | `function(V, label)`                                                                                                         |
| Takes                | `V: AbstractFunctionSpace`, `label: str \| tuple` |
| Properties           | `space` $V$, an `AbstractFunctionSpace`, and `label`, a `str` or a `tuple` with the rules of the label of an [`Index`](expressions.md#index) |
| Static type          | The value shape, scalar type and argument spaces of $V$, the domain of $V$, and no free indices |

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
load = action(A, 1, u)
value = action(load, 0, conjugate(v))
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
binds indices; an action supplies an argument as [`Indexed`](expressions.md#indexed) supplies
an index to an axis.

### Argument spaces

Every expression derives a `tuple` of argument spaces in position order, or
`()` for a non-form. Each entry is an `AbstractFunctionSpace` for one function,
or a `tuple` of at least two spaces for a product-space argument.

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

`ConjugateSpace` is the complex conjugate space $\bar V$ of a complex space $V$.

| `ConjugateSpace(V)` |                                                                                                   |
| ------------------- | ------------------------------------------------------------------------------------------------- |
| Extends             | `AbstractFunctionSpace` |
| Free function       | `conjugate_space(V)` |
| Takes               | `V: AbstractFunctionSpace` |
| Returns             | `V` for a real `V`, `V.space` for a `ConjugateSpace`, and `ConjugateSpace(V)` otherwise |
| Properties          | `space` $V$, an `AbstractFunctionSpace` with the complex scalar type, not a `ConjugateSpace`       |
| Derived             | `domain`, `value_shape` and `scalar_type`, those of $V$; `argument_spaces`, those of $V$ with each space conjugated, member by member for a `tuple` |
| Raises              | `ValueError` for a real $V$ or a `ConjugateSpace`  |

$\bar V$ has the functions of $V$ with the conjugate scalar multiplication, so
a function in $\bar V$ is the conjugate of a function in $V$, and a linear
form on $\bar V$ is an antilinear form on $V$. The core uses $\bar V$ only as
an argument space, to record that a form is antilinear in an argument. The conjugate
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

`DualSpace` is the space of bounded antilinear functionals on a space or a
product of spaces. All spaces in the language are reflexive: the core
identifies the double dual with the original space, but keeps a space
and its dual distinct. A Riesz map between them requires a chosen inner
product and belongs to the supplying library.

| `DualSpace(S)` | |
| -------------- | - |
| Extends | `AbstractFunctionSpace` |
| Takes | `S: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]` |
| Properties | `space` $S$, a space other than a `DualSpace`, or a `tuple` of at least two spaces with the same scalar type |
| Derived | `domain`, `value_shape`, `scalar_type` and `argument_spaces`, below |
| Raises | `ValueError` for a `DualSpace` or an invalid tuple |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `dual_space(S)` | | `S: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]` | `S.space` for a `DualSpace`, otherwise `DualSpace(S)` |

`DualSpace(S)` provides:

- No domain, and value shape `()`.
- The argument spaces `(S',)`, with $S'$ = `conjugate_space(S)` for a
  space, or the tuple of the conjugate spaces for a `tuple` $S$.
- The scalar type of $S$, common to the members of a `tuple`.

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
`action(l, 0, conjugate(w))`: the argument space of `l` is conjugated, and
[actions](#action) are linear in their arguments.
The library supplying `l` evaluates the pairing.

#### The double dual

`dual_space(dual_space(V))` returns `V`. The canonical identification is

$$
\iota_V(u)(\ell) = \overline{\ell(u)},
$$

with the conjugation absent over real scalars. Thus
`form(conjugate(action(l, 0, conjugate(u))), (conjugate(l),))`
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
pairing = action(l, 0, conjugate(v))
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
pairing_mixed = action(m, 0, (conjugate(v), conjugate(q)))
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

`Form` binds some functions of an expression as arguments of a variational form.

| `Form(expression, arguments)` |                                                                                                                                                                                                    |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                       | `AbstractExpression` |
| Free function                 | `form(expression, arguments)`                                                                                                                                                                      |
| Takes                         | `expression: AbstractExpression`, `arguments: tuple[Function \| Conjugate \| tuple[Function \| Conjugate, ...], ...]` |
| Properties                    | `expression`, an expression, and `arguments`, a `tuple` of entries, each a [function](#function) `v`, the [conjugate](expressions.md#conjugate) `conjugate(v)` of a function, or a `tuple` of at least two of these |
| Requires                      | Every function in `arguments` distinct and [bindable](#bindable-functions) |
| Static type                   | That of the expression, with the argument spaces of the expression followed by the space $V$ of each entry `v`, `conjugate_space(V)` of each entry `conjugate(v)`, or the `tuple` of the spaces of the members of a `tuple` entry, in order |
| Raises                        | `ValueError` for a conjugated expression that is not a function, a repeated or non-bindable function, or a `tuple` of fewer than two members |

An entry's binder index is its index in `arguments`; its argument position
is the arity of `expression` plus that index. An entry `v` declares
linearity in `v`; `conjugate(v)` declares linearity in $\bar v$, hence
antilinearity in `v`. They give the same argument space for a real space.
A `tuple` entry declares joint linearity in its members, so a term may depend
on only one member. The user is responsible for these declarations. Construction checks the
structural requirements above; an optional [linearity check](#checking-linearity)
may help diagnose mistakes.

Binding is lexical: an occurrence refers to the innermost enclosing `Form`
that binds the function, and is fixed otherwise. Thus a function may be an
argument in one form and fixed in another. Substitution renames bound functions
to avoid capture (see [Substitution and binding](evaluation.md#substitution-and-binding)).

The expression may itself be a form; new arguments follow its existing
arguments, as in a [Gateaux derivative](#gateaux-derivative). No core pass
flattens `Form(Form(e, a), b)` into `Form(e, a + b)`, because products of forms
would require argument reordering. Assembly reads argument order from the
static type.

A dual argument accepts a terminal load or a symbolic form with the same
functional signature. Substituting a symbolic form into a pairing creates
an action on that form (see [Applying actions](evaluation.md#applying-actions)).

Tuple arguments express mixed problems without mixed spaces. Their members may
have different value shapes and domains, such as a Raviart–Thomas field and a
discontinuous field, or a boundary Lagrange multiplier. A library may assemble
one system or separate blocks. A block retains one member of each argument
and replaces the other members by `Zero`.

The expression may have any value shape and domain. A form of a field is a
linear operator; a form of integrals has no domain, even if its arguments do.

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

### `Action`

`Action` gives an expression to one argument of a form.

| `Action(f, position, g)` |                                                                                                                                                        |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Extends                  | `AbstractExpression` |
| Free function            | `action(f, position, g)`                                                                                                                               |
| Takes                    | `f: AbstractExpression`, `position: int`, `g: AbstractExpression \| tuple[AbstractExpression, ...]` |
| Properties               | `form` $f$, `position`, an `int`, and `argument`, the expression $g$, or a `tuple` of expressions |
| Requires                 | `0 <= position < len(f.argument_spaces)` and the conditions below |
| Static type              | That of $f$, with the free indices of $f$ followed by those of $g$, each `Index` taken once, and the argument spaces of $f$ without the entry at `position` |

| Argument space at `position` | Admissible $g$ |
| ---------------------------- | -------------- |
| $V$ | Same value shape and argument spaces as a function in $V$; its scalar type or real; its domain or no domain |
| `ConjugateSpace(V)` | The conjugate type of an admissible value for $V$, including argument spaces; supplies the value of `conjugate(v)` |
| Tuple of spaces | A tuple with one admissible member per space |

An index free in both $f$ and $g$ must have the same extent in both.

Free indices in $g$ give a family of actions, one per index value. For a
`tuple` $g$, its free indices are those of its members.

An admissible real value for a complex argument is interpreted in that
argument's complex scalar type. Evaluation makes this embedding explicit
before substitution (see [Type compatibility of rewrites](evaluation.md#type-compatibility-of-rewrites)).

Let $p$ be `position` and $n$ the arity of $e$. The value of
`Action(Form(e, arguments), p, g)` is determined by the argument's position:

| Position | Value |
| -------- | ----- |
| $p < n$ | `Form(Action(e, p, g), arguments)` |
| $p \ge n$ | Remove argument entry $p-n$ and substitute $g$ for it in $e$; bind the remaining entries in `Form`, or return $e$ if none remain |

Substitution renames a remaining bound function if it is free in $g$ (see
[Substitution and binding](evaluation.md#substitution-and-binding)). For `conjugate(v)`,
substitute `conjugate(g)` for `v`: $g$ is the value of the entry, making the
action linear in $g$ at every position. For a tuple entry, substitute the
members of $g$ simultaneously for the matching members.

Actions on other forms follow by linearity. An action on a terminal is
defined by its library, for each value of the free indices of $g$. It may
leave argument slots open; further actions remain nested, with positions
counted among the remaining slots. The library evaluates these partial or
complete applications. A dual pairing is the one-slot case.

Spatial differentiation and pull-back transform the output of a form, keeping
the supplied argument fixed:

| Expression | Value |
| ---------- | ----- |
| `gradient(action(A, p, g), P)` | `action(gradient(A, P), p, g)` |
| `pull_back(action(A, p, g), psi)` | `action(pull_back(A, psi), p, g)` |

The transformed form retains its argument spaces. For $A(g) = a g$,
$(\nabla A)(g) = \nabla(a g)$, including the derivative of $g$.

A dual argument may receive a symbolic form. With `P` from the
[`DualSpace`](#dualspace) example and scalar fields `f`, `w` in `V` on `D`:

```python
L = form(integral(f * conjugate(w), D), (conjugate(w),))
transferred = action(P, 0, L)
```

Its value is `form(integral(f * conjugate(v), D), (conjugate(v),))`.
For a conjugated dual binder, the same substitution rule applies: supplying
`conjugate(L)` for `conjugate(l)` substitutes `L` for `l`, in value.

#### Examples

For scalar functions `u`, `v` and `w` in `V` on `D`, set:

```python
F = form(integral(u * v, D), (v, u))
M = form(integral(u * conjugate(v), D), (conjugate(v), u))
```

| Action | Value |
| ------ | ----- |
| `action(F, 1, w)` | `form(integral(w * v, D), (v,))` |
| `action(F, 1, v)` | `form(integral(v * v_fresh, D), (v_fresh,))`, renaming the binder to keep the supplied `v` free |
| `action(M, 0, conjugate(w))` | `form(integral(u * conjugate(w), D), (u,))` |

Actions distribute over compatible sums. `action(F, 1, directions[k])`
retains the free index `k`. For `stokes = form(e, ((v, q), (u, p)))`,
`action(stokes, 1, (w, r))` supplies both members of the second argument
and leaves a linear form in `(v, q)`.

```{dropdown} Legacy UFL
Legacy UFL operations are actions and rebindings. The action on the last argument of a bilinear form is `action(f, 1, w)`. The adjoint of `form(e, (conjugate(v), u))` is `form(conjugate(e), (conjugate(u), v))`, because `conjugate(e)` is linear in `v` and antilinear in `u`. The Gateaux derivative with respect to a function adds an argument (see [Gateaux derivative](#gateaux-derivative)).
```

Worked examples for this chapter: [the Laplace-Beltrami operator](../examples/surfaces.md#the-laplace-beltrami-operator), [the flux through a surface](../examples/surfaces.md#the-flux-through-a-surface), [Helmholtz with a perfectly matched layer](../examples/forms.md#helmholtz-with-a-perfectly-matched-layer), [mixed Poisson with Raviart-Thomas elements](../examples/forms.md#mixed-poisson-with-raviart-thomas-elements), [dual-weighted residual error estimation](../examples/estimation.md#dual-weighted-residual-error-estimation), [a classic FEniCS layer for the Poisson problem](../examples/forms.md#a-classic-fenics-layer-for-the-poisson-problem).

### Gateaux derivative

`gateaux_derivative(F, x, a, w)` constructs a derivative form, adding an
argument $w$ after the arguments of $F$. It composes `Form` and
`DirectionalDerivative` with tensor-aware arithmetic; it introduces no node class.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `gateaux_derivative(F, x, a, w)` | | `F: AbstractExpression`, `x: Function \| tuple`, `a: AbstractExpression \| tuple`, `w: Function \| Conjugate \| tuple` | The `Form` constructed below |

The variable $x$ is a [bindable function](#bindable-functions), or a tuple of
at least two distinct bindable functions. The point $a$ and entry $w$ match
its tuple structure and length. Member by member, $a$ is an admissible value
for the space of $x$ as in [`Action`](#action), and the function underlying
$w$ belongs to that space. The functions of $w$ are distinct, none is free
in $F$ or $a$, and tuple members of $w$ are either all plain functions or all
conjugates. Shared free indices of $F$ and the points have matching extents.

The types of $F$, $x$, $a$ and $w$ are checked statically. Construction raises
`ValueError` if a remaining structural requirement fails.
It checks these requirements, not existence of the differential
or holomorphicity. The resulting form has the value shape of $F$, the free
indices of $F$ and $a$, and the argument spaces of $F$ followed by the space
of entry $w$. Its scalar type is that of $F$ for real $x$, and complex for
complex $x$.

The members of a `tuple` $x$ have one scalar type, the scalar type of $x$.
The linear and the antilinear parts below are defined only on a complex
space, and the product of a real and a complex space is a real space, so a
`tuple` that mixes the two is rejected.

```{dropdown} Open: mixed real and complex tuples
A complex state with a real unknown, such as a Lagrange multiplier, needs a
mixed `tuple`. A convention could put the part in the real members into
$J_+$. Until then, take the derivatives separately.
```

The real differential in a direction $v$ is
$\delta F[v]$ = `directional_derivative(F, x, a, v)`, with the members of a
`tuple` variable substituted at once (see
[`DirectionalDerivative`](geometry.md#directionalderivative)). $\delta F$ is
real-linear in $v$. Write $h$ for the function underlying $w$, member by
member for a tuple; for a plain entry, $h=w$. Construction uses
`i_unit = Complex(Integer(0), Integer(1))` and binds the entry as follows:

| Space and entry | Construction |
| --------------- | ------------ |
| A real space, either entry convention | `Form(δF[h], (w,))` |
| A complex space, plain entry `w` | `Form(Integer(1) / Integer(2) * (δF[h] - i_unit * δF[i_unit * h]), (w,))` |
| A complex space, conjugate entry `w` | `Form(Integer(1) / Integer(2) * (δF[h] + i_unit * δF[i_unit * h]), (w,))` |

For a tuple direction, multiplication by `i_unit` is member by member. The
complex bodies are respectively the linear part $J_+[h]$ and the antilinear
part $J_-[h]$ of the real differential. Binding follows the ordinary scopes
of `Form` and `DirectionalDerivative`; no separate derivative-form binding
rule or evaluation pass is needed.

The decomposition $\delta F = J_+ + J_-$ is unique. If the differential
is complex-linear, $J_- = 0$; if complex-antilinear, $J_+ = 0$.
Recognising these cases during lowering is optional. Users who know the
dependence may construct the corresponding forms directly.

The derivative form binds $w$ after the arguments of $F$. An action on that
new argument gives the selected differential part. For complex variables,
recover the full real differential by acting on $J_+$ with $h$ and on $J_-$
with $\overline h$, then adding the results. The two forms themselves have
different argument spaces and cannot be added.

See [Complex residuals and Hessians](../examples/forms.md#complex-residuals-and-hessians)
for Wirtinger derivatives and Newton's method.

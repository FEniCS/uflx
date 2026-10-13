# Core UFLx language definition

This note describes the smallest core language behind UFLx. It intentionally
avoids implementation details, compiler algorithms, and broad redesign ideas.
It assumes knowledge of differentiable geometry and functional analysis. Some
advanced differential geometry concepts are explained in more detail.

This document is aimed at users and developers with extensive knowledge of
domain specific languages for variational forms (e.g. UFL, ngsolve).

## What UFLx is and isn't

In one sentence:

> UFLx is a typed symbolic language writing variational forms.

UFLx is a symbolic language embedded in Python for writing variational forms.

The UFLx core language (henceforth core) should not use finite element or
numerical analysis terminology - differential geometry and functional analysis
(and their fascinating intersection), without finite element (numerical
analysis) concepts, should be all that should is required.

Finite element and numerical analysis concepts, when necessary, will be added
on top of the core, or provided by external libraries (Basix being the
'in-house' example). Basix and DOLFINx and their assumptions should not be
visible in the core.

The core requires an understanding of modern differential geometry and
functional analysis. UFLx will then provide higher-level functionality based on
this core which mimics the 'easier' interface presented by other competing
packages.

UFLx is not a symbolic language for reasoning about problems in functional
analysis or differential geometry. Hence, the core should only contain the
minimal concepts and relations necessary for writing a wide range of
variational forms as code for subsequent compilation.

A UFLx program is ordinary Python code that constructs symbolic objects and
connects them through a graph. The primary result is usually a Form, such as a
bilinear form, linear form, or functional. Python code will be used to motivate
examples but concrete implementation is out of scope.

UFLx does not handle concrete data or numerical algorithms; e.g. actual meshes,
boundary markers, numerical values of coefficients, or assembly procedures.
Those belong to the surrounding problem-solving environment, such as DOLFINx,
Firedrake, FIAT and Basix. UFLx instead provides a specified set of leaf nodes
which the external library must provide.

## Basic principles

- All core object definitions should be closely related to their strict
  mathematical definition.
- Extra definitions, properties or terminology outside of the established
  mathematical fields noted must be avoided, in particular ideas from more
  applied fields such as mechanics and engineering-oriented computational
  sciences.
- Strictness and expressivity, not brevity, is key in the core - brevity can be
  added later with higher-level APIs.
- The mathematical fields provide a rich foundation for writing variational
  problems in a uniform way across geometric dimensions, i.e. 'special cases'
  for particular dimensions do not exist, and so should not exist in the core.
- The core should not contain any algorithm that is designed to produce
  convenient output for form compilers.

## Basic implementation principles

- UFLx should be deterministic on the same machine using the same byte-code.
- Performance is *important* and so any data structures and algorithms in the
  core should be carefully chosen and tested for scalability.
- The minimal language should only need a relatively small number of object
  kinds. When designing the language, consider whether objects can be composed,
  rather than extended.
- Explicit is better than implicit.
- Boilerplate and error-prone patterns should be minimised.
- Free functions should be preferred over class members. Class members should
  immediately call free functions wherever possible. Free functions maximise
  the chance of code re-use and are easier to reason about.
- A parent class should avoid returning a subclass.
- Simplifications and rewriting should be avoided at all costs during user
  graph construction, maintaining the user's intent. Only later algorithmic
  passes should re-write. A good example is `(a + b) + c` in floating point
  arithmetic. The user may have very good reason to have `(a + b)` performed
  before adding the answer to `c`!

## The minimal object model

### The base object `Base`

The base object in UFLx must be extended by all objects that are included in
the graph. The graph is a directed acyclic graph (DAG).

Everything needed to define the variational problem must be recoverable from
the graph. No information is held outside of the graph.

Objects must satisfy the following rules:

1. All construction arguments must be deeply immutable values. For example, a
   `tuple` of `list`s is not allowed.
1. All construction arguments must be compared by type and value, so that `1`,
   `1.0` and `True` are distinct, `0.0` and `-0.0` are also distinct, and every
   value is equal to itself. Every construction argument must have a hash that
   is deterministic across processes.
1. Construction must depend only on the construction arguments, not on global
   state.
1. An object must be reconstructible from its class and construction arguments,
   and the reconstructed object must be equal to the original.
1. Any results produced by the class after construction must be a result of
   only the class and the construction arguments.
1. The state of the class must be immutable after construction.
1. Any ordering of objects must come from a total order that is deterministic
   across processes.
1. The edges of an object must be exactly the base objects among its
   construction arguments, in order.
1. Equality (`==`) must be structural and return a `bool`. Equality does not
   decide mathematical equivalence, so $a + b$ and $b + a$ are equal only if
   construction normalises $a + b$ and $b + a$ to equal construction arguments.

Hence two objects are equal if and only if the two objects have the same class
and equal construction arguments.

To satisfy the above, the base object must only accept a small curated set of
types.

Because objects are immutable, an object need only be validated once, at
construction. Every object that exists is therefore valid.

Terminal names are construction arguments, so two terminals with the same name
and function space are equal (see [Terminals](#terminals)).

The base object should help consumers follow these rules, but should not
prevent downstream users from working around these rules.

An object cannot be changed. To obtain a modified object, construct a new
object and reconstruct every object that depends on the original object, up to
the root of the graph.

The above design implies that construction arguments should be small and cheap
to hash and compare, and the graph implementation must be performant.

### Class hierarchy

Every object in the graph extends `Base`. The classes defined in this document
form the following hierarchy:

```
Base
├── AbstractExpression
│   ├── AbstractScalar                    literals, no set
│   │   ├── Integer
│   │   ├── Rational
│   │   ├── RealScalar
│   │   └── ComplexScalar
│   ├── Coordinate                        terminal, X on P
│   ├── Tensor                            built from scalar entries
│   ├── Indexed                           full indexing only
│   ├── UnaryOperator
│   │   ├── Gradient                      grad, on coordinate domains only
│   │   ├── Transpose
│   │   ├── Determinant
│   │   └── SquareRoot
│   └── BinaryOperator
│       ├── Sum                           +
│       ├── Product                       *
│       ├── Power                         **
│       └── MatrixProduct                 @
│
└── AbstractDomain                        parameter_domain(D), mapping(D)
    ├── AbstractCoordinateDomain          stores n, mapping is X
    │   └── RD                            R^n
    └── ParameterisedDomain               stores P and ψ, mapping is ψ
```

The edges of the graph point from an object to the objects that the object
stores, and never back:

```
ParameterisedDomain ──► P (AbstractCoordinateDomain)
        │
        └─────────────► ψ (AbstractExpression) ──► ... ──► Coordinate ──► P

Gradient ──► f        Indexed ──► expression        Tensor ──► entries
```

### Free functions

Free functions are not objects in the graph. A free function builds or returns
objects in the graph. Operator overloads call the free functions.

| Group                | Free function                  | Overload  | Returns                                |
| -------------------- | ------------------------------ | --------- | -------------------------------------- |
| Conversion           | `as_literal(value)`            |           | `AbstractScalar`                       |
|                      | `as_expression(value)`         |           | `AbstractExpression`                   |
| Tensors and indexing | `tensor(entries)`              |           | `Tensor`                               |
|                      | `vector(entries)`              |           | `Tensor` of rank 1                     |
|                      | `matrix(entries)`              |           | `Tensor` of rank 2                     |
|                      | `indexed(f, indices)`          | `f[i]`    | `Indexed`                              |
| Arithmetic           | `add(f, g)`                    | `+`       | `Sum`                                  |
|                      | `subtract(f, g)`               | `-`       | not yet specified                      |
|                      | `multiply(f, g)`               | `*`       | `Product`                              |
|                      | `divide(f, g)`                 | `/`       | not yet specified                      |
|                      | `negate(f)`                    | unary `-` | not yet specified                      |
|                      | `absolute(f)`                  | `abs`     | not yet specified                      |
|                      | `power(f, g)`                  | `**`      | `Power`                                |
|                      | `matrix_product(f, g)`         | `@`       | `MatrixProduct`                        |
| Tensor algebra       | `transpose(f)`                 |           | `Transpose`                            |
|                      | `determinant(f)`               |           | `Determinant`                          |
| Elementary functions | `square_root(f)`               |           | `SquareRoot`                           |
| Differentiation      | `grad(f)`                      |           | `Gradient`                             |
|                      | `derivative(a, u, du)`         |           | Gateaux derivative, see section 4      |
| Domains              | `parameterised_domain(P, psi)` |           | `ParameterisedDomain`                  |
|                      | `ambient_domain(D)`            |           | `RD(n)`                                |
|                      | `parameter_domain(D)`          |           | `AbstractCoordinateDomain`             |
|                      | `mapping(D)`                   |           | `AbstractExpression` with shape `(n,)` |
|                      | `coordinate(P)`                |           | `Coordinate`                           |
| Geometry             | `jacobian(D)`                  |           | `Gradient` of `mapping(D)`             |
|                      | `metric(D)`                    |           | `MatrixProduct` of $J^T$ and $J$       |
|                      | `measure_density(D)`           |           | `SquareRoot` of the determinant of $g$ |

## The lowest objects

### Expressions

An expression denotes an array-valued function defined on a set of the
expression. An expression with no set denotes a single array. An
expression can be real or complex valued.

#### `AbstractExpression`

This base class is sketched in pseudo-code according to the principles outlined
above. Later classes will not be sketched in detail.

```
# expressions.py
class AbstractExpression(Base, ABC):
    """Abstract base class for expressions."""
    @abstractmethod
    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression"

    # Operator overloading must be in the class.
    def __add__(self, other: object) -> AbstractExpression:
        """Add."""
        # Calls the free function arithmetic.add.
        ...

    def __radd__(self, other: object) -> AbstractExpression:
        """Add."""
        # Calls the free function arithmetic.add.
        ...

    def __sub__(self, other: object) -> AbstractExpression:
        """Subtract."""
        # Calls the free function arithmetic.subtract.
        ...

    def __rsub__(self, other: object) -> AbstractExpression:
        """Subtract."""
        # Calls the free function arithmetic.subtract.
        ...

    def __mul__(self, other: object) -> AbstractExpression:
        """Multiply."""
        # Calls the free function arithmetic.multiply.
        ...

    def __rmul__(self, other: object) -> AbstractExpression:
        """Multiply."""
        # Calls the free function arithmetic.multiply.
        ...

    def __matmul__(self, other: object) -> AbstractExpression:
        """Matrix multiply."""
        # Calls the free function arithmetic.matrix_product.
        ...

    def __truediv__(self, other: object) -> AbstractExpression:
        """Divide."""
        # Calls the free function arithmetic.divide.
        ...

    def __rtruediv__(self, other: object) -> AbstractExpression:
        """Divide."""
        # Calls the free function arithmetic.divide.
        ...

    def __neg__(self) -> AbstractExpression:
        """Negate."""
        # Calls the free function arithmetic.negate.
        ...

    def __abs__(self) -> AbstractExpression:
        """Absolute value."""
        # Calls the free function arithmetic.absolute.
        ...

    def __pow__(self, power: object) -> AbstractExpression:
        """Raise to an integer power."""
        # Calls the free function arithmetic.power.
        ...

    def __getitem__(self, indices: object) -> AbstractExpression:
        """Index."""
        # Calls the free function indexed.
        ...

    def __iter__(self) -> NoReturn:
        """Raise TypeError, because expressions are not iterable."""
        # Without __iter__, Python iterates through __getitem__.
        ...

# arithmetic.py

# A free function
def add(first: AbstractExpression, second: AbstractExpression) -> AbstractExpression:
    """Add two expressions.

    Args:
        first: The first summand.
        second: The second summand.

    Returns:
        The sum of the two expressions.

    Raises:
        ValueError: If the value shapes of the summands differ.
    """
    # Checks on value_shape
    return Sum((first, second))
```

### Literals

A literal is a number written down in an expression, e.g. `Integer(2) * u`. A
literal's value is part of the expression itself.

> Literals are Expressions, but they do not depend on a set (no domain).

Base:

```
AbstractScalar(AbstractExpression)
```

Concrete:

```
Integer(2)
Rational(2, 4)
RealScalar(3.5)
ComplexScalar(RealScalar(3.0), RealScalar(4.0))
```

Literals should not fold on construction, e.g.

```
Integer(2) + Integer(3) != Integer(5)
Integer(2) + Integer(3) == Sum(Integer(2), Integer(3))
```

#### Conversion of Python values with `as_literal` and `as_expression`

>>> LLM generated content.

Two free functions convert Python values into expressions. The free function
`as_literal` converts a Python number into a literal. The free function
`as_expression` accepts either an `AbstractExpression` or a Python number.

Every user-facing free function and every operator overload applies
`as_expression` to each operand or entry that the function accepts as an
expression. Class constructors apply neither `as_literal` nor `as_expression`.
Class constructors accept only the types that the class stores.

##### `as_literal`

>>> Syntax nicety, but quite nice for users...

The free function `as_literal(value)` converts a Python number into an
`AbstractScalar` literal. The function `as_literal` dispatches on the exact
type of `value`, not on subclasses of the type:

| `type(value)`        | Result                                                          |
| -------------------- | --------------------------------------------------------------- |
| `int`                | `Integer(value)`                                                |
| `fractions.Fraction` | `Rational(value.numerator, value.denominator)`                  |
| `float`              | `RealScalar(value)`                                             |
| `complex`            | `ComplexScalar(RealScalar(value.real), RealScalar(value.imag))` |

The function `as_literal` raises:

- `TypeError` for every other type. The rejected types include:
  - every `AbstractExpression`, including literals;
  - `bool`, although `bool` is a subclass of `int`;
  - every other subclass of `int`, `float`, `complex` or `Fraction`;
  - NumPy scalars;
  - `str` and `tuple`.
- `ValueError` for a `float` that is not finite, and for a `complex` with a
  real or imaginary part that is not finite.

Exact type dispatch keeps the result of `as_literal` deterministic. Base rule 2
compares construction arguments by type and value, so `True` and `1` must not
convert to the same literal. A NaN is not equal to itself, so a NaN cannot be a
construction argument.

The function `as_literal` preserves the Python type and the value of `value`
exactly:

- `as_literal(1)` is `Integer(1)` and `as_literal(1.0)` is `RealScalar(1.0)`,
  so `as_literal(1) != as_literal(1.0)`.
- `as_literal(-0.0)` is `RealScalar(-0.0)`, and `RealScalar(-0.0)` is not equal
  to `RealScalar(0.0)`.
- `as_literal(1j)` is `ComplexScalar(RealScalar(0.0), RealScalar(1.0))`. The
  zero real part is kept.
- `as_literal(Fraction(2, 4))` is `Rational(1, 2)`, because Python normalises
  `Fraction(2, 4)` to `Fraction(1, 2)` before `as_literal` sees the value.

The function `as_literal` performs no folding and no simplification.

##### `as_expression`

>> Syntax nicety, but quite nice for users...

The free function `as_expression(value)` returns `value` unchanged if `value`
is an `AbstractExpression`. Otherwise the function `as_expression` returns
`as_literal(value)`, and raises the errors of `as_literal`.

The function `as_expression` is the only place where a value is checked for
being an `AbstractExpression` before conversion. The function `as_expression`
is idempotent, so `as_expression(as_expression(v)) == as_expression(v)` for
every accepted `v`.

The functions `as_literal` and `as_expression` convert a single value. Neither
function converts a `tuple`. Free functions that accept nested `tuple`s, such
as `tensor`, apply `as_expression` to each entry.

### Operators

```UnaryOperator(AbstractExpression)```
```BinaryOperator(AbstractExpression)```

Downstream basic ops are relatively obvious.

### Basic operators

Concrete subclasses of `UnaryOperator` needed for the first example. Each class
stores one expression, which is the only edge, and derives the set of that
expression. Each free function applies `as_expression`, raises `ValueError` for
a wrong value shape and performs no simplification.

| Class         | Free function | Argument shape | Result shape |
| ------------- | ------------- | -------------- | ------------ |
| `Transpose`   | `transpose`   | `(m, n)`       | `(n, m)`     |
| `Determinant` | `determinant` | `(n, n)`       | `()`         |
| `SquareRoot`  | `square_root` | `()`           | `()`         |

Concrete subclasses of `BinaryOperator` needed for the first example. Each
class stores two expressions, which are the two edges in order. The set is the
set of the operands. Operands with a set must have the same set. Each free
function applies `as_expression` to both operands, raises `ValueError` for
wrong value shapes and performs no simplification.

| Class           | Free function    | Overload | Operand shapes     | Result shape |
| --------------- | ---------------- | -------- | ------------------ | ------------ |
| `Sum`           | `add`            | `+`      | `s`, `s`           | `s`          |
| `Product`       | `multiply`       | `*`      | `()`, `()`         | `()`         |
| `Power`         | `power`          | `**`     | `()`, `()`         | `()`         |
| `MatrixProduct` | `matrix_product` | `@`      | `(m, k)`, `(k, n)` | `(m, n)`     |

## Differential geometry

We begin with differential geometry because the function spaces and operators
of UFLx are defined on a domain. The function spaces and operators require two
kinds of structure on the domain:

- topological structure, which is the dimension, boundary, orientation and
  decomposition into cells of the domain,
- geometric structure, which is the metric of the domain and the measure
  induced by the metric.

### AbstractDomain

An abstract base class for domains `AbstractDomain`.

Defined by:

- The ambient dimension $n$ of the Euclidean space $R^n$ whose coordinates name the
  points of this domain (geometric dimension). 
- The topological dimension, which is the maximum dimension of the domains
  cells (standard Lebesgue definition). 

Examples:
- A domain consisting of triangle and interval domain has topological dimension
  2.
- A sphere in $R^3$ has a geometric dimension 3 and topological dimension 2.

Every `AbstractDomain` $D$ with topological dimension $k$ and geometric
dimension $n$ provides two free functions:

- `parameter_domain(D)`, an `AbstractCoordinateDomain` $P$ of dimension $k$.
- `mapping(D)`, an `AbstractExpression` $\psi$ on $P$ with value shape `(n,)`.
  The mapping $\psi$ maps each point of $P$ to a point of $D$, and the image
  $\psi(P)$ is $D$.

The mapping of $D$ gives the coordinates in $R^n$ of the points of $D$. The
geometry of $D$ is derived from the mapping (see
[Geometry of a domain](#geometry-of-a-domain)).

### AbstractCoordinateDomain

An abstract base class for coordinate domains `AbstractCoordinateDomain`, and a
subclass of `AbstractDomain`.

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$.
Hence the topological dimension of a coordinate domain equals the geometric
dimension $n$. The main purpose of this class over `AbstractDomain` is to
enforce that the topological dimension equals the geometric dimension.

Every concrete `AbstractCoordinateDomain` $P$ stores:

- The dimension $n$, provided by the property `dimension`.

An `AbstractCoordinateDomain` $P$ derives:

- The geometric dimension and the topological dimension, which are both $n$.
- The ambient domain `ambient_domain(P)`, which is `RD(n)`.
- The coordinate function `coordinate(P)` (see [Coordinate](#coordinate)).
- The parameter domain `parameter_domain(P)`, which is $P$.
- The mapping `mapping(P)`, which is the identity map `coordinate(P)`.

Examples:

- A reference cell in $R^n$ and $R^n$ itself are coordinate domains.
- A facet of a reference cell in $R^n$ is *not* a coordinate domain, because the
  facet has empty interior in $R^n$. The facet is the image of a reference cell
  in $R^{n-1}$ under a map, so the facet is a parameterised domain (defined
  below).

### RD

A concrete subclass of `AbstractCoordinateDomain` for the whole of $R^n$.

`RD(n)` stores the dimension $n$, an `int` with $n \ge 0$. The integer $n$ is
not a base object, so `RD(n)` has no edges in the graph. `RD(m) == RD(n)` if
and only if $m = n$. Everything else comes from `AbstractCoordinateDomain`, so
`ambient_domain(RD(n)) == RD(n)`.

Example:

- `RD(2)` is the parameter domain of the paraboloid example below.

### Coordinate

A concrete subclass of `AbstractExpression` for the coordinate function of a
coordinate domain.

Every concrete coordinate domain $P$ provides a coordinate function
`Coordinate(P)`, obtained with the free function `coordinate(P)`. The
coordinate function $X$ of $P$ is the terminal expression that maps each point
of $P$ to the $n$-tuple of that point:

$X : P \to R^n$, $X(p) = p$.

The coordinate function $X$ is a function on $P$, not a point of $P$.
Evaluating $X$ at a point $p$ of $P$ gives the $n$-tuple of $p$.

`Coordinate(P)` stores the coordinate domain $P$, which is the only edge.
`Coordinate(P)` derives the value shape `(n,)` and the set $P$. Two coordinate
functions are equal if and only if the two coordinate domains are equal.

Only coordinate domains provide a coordinate function. The coordinates of the
points of any other domain are given by `mapping` (see
[AbstractDomain](#abstractdomain)), which is built from the coordinate function
of the parameter domain.

Examples:

- `coordinate(RD(3))` is `Coordinate(RD(3))`, with value shape `(3,)` and the
  set `RD(3)`.
- `coordinate(RD(3)) == coordinate(RD(3))` is `True`, and
  `coordinate(RD(2)) == coordinate(RD(3))` is `False`.
- In the paraboloid example, `X = coordinate(RD(2))` has value shape `(2,)`.

### Tensor

A concrete subclass of `AbstractExpression` that builds a tensor-valued
expression from scalar expressions.

`Tensor(entries)` stores:

- The entries, a nested `tuple` of scalar `AbstractExpression`s. The nesting
  depth $r \ge 1$ is the rank of the tensor.

The nested `tuple` must be rectangular. Every `tuple` at the same depth must
have the same length, and every length must be at least 1. Every entry must
have value shape `()`.

`Tensor(entries)` derives:

- The value shape `(n_0, ..., n_{r-1})`, where $n_i$ is the length of every
  `tuple` at depth $i$.
- The set of the tensor. The set of the tensor is the set of the entries. All
  entries with a set must have the same set. If no entry has a set, the tensor
  has no set.

The edges of `Tensor(entries)` are the entries in row-major order, so the
entry `entries[i_0]...[i_{r-1}]` comes before every later entry in
lexicographic order of the indices $(i_0, \ldots, i_{r-1})$.

Construction performs no simplification. For example,
`vector((Integer(1) + Integer(1),)) != vector((2,))`.

The free function `tensor(entries)` applies `as_expression` to every entry and
returns `Tensor(entries)`. The free functions `vector` and `matrix` are
`tensor` restricted to nesting depth 1 and nesting depth 2, and raise
`ValueError` for any other nesting depth.

Examples:

- `vector((1, 2))` has value shape `(2,)` and no set.
- `matrix(((1, 0), (0, 1)))` has value shape `(2, 2)` and no set. The entries
  are `Integer`s.
- `matrix(((1.0, 0.0), (0.0, 1.0)))` has `RealScalar` entries, so
  `matrix(((1, 0), (0, 1))) != matrix(((1.0, 0.0), (0.0, 1.0)))`.
- `matrix(((1, 2),))` has value shape `(1, 2)`.
- `tensor((((1,),),))` has value shape `(1, 1, 1)`.

### Indexed

>>> LLM generated content.

A concrete subclass of `AbstractExpression` that selects one component of a
tensor-valued expression.

`Indexed(expression, indices)` stores:

- The expression, an `AbstractExpression` of rank $r \ge 1$. The expression is
  the only edge of `Indexed(expression, indices)`.
- The indices, a `tuple` of $r$ `int`s.

Indexing is full indexing only. The number of indices must equal the rank $r$
of the expression, and each index $i_j$ must satisfy $0 \le i_j < n_j$, where
`(n_0, ..., n_{r-1})` is the value shape of the expression.

`Indexed(expression, indices)` derives the value shape `()` and the set of the
expression.

The free function `indexed(expression, indices)` accepts a single `int` $i$,
which the free function converts into the `tuple` `(i,)`, or a `tuple` of
`int`s, and returns `Indexed(expression, indices)`. The overload
`expression[indices]` calls `indexed`, so `X[0]` and `A[0, 1]` are
`indexed(X, 0)` and `indexed(A, (0, 1))`. The free function `indexed` raises:

| Error        | Cause                                                        |
| ------------ | ------------------------------------------------------------ |
| `TypeError`  | An index whose type is not exactly `int`, such as `bool`.    |
| `TypeError`  | An expression with value shape `()`.                         |
| `IndexError` | A number of indices that differs from the rank.              |
| `IndexError` | A negative index, or an index $i_j \ge n_j$.                 |

A negative index is rejected rather than converted, so the stored indices are
always the indices that the user wrote. Expressions are not iterable (see
`AbstractExpression.__iter__`).

Construction performs no simplification. For example, `vector((a, b))[0]` is
`Indexed(Tensor((a, b)), (0,))` and is not equal to `a`.

Partial indexing, slicing and free indices are not part of `Indexed`.

Examples:

- For the coordinate function $X$ of `RD(3)`, `X[2]` has value shape `()`.
- For `A = matrix(((1, 2), (3, 4)))`, `A[1, 0]` has value shape `()`, and
  `A[1]` raises `IndexError`.
- `X[-1]` and `X[True]` raise errors.
- For the coordinate function $X$ of `RD(2)`, `Tensor((X[0], X[1])) != X`,
  because construction performs no simplification.
- `vector((1 - X[0], X[0]))` has value shape `(2,)` and the set `RD(2)`.
- `matrix(((X[0], X[1]),))` has value shape `(1, 2)`.
- In the paraboloid example, the height `X[0]**2 + X[1]**2` is an expression
  on `RD(2)` with value shape `()`.

### Gradient

> LLM generated content.

A concrete subclass of `UnaryOperator` for the gradient of an expression $f$
on a coordinate domain $P$ of dimension $n$. The gradient of $f$ is the
derivative of $f$ with respect to the coordinate function $X$ of $P$.

`Gradient(f)` stores $f$, which is the only edge. The set of $f$ must be an
`AbstractCoordinateDomain` $P$. `Gradient(f)` derives the value shape
`f.value_shape + (n,)` and the set $P$.

The free function `grad(f)` applies `as_expression` to $f$ and returns
`Gradient(f)`. The free function raises `ValueError` if $f$ has no set, or if
the set of $f$ is not a coordinate domain.

`Gradient` is lazy. Construction does not evaluate the derivative. A later pass
evaluates the derivative with the standard rules, including `grad(X)` $= I$
and a zero derivative for every literal.

The gradient on a parameterised domain is not a primitive. The gradient on a
parameterised domain is a composite of `grad` on the parameter domain and the
geometry of the parameterised domain (see
[Geometry of a domain](#geometry-of-a-domain)).

Examples, with `X = coordinate(RD(2))`:

- `grad(X)` has value shape `(2, 2)`. A later pass evaluates `grad(X)` to the
  identity.
- `grad(X[0]**2 + X[1]**2)` has value shape `(2,)`. A later pass evaluates the
  gradient to $(2 X_0, 2 X_1)$.
- `grad(1)` raises `ValueError`, because the literal `Integer(1)` has no set.

### ParameterisedDomain

A concrete subclass of `AbstractDomain` for parameterised domains.

A parameterised domain $D$ is the image $\psi(P)$ of a coordinate domain
$P \subset R^k$ under a map $\psi : P \to R^n$.

`ParameterisedDomain(P, psi)` stores:

- The parameter domain $P$, an `AbstractCoordinateDomain` of dimension $k$,
  returned by `parameter_domain(D)`.
- The mapping $\psi$, an `AbstractExpression` on $P$ with value shape `(n,)`,
  returned by `mapping(D)`. The mapping $\psi$ is built from the coordinate
  function $X$ of $P$.

$P$ and $\psi$ are the two edges, in that order.

`ParameterisedDomain(P, psi)` derives:

- The topological dimension $k$, which is the dimension of $P$.
- The geometric dimension $n$, which is taken from the value shape `(n,)` of
  $\psi$. The geometric dimension is derived, not stored, so the geometric
  dimension cannot disagree with the value shape of $\psi$.
- The ambient domain `ambient_domain(D)`, which is `RD(n)`.

A parameterised domain does not provide a coordinate function.

The free function `parameterised_domain(P, psi)` applies `as_expression` to
$\psi$ and returns `ParameterisedDomain(P, psi)`. The free function raises
`ValueError` if $\psi$ does not have rank 1, if $\psi$ has a set other than
$P$, or if $k > n$.

A coordinate domain is *not* a parameterised domain. `AbstractCoordinateDomain`
and `ParameterisedDomain` are both subclasses of `AbstractDomain`. Both classes
provide `parameter_domain` and `mapping`.

Examples:

- The hypotenuse of the reference triangle is a parameterised domain with
  $k = 1$ and $n = 2$. $P$ is the reference interval in $R^1$, and
  $\psi(X) = (1 - X_0, X_0)$.
- `parameterised_domain(RD(2), vector((X[0], X[1], X[0]**2 + X[1]**2)))`, with
  `X = coordinate(RD(2))`, is the paraboloid, with $k = 2$ and $n = 3$.
  `ambient_domain` of the paraboloid is `RD(3)`.
- A physical triangle in $R^2$ is a parameterised domain with $k = n = 2$. $P$
  is a reference triangle, and $\psi$ is the affine map from $P$ onto the
  physical triangle. See the affinely mapped triangle example.

### Geometry of a domain

The geometry of every domain $D$, coordinate or parameterised, comes from the
mapping $\psi$ = `mapping(D)` on the parameter domain $P$ =
`parameter_domain(D)` of dimension $k$:

| Free function        | Definition                       | Value shape |
| -------------------- | -------------------------------- | ----------- |
| `jacobian(D)`        | $J$ = `grad(mapping(D))`         | `(n, k)`    |
| `metric(D)`          | $g$ = `transpose(J) @ J`         | `(k, k)`    |
| `measure_density(D)` | `square_root(determinant(g))`    | `()`        |

The Jacobian, metric and measure density are lazy expressions. Construction
performs no simplification.

For a coordinate domain, the mapping is the coordinate function $X$. The
general rules give $J$ = `grad(X)` $= I$, $g = I$ and a measure density of 1,
so a coordinate domain needs no special case.

### Example: a paraboloid

The paraboloid over `RD(2)` collects the examples of the sections above.

```python
# AbstractCoordinateDomain: RD(2), no edges.
P = RD(2)

# AbstractExpression: Coordinate(RD(2)), value shape (2,).
X = coordinate(P)

# AbstractExpression: Tensor((Indexed(X, (0,)), Indexed(X, (1,)), Sum(...))),
# value shape (3,).
psi = vector((X[0], X[1], X[0]**2 + X[1]**2))

# AbstractDomain: ParameterisedDomain(RD(2), psi), k = 2, n = 3.
D = parameterised_domain(P, psi)

# AbstractExpression: Gradient(psi), value shape (3, 2).
J = jacobian(D)

# AbstractExpression: MatrixProduct(Transpose(J), J), value shape (2, 2).
g = metric(D)

# AbstractExpression: SquareRoot(Determinant(g)), value shape ().
rho = measure_density(D)
```

Every object in the example is lazy.

Note that at this stage we do not need the concept of a `Point` as we have not
evaluated anything at a `Point` so far.

### Example: an affinely mapped triangle

The physical triangle with vertices $v_0 = (1, 0)$, $v_1 = (3, 1)$ and
$v_2 = (2, 3)$ is the image of a reference triangle under an affine map.

Important point is that there is no topology here, just geometry.

$\{p \in R^n : p_i \ge 0 \text{ for every } i, \ \sum_i p_i \le 1\}$

as a concrete `AbstractCoordinateDomain`. A concrete coordinate domain provides
only the dimension $n$. Everything else comes from `AbstractCoordinateDomain`.

```python
class UnitSimplex(AbstractCoordinateDomain):
    def __init__(self, n: int):
        if n < 0:
            raise ValueError("Dimension must be non-negative.")
        super().__init__(n)
        self._n = n

    @property
    def dimension(self) -> int:
        return self._n
```

`UnitSimplex(2)` is the reference triangle with vertices $(0, 0)$, $(1, 0)$ and
$(0, 1)$. The affine map from `UnitSimplex(2)` onto the physical triangle is

$$\psi(X) = v_0 + X_0 (v_1 - v_0) + X_1 (v_2 - v_0)
= (1 + 2 X_0 + X_1, X_0 + 3 X_1).$$

```python
# AbstractCoordinateDomain: UnitSimplex(2), no edges.
T = UnitSimplex(2)

# AbstractExpression: Coordinate(UnitSimplex(2)), value shape (2,).
X = coordinate(T)

# AbstractExpression: Tensor(..., ...), value shape (2,).
psi = vector((1 + 2 * X[0] + X[1], X[0] + 3 * X[1]))

# AbstractDomain: ParameterisedDomain(UnitSimplex(2), psi), k = n = 2.
K = parameterised_domain(T, psi)

# AbstractExpression: Gradient(psi), value shape (2, 2).
J = jacobian(K)

# AbstractExpression: SquareRoot(Determinant(MatrixProduct(Transpose(J), J))),
# value shape ().
rho = measure_density(K)
```

A later pass evaluates the Jacobian to the constant matrix

$$J = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix},$$

so $\det J = 5$ and $\sqrt{\det g} = |\det J| = 5$. The area of the physical
triangle is $5$ times the area $1/2$ of the reference triangle.

The physical triangle is a parameterised domain with $k = n$, not a coordinate
domain, although the physical triangle has non-empty interior in $R^2$. The
class records how the domain is built.

### Cells

A domain can be divided into finitely many cells. A cell is a piece of the
domain with the same dimension as the domain. Two cells meet only on the
boundaries of the two cells.

A facet is a piece of the boundary of a cell, one dimension lower than the
cell. A facet in the boundary of one cell is an exterior facet. A facet shared
by two cells is an interior facet.

Cells and facets are topological. Defining cells and facets does not require a
metric.

# JSH worked to here.

### Finite elements

A finite element describes the local finite-dimensional field used on a cell.

At minimum, an element is defined by:

- the cell type,
- the element family,
- the polynomial degree or equivalent order information,
- the value shape,
- and, for mapped elements, the required pullback or mapping behaviour.

Examples:

```python
P1 = FiniteElement(family="Lagrange", cell=triangle, degree=1, shape=())
vP2 = FiniteElement(family="Lagrange", cell=triangle, degree=2, shape=(2,))
```

The element describes what kind of discrete field is being represented. It does
not store actual coefficient values.

### Domains

A domain is the integration domain of a form. For a minimal language
definition, a domain only needs to provide:

- topological dimension,
- geometric dimension,
- and a coordinate element, if geometry is represented using a finite element.

The actual mesh is still external to UFLx.

### Function spaces

A function space defines a finite-dimensional space of functions. There are
three types of function spaces:

1. standard finite-element spaces, defined by a domain and an element,
1. constant function spaces, defined by a shape and a scalar type,
1. non-finite-element function spaces, defined by domain, shape and scalar
   type, but no element.

Examples:

```python
V = FunctionSpace(domain=domain, element=element)  # standard finite-element space
C = FunctionSpace(shape=(2,), scalar_type="real")  # constant function space
N = FunctionSpace(
    domain=domain, shape=(2,), scalar_type="real"
)  # non-finite-element function space
```

In older UFL examples, `Argument` and `Coefficient` are sometimes created
directly from an element. This is not allowed.

### Terminals

Terminals are the leaves of UFLx expression trees.

The essential terminals are:

- `Argument`, representing an arbitrary basis function of a function space,

  It has a `number` attribute to represent the position of the `Argument` in
  the form. For example, consider a bilinear form `a(u, v)`. The `Argument`
  representing `u` has `number=1`, and the `Argument` representing `v` has
  `number=0`.

  > Experimental (API breaking) note: Existence of an `Argument` and a
  > `Coefficient` is rooted deep into the legacy UFL. However, more useful
  > concept would be to only allow for `Function` terminals, which can
  > represent both known and unknown functions. It would be a responsibility of
  > the user, when creating a form, to define which of the `Function` terminals
  > are to be considered as "arguments" and which not. This would be much
  > closer to a mathematical notation and would solve many user error of having
  > to treat `Argument` and `Coefficient` differently.

- `Coefficient`, representing a known function in a function space,

- literals such as numbers,

- and geometric quantities such as coordinates, normals, Jacobians, cell
  diameter, facet area, etc.

Example:

```python
v = Argument(function_space=V, name="v", number=0)
u = Argument(function_space=V, name="u", number=1)
f = Coefficient(function_space=V, name="f")
```

> Legacy note: There is no `Constant` in UFLx as a means to represent literal
> value which is not known at compile time.

UFLx supports named symbolic terminals. This allows for better error messages
and better postprocessing of expressions. In addition, it improves the
stability of the expression tree.

> Legacy note: In legacy UFL, symbolic objects do not have names. It requires
> expensive and error prone renumbering of terminals in expressions before
> lowering and code generation.

### Expressions

An expression is either a terminal or an operator applied to other expressions.

Examples of operators are addition, multiplication, indexing, `grad`, `div`,
`dot`, `inner`, `jump`, and `avg`.

Expressions are symbolic. They are not evaluated when they are created.

A minimal expression carries these static properties:

- tensor shape,
- free indices,
- associated domain or domains,
- scalar type (real or complex),
- continuity information,
- (possibly) physical dimension.

These properties are enough to reject many meaningless expressions before
lowering or code generation.

### Measures

A measure describes where and how an integrand is integrated.

TODO.

The standard measures are:

- `dx` for cell integrals,
- `ds` for exterior-facet integrals,
- `dS` for interior-facet integrals.

Examples:

```python
dx
ds
dS
```

### Integrals

An integral is made by multiplying a scalar expression by a measure.

Examples:

```python
v * f * dx
inner(grad(v), grad(u)) * dx
jump(v) * dS
```

A valid integral has a scalar-valued integrand with no unresolved free indices.

### Forms

A form is a multi-linear functional
$a(u, v, w, \ldots) : U \times V \times W \times \ldots \to \mathbb{R}$ or
$\mathbb{C}$.

Arguments determine the arity of a form.

> Legacy note: In contrast to the legacy UFL, in UFLx a form does not need to
> be defined as a sum of integrals. For example, the following is a valid UFLx
> form: `a = (u * dx) * (v * dx)`, or
> `a = PointEvaluation(u, (0.5, 0.5)) * PointEvaluation(v, (0.5, 0.5))`.

Examples:

```python
a = inner(grad(v), grad(u)) * dx
L = v * f * dx
F = a - L
```

> Experimental (API breaking) note: Similar to the `Argument` and `Coefficient`
> terminals, the `Form` in the legacy UFL is defined implicitly, by the
> presence of `Argument` terminals. In UFLx, a `Form` could be an explicit
> object that can be created by the user, e.g.

```python
a = Form(inner(grad(v), grad(u)) * dx, arguments=[u, v])
```

## 3. Minimal validity rules

The language is small, but it is typed.

- Addition and subtraction require compatible shapes and compatible free
  indices.
- Multiplication follows scalar, tensor, matrix-vector, matrix-matrix, and
  indexed-expression rules.
- Tensor contractions such as `dot` and `inner` require compatible tensor
  shapes.
- Spatial derivatives require a spatial domain.
- Facet quantities such as `FacetNormal` are only meaningful in facet
  integration contexts.
- Interior-facet expressions may use side restrictions:
  ```python
  u("+")
  u("-")
  ```
- Operators such as `jump` and `avg` are shorthand for combinations of such
  restrictions.
- An integral is valid only when the final integrand is scalar and has no free
  indices.
- A form is valid when all of its integrals are valid and their domains,
  arguments, and coefficients are mutually consistent.

## 4. What transformation procedures does the UFLx core provide

There are two key transformation procedures that the UFLx core provides:

1. Automatic/symbolic differentiation.

   Similar to UFL, UFLx provides a symbolic differentiation via a directional
   (Gateaux) derivative. The derivative of an Expression, Integral or a Form is
   lazy/eagerly evaluated. That means the derivative is stored as an expression
   tree operator, and the actual differentiation is performed when the
   derivative is used in a context where it needs to be evaluated, such as
   lowering or code generation.

   Example:

   ```python
    a = inner(grad(v), grad(u)) * dx
    da = derivative(a, u, du)
   ```

1. Expression and integral transformation between different configurations.

   > Legacy note: In legacy UFL, there is a very central perspective on the
   > configuration of where objects are defined. Usually, all expressions are
   > defined in the physical configuration, and UFL applies pullbacks to the
   > entire expression tree.

   In UFLx, configuration where expressions are defined are a consequence of
   the domain. Since the domain is part of the Expression's static attributes,
   it is propagated through the expression tree based on the language rules.

   TODO.

## 5. What belongs to the language

The minimal language includes:

- symbolic finite-element objects,
- symbolic domains/function spaces,
- symbolic terminals,
- tensor-valued expressions,
- index notation,
- differential operators,
- DG side restrictions,
- integration measures,
- integrals,
- and forms.

This is enough to describe ordinary weak forms such as mass, stiffness,
elasticity, Stokes, mixed Poisson, and DG forms.

## 6. What does not belong to the minimal language

The following are not primitive language concepts:

- actual mesh connectivity,
- actual coefficient vectors,
- basis tabulation,
- quadrature-loop generation,
- form compilation,
- assembly,
- solver setup,
- boundary-condition application,
- expression simplification strategy,
- code generation.

These are consumers or transformations of the symbolic language, not the core
language itself.

## 7. Extensibility of the attributes

> Legacy note: In legacy UFL, the set of static attributes is fixed and
> hardcoded. This led to many language features being hacked in (e.g. complex
> numbers support, dual spaces, etc.).

In UFLx, the set of static attributes is extensible. This allows for cleaner
language extensions in the future, such as support for physical dimensions,
dual spaces, error estimates, custom quadrature rules, custom
integrability/continuity rules, etc.

## 8. Examples

### Poisson equation

```python
coordinate_element = FiniteElement("Lagrange", triangle, 1, shape=(2,))
element = FiniteElement("Lagrange", triangle, 1)

mesh = Mesh(coordinate_element)
V = FunctionSpace(mesh, element)

u = TrialFunction(V, name="u")
v = TestFunction(V, name="v")
f = Coefficient(V, name="f")

a = inner(grad(v), grad(u)) * dx
L = v * f * dx
```

# Expressions

## Expressions and static type

### `AbstractExpression`

| `AbstractExpression` |                                                                                      |
| -------------------- | ------------------------------------------------------------------------------------ |
| Extends              | `Base`                                                                               |
| Provides             | `value_shape`, `free_indices`, `scalar_type`, `domain`, `argument_spaces` |

These five derived properties form the static type (see
[Argument spaces](forms.md#argument-spaces)). The value shape is a `tuple` of
[extents](#extents), one per axis; `()` denotes a scalar. The array rank,
written simply as rank, is the number of axes.

### Extents

An extent is the length of an axis of a value shape, or the extent of a free
index. An extent is an `int` $n \ge 0$ or an `AbstractExtent`.

`AbstractExtent` represents an extent whose value is not known during
construction. Its supplied value must be a nonnegative integer.
An `AbstractExtent` is not an expression. A concrete
`AbstractExtent` is a subclass outside the core, such as `RuntimeExtent` in the
[cut cell example](../examples/integration.md#integrating-over-a-cut-cell).

| `AbstractExtent` |                                      |
| ---------------- | ------------------------------------ |
| Extends          | `Base`                               |
| Defined by       | A subclass outside the core          |

An `AbstractExtent` has no arithmetic and no ordering, so the rules of the
core, which compare extents by equality, accept an `AbstractExtent` unchanged.
A concrete `AbstractExtent` stores the arguments that identify the extent, such
as a label.

Example:

- With an `AbstractExtent` `q`, an expression `points` of value shape `(q, 2)`
  and an expression `weights` of value shape `(q,)` share the extent `q`, so
  `weights[i] * points[i, 0]` is accepted. With a second `AbstractExtent` `r`,
  an expression `weights` of value shape `(r,)` raises `ValueError` in
  `weights[i] * points[i, 0]`, because the index `i` has two extents.

```{dropdown} Open: ragged indexing
An axis whose length depends on another index, such as the quadrature points
of each cut cell, cannot be written. A dependent extent, an `AbstractExtent`
whose label contains the other index, would cover the case. An index value
read from a table is a [`GatheredIndex`](#gatheredindex).
```

### `Index`

`Index` is a free index, not an expression. An `Index` names an axis of a
tensor, as `i` and `j` in `A[i, j]` (see [Index notation](#index-notation)).

| `Index(label)` |                                                                                                 |
| -------------- | ----------------------------------------------------------------------------------------------- |
| Extends        | `Base` |
| Free function  | `index(label)`                                                                                  |
| Takes          | `label: str \| tuple` |
| Properties     | `label`, a `str` for a user index, or a `tuple` that contains graph objects for a derived index |
| Raises         | `TypeError` if the label is neither a `str` nor a `tuple`                                       |

An `Index` takes its extent from the indexed axis; it has no extent of its own.

A user index has a `str` label, which the user chooses. Two free occurrences
of one `Index` in one scope are the same index, so an index that the user
means to be new needs a new name. The scope of an index is given in
[`ComponentTensor` and `IndexSum`](#componenttensor-and-indexsum).

A derived index has a `tuple` label, which a free function or a pass builds
from the operands, so that the index is fresh without a global counter (see
[Derived indices](#derived-indices)).

### Free indices

Every `AbstractExpression` derives the free indices, a `tuple` of pairs. Each
pair is an `Index` and an extent. The order is the order of first occurrence, never the
order of a `set`.

Unless a node class states otherwise, the free indices of a node are derived
from the edges of the node that are expressions:

- The free indices of the edges are taken in edge order, with each `Index`
  taken once, at the first occurrence.
- An `Index` that occurs in two edges must have the same extent in both edges,
  or construction raises `ValueError`.

Summation is never implicit. `a[i] * b[i]` has the free index `i` and denotes
the component-wise product. The sum over `i` is `index_sum(a[i] * b[i], (i,))`.
The operands of a node may have different free indices, so `a[i] + b[j]` has
the free indices `i` and `j` and denotes the outer sum.

### Scalar type

The scalar type is real or complex. Unless the [static type
table](reference.md#appendix-static-type-table) states otherwise, an expression is complex
if and only if at least one edge is a complex expression.

The scalar type is the arithmetic in which the expression is evaluated, not a
statement about the value. A real `SquareRoot`, `Logarithm` or `Power` whose
value is not a real number, such as `square_root(Integer(-1))`, is undefined.

### Domain

An expression's `domain` is where it is defined: an
[`AbstractDomain`](geometry.md#abstractdomain), or an [`AbstractTopology`](geometry.md#abstracttopology)
whose cell types are reference cells. It is `None` for an expression independent
of a domain, such as a literal.

Unless the [static type table](reference.md#appendix-static-type-table) states otherwise,
a node inherits the domain of its expression operands with a domain. All such
operands must have the same domain, or construction raises `ValueError`. A node
with no such operand has domain `None`.

The property always describes the resulting expression. An `Integral` stores
its region as `integration_domain` but has domain `None`; a `PullBack` stores
its `map` and has that map's `parameter_domain`. A `PushForward` stores its
`mesh`, which is also its domain. For `Coordinate`, `Gradient`, `Divergence`
and `Curl`, the stored domain is also the result's domain.

### Expressions as functions

An expression denotes an array-valued function on the domain of the
expression. An expression on a topology denotes a function on the reference
cell of each cell of the topology (see [`CellMap`](geometry.md#cellmap)). An expression
with no domain denotes a single array. An expression with argument spaces other than `()` denotes a multilinear map from
the argument spaces to such functions (see [Variational
forms](forms.md#variational-forms)).

### Node class or free function

An operation uses a construction free function when it has a fixed
composition under the [construction policy](rules.md#construction-policy). A node
class is needed when choosing the result requires expression structure, when
the compiler must choose a numerical algorithm, or when a pass must retain
the operation as written. In the last case, the node's section states its
composition; `Integral` over a mesh is an example.

For example, the transpose is always `component_tensor(A[k, l], (l, k))`,
so it is a free function. `Inverse` and `Determinant` keep the algorithm
choice for the compiler.

A construction free function is part of the core only if a definition or a
rule of the core uses it, such as `density` in [`Integral` over a
mesh](geometry.md#integral-over-a-mesh) and `cell_map` in
[Restriction](geometry.md#restriction), or if a Python operator calls the free function
as an overload, such as `matrix_product` for `@` (see [Operator
overloads](rules.md#operator-overloads)). `transpose` is in the core because the rule
for `Inverse` of [Applying derivatives](evaluation.md#applying-derivatives) uses it. Every
other construction free function belongs to a library, such as the operations
of [Tensor algebra](../examples/tensors.md#tensor-algebra).

A node class is part of the core only if a definition or a rule of the core
needs the node class, or if the node class is an operation of the scalar
field, such as conjugation. A rule for the node class itself does not count.
Otherwise the node class belongs to a
[language extension](extensions.md#language-extensions). For example, the core writes a
rational number as the `Division` of two `Integer`s, so `Division` is part of
the core.

The [node-class rationale](reference.md#appendix-node-class-rationale) explains the
individual choices.

## Numbers and scalar operations

### Literals

A literal is a number written in an expression, such as `Integer(2) * u`.
`Integer` and `Float` are real terminals with value shape `()`, no free
indices, no domain and argument spaces `()`.

| `Integer(v)` |                       |
| ------------ | --------------------- |
| Extends      | `AbstractExpression` |
| Properties   | `value` $v$, an `int` |

| `Float(v)` |                               |
| ---------- | ----------------------------- |
| Extends    | `AbstractExpression`         |
| Properties | `value` $v$, a finite `float` |

`Integer` and `Float` accept only exact Python types, as in rule 1 of
[the base object](rules.md#the-base-object-base). So `Integer(True)` raises
`TypeError`, and `Float(-0.0) != Float(0.0)`.

A rational number is the `Division` of two `Integer`s, such as
`Integer(1) / Integer(3)`. The compiler chooses its numerical evaluation;
the core has no rational literal.

### `Complex`

`Complex` is a scalar complex constant built from two real literals. It is a
terminal, and its value is
`real` $+\,i\,$ `imag`.

| `Complex(real, imag)` | |
| --------------------- | - |
| Extends               | `AbstractExpression` |
| Takes                 | `real, imag: Integer \| Float` |
| Properties            | `real` and `imag`, the real and imaginary parts |
| Static type           | Value shape `()`, complex, no free indices, no domain and argument spaces `()` |
| Raises                | `TypeError` for a part of another class |

```python
Complex(Integer(0), Integer(1))  # The imaginary unit.
Complex(Integer(-1), Integer(0))  # Minus one with complex scalar type.
```

`Complex` retains its parts and stays complex even when the imaginary part
is zero. As a terminal, it has no evaluation operands (see
[Evaluation traversal](evaluation.md#evaluation-traversal)).

General real expressions `a` and `b` form a complex expression through
`a + Complex(Integer(0), Integer(1)) * b`. A real expression `r` is embedded
in complex arithmetic by `Complex(Integer(1), Integer(0)) * r`, preserving
its value shape, free indices, domain and argument spaces.

### `Zero`

`Zero` represents zero with a specified value shape, free indices, argument
spaces and scalar type (real or complex).

| `Zero(value_shape, free_indices, argument_spaces, scalar_type)` |                                                                                                                                |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| Extends                                            | `AbstractExpression` |
| Properties                                         | `value_shape`, `free_indices`, `argument_spaces` and `scalar_type`, as in [Extents](#extents), [Free indices](#free-indices) [Argument spaces](forms.md#argument-spaces) and [Scalar type](#scalar-type) |
| Static type                                        | The stored properties, with no domain                                                                                         |
| Raises                                             | `TypeError` for an argument of another type, and `ValueError` for a negative extent or an index that occurs twice |

For example, `Zero((3,), (), (), real)` is a real zero vector and
`Zero((), (), (V,), complex)` is a complex zero form on `V`. Here `real` and
`complex` denote the scalar types, not Python numeric constructors.

### Arithmetic

`Sum`, `Product`, `Division` and `Power` denote the sum, product, quotient and
power of two scalar expressions. Both operands must have value shape `()`, or
construction raises `ValueError`; the result has value shape `()`. Free indices
follow the [default rule](#free-indices).

The free functions below construct these nodes for scalar operands. `add`,
`multiply` and `divide` also support [tensor arithmetic](#tensor-arithmetic).

| Class      | Free function | Overload |
| ---------- | ------------- | -------- |
| `Sum`      | `add`         | `+`      |
| `Product`  | `multiply`    | `*`      |
| `Division` | `divide`      | `/`      |
| `Power`    | `power`       | `**`     |

### `Sum`

| `Sum(first, second)` |                                             |
| -------------------- | ------------------------------------------- |
| Extends              | `AbstractExpression` |
| Free function        | [`add`](#tensor-arithmetic), with the overload `+` |
| Takes                | `first, second: AbstractExpression` |
| Properties           | `first` and `second` |
| Requires             | The operands of [Arithmetic](#arithmetic)   |

### `Product`

| `Product(first, second)` |                                                  |
| ------------------------ | ------------------------------------------------ |
| Extends                  | `AbstractExpression` |
| Free function            | [`multiply`](#tensor-arithmetic), with the overload `*` |
| Takes                    | `first, second: AbstractExpression` |
| Properties               | `first` and `second` |
| Requires                 | The operands of [Arithmetic](#arithmetic)        |

### `Power`

| `Power(first, second)` |                                                                             |
| ---------------------- | --------------------------------------------------------------------------- |
| Extends                | `AbstractExpression` |
| Free function          | `power(first, second)`, with the overload `**`                              |
| Takes                  | `first, second: AbstractExpression` |
| Properties             | `first` and `second` |
| Requires               | The operands of [Arithmetic](#arithmetic)                                   |

Complex powers use the principal value, consistent with `Logarithm`.
The exponent may have domain. The [derivative rules](evaluation.md#derivatives-of-the-node-classes)
distinguish constant and variable exponents.

### Negation and subtraction

Negation and subtraction are free functions, not node classes:

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `negate(f)` | `-f` | `f: AbstractExpression` | `multiply(Integer(-1), f)` |
| `subtract(f, g)` | `f - g` | `f, g: AbstractExpression` | `add(f, negate(g))` |

These compositions inherit the tensor
behaviour of `multiply` and `add` (see [Tensor arithmetic](#tensor-arithmetic)).

### `Division`

| `Division(first, second)` |                                                                                        |
| ------------------------- | -------------------------------------------------------------------------------------- |
| Extends                   | `AbstractExpression` |
| Free function             | [`divide`](#tensor-arithmetic), with the overload `/` |
| Takes                     | `first, second: AbstractExpression` |
| Properties                | `first` and `second` |
| Requires                  | The operands of [Arithmetic](#arithmetic)                                              |

### `Absolute`

`Absolute` is the absolute value, or the modulus of a complex scalar.

| `Absolute(f)` |                                           |
| ------------- | ----------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `absolute(f)`, with the overload `abs(f)` |
| Takes         | `f: AbstractExpression` |
| Properties    | `operand` $f$ |
| Requires      | A scalar $f$                              |
| Static type   | Value shape `()`, real                    |

### `Conjugate`

`Conjugate` is the complex conjugate of an expression of any value shape, entry
by entry.

| `Conjugate(f)` |                                                                                             |
| -------------- | ------------------------------------------------------------------------------------------- |
| Extends        | `AbstractExpression` |
| Free function  | `conjugate(f)`                                                                              |
| Takes          | `f: AbstractExpression` |
| Properties     | `operand` $f$ |
| Static type    | The value shape and the free indices of $f$, and the argument spaces of $f$ with each space conjugated (see [`ConjugateSpace`](forms.md#conjugatespace)) |

`Conjugate` is the only scalar-field node that accepts any value shape, because
[`Form`](forms.md#form) binds `conjugate(v)` even for a vector argument. Its entrywise
meaning is `conjugate(f)[e]` = `conjugate(f[e])` in value.

### `RealPart` and `ImaginaryPart`

`RealPart` and `ImaginaryPart` give the real and imaginary parts of a scalar.

| `RealPart(f)`, `ImaginaryPart(f)` |                                       |
| --------------------------------- | ------------------------------------- |
| Extends                           | `AbstractExpression` |
| Free function                     | [`real_part`](#tensor-arithmetic), [`imaginary_part`](#tensor-arithmetic) |
| Takes                             | `f: AbstractExpression` |
| Properties                        | `operand` $f$ |
| Requires                          | A scalar $f$                          |
| Static type                       | Value shape `()`, real                |

These nodes and `Absolute` are the only nodes that take a complex operand to
a real expression. The composition `(f + conjugate(f)) / Integer(2)` has the
value of the real part but remains complex by the [scalar-type rule](#scalar-type).

Their free functions apply entrywise to tensors (see
[Tensor arithmetic](#tensor-arithmetic)).

The operand may be real. Construction does not simplify, so `real_part(f)`
of a real $f$ is not $f$, and `imaginary_part(f)` of a real $f$ is not
`Zero`.

For `z = Complex(Integer(1), Integer(2))`, `real_part(z)` and
`imaginary_part(z)` are real expressions with the values 1 and 2, respectively.

### Elementary functions

`SquareRoot` and `Logarithm` give the square root and natural logarithm of a
scalar. The core provides these two; other layers are expected to provide
a wider range of elementary functions.

| `SquareRoot(f)`, `Logarithm(f)` |                                     |
| ------------------------------- | ----------------------------------- |
| Extends         | `AbstractExpression` |
| Free function   | `square_root(f)` and `logarithm(f)` |
| Takes           | `f: AbstractExpression` |
| Properties      | `operand` $f$ |
| Requires        | A scalar $f$                        |
| Static type     | Value shape `()`                    |

For a complex operand, `Logarithm` is the principal branch, with the imaginary
part in $(-\pi, \pi]$, and `SquareRoot` is the principal square root.

~~~{dropdown} Why `SquareRoot` is a node class
`SquareRoot` is not the composition `power(f, Integer(1) / Integer(2))`,
because IEEE 754 requires a correctly rounded square root and not a correctly
rounded power, so the result of `power` depends on the math library.
~~~

## Tensors and indexing

### `Tensor`

`Tensor` stacks expressions along a new first axis.

| `Tensor(entries)` |                                                                                                           |
| ----------------- | --------------------------------------------------------------------------------------------------------- |
| Extends           | `AbstractExpression` |
| Free function     | `tensor(entries)`                                                                                         |
| Takes             | `entries: tuple[AbstractExpression, ...]` |
| Properties        | `entries`, a non-empty `tuple` of `AbstractExpression`s with the same value shape `s`                     |
| Static type       | Value shape `(n,) + s`, with $n$ the number of entries                                                    |

A matrix is a `Tensor` of `Tensor`s.

Value shapes may contain extent 0, as in `coordinate(RD(0))`. An empty
`Tensor` cannot provide the entry shape `s`; use `Zero` for an array with an
extent 0, which has no entries.

Examples, with `one, two = Integer(1), Integer(2)`:

- `tensor((one, two))` is a vector of value shape `(2,)`.
- `tensor((one, tensor((one, two))))` raises `ValueError`, because the entries
  have different value shapes.

### Index notation

Index notation writes tensor algebra with four node classes. `Tensor` builds a
tensor from scalars, `Indexed` selects one component, and `ComponentTensor`
and `IndexSum` bind the free indices of a scalar, as axes or as a sum. Tensor
rearrangements and contractions are compositions of these, so `transpose`,
`dot` and `inner` are free functions, not node classes (see [Node class or
free function](#node-class-or-free-function)).

With `i, j = index("i"), index("j")` and a matrix `A` of value shape `(2, 3)`:

- `A[i, j]` is a scalar with the free indices `i`, of extent 2, and `j`, of
  extent 3. The scalar denotes one value for each value of `i` and `j`.
- `component_tensor(A[i, j], (j, i))` binds `j` and `i` as axes, in that
  order, and is the transpose of `A`, of value shape `(3, 2)`.
- `index_sum(A[i, j] * A[i, j], (i, j))` sums over `i` and `j`, and is the sum
  of the squares of the entries of `A`, with no free index.

The sections below give the rules. Tuples of ordinary indices, including
empty tuples for scalars, serve tensors of every rank. A
[`GatheredIndex`](#gatheredindex) reads an index value from a table that the
solver gives, such as the cell of each facet. A free function or a pass that
needs a fresh index derives the index from the operands (see [Derived
indices](#derived-indices)).

### `Indexed`

`Indexed` selects one component of a tensor. The overload `[]` calls
[`component`](#partial-indexing), which also supports partial indexing.

| `Indexed(expression, indices)` |                                                                                   |
| ------------------------------ | --------------------------------------------------------------------------------- |
| Extends                        | `AbstractExpression` |
| Free function                  | `indexed(expression, indices)`                                                    |
| Takes                          | `expression: AbstractExpression`, `indices: tuple` |
| Properties                     | `expression`, and `indices`, a non-empty `tuple` of entries of the kinds below    |
| Static type                    | Value shape `()`, and the free indices below                                      |
| Raises                         | `TypeError`, `ValueError` and `IndexError`, below                                 |

The factory `indexed(f, ())` returns `f` if `f` is scalar, and raises
`IndexError` otherwise: full indexing must cover every axis. `Indexed(f, ())`
still raises `ValueError`; the factory creates no node for the scalar case.

Indexing is full indexing. Each entry covers one axis of the
expression, and the entries together cover every axis exactly once. With
`(n_0, ..., n_{r-1})` the value shape of the expression, an entry is one of
the following kinds:

| Entry                                | Covers                                | Rule                                                                                               | Raises       |
| ------------------------------------ | ------------------------------------- | -------------------------------------------------------------------------------------------------- | ------------ |
| An `int` $k$                         | One axis, of extent $n_j$             | $n_j$ is an `int`, and $0 \le k < n_j$                                                             | `IndexError` |
| An [`Index`](#index)                 | One axis, of extent $n_j$             | The `Index` takes the extent $n_j$                                                                 |              |
| A [`GatheredIndex`](#gatheredindex)  | One axis, of extent $n_j$             | The extent of the `GatheredIndex` equals $n_j$                                                     | `IndexError` |

The free indices of the `Indexed` are the free indices of the expression,
followed by each `Index` entry with the extent that the
entry takes, followed by the free indices of each `GatheredIndex` entry, with
each index taken once.

An `Index` may occur in two or more entries, or may already be free in the
expression. All the occurrences are the same index and must have the same
extent, so the indexed axes take the same value. `A[i, i]` is the diagonal
entry of a matrix `A`.

Construction also raises `TypeError` for an entry of another kind,
`ValueError` for no entries or for an index with two extents, in two entries
or in an entry and the expression, and `IndexError` if the entries do not
cover exactly the rank. An expression with value shape `()` has rank 0.

For a matrix `A` of value shape `(2, 2)`, `A[1, 0]` has value shape `()`,
whereas `indexed(A, (1,))` raises `IndexError` because one axis is unindexed.

#### Partial indexing

The overload `[]` calls the free function `component`, which also takes fewer
entries than the rank and indexes the leading axes:

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `component(f, entries)` | `[]` | `f: AbstractExpression`, `entries: tuple` | `f` if `entries == ()`; otherwise `component_tensor(indexed(f, entries + K), K)` if the entries are fewer than the rank of $f$, otherwise `indexed(f, entries)` |

`K` is the `tuple` of the derived indices
`Index(("uflx.component", f, entries, m))`, one for each axis that the entries
do not cover (see [Derived indices](#derived-indices)). `component` reads
the rank of $f$ and the number of entries.
In particular, `component(f, ())`, also written `f[()]`, returns `f` at any rank.
`Indexed` validates the entries, so too many entries raise `IndexError`.

A free function of the core calls `indexed`, never `component`, so a
composition in this document that writes `f[K]` is `indexed(f, K)`, and an
operand of the wrong rank raises from `Indexed`.

For the same matrix `A`, `A[1]` is row 1, of value shape `(2,)`. With `k`
the derived index, its construction is:

```python
entry = indexed(A, (1, k))
row = component_tensor(entry, (k,))
```

### `ComponentTensor` and `IndexSum`

`ComponentTensor` binds free indices as tensor axes; `IndexSum` sums over them.

If any summed index has extent zero, `IndexSum` has value zero with the
result's static type, without evaluating its summand. This also applies when
a symbolic extent is supplied as zero; it is an evaluation rule, not a
construction simplification.

| `ComponentTensor(expression, indices)`, `IndexSum(expression, indices)` |                                                                                                                                                                                                                   |
| ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                                                                 | `AbstractExpression` |
| Free function                                                           | `component_tensor(expression, indices)`, `index_sum(expression, indices)`                                                                                                                                         |
| Takes                                                                   | `expression: AbstractExpression`, `indices: tuple[Index, ...]` |
| Properties                                                              | `expression`, a scalar, and `indices`, a non-empty `tuple` of distinct `Index` entries free in the expression                                                                                 |
| Static type                                                             | The free indices of the expression without the indices. `ComponentTensor` has the value shape of the extents of the indices in order, and `IndexSum` the value shape `()` |

The factories `component_tensor(f, ())` and `index_sum(f, ())` return the
scalar expression `f` directly, preserving its entire static type. They raise
`ValueError` for a tensor-valued `f`. The node constructors still require
nonempty index tuples. Summing over no indices is the identity; summing over
an index of extent zero gives zero.

Scope is lexical. An index occurrence refers to the innermost enclosing
`ComponentTensor` or `IndexSum` that binds it; otherwise it is free. A bound
index is not free in the binding node. Nested binders may bind the same index.
For example, with a vector `a`:

```python
i = index("i")
p = component_tensor(a[i], (i,))
entry = p[i]
```

`p` has no free index. In `entry`, the outer `i` is free and does not refer
to the bound `i` inside `p`.

A pass that substitutes into the expression of a binder renames the bound
indices that would capture a free index of the substituted expression (see
[Substitution and binding](evaluation.md#substitution-and-binding)).

For a square matrix `A`, `index_sum(A[i, i], (i,))` is its trace, and
`component_tensor(A[i, i], (i,))` is its diagonal. Both work when the extent is
an `AbstractExtent`, without knowing its value during construction.

### `AbstractIndexTable`

`AbstractIndexTable` represents a table of index values supplied by the solver,
such as the cell on one side of each facet, or the local facet number in a facet
integral. The values are fixed but not known during construction. An
`AbstractIndexTable` is not an expression and not a free index.

| `AbstractIndexTable` |                                                        |
| -------------------- | ------------------------------------------------------ |
| Extends              | `Base` |
| Provides             | `extent`, an extent, and `shape`, a `tuple` of extents |
| Defined by           | A subclass outside the core                            |

The table has the shape `shape`, and each value is an `int` $v$ with
$0 \le v <$ `extent`. An `AbstractIndexTable` with the shape `()` is one
value. An expression reads the table only through a
[`GatheredIndex`](#gatheredindex), so one value is read as
`gathered_index(v, ())`.

### `GatheredIndex`

`GatheredIndex` reads an index value from a table. It is neither an expression
nor a free index; its table indices determine its free indices.

| `GatheredIndex(table, indices)` |                                                                                                                                    |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| Extends                         | `Base` |
| Free function                   | `gathered_index(table, indices)`                                                                                                   |
| Takes                           | `table: AbstractIndexTable`, `indices: tuple` |
| Properties                      | `table`, an `AbstractIndexTable`, and `indices`, an `int`, an `Index` or a `GatheredIndex` for each axis of the shape of the table |
| Derived                         | `extent` and `free_indices`, below |
| Raises                          | `TypeError` for an invalid table or entry type, `IndexError` for an incorrect number of entries or invalid bounds or gather extent, and `ValueError` for inconsistent extents of a repeated index |

`indices` contains one `int`, `Index` or `GatheredIndex` per table axis, with
the same bounds and extent checks as [`Indexed`](#indexed). A scalar table
therefore takes `()`. Both constructor and factory follow this rule.
Gathers may nest. The result has the table's
extent, and its free indices are the `Index` entries with their axis extents,
followed by the free indices of the `GatheredIndex` entries, taking each index
once.

#### Examples

The examples use the following objects:

- A table `cells` of the shape `(4,)`, the extent 2 and the values
  `(0, 1, 1, 0)`, such as the cell of each of the four boundary edges of the
  square divided into two triangles.
- An array `a` of value shape `(2,)`, with one value for each triangle.

```python
e = index("e")
g = gathered_index(cells, (e,))

# The value of a on the triangle of each edge, with the free index e of extent 4.
a_e = a[g]

# The array (a[0], a[1], a[1], a[0]) of value shape (4,).
b = component_tensor(a_e, (e,))
```

- Let `cell_vertices` be a table of the vertices of each triangle, of the shape
  `(2, 3)` and the extent 4. `gathered_index(cell_vertices, (g, w))` is the
  vertex `w` of the triangle of edge `e`, with the free indices `e` and `w`.
- For an array `c` of value shape `(3,)`, `c[g]` raises `IndexError`, because
  the extent 2 of `g` differs from the extent 3 of the axis.

### Derived indices

A free function that needs fresh indices derives the indices from the
operands. The label of a derived index is a `tuple` of the identifier of the
free function, the arguments of the free function and a number. The
identifier is the public name of the free function, qualified by the package
name, such as `"uflx.transpose"` for the core or `"library.outer"` for a
library. For example, `transpose(A)` (see [Tensor
arithmetic](#tensor-arithmetic)) uses `Index(("uflx.transpose", A, 0))` and
`Index(("uflx.transpose", A, 1))`.
A derived index is fresh with respect to the operands, so the free function
does not inspect the operands to find a fresh index. A global counter is not
allowed, because a counter makes construction depend on global state.

A pass that needs fresh indices derives the indices in the same way, with the
identifier of the pass, the arguments of the rewritten node and a number. The
renamed index of [Substitution and binding](evaluation.md#substitution-and-binding) has the
label `("uflx.rename", i, body, (s_1, ..., s_n))`.

A label refers only to objects that the index is fresh for, and never to an
object that contains the index. A replacement rebuilds a derived index whose
label contains a replaced object, so that the replacement commutes with the
free functions. A renamed name keeps its label. The substitution that renames
a name continues into the operands that contain the new name, and rebuilding
the label there would substitute into $s_1, \dots, s_n$ again, which may
contain the binding node, so the substitution would not terminate.

The edges of an index to the objects in the label of the index identify the
index and are not operands (see [Evaluation traversal](evaluation.md#evaluation-traversal)).

Label contents do not introduce free occurrences of functions or indices.
This also applies to function labels: a function whose label mentions `u`
is not thereby a function of `u`. Labels still determine identity and are
rebuilt according to the substitution rules above.

A free function or a pass derives one `Index` per axis.
`component_tensor(e, K)[J]` is $e$ with each index of `K` replaced by the
corresponding entry of `J`.

## Tensor operations

### Tensor arithmetic

Tensor arithmetic uses fixed compositions of scalar nodes and index notation.
The free functions read operand ranks to choose the compositions below.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `add(first, second)` | `+` | `first, second: AbstractExpression` | `component_tensor(Sum(first[K], second[K]), K)` |
| `multiply(first, second)` | `*` | `first, second: AbstractExpression` | `component_tensor(Product(first, second[K]), K)` if `first` is scalar, `component_tensor(Product(first[K], second), K)` if only `second` is scalar, otherwise `Product(first, second)` |
| `divide(first, second)` | `/` | `first, second: AbstractExpression` | `component_tensor(Division(first[K], second), K)` |
| `real_part(f)` | | `f: AbstractExpression` | `component_tensor(RealPart(f[K]), K)` |
| `imaginary_part(f)` | | `f: AbstractExpression` | `component_tensor(ImaginaryPart(f[K]), K)` |
| `matrix_product(first, second)` | `@` | `first, second: AbstractExpression` | The contraction of the last axis of `first` with the first axis of `second`, below |
| `transpose(A)` | | `A: AbstractExpression` | `component_tensor(A[k, l], (l, k))`, with `(k, l)` = `K`, for $A$ of rank 2 |

In each composition, `K` is the `tuple` of the derived indices
`Index((name, first, second, m))`, or `Index((name, f, m))`, for
$m = 0, \dots, n - 1$, with `name` the identifier of the free function, such
as `"uflx.add"` (see [Derived indices](#derived-indices)). For addition and
division, $n$ is the rank of `first`; for real and imaginary parts, the rank
of `f`; for multiplication, the rank of the operand indexed in the selected
composition. Thus `K == ()` for scalar operands, and the empty-index factory
identities give the scalar operation directly. Explicit scalar constructors
inside the compositions avoid recursive calls to the tensor factories. The
indices of `K` are fresh with respect to the operands, so the free indices
of the operands remain free in the result.

The free functions read only the ranks; constructors reject incompatible
operands (see [The base object `Base`](rules.md#the-base-object-base)). There is no
broadcasting, and no product of two tensors other than `@`:

- `A + B` for ranks that differ raises `IndexError` from `Indexed`, and for
  extents that differ raises `ValueError` from `Sum`, because an index of `K`
  has two extents.
- `a + T` for a scalar `a` and a tensor `T` also raises `IndexError`: there is no scalar broadcasting.
- `A * B` for two tensors raises `ValueError` from `Product`, and `a / T`
  raises `ValueError` from `Division`.

A scalar operation's free function accepts tensors only when the operation has
no competing tensor meaning. Thus `conjugate`, `real_part` and `imaginary_part`
apply entrywise. The following remain scalar and raise `ValueError` for tensors:

| Free function | Second meaning on a tensor |
| ------------- | -------------------------- |
| `power`, `**` | The matrix power |
| `absolute`, `abs` | The norm of a vector |
| `square_root` | The matrix square root |
| `logarithm` | The matrix logarithm |

An elementary function added later, such as the exponential, stays scalar for
the same reason, because the matrix exponential is a second meaning. The
entry-by-entry application of a scalar operation `op` is written as
`component_tensor(op(T[K]), K)`, with one index in `K` for each axis of `T`. This
also applies to scalar `T`, with `K == ()`. A full contraction similarly uses
`index_sum(indexed(A, K) * indexed(B, K), K)` for equal-shaped operands;
for scalars it is their product.

The graph of `A + B` is the composition, with no node for the sum of tensors.
A language extension adds such a node if a pass must meet `A + B` as written
(see [Language extensions](extensions.md#language-extensions)).

`matrix_product(A, B)` takes operands `A` and `B` each of rank 1 or 2. With
the derived indices $i$, $k$ and $j$ =
`Index(("uflx.matrix_product", A, B, m))` for $m = 0, 1, 2$, the free
function reads the ranks and builds:

| Value shapes of `A` and `B` | `matrix_product(A, B)` |
| --------------------------- | ---------------------- |
| `(n,)` and `(n,)` | `index_sum(A[k] * B[k], (k,))` |
| `(m, n)` and `(n,)` | `component_tensor(index_sum(A[i, k] * B[k], (k,)), (i,))` |
| `(n,)` and `(n, p)` | `component_tensor(index_sum(A[k] * B[k, j], (k,)), (j,))` |
| `(m, n)` and `(n, p)` | `component_tensor(index_sum(A[i, k] * B[k, j], (k,)), (i, j))` |

`matrix_product` indexes an operand of rank 2 with two indices and every other
operand with one index, so `Indexed` raises `IndexError` for an operand of
rank 0 or of rank 3 or more. Unequal inner extents raise `ValueError` from
`Product`, because $k$ has two extents. Neither operand is conjugated, so
`x @ y` of two complex vectors is not the inner product.

For ranks 1 and 2, `matrix_product` agrees with `@` of NumPy. NumPy treats an
operand of rank 3 or more as a stack of matrices, so the core rejects such an
operand and does not give `@` a second meaning. A contraction of higher rank
is written in index notation, or with a free function of a library, such as
`dot` (see [Tensor algebra](../examples/tensors.md#tensor-algebra)).

Examples, with a scalar `a`, vectors `x` and `y` of value shape `(3,)`,
a matrix `A` of value shape `(2, 3)` and `i = index("i")`:

- `a * x` scales the vector and has value shape `(3,)`.
- `x[i] * y` also has value shape `(3,)`, but retains the free index `i`.
- `A @ x` contracts the last axis of `A`, giving value shape `(2,)`.

### `Inverse`

`Inverse` is the inverse of a square matrix, or the pseudo-inverse of a matrix
with more rows than columns.

| `Inverse(A)`  |                                                                                         |
| ------------- | --------------------------------------------------------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `inverse(A)`                                                                            |
| Takes         | `A: AbstractExpression` |
| Properties    | `operand` $A$ |
| Requires      | $A$ with value shape `(n, k)`; either `n == k`, or both extents are integers with `k < n` and $A$ is real |
| Static type   | Value shape `(k, n)`                                                                    |

Equal extents establish the square case even when they are symbolic.
Unequal extents must both be integers; symbolic extents cannot be ordered.

`Inverse(A)` is the pseudo-inverse $A^+ = (A^T A)^{-1} A^T$ of a matrix $A$
of matrix rank $k$, the left inverse of $A$ with $A^+ A$ the identity. For
$k = n$, $A^+$ is the inverse matrix $A^{-1}$.

A pass chooses how to evaluate `Inverse`, in closed form for a small matrix or
by a factorisation. The definition fixes the value and not the algorithm: the
normal equations square the condition number of $A$, and a factorisation such
as QR does not.

Examples:

```python
one, two, three = Integer(1), Integer(2), Integer(3)
A = tensor((tensor((two, one)), tensor((one, three))))
J = tensor((tensor((one,)), tensor((one,))))
```

- `inverse(A)`
  has value shape `(2, 2)` and value
  $\frac{1}{5} \begin{pmatrix} 3 & -1 \\ -1 & 2 \end{pmatrix}$.
- `J` has value shape `(2, 1)`. `inverse(J)` has value shape `(1, 2)` and value
  $\begin{pmatrix} 1/2 & 1/2 \end{pmatrix}$. This is a left inverse.
- `inverse(transpose(J))` raises `ValueError`: the operand has more columns
  than rows.

### `Determinant`

`Determinant(A)` is the determinant of a square matrix $A$.

| `Determinant(A)` | |
| ---------------- | - |
| Extends | `AbstractExpression` |
| Free function | `determinant(A)` |
| Takes | `A: AbstractExpression` |
| Properties | `operand` $A$ |
| Requires | Value shape `(n, n)` and argument spaces `()` |
| Static type | Value shape `()` |
| Raises | `TypeError` for an $A$ that is not an expression, and `ValueError` if a requirement fails |

The extent $n$ may be an `int` or an `AbstractExtent`; the two axes must
have equal extents. The empty determinant, of shape `(0, 0)`, has value 1.
Domain, free indices and scalar type follow the operand.

The node remains in the core normal form. The compiler chooses its evaluation,
such as an explicit formula for a small matrix or a factorisation. The
[derivative rule](evaluation.md#derivatives-of-the-node-classes) uses `Inverse`.

Worked examples for this chapter: [transpose](../examples/tensors.md#transpose), [tensor algebra](../examples/tensors.md#tensor-algebra).

## Calculus on a domain

### `Coordinate`

`Coordinate` maps each point of a domain to its coordinates in the ambient
coordinate domain. It is a terminal.

| `Coordinate(P)` |                                                    |
| --------------- | -------------------------------------------------- |
| Extends         | `AbstractExpression` |
| Free function   | `coordinate(P)`                                    |
| Takes           | `P: AbstractDomain` |
| Properties      | `domain`, an `AbstractDomain` $P$                  |
| Static type     | Value shape `(n,)`, real, the domain $P$          |
| Raises          | `TypeError` if $P$ is not an `AbstractDomain`      |

$n$ is the dimension of the ambient coordinate domain of $P$. For a
coordinate domain $P$, the coordinate function is the identity,
$X : P \to R^n$, $X(p) = p$, and the gradient of $X$ is the identity matrix
(see [`Gradient`](#gradient)).

For a mesh $D$, `coordinate(D)` is the physical coordinate field on $D$.
Its pull-back by a cell map is that cell's parameterisation (see
[`Mesh`](geometry.md#mesh)). This differs from `D.coordinate_function`, which is
the stored coordinate expression on the topology.
For a mesh `image(x, Q)`, the pull-back by `cell_map(D, j, e)` is
`pull_back(x, cell_map(T, j, e))`.

For example, `coordinate(RD(3))` has value shape `(3,)` and domain `RD(3)`.

### Valid pair

An expression and a domain $(f, P)$ are a valid pair if $P$ is an
[`AbstractDomain`](geometry.md#abstractdomain) and $f$ has the domain $P$ or no domain.
Otherwise construction raises `TypeError` for $P$ or `ValueError` for $f$. An
$f$ with no domain is a constant on $P$.

### `Gradient`

`Gradient` is the gradient of an expression on a domain $P$. The gradient adds a
last axis in the coordinate frame of the ambient coordinate domain of $P$, of
dimension $n$.

| `Gradient(f, P)` |                                                                                                   |
| ---------------- | ------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `gradient(f, P)`                                                                                  |
| Takes            | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties       | `expression` $f$ and `domain` $P$, with $n$ the dimension of the ambient coordinate domain of $P$ |
| Requires         | A [valid pair](#valid-pair) $(f, P)$                                                              |
| Static type      | Value shape `f.value_shape + (n,)`                                                                |
| Raises           | `TypeError` for $P$ and `ValueError` for $f$ if a requirement fails                               |

The gradient acts on each component of an array-valued $f$. Entry `j` of the
last axis of `Gradient(f, P)` is the partial derivative $\partial_j f$ in the
coordinate $X_j$ of the frame.

On a mesh, `Gradient`, `Divergence` and `Curl` act cellwise.
They do not include distributional interface terms for discontinuous fields;
traces and interface integrals are written explicitly.

#### The frame

For a domain $P$ of dimension $k$, the frame is the coordinate frame of the
ambient coordinate domain of $P$, of dimension $n$. A coordinate domain has
$k = n$.

- For $k = n$, $\partial_j$ is the partial derivative in the coordinate
  $X_j$ of $P$.
- For $k < n$, such as a mesh that is a surface in $R^3$, the
  coordinate directions of the ambient coordinate domain are not tangent to
  $P$. `Gradient(f, P)` is then the tangential gradient, the gradient of an
  extension of $f$ off $P$, projected orthogonally onto the tangent space of
  $P$ with the Euclidean metric of the ambient coordinate domain. The
  tangential gradient does not depend on the extension.

The divergence and the curl are [`Divergence` and `Curl`](#divergence-and-curl),
defined by `Gradient`. Other operators of vector calculus are fixed
compositions of `Gradient`, which a library provides as free functions (see
[New operations](extensions.md#new-operations)).

#### Pull-backs of `Gradient`

The gradient depends on the coordinates of $P$, so a pull-back does not
commute with `Gradient`. `PullBack(Gradient(f, Q), D)` is the gradient of $f$
on $Q$, evaluated at $\psi(p)$.

#### Examples

Start with a scalar expression on `RD(3)`:

```python
P = RD(3)
X = coordinate(P)
x, y, z = X[0], X[1], X[2]
f = x * y * z
```

| Expression                              | Value shape | Meaning                                       |
| --------------------------------------- | ----------- | --------------------------------------------- |
| `gradient(f, P)`                        | `(3,)`      | The gradient $(yz, xz, xy)$, after evaluation |
| `gradient(gradient(f, P), P)`           | `(3, 3)`    | The Hessian of $f$                            |
| `gradient(X, P)`                        | `(3, 3)`    | The identity matrix, after evaluation         |

- Let `M` be a mesh of dimension 2 in `RD(3)`, and `w` a field
  on `M`. `gradient(w, M)` is the surface gradient of `w`, with value shape
  `(3,)`. `gradient(coordinate(M), M)` is the orthogonal projection onto the
  tangent plane of `M`, with value shape `(3, 3)`.
- `gradient(X, RD(2))` raises `ValueError`, because `X` has the domain
  `RD(3)`.

### `Divergence` and `Curl`

`Divergence` and `Curl` are the divergence and curl of an expression $u$ on a
domain $P$, with $n$ the dimension of its ambient coordinate domain.

| `Divergence(u, P)` |                                                                     |
| ------------------ | ------------------------------------------------------------------- |
| Extends            | `AbstractExpression` |
| Free function      | `divergence(u, P)`                                                  |
| Takes              | `u: AbstractExpression`, `P: AbstractDomain` |
| Properties         | `expression` $u$ and `domain` $P$                                   |
| Requires           | A [valid pair](#valid-pair) $(u, P)$, and $u$ of rank at least 1 with the last extent $n$ |
| Static type        | Value shape `u.value_shape[:-1]`                                    |
| Raises             | `TypeError` for $P$ and `ValueError` for $u$ if a requirement fails |

| `Curl(u, P)`  |                                                                     |
| ------------- | ------------------------------------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `curl(u, P)`                                                        |
| Takes         | `u: AbstractExpression`, `P: AbstractDomain` |
| Properties    | `expression` $u$ and `domain` $P$                                   |
| Requires      | A [valid pair](#valid-pair) $(u, P)$, $n$ 2 or 3, and $u$ of value shape `(n,)` |
| Static type   | Value shape `(3,)` for $n = 3$ and `()` for $n = 2$                 |
| Raises        | `TypeError` for $P$ and `ValueError` for $u$ if a requirement fails |

The value of `Divergence(u, P)` is the trace of `gradient(u, P)` over the last
axis of $u$ and the axis of the gradient. With `G = gradient(u, P)`, so
$G_{ij}=\partial_j u_i$, `Curl(u, P)` is the antisymmetric contraction:

| Ambient dimension | Value |
| ----------------- | ----- |
| 2 | `G[1, 0] - G[0, 1]` |
| 3 | `tensor((G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]))` |

`Gradient` stays the one derivative that defines the others. The two node
classes are in the core so that the identity of Piola
holds structurally (see [`PushForward`](geometry.md#pushforward)).

### `DirectionalDerivative`

`DirectionalDerivative` differentiates an expression $f$ with respect to a
variable $x$ at a point $a$ in a direction $h$,

$$
\mathrm{D}(f, x, a, h) = \frac{d}{d\varepsilon} f[x := a + \varepsilon h] \Big|_{\varepsilon = 0},
$$

with $\varepsilon$ real and $f[x := s]$ the substitution of $s$ for $x$ in $f$
that does not enter domains, domain maps or the spaces of functions. The
variable $x$ is a function or a `tuple` of functions. For a `tuple`, $a$ and $h$
are `tuple`s with one member for each member of $x$, and the substitution
replaces every member $x_k$ by $a_k + \varepsilon h_k$ at once. The functions of
$x$ are bound in $f$ (see [Substitution and
binding](evaluation.md#substitution-and-binding)).

| `DirectionalDerivative(f, x, a, h)` |                                                                                                       |
| -------------------------------- | ----------------------------------------------------------------------------------------------------- |
| Extends                          | `AbstractExpression` |
| Takes                            | `f: AbstractExpression`, `x: Function \| tuple[Function, ...]`, `a, h: AbstractExpression \| tuple[AbstractExpression, ...]` |
| Properties                       | `expression` $f$, `variable` $x$, `point` $a$ and `direction` $h$ |
| Requires                         | $x$ a [bindable function](forms.md#bindable-functions), or a `tuple` of at least two distinct bindable functions; member by member, $a$ and $h$ admissible values for the space of $x$ as in [`Action`](forms.md#action) |
| Static type                      | The value shape of $f$, the free indices of $f$ followed by those of the members of $a$ and $h$, each `Index` taken once, and the argument spaces and the scalar type of $f$ |
| Raises                           | `TypeError` for an $x$ that is not a `Function` or a `tuple` of functions, and `ValueError` if another requirement fails |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `directional_derivative(f, x, a, h)` | | `f: AbstractExpression`, `x: Function \| tuple`, `a, h: AbstractExpression \| tuple` | `DirectionalDerivative(f, x, a, h)` |

This is a partial derivative: substitution replaces only occurrences of $x$.
The point $a$ need not be a function. A `tuple` variable differentiates jointly
in its members, as for the unknowns of a mixed problem. Real points and
directions for a complex variable use the complex embedding of
[type compatibility](evaluation.md#type-compatibility-of-rewrites) before substitution.

A dual function may be a variable. Its point and direction may be symbolic
forms with the same functional signature; differentiation preserves the
argument spaces of the expression being differentiated.

The substitution does not enter domains, so the derivative holds the geometry
fixed. A shape derivative, with respect to the coordinates of a mesh, is a
language extension (see [Shape derivatives](extensions.md#shape-derivatives)).

A real-linear differential is assumed. The gradient of a direction with no
domain is zero, so the derivative with respect to a function `u` in a constant
direction holds `gradient(u, P)` fixed. A direction that is a function in the
space of `u` gives the Gateaux derivative (see [Gateaux
derivative](forms.md#gateaux-derivative)).

For `u`, `a` and `h` in a real scalar space of constants,
`directional_derivative(u * u, u, a, h)` has the value $2 a h$. The variable
is `u`, the evaluation point is `a`, and the direction is `h`.

An indexed direction gives a family of derivatives. For scalar constants
`h0` and `h1` in the same space:

```python
i = index("i")
directions = tensor((h0, h1))
d = directional_derivative(u * u, u, a, directions[i])
```

`d` has free index `i` and value `Integer(2) * a * directions[i]`. Tensor derivatives
use families of unit directions; the directions specify any symmetry. See
[Hyperelasticity and the stress](../examples/forms.md#hyperelasticity-and-the-stress).

### `Integral`

`Integral` integrates an expression $f$ over a domain $P$. A basic integral
needs no topological notion of the domain.

| `Integral(f, P)` |                                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `integral(f, P)`                                                                            |
| Takes            | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties       | `expression` $f$ and `integration_domain` $P$                                                           |
| Requires         | A [valid pair](#valid-pair) $(f, P)$, or $P$ a mesh and $f$ with the domain `P.ambient_domain` |
| Static type      | The value shape and free indices of $f$, with no domain |
| Raises           | `TypeError` for $P$ and `ValueError` for $f$ if a requirement fails                         |

For a coordinate domain $P$ of dimension $n$, `Integral(f, P)` is the
Lebesgue integral in the coordinates $X$ of $P$,

$$
\int_P f \,\mathrm{d}X.
$$

For $n = 0$, the integral is the value of $f$ at the single point of $R^0$.

On a mesh, the integral sums the integrals over its pieces with
their induced densities (see [the expansion](geometry.md#integral-over-a-mesh)).

An array-valued $f$ integrates component by component. A later pass might
approximate an `Integral` over a coordinate domain by
[quadrature](geometry.md#quadrature), or evaluate it exactly when the rule is exact.

#### Examples

With `T = UnitSimplex(2)` and `X = coordinate(T)`:

- `integral(X[0], T)` is $\int_T X_0 \,\mathrm{d}X$. A later pass evaluates
  the integral to $1/6$.
- `integral(Integer(1), T)` is the area $1/2$ of `T`.
- `integral(X[0], T) * integral(X[1], T)` is a product of two integrals, with
  value shape `()` and no domain.

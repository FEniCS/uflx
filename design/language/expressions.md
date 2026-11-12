# Expressions

API signatures use descriptive parameter names. A stored property has the same
name as its constructor parameter. API names, attribute access and literal
values use code typography; mathematical notation describes values and
relationships. Examples may use short local names such as `A` or `i`.

## Expressions and static type

### `AbstractExpression`

`AbstractExpression` is the base class for symbolic values. Its static type
describes the result of an expression without evaluating it.

| `AbstractExpression` |                                                                                      |
| -------------------- | ------------------------------------------------------------------------------------ |
| Extends              | `Base`                                                                               |
| Provides             | `value_shape`, `free_indices`, `scalar_type`, `domain`, `argument_spaces` |

| Property | Meaning |
| -------- | ------- |
| `value_shape` | Length of each array axis, in order; `()` for a scalar. |
| `free_indices` | Unbound indices on which the value depends, paired with their extents in first-occurrence order. |
| `scalar_type` | Real or complex arithmetic used to evaluate the expression. |
| `domain` | Where the expression is defined; `None` when it has no spatial argument. |
| `argument_spaces` | Spaces of the multilinear arguments on which the expression depends; `()` when there are none. |

These five derived properties form the static type (see
[Argument spaces](forms.md#argument-spaces)). The value shape is a `tuple` of
[extents](#extents), one per axis; `()` denotes a scalar. The array rank,
written simply as rank, is the number of axes.

### Extents

An extent is the length of an axis of a value shape, or the extent of a free
index. An extent is a nonnegative `int` or an `AbstractExtent`.

`AbstractExtent` represents an extent whose value is not known during
construction. Its supplied value must be a nonnegative integer.
An `AbstractExtent` is not an expression. A concrete
`AbstractExtent` is a subclass outside the core, such as `RuntimeExtent` in the
[cut cell example](../examples/integration.md#integrating-over-a-cut-cell).

| `AbstractExtent` |                                      |
| ---------------- | ------------------------------------ |
| Extends          | `Base`                               |
| Defined by       | A subclass outside the core          |

An `AbstractExtent` has no arithmetic and no ordering. Rules that compare
extents by equality accept an `AbstractExtent` unchanged.
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

The Python representation of the scalar types is implementation-defined. When
stored as a construction argument, it must satisfy rule 1.

The scalar type is the arithmetic in which the expression is evaluated, not a
statement about the value. A real `SquareRoot`, `Logarithm` or `Power` whose
value is not a real number, such as `square_root(Integer(-1))`, is undefined.

### Domain

An expression's `domain` identifies where it is defined:

```python
domain: AbstractTopology | AbstractDomain | None
```

An [`AbstractTopology`](topology.md#abstracttopology) supplies labelled reference-cell
slots; its cell types must be reference cells when used as a field domain.
An [`AbstractDomain`](topology.md#abstractdomain) supplies the ambient-coordinate
contract. `None` means the result has no spatial argument, as for a literal or
an integral.

Unless the [static type table](reference.md#appendix-static-type-table) states otherwise,
a node inherits the domain of its expression operands with a domain. All such
operands must have the same domain, or construction raises `ValueError`. A node
with no such operand has domain `None`.

The property always describes the resulting expression. An `Integral` stores
its region as `integration_domain` but has domain `None`; a `PullBack` stores
its `map` and has the domain from which that map starts: `parameter_domain`
for a `DomainMap`, `reference_cell` for a `CellMap`, or `source`
for a `SideMap`. A `PushForward` stores its
`mesh`, which is also its domain. For `Coordinate`, `Gradient`, `Divergence`
and `Curl`, the stored domain is also the result's domain.

### Domain-free expressions

An expression with `domain is None` has no spatial argument. For fixed values
of its free indices and form arguments, it denotes one array. Where an operation
accepts it on a domain $P$, use that same array at every point of $P$.

| Property or operation | Meaning |
| --------------------- | ------- |
| Domain-free value on $P$ | A spatially constant field; the expression itself keeps domain `None` |
| Free indices and form arguments | May vary the value; the constant-field interpretation applies separately to each assignment |
| Functions and stored geometry | May affect the value. Having no spatial domain does not imply independence from these inputs. |
| Spatial derivative | Zero, with the required type and derivative index |
| Pull-back | The same value at every parameter point, with type completion for any added free indices |
| Integration | Integrate that constant field with the domain's measure |

For example, `integral(u, P)` has no spatial argument but may depend on `u`
and the geometry of $P$. A function in a constant space is spatially constant
but still varies when differentiated with respect to that function.

### Expressions as functions

An expression denotes an array-valued function on the domain of the
expression. An expression on a topology denotes a function on the reference
cell of each cell of the topology (see [`CellMap`](topology.md#cellmap)). An expression
with no domain has the [constant-field interpretation](#domain-free-expressions).
An expression with argument spaces other than `()` denotes a multilinear map from
the argument spaces to such functions (see [Variational
forms](forms.md#variational-forms)).

### Node class or free function

Use a construction free function for a fixed composition permitted by the
[construction policy](rules.md#construction-policy). Use a node class when
choosing the result requires inspecting expression structure, the compiler
must choose an algorithm, or a pass must retain the operation as written.
In the last case, the node's section states its composition.
For example, transpose is the fixed composition
`component_tensor(A[k, l], (l, k))`; `Inverse` retains the algorithm choice.

A construction free function belongs to the core only if a core definition,
rule or Python overload uses it. A node class belongs to the core only if
another core definition or rule needs it, or it is an operation of the scalar
field. Its own rules do not establish that need.

Other compositions belong to libraries; other node classes belong to
[language extensions](extensions.md#language-extensions). The
[node-class rationale](reference.md#appendix-node-class-rationale) explains
individual choices.

## Numbers and scalar operations

### Literals

A literal is a number written in an expression, such as `Integer(2) * u`.
`Integer` and `Float` are real terminals with value shape `()`, no free
indices, no domain and argument spaces `()`.

| `Integer(value)` |                       |
| ------------ | --------------------- |
| Extends      | `AbstractExpression` |
| Properties   | `value`, the stored integer, of exact Python type `int` |

| `Float(value)` |                               |
| ---------- | ----------------------------- |
| Extends    | `AbstractExpression`         |
| Properties | `value`, the stored finite number, of exact Python type `float` |

`Integer` and `Float` accept only exact Python types, as in rule 1 of
[the base object](rules.md#the-base-object-base). So `Integer(True)` raises
`TypeError`, and `Float(-0.0) != Float(0.0)`.

A rational number is the `Division` of two `Integer`s, such as
`Integer(1) / Integer(3)`. The compiler chooses its numerical evaluation;
the core has no rational literal.

### `Complex`

`Complex` is a scalar complex constant built from two real literals. It is a
terminal. The `real` property supplies its real part and `imag` supplies its
imaginary part.

| `Complex(real, imag)` | |
| --------------------- | - |
| Extends               | `AbstractExpression` |
| Takes                 | `real, imag: Integer \| Float` |
| Properties            | `real` and `imag`, the real and imaginary parts |
| Static type           | Value shape `()`, complex, no free indices, no domain and argument spaces `()` |

```python
Complex(Integer(0), Integer(1))  # The imaginary unit.
Complex(Integer(-1), Integer(0))  # Minus one with complex scalar type.
```

`Complex` retains its parts and stays complex even when the imaginary part
is zero. As a terminal, it has no evaluation operands (see
[Evaluation traversal](rewriting.md#evaluation-traversal)).

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
| Raises                                             | `ValueError` for a negative extent or an index that occurs twice |

For example, `Zero((3,), (), (), real)` is a real zero vector and
`Zero((), (), (V,), complex)` is a complex zero form on `V`. Here `real` and
`complex` denote the scalar types, not Python numeric constructors.

### Arithmetic

Arithmetic nodes retain scalar and tensor operations as written. Their free
functions construct the corresponding nodes directly. Invalid operand shapes
raise `ValueError`. Free indices follow the [default rule](#free-indices).

| Node | Operand shapes | Value shape | Meaning |
| ---- | -------------- | ----------- | ------- |
| `Sum(first, second)` | Equal | Shape of `first` | Entrywise sum |
| `Product(first, second)` | At least one scalar | Shape of the other operand, or `()` if both scalar | Scalar multiplication, preserving operand order |
| `Division(numerator, denominator)` | Scalar `denominator` | Shape of `numerator` | Entrywise division by `denominator` |
| `Power(base, exponent)` | Both scalar | `()` | Scalar power |

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
| Properties | `first` and `second`, the expressions to add |
| Requires | `first` and `second` have equal value shapes |

### `Product`

| `Product(first, second)` |                                                  |
| ------------------------ | ------------------------------------------------ |
| Extends                  | `AbstractExpression` |
| Free function            | [`multiply`](#tensor-arithmetic), with the overload `*` |
| Takes                    | `first, second: AbstractExpression` |
| Properties | `first` and `second`, the factors in their stored order |
| Requires | At least one of `first` and `second` is scalar |

### `Power`

| `Power(base, exponent)` |                                                                             |
| ---------------------- | --------------------------------------------------------------------------- |
| Extends                | `AbstractExpression` |
| Free function          | `power(base, exponent)`, with the overload `**`                              |
| Takes                  | `base, exponent: AbstractExpression` |
| Properties | `base`, the value to raise, and `exponent`, the power to raise it to |
| Requires | `base` and `exponent` are scalar |

| Operands | Value |
| -------- | ----- |
| Defined base, exponent zero | One, including `0**0`, with the scalar type of the power |
| Defined base, exponent one | The base value, with the scalar type of the power |
| Other defined powers | The real power for real operands, or the principal complex value consistent with `Logarithm` |
| Either operand undefined | Undefined, including an undefined base raised to zero |

For example, constructor keywords and property access use the same names:

```python
squared = Power(base=expression, exponent=Integer(2))
squared.base == expression
squared.exponent == Integer(2)
```

These are value rules; construction retains `Power` and its static type.
The exponent may have domain. The [derivative rules](differentiation.md#derivatives-of-the-node-classes)
distinguish a fixed zero or one exponent from an exponent that only takes that
value at the differentiation point.

### Negation and subtraction

Negation and subtraction are free functions, not node classes:

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `negate(operand)` | `-operand` | `operand: AbstractExpression` | `multiply(Integer(-1), operand)` |
| `subtract(first, second)` | `first - second` | `first, second: AbstractExpression` | `add(first, negate(second))` |

These compositions inherit the tensor
behaviour of `multiply` and `add` (see [Tensor arithmetic](#tensor-arithmetic)).

### `Division`

| `Division(numerator, denominator)` |                                                                                        |
| ------------------------- | -------------------------------------------------------------------------------------- |
| Extends                   | `AbstractExpression` |
| Free function             | [`divide`](#tensor-arithmetic), with the overload `/` |
| Takes                     | `numerator, denominator: AbstractExpression` |
| Properties | `numerator`, the expression to divide, and `denominator`, the scalar divisor |
| Requires | `denominator` is scalar; `numerator` may have any value shape |

### `Absolute`

`Absolute` is the absolute value, or the modulus of a complex scalar.

| `Absolute(operand)` |                                           |
| ------------- | ----------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `absolute(operand)`, with the overload `abs(operand)` |
| Takes         | `operand: AbstractExpression` |
| Properties    | `operand`, the scalar whose absolute value or complex modulus is taken |
| Requires      | A scalar `operand`                              |
| Static type   | Value shape `()`, real                    |

### `Conjugate`

`Conjugate` is the complex conjugate of an expression of any value shape, entry
by entry.

| `Conjugate(operand)` |                                                                                             |
| -------------- | ------------------------------------------------------------------------------------------- |
| Extends        | `AbstractExpression` |
| Free function  | `conjugate(operand)`                                                                              |
| Takes          | `operand: AbstractExpression` |
| Properties     | `operand`, the expression to conjugate entry by entry |
| Static type    | The value shape and the free indices of `operand`, and the argument spaces of `operand` with each space conjugated (see [`ConjugateSpace`](forms.md#conjugatespace)) |

`Conjugate` accepts any value shape because
[`Form`](forms.md#form) binds `conjugate(v)` even for a vector argument. Its entrywise
meaning is `conjugate(f)[e]` = `conjugate(f[e])` in value.

### `RealPart` and `ImaginaryPart`

`RealPart` and `ImaginaryPart` give the real and imaginary parts of a scalar.

| `RealPart(operand)`, `ImaginaryPart(operand)` |                                       |
| --------------------------------- | ------------------------------------- |
| Extends                           | `AbstractExpression` |
| Free function                     | [`real_part`](#tensor-arithmetic), [`imaginary_part`](#tensor-arithmetic) |
| Takes                             | `operand: AbstractExpression` |
| Properties                        | `operand`, the scalar whose real or imaginary part is taken |
| Requires                          | A scalar `operand`                          |
| Static type                       | Value shape `()`, real                |

These nodes and `Absolute` are the only nodes that take a complex operand to
a real expression. The composition `(f + conjugate(f)) / Integer(2)` has the
value of the real part but remains complex by the [scalar-type rule](#scalar-type).

Their free functions apply entrywise to tensors (see
[Tensor arithmetic](#tensor-arithmetic)).

The operand may be real. Construction does not simplify, so `real_part(f)`
of a real `f` is not `f`, and `imaginary_part(f)` of a real `f` is not
`Zero`.

For `z = Complex(Integer(1), Integer(2))`, `real_part(z)` and
`imaginary_part(z)` are real expressions with the values 1 and 2, respectively.

### Elementary functions

`SquareRoot` and `Logarithm` give the square root and natural logarithm of a
scalar. The core provides these two; other layers are expected to provide
a wider range of elementary functions.

| `SquareRoot(operand)`, `Logarithm(operand)` |                                     |
| ------------------------------- | ----------------------------------- |
| Extends         | `AbstractExpression` |
| Free function   | `square_root(operand)` and `logarithm(operand)` |
| Takes           | `operand: AbstractExpression` |
| Properties      | `operand`, the scalar to which the square root or logarithm is applied |
| Requires        | A scalar `operand`                        |
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
| Static type       | Value shape `(n,) + s`, with `n` the number of entries                                                    |

A matrix is a `Tensor` of `Tensor`s.

Value shapes may contain extent 0, as in `coordinate(EuclideanSpace(0))`. An empty
`Tensor` cannot provide the entry shape `s`. Use `Zero` for a typed zero or
[`broadcast_to`](#broadcastto) with an operand supplying the member type.

Examples, with `one, two = Integer(1), Integer(2)`:

- `tensor((one, two))` is a vector of value shape `(2,)`.
- `tensor((one, tensor((one, two))))` raises `ValueError`, because the entries
  have different value shapes.

### `BroadcastTo`

`BroadcastTo` repeats an expression over a complete target shape. Source axes
align with the trailing target axes. It binds no indices.

| `BroadcastTo(value, value_shape)` | |
| -------------------------------- | - |
| Extends | `AbstractExpression` |
| Free function | `broadcast_to(value, value_shape)` |
| Takes | `value: AbstractExpression`, `value_shape: tuple[int \| AbstractExtent, ...]` |
| Properties | `value` and the complete target `value_shape` |
| Requires | Target entries are [extents](#extents); target rank is at least the source rank; every aligned source extent equals the target extent or is the integer 1 |
| Static type | Target value shape; free indices, domain, scalar type and argument spaces of `value` |
| Raises | `ValueError` for an invalid extent, lower target rank or incompatible aligned extents |

| Case | Meaning |
| ---- | ------- |
| Leading target axis | Repeat the source along that axis |
| Aligned source extent 1 | Use source entry 0 for every target entry |
| Equal aligned extents | Use the same entry on both axes; symbolic extents must be structurally equal |
| Zero target extent | No member value is evaluated. The operand still supplies the static type and remains a graph dependency; construction validates both operand and target shape. |
| Target `()` | Valid only for a scalar source |
| Target extent `-1` | Invalid; extents are never inferred |

```python
broadcast_to(Integer(1), (n,))[i]  # One for each i; i remains free.
broadcast_to(v, (n,) + v.value_shape)  # n copies of v.
broadcast_to(v, (0,) + v.value_shape)  # An empty array with v's member type.
```

Selecting a broadcast uses the ordinary indexing rules, including free-index
extent checks. A free selector remains free even when every selected value is
the same. Arithmetic does not broadcast implicitly.

### Index notation

Index notation uses `Tensor` to build a tensor from scalars, `BroadcastTo` to
repeat values, and `Indexed` to select one component. `ComponentTensor` and
`IndexSum` bind the free indices of a scalar, as axes or as a sum. Tensor
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
| Properties                     | `expression`, the tensor to read, and `indices`, a non-empty `tuple` selecting one entry per axis as described below    |
| Static type                    | Value shape `()`, and the free indices below                                      |
| Raises                         | `ValueError` and `IndexError`, below                                 |

The factory `indexed(f, ())` returns `f` if `f` is scalar, and raises
`IndexError` otherwise: full indexing must cover every axis. `Indexed(f, ())`
still raises `ValueError`; the factory creates no node for the scalar case.

Indexing is full indexing. Each entry covers one axis of the
expression, and the entries together cover every axis exactly once. With
`(n_0, ..., n_{r-1})` the value shape of the expression, an entry is one of
the following kinds:

| Entry                                | Covers                                | Rule                                                                                               | Raises       |
| ------------------------------------ | ------------------------------------- | -------------------------------------------------------------------------------------------------- | ------------ |
| An `int` $k$                         | One axis, of extent $n_j$             | $k \ge 0$; if $n_j$ is an `int`, $k < n_j$; otherwise the caller ensures the upper bound, as below | `IndexError` |
| An [`Index`](#index)                 | One axis, of extent $n_j$             | The `Index` takes the extent $n_j$                                                                 |              |
| A [`GatheredIndex`](#gatheredindex)  | One axis, of extent $n_j$             | The extent of the `GatheredIndex` equals $n_j$                                                     | `IndexError` |

For a symbolic extent, construction cannot check the upper bound. The caller
must ensure that an integer entry is smaller than the supplied extent.
Construction checks nonnegativity and any known upper bound.

The free indices of the `Indexed` are the free indices of the expression,
followed by each `Index` entry with the extent that the
entry takes, followed by the free indices of each `GatheredIndex` entry, with
each index taken once.

An `Index` may occur in two or more entries, or may already be free in the
expression. All the occurrences are the same index and must have the same
extent, so the indexed axes take the same value. `A[i, i]` is the diagonal
entry of a matrix `A`.

Static checking rejects entries of another kind. Construction raises
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
| `component(expression, entries)` | `[]` | `expression: AbstractExpression`, `entries: tuple` | `expression` if `entries == ()`; otherwise `component_tensor(indexed(expression, entries + K), K)` if the entries are fewer than the rank of `expression`, otherwise `indexed(expression, entries)` |

`K` is the `tuple` of the derived indices
`Index(("uflx.component", expression, entries, m))`, one for each axis that
`entries` does not cover (see [Derived indices](#derived-indices)). `component`
reads the rank of `expression` and the number of entries.
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
| Properties                                                              | `expression`, the scalar body, and `indices`, a non-empty `tuple` of distinct free `Index` entries to bind as axes or sum over, in order                                                                                 |
| Static type                                                             | The free indices of `expression` with the bound `indices` removed. `ComponentTensor` has the value shape of the extents of the indices in order, and `IndexSum` the value shape `()` |

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
[Substitution and binding](substitution.md#substitution-and-binding)).

For a square matrix `A`, `index_sum(A[i, i], (i,))` is its trace, and
`component_tensor(A[i, i], (i,))` is its diagonal. Both work when the extent is
an `AbstractExtent`, without knowing its value during construction.

### `AbstractIndexTable`

`AbstractIndexTable` represents a table of index values supplied by the solver,
such as a local basis permutation, a local-to-patch degree-of-freedom map,
or the local facet number selecting a pretabulated reference map. The values are fixed but not known during construction. An
`AbstractIndexTable` is not an expression and not a free index.

| `AbstractIndexTable` |                                                        |
| -------------------- | ------------------------------------------------------ |
| Extends              | `Base` |
| Provides             | `extent`, an extent giving the exclusive upper bound on each stored index value, and `shape`, a `tuple` of extents giving the table dimensions |
| Defined by           | A subclass outside the core                            |

The table has the shape `shape`. Each stored value is a nonnegative `int`
smaller than `extent`. An `AbstractIndexTable` with the shape `()` is one
value. An expression reads the table only through a
[`GatheredIndex`](#gatheredindex), so one value is read as
`gathered_index(table, ())`.

### `GatheredIndex`

`GatheredIndex` reads an index value from a table. It is neither an expression
nor a free index; its table indices determine its free indices.

| `GatheredIndex(table, indices)` |                                                                                                                                    |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| Extends                         | `Base` |
| Free function                   | `gathered_index(table, indices)`                                                                                                   |
| Takes                           | `table: AbstractIndexTable`, `indices: tuple` |
| Properties                      | `table`, the index table to read, and `indices`, a `tuple` selecting one entry per table axis using an `int`, `Index` or `GatheredIndex` |
| Derived                         | `extent` and `free_indices`, below |
| Raises                          | `IndexError` for an incorrect number of entries or invalid bounds or gather extent, and `ValueError` for inconsistent extents of a repeated index |

`indices` contains one `int`, `Index` or `GatheredIndex` per table axis, with
the same bounds and extent checks as [`Indexed`](#indexed). A scalar table
therefore takes `()`. Both constructor and factory follow this rule.
Gathers may nest. The result has the table's
extent, and its free indices are the `Index` entries with their axis extents,
followed by the free indices of the `GatheredIndex` entries, taking each index
once.

#### Examples

A local permutation table `permutation` has shape `(3,)`, extent 3 and
values `(2, 0, 1)`. The package supplies it to account for a local orientation.
For a local coefficient array `a` of shape `(3,)`:

```python
i = index("i")
g = gathered_index(permutation, (i,))
b = component_tensor(a[g], (i,))
```

The result is `(a[2], a[0], a[1])`. For an array `c` of shape `(4,)`, `c[g]`
raises `IndexError`, because the gather extent 3 differs from the axis extent 4.
A second local permutation `p` can be composed as
`gathered_index(p, (g,))`. These operations index local arrays; they do not
select cells from a mesh.

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
renamed index of [Substitution and binding](substitution.md#substitution-and-binding) has the
label `("uflx.rename", i, body, (s_1, ..., s_n))`.

A label refers only to objects that the index is fresh for, and never to an
object that contains the index. A replacement rebuilds a derived index whose
label contains a replaced object, so that the replacement commutes with the
free functions. A renamed name keeps its label.

```{dropdown} Why renamed labels are fixed
Renaming continues into operands containing the new name. Rebuilding its
label there would substitute into $s_1, \dots, s_n$ again; these may contain
the binding node, preventing termination.
```

The edges of an index to the objects in the label of the index identify the
index and are not operands (see [Evaluation traversal](rewriting.md#evaluation-traversal)).

Label contents add no [free occurrences](substitution.md#free-occurrences),
including in function labels. They still determine identity and are rebuilt
according to the substitution rules above.

A free function or a pass derives one `Index` per axis.
`component_tensor(e, K)[J]` is $e$ with each index of `K` replaced by the
corresponding entry of `J`.

## Tensor operations

### Tensor arithmetic

Addition, scalar multiplication and division retain their arithmetic nodes.
[Component expansion](lowering.md#expanding-tensor-arithmetic) introduces index
notation during lowering. The other operations below use fixed compositions.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `add(first, second)` | `+` | `first, second: AbstractExpression` | `Sum(first, second)` |
| `multiply(first, second)` | `*` | `first, second: AbstractExpression` | `Product(first, second)` |
| `divide(numerator, denominator)` | `/` | `numerator, denominator: AbstractExpression` | `Division(numerator, denominator)` |
| `real_part(operand)` | | `operand: AbstractExpression` | `component_tensor(RealPart(operand[K]), K)` |
| `imaginary_part(operand)` | | `operand: AbstractExpression` | `component_tensor(ImaginaryPart(operand[K]), K)` |
| `matrix_product(first, second)` | `@` | `first, second: AbstractExpression` | The contraction of the last axis of `first` with the first axis of `second`, below |
| `transpose(matrix)` | | `matrix: AbstractExpression` | `component_tensor(matrix[k, l], (l, k))`, with `(k, l)` = `K`, for `matrix` of rank 2 |

In each composition, `K` contains one [derived index](#derived-indices)
per axis of the indexed operand. Its label contains the qualified function
identifier, the function arguments in signature order, and the axis number.
For example, `real_part(operand)` uses `Index(("uflx.real_part", operand, m))`,
where `m` runs from zero to rank minus one.
The empty-index factory identities give scalar operations when `K == ()`.
Scalar constructors inside the compositions avoid recursive factory calls.
Operand free indices remain free.

There is no implicit broadcasting. Use [`broadcast_to`](#broadcastto) explicitly.

| Invalid operands | Exception |
| ---------------- | --------- |
| Addition: unequal shapes | `ValueError` from `Sum` |
| Multiplication: two tensors | `ValueError` from `Product` |
| Division: tensor denominator | `ValueError` from `Division` |

A scalar operation's free function accepts tensors only when the operation has
no competing tensor meaning. Thus `conjugate`, `real_part` and `imaginary_part`
apply entrywise. The following remain scalar and raise `ValueError` for tensors:

| Free function | Second meaning on a tensor |
| ------------- | -------------------------- |
| `power`, `**` | The matrix power |
| `absolute`, `abs` | The norm of a vector |
| `square_root` | The matrix square root |
| `logarithm` | The matrix logarithm |

The same policy applies to future elementary functions, such as the
exponential. Entrywise application is
`component_tensor(op(indexed(T, K)), K)`, with one index per axis of `T`.
A full contraction is `index_sum(indexed(A, K) * indexed(B, K), K)` for
equal-shaped operands. Both also work for scalars, with `K == ()`.

These operations construct their compositions. Retaining a tensor operation
as a node requires a [language extension](extensions.md#language-extensions).

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

For ranks 1 and 2, `@` agrees with NumPy. Higher-rank contractions use index
notation or a library function such as [dot](../examples/tensors.md#tensor-algebra);
NumPy's stacked-matrix convention is not supported.

Examples, with a scalar `a`, vectors `x` and `y` of value shape `(3,)`,
a matrix `A` of value shape `(2, 3)` and `i = index("i")`:

- `a * x` scales the vector and has value shape `(3,)`.
- `x[i] * y` also has value shape `(3,)`, but retains the free index `i`.
- `A @ x` contracts the last axis of `A`, giving value shape `(2,)`.

### `Inverse`

`Inverse` is the inverse of a square matrix, or the pseudo-inverse of a matrix
with more rows than columns.

| `Inverse(matrix)`  |                                                                                         |
| ------------- | --------------------------------------------------------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `inverse(matrix)`                                                                            |
| Takes         | `matrix: AbstractExpression` |
| Properties    | `matrix`, the matrix to invert |
| Requires      | `matrix` with value shape `(n, k)` for fixed nonnegative integers `n` and `k`, and argument spaces `()`; either `n == k`, or `k < n` and `matrix` is real |
| Static type   | Value shape `(k, n)`                                                                    |

Let $A$ denote the value of `matrix`, with matrix rank $k$. The value of
`Inverse(matrix)` is the pseudo-inverse $A^+ = (A^T A)^{-1} A^T$,
the left inverse of $A$ with $A^+ A$ the identity. For
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

`Determinant(matrix)` is the determinant of the square matrix `matrix`.

| `Determinant(matrix)` | |
| ---------------- | - |
| Extends | `AbstractExpression` |
| Free function | `determinant(matrix)` |
| Takes | `matrix: AbstractExpression` |
| Properties | `matrix`, the matrix whose determinant is taken |
| Requires | Value shape `(n, n)` for a fixed nonnegative integer `n`, and argument spaces `()` |
| Static type | Value shape `()` |
| Raises | `ValueError` if a requirement fails |

The empty determinant, of shape `(0, 0)`, has value 1.
Domain, free indices and scalar type follow the operand.

The node remains in the core normal form. The compiler chooses its evaluation,
such as an explicit formula for a small matrix or a factorisation. The
[derivative rule](differentiation.md#derivatives-of-the-node-classes) uses cofactors.

Worked examples for this chapter: [transpose](../examples/tensors.md#transpose), [tensor algebra](../examples/tensors.md#tensor-algebra).

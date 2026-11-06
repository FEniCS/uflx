# Reference

## Appendix: Node-class rationale

This table explains why each concrete expression node is needed. The rules
are in the class sections.

| Node class                                       | Why a node class                                                                                              |
| ------------------------------------------------ | ------------------------------------------------------------------------------------------------------------- |
| `Integer`, `Float`                                | A number written down is a terminal |
| `Complex`                                        | A scalar complex constant retains its two real literals without arithmetic construction |
| `Zero`                                           | A zero carries its shape, free indices, argument spaces and scalar type without edges to other expressions                              |
| `Sum`, `Product` | Retain addition and scalar multiplication at their original value shape |
| `Power`, `Absolute` | Primitive scalar operations |
| `Conjugate`                                      | An operation of the scalar field, and `Form` meets `conjugate(v)` of a function of any value shape as written |
| `RealPart`                                       | `(f + conjugate(f)) / Integer(2)` has the value of the real part but is complex, and the real part is real    |
| `ImaginaryPart`                                  | `real_part(Complex(Integer(0), Integer(1)) * conjugate(f))` multiplies the real part by zero, which is NaN for an infinite real part |
| `Division`                                       | Preserves the written quotient rather than replacing it with multiplication by a reciprocal |
| `SquareRoot`                                     | IEEE 754 rounds the square root correctly and the power not, so `power(f, Integer(1) / Integer(2))` differs   |
| `Logarithm`                                      | An elementary function, which nothing else composes                                                           |
| `Inverse`, `Determinant`                         | The compiler chooses the numerical algorithm; construction retains the operation |
| `BroadcastTo` | Repeats values over explicit axes, including symbolic and empty shapes, without artificial arithmetic |
| `Tensor`, `Indexed`, `ComponentTensor`, `IndexSum` | The primitives of index notation                                                                            |
| `Coordinate`, `Function`, `Normal`               | Terminals: the coordinate function of a domain, a member of a space, the unit side conormal |
| `Gradient`                                       | Nothing composes to differentiate                                                                             |
| `Divergence`, `Curl` | Retain the operators so lowering can apply Piola identities before expanding them into gradients |
| `DirectionalDerivative`                          | Nothing composes to differentiate with respect to an expression                                               |
| `Integral`                                       | Nothing composes to integrate over a coordinate domain, and over a mesh a pass must meet the integral as written |
| `PullBack`                                       | Nothing composes to move an expression between domains                                                        |
| `PushForward` | Represents a topology field on a mesh and retains its transformation law for cellwise lowering and Piola identities |
| `Form`, `Apply`                                 | Binding arguments and fully applying a map are not compositions                                |

## Appendix: Class hierarchy

Every object extends `Base`. The diagram shows the abstract bases; the lists
below give concrete classes and their parents. Arrows point to parents.

<!-- generated: class-hierarchy, run design/appendices.py to update -->
```mermaid
---
config:
  theme: neutral
  look: classic
  class:
    hideEmptyMembersBox: true
---
classDiagram
    direction LR
    accTitle: Abstract class hierarchy
    Base <|-- AbstractExpression
    Base <|-- AbstractExtent
    Base <|-- AbstractIndexTable
    Base <|-- AbstractCellType
    Base <|-- AbstractTopology
    Base <|-- AbstractDomain
    AbstractDomain <|-- AbstractCoordinateDomain
    AbstractCoordinateDomain <|-- AbstractReferenceCell
    AbstractCellType <|-- AbstractReferenceCell
    Base <|-- AbstractFunctionSpace
```

### Expressions and static type

Extend `Base`:

[`Index`](expressions.md#index).

### Numbers and scalar operations

Extend `AbstractExpression`:

[`Integer`](expressions.md#literals), [`Float`](expressions.md#literals), [`Complex`](expressions.md#complex), [`Zero`](expressions.md#zero), [`Sum`](expressions.md#sum), [`Product`](expressions.md#product), [`Power`](expressions.md#power), [`Division`](expressions.md#division), [`Absolute`](expressions.md#absolute), [`Conjugate`](expressions.md#conjugate), [`RealPart`](expressions.md#realpart-and-imaginarypart), [`ImaginaryPart`](expressions.md#realpart-and-imaginarypart), [`SquareRoot`](expressions.md#elementary-functions), [`Logarithm`](expressions.md#elementary-functions).

### Tensors and indexing

Extend `AbstractExpression`:

[`Tensor`](expressions.md#tensor), [`BroadcastTo`](expressions.md#broadcastto), [`Indexed`](expressions.md#indexed), [`ComponentTensor`](expressions.md#componenttensor-and-indexsum), [`IndexSum`](expressions.md#componenttensor-and-indexsum).

Extend `Base`:

[`GatheredIndex`](expressions.md#gatheredindex).

### Tensor operations

Extend `AbstractExpression`:

[`Inverse`](expressions.md#inverse), [`Determinant`](expressions.md#determinant).

### Domains

Extend `AbstractCoordinateDomain`:

[`RD`](topology.md#rd).

### Reference cells

Extend `Base`:

[`CellMap`](topology.md#cellmap).

### Restriction

Extend `Base`:

[`SideMap`](topology.md#sidemap).

### Meshes

Extend `AbstractDomain`:

[`Mesh`](geometry.md#mesh).

### Calculus on a domain

Extend `AbstractExpression`:

[`Coordinate`](geometry.md#coordinate), [`Gradient`](geometry.md#gradient), [`Divergence`](geometry.md#divergence-and-curl), [`Curl`](geometry.md#divergence-and-curl), [`DirectionalDerivative`](geometry.md#directionalderivative), [`Integral`](geometry.md#integral).

### Maps between domains

Extend `Base`:

[`DomainMap`](geometry.md#domainmap).

Extend `AbstractExpression`:

[`PullBack`](geometry.md#pullback).

### Fields and integration on a mesh

Extend `Base`:

[`TransformationLaw`](geometry.md#transformationlaw).

Extend `AbstractExpression`:

[`PushForward`](geometry.md#pushforward), [`Normal`](geometry.md#normal).

### Function spaces

Extend `AbstractExpression`:

[`Function`](forms.md#function).

### Variational forms

Extend `AbstractFunctionSpace`:

[`ConjugateSpace`](forms.md#conjugatespace), [`DualSpace`](forms.md#dualspace).

Extend `AbstractExpression`:

[`Form`](forms.md#form), [`Apply`](forms.md#apply).
<!-- end generated: class-hierarchy -->

## Appendix: Free functions

<!-- generated: free-functions, run design/appendices.py to update -->
### Expressions and static type

| Function                               | Overload | Takes                 | Returns |
| -------------------------------------- | -------- | --------------------- | ------- |
| [`index(label)`](expressions.md#index) |          | `label: str \| tuple` | `Index` |

### Numbers and scalar operations

| Function                                                    | Overload | Takes                               | Returns                    |
| ----------------------------------------------------------- | -------- | ----------------------------------- | -------------------------- |
| [`power(first, second)`](expressions.md#power)              | `**`     | `first, second: AbstractExpression` | `Power`                    |
| [`negate(f)`](expressions.md#negation-and-subtraction)      | `-f`     | `f: AbstractExpression`             | `multiply(Integer(-1), f)` |
| [`subtract(f, g)`](expressions.md#negation-and-subtraction) | `f - g`  | `f, g: AbstractExpression`          | `add(f, negate(g))`        |
| [`absolute(f)`](expressions.md#absolute)                    | `abs(f)` | `f: AbstractExpression`             | `Absolute`                 |
| [`conjugate(f)`](expressions.md#conjugate)                  |          | `f: AbstractExpression`             | `Conjugate`                |
| [`square_root(f)`](expressions.md#elementary-functions)     |          | `f: AbstractExpression`             | `SquareRoot`               |
| [`logarithm(f)`](expressions.md#elementary-functions)       |          | `f: AbstractExpression`             | `Logarithm`                |

### Tensors and indexing

| Function                                                                               | Overload | Takes                                                                         | Returns                                                                                                                                                         |
| -------------------------------------------------------------------------------------- | -------- | ----------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`tensor(entries)`](expressions.md#tensor)                                             |          | `entries: tuple[AbstractExpression, ...]`                                     | `Tensor`                                                                                                                                                        |
| [`broadcast_to(value, value_shape)`](expressions.md#broadcastto)                       |          | `value: AbstractExpression`, `value_shape: tuple[int \| AbstractExtent, ...]` | `BroadcastTo`                                                                                                                                                   |
| [`indexed(expression, indices)`](expressions.md#indexed)                               |          | `expression: AbstractExpression`, `indices: tuple`                            | `Indexed`                                                                                                                                                       |
| [`component(f, entries)`](expressions.md#partial-indexing)                             | `[]`     | `f: AbstractExpression`, `entries: tuple`                                     | `f` if `entries == ()`; otherwise `component_tensor(indexed(f, entries + K), K)` if the entries are fewer than the rank of $f$, otherwise `indexed(f, entries)` |
| [`component_tensor(expression, indices)`](expressions.md#componenttensor-and-indexsum) |          | `expression: AbstractExpression`, `indices: tuple[Index, ...]`                | `ComponentTensor`                                                                                                                                               |
| [`index_sum(expression, indices)`](expressions.md#componenttensor-and-indexsum)        |          | `expression: AbstractExpression`, `indices: tuple[Index, ...]`                | `IndexSum`                                                                                                                                                      |
| [`gathered_index(table, indices)`](expressions.md#gatheredindex)                       |          | `table: AbstractIndexTable`, `indices: tuple`                                 | `GatheredIndex`                                                                                                                                                 |

### Tensor operations

| Function                                                            | Overload | Takes                               | Returns                                                                            |
| ------------------------------------------------------------------- | -------- | ----------------------------------- | ---------------------------------------------------------------------------------- |
| [`add(first, second)`](expressions.md#tensor-arithmetic)            | `+`      | `first, second: AbstractExpression` | `Sum(first, second)`                                                               |
| [`multiply(first, second)`](expressions.md#tensor-arithmetic)       | `*`      | `first, second: AbstractExpression` | `Product(first, second)`                                                           |
| [`divide(first, second)`](expressions.md#tensor-arithmetic)         | `/`      | `first, second: AbstractExpression` | `Division(first, second)`                                                          |
| [`real_part(f)`](expressions.md#tensor-arithmetic)                  |          | `f: AbstractExpression`             | `component_tensor(RealPart(f[K]), K)`                                              |
| [`imaginary_part(f)`](expressions.md#tensor-arithmetic)             |          | `f: AbstractExpression`             | `component_tensor(ImaginaryPart(f[K]), K)`                                         |
| [`matrix_product(first, second)`](expressions.md#tensor-arithmetic) | `@`      | `first, second: AbstractExpression` | The contraction of the last axis of `first` with the first axis of `second`, below |
| [`transpose(A)`](expressions.md#tensor-arithmetic)                  |          | `A: AbstractExpression`             | `component_tensor(A[k, l], (l, k))`, with `(k, l)` = `K`, for $A$ of rank 2        |
| [`inverse(A)`](expressions.md#inverse)                              |          | `A: AbstractExpression`             | `Inverse`                                                                          |
| [`determinant(A)`](expressions.md#determinant)                      |          | `A: AbstractExpression`             | `Determinant`                                                                      |

### Topology

| Function                                           | Overload | Takes                           | Returns |
| -------------------------------------------------- | -------- | ------------------------------- | ------- |
| [`cell_index(T, j)`](topology.md#abstracttopology) |          | `T: AbstractTopology`, `j: int` | `Index` |

### Reference cells

| Function                                   | Overload | Takes                                                                       | Returns                                                         |
| ------------------------------------------ | -------- | --------------------------------------------------------------------------- | --------------------------------------------------------------- |
| [`cell_map(X, j, e)`](topology.md#cellmap) |          | `X: AbstractTopology \| Mesh`, `j: int`, `e: int \| Index \| GatheredIndex` | `CellMap(X, j, e)`                                              |
| [`cell_map(X, j)`](topology.md#cellmap)    |          | `X: AbstractTopology \| Mesh`, `j: int`                                     | `CellMap(X, j, cell_index(T, j))`, with $T$ the topology of $X$ |

### Restriction

| Function                                   | Overload | Takes                                      | Returns            |
| ------------------------------------------ | -------- | ------------------------------------------ | ------------------ |
| [`side_map(X, P, s)`](topology.md#sidemap) |          | `X, P: AbstractTopology \| Mesh`, `s: str` | `SideMap(X, P, s)` |

### Meshes

| Function                         | Overload | Takes                                                  | Returns |
| -------------------------------- | -------- | ------------------------------------------------------ | ------- |
| [`mesh(x, Q)`](geometry.md#mesh) |          | `x: AbstractExpression`, `Q: AbstractCoordinateDomain` | `Mesh`  |

### Calculus on a domain

| Function                                                                  | Overload | Takes                                                                                | Returns                             |
| ------------------------------------------------------------------------- | -------- | ------------------------------------------------------------------------------------ | ----------------------------------- |
| [`coordinate(P)`](geometry.md#coordinate)                                 |          | `P: AbstractDomain`                                                                  | `Coordinate`                        |
| [`gradient(f, P)`](geometry.md#gradient)                                  |          | `f: AbstractExpression`, `P: AbstractDomain`                                         | `Gradient`                          |
| [`divergence(u, P)`](geometry.md#divergence-and-curl)                     |          | `u: AbstractExpression`, `P: AbstractDomain`                                         | `Divergence`                        |
| [`curl(u, P)`](geometry.md#divergence-and-curl)                           |          | `u: AbstractExpression`, `P: AbstractDomain`                                         | `Curl`                              |
| [`directional_derivative(f, x, a, h)`](geometry.md#directionalderivative) |          | `f: AbstractExpression`, `x: Function \| tuple`, `a, h: AbstractExpression \| tuple` | `DirectionalDerivative(f, x, a, h)` |
| [`integral(f, P)`](geometry.md#integral)                                  |          | `f: AbstractExpression`, `P: AbstractDomain`                                         | `Integral`                          |

### Maps between domains

| Function                                                  | Overload | Takes                                                                                                       | Returns                                            |
| --------------------------------------------------------- | -------- | ----------------------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| [`domain_map(P, psi, Q)`](geometry.md#domainmap)          |          | `P: AbstractCoordinateDomain \| AbstractTopology`, `psi: AbstractExpression`, `Q: AbstractCoordinateDomain` | `DomainMap`                                        |
| [`pull_back(f, D)`](geometry.md#pullback)                 |          | `f: AbstractExpression`, `D: DomainMap \| SideMap \| CellMap`                                               | `PullBack`                                         |
| [`jacobian(D)`](geometry.md#jacobian-and-density)         |          | `D: DomainMap \| CellMap`                                                                                   | `gradient(D.parameterisation, D.parameter_domain)` |
| [`inverse_jacobian(D)`](geometry.md#jacobian-and-density) |          | `D: DomainMap \| CellMap`                                                                                   | `inverse(jacobian(D))`                             |
| [`density(D)`](geometry.md#jacobian-and-density)          |          | `D: DomainMap \| CellMap`                                                                                   | The density $\rho$ below                           |

### Fields and integration on a mesh

| Function                                                                 | Overload | Takes                                                                         | Returns                                    |
| ------------------------------------------------------------------------ | -------- | ----------------------------------------------------------------------------- | ------------------------------------------ |
| [`transformation_law(variances, weight)`](geometry.md#transformationlaw) |          | `variances: tuple[Literal["covariant", "contravariant"], ...]`, `weight: int` | `TransformationLaw`                        |
| [`push_forward(f, D, m)`](geometry.md#pushforward)                       |          | `f: AbstractExpression`, `D: Mesh`, `m: TransformationLaw`                    | `PushForward(f, D, m)`                     |
| [`side_domain(D, S, s)`](geometry.md#side-domains)                       |          | `D: Mesh`, `S: AbstractTopology`, `s: str`                                    | `mesh(pull_back(x, side_map(T, S, s)), Q)` |
| [`normal(sm)`](geometry.md#normal)                                       |          | `sm: SideMap`                                                                 | `Normal`                                   |

### Function spaces

| Function                                  | Overload | Takes                                             | Returns    |
| ----------------------------------------- | -------- | ------------------------------------------------- | ---------- |
| [`function(V, label)`](forms.md#function) |          | `V: AbstractFunctionSpace`, `label: str \| tuple` | `Function` |

### Variational forms

| Function                                              | Overload | Takes                                                                                                                 | Returns                                                                                 |
| ----------------------------------------------------- | -------- | --------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| [`conjugate_space(V)`](forms.md#conjugatespace)       |          | `V: AbstractFunctionSpace`                                                                                            | `V` for a real `V`, `V.space` for a `ConjugateSpace`, and `ConjugateSpace(V)` otherwise |
| [`dual_space(S)`](forms.md#dualspace)                 |          | `S: AbstractFunctionSpace \| tuple[AbstractFunctionSpace, ...]`                                                       | `S.space` for a `DualSpace`, otherwise `DualSpace(S)`                                   |
| [`form(expression, arguments)`](forms.md#form)        |          | `expression: AbstractExpression`, `arguments: tuple[Function \| Conjugate \| tuple[Function \| Conjugate, ...], ...]` | `Form`                                                                                  |
| [`apply(f, arguments)`](forms.md#apply)               |          | `f: AbstractExpression`, `arguments: tuple[AbstractExpression \| tuple[AbstractExpression, ...], ...]`                | `f` for `apply(f, ())` when `f.argument_spaces == ()`, otherwise `Apply(f, arguments)`  |
| [`partial(F, mapping)`](forms.md#partial-application) |          | `F: Form`, `mapping: dict`                                                                                            | The form with selected outer entries supplied, or its body when none remain             |
| [`rebind(b, mapping)`](forms.md#rebinding)            |          | `b: Form \| ComponentTensor \| IndexSum`, `mapping: dict`                                                             | The binder with selected declarations renamed                                           |

### Pass recipes

| Function                                                            | Overload | Takes                   | Returns                                                  |
| ------------------------------------------------------------------- | -------- | ----------------------- | -------------------------------------------------------- |
| [`expand_arithmetic(f)`](evaluation.md#expanding-tensor-arithmetic) |          | `f: AbstractExpression` | An equivalent expression with tensor arithmetic expanded |

### Standard lowering recipe

| Function                                             | Overload | Takes                   | Returns                                      |
| ---------------------------------------------------- | -------- | ----------------------- | -------------------------------------------- |
| [`lower(f)`](evaluation.md#standard-lowering-recipe) |          | `f: AbstractExpression` | An equivalent expression in core normal form |
<!-- end generated: free-functions -->

## Appendix: Static type table

This normative table gives the argument spaces of each node class. An operand
not named in the argument spaces column must have the argument spaces `()`,
except binders and functional values accepted by `Form`, `Apply` and the
derivatives. Their signatures are checked by the respective class contracts
and are not added to the enclosing expression's argument spaces.

| Class                                    | Argument spaces of the node                                    |
| ---------------------------------------- | -------------------------------------------------------------- |
| `Integer`, `Float`                        | `()` |
| `Complex`                                | Those of the parts, which must be equal |
| `Zero`                                   | The stored argument spaces                                     |
| `Sum`                                    | Those of the operands, which must be equal                     |
| `Product`                                | Those of the operands, concatenated                            |
| `Power`, `Absolute`                      | `()`                                                           |
| `RealPart`, `ImaginaryPart`              | `()`                                                           |
| `Division`                               | Those of the first operand                                     |
| `Conjugate`                              | Those of the operand, each conjugated                          |
| `SquareRoot`, `Logarithm`, `Inverse`, `Determinant` | `()`                                                           |
| `Tensor`                                 | Those of the entries, which must be equal                      |
| `BroadcastTo` | Those of the value |
| `Indexed`, `ComponentTensor`, `IndexSum` | Those of the expression                                        |
| `Coordinate`                             | `()`                                                           |
| `Normal`                                 | `()`                                                           |
| `Function`                               | Those of the space                                             |
| `Form`                                   | Those of the expression, followed by the space of each entry, conjugated for a `conjugate(v)` entry |
| `Apply`                                 | `()`              |
| `Gradient`, `Divergence`, `Curl`         | Those of the operand                                           |
| `PushForward`                            | Those of the expression                                        |
| `DirectionalDerivative`                  | Those of the expression                                        |
| `Integral`                               | Those of the integrand                                         |
| `PullBack`                               | Those of the expression                                        |

The value shape and the free indices of each node class are in the section
of the class. The domain and the scalar type of the node classes that differ
from the default rules of [Domain](expressions.md#domain) and [Scalar type](expressions.md#scalar-type):

| Class                    | Domain of the node                | Scalar type of the node           |
| ------------------------ | ---------------------------------- | --------------------------------- |
| `Complex`          | None                               | Complex                           |
| `Zero`                   | None                               | The stored scalar type            |
| `Absolute`               | The domain of the operand         | Real                              |
| `RealPart`, `ImaginaryPart` | The domain of the operand      | Real                              |
| `Coordinate`             | The domain $P$                     | Real                              |
| `Gradient`, `Divergence`, `Curl` | The domain $P$               | The scalar type of the expression |
| `DirectionalDerivative`  | The domain of $f$                 | The scalar type of $f$            |
| `Integral`               | None                               | The scalar type of the integrand  |
| `PullBack`               | The parameter domain of the map | The scalar type of the expression |
| `Normal`                 | $E$                                | Real                              |
| `Function`               | The domain of the space, or none   | The scalar type of the space      |
| `PushForward`            | The mesh $D$                       | The scalar type of the expression |
| `Form`                   | The domain of the expression      | The scalar type of the expression |
| `Apply`                 | The domain of the form            | The scalar type of the form       |

## Appendix: Summary

Counts of the core, from the sections.

<!-- generated: summary, run design/appendices.py to update -->
| Quantity                        | Count |
| ------------------------------- | ----- |
| Classes                         | 53    |
| Abstract classes                | 10    |
| Concrete classes                | 43    |
| Expression node classes         | 33    |
| Free functions                  | 53    |
| Free functions with an overload | 9     |
<!-- end generated: summary -->

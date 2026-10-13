# Core UFLx language definition

This note describes the smallest core language behind UFLx. It intentionally
avoids implementation details, compiler algorithms, and broad redesign ideas.
It assumes knowledge of differentiable geometry and functional analysis. Some
advanced differential geometry concepts are explained in more detail.

This document is aimed at users and developers with extensive knowledge of
domain specific languages for variational forms (e.g. UFL, ).

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

UFLx does not handle concrete data or algorithms; e.g. actual meshes, boundary
markers, numerical values of coefficients, or assembly procedures. Those belong
to the surrounding problem-solving environment, such as DOLFINx, Firedrake,
FIAT and Basix. UFLx instead provides a specified set of leaf nodes which the
external library must provide.

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

May be necessary:

1. No attribute may be assigned, reassigned or deleted after construction,
   except to cache state derived according to the rule beginning (any state
   beyond...).

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

#### Members vs free functions principles

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

### Operators

```UnaryOperator(AbstractExpression)```
```BinaryOperator(AbstractExpression)```

Obvious.

## Differential geometry

We begin with differential geometry because the function spaces and operators
of UFLx are defined on a domain. The function spaces and operators require two
kinds of structure on the domain:

- topological structure, which is the dimension, boundary, orientation and
  decomposition into cells of the domain,
- geometric structure, which is the metric of the domain and the measure
  induced by the metric.

### Cells

A domain is divided into finitely many cells. A cell is a piece of the domain
with the same dimension as the domain, and with the shape of a ball up to
continuous deformation. Two cells meet only on the boundaries of the two cells.

A facet is a piece of the boundary of a cell, one dimension lower than the
cell. A facet in the boundary of one cell is an exterior facet. A facet shared
by two cells is an interior facet.

Cells and facets are topological. Defining cells and facets does not require a
metric.

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

A form is a multi-linear functional $a(u, v, w, \\ldots): U \\times V \\times W
\\times \\ldots \\to {\\mathbb{R}, \\mathbb C}$.

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

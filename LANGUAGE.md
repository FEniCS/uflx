# Minimal UFLx language definition

This note describes the smallest useful language model behind UFLx. It intentionally
avoids implementation details, compiler algorithms, and broad redesign ideas.

## 1. What UFLx is

UFLx is a symbolic language embedded in Python for writing finite-element variational
forms.

A UFLx program is ordinary Python code that constructs symbolic objects. The important
result is usually a `Form`, such as a bilinear form, linear form, or functional.

UFLx is not a full PDE language. It does not own the actual mesh, boundary markers,
numerical values of coefficients, or assembly procedure. Those belong to the surrounding
problem-solving environment, such as DOLFIN or DOLFINx.

In one sentence:

> UFLx is a typed symbolic expression language for finite-element forms.

## 2. The minimal object model

The minimal language needs only a small number of object kinds.

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

The element describes what kind of discrete field is being represented. It does not
store actual coefficient values.

### Domains

A domain is the integration domain of a form. For a minimal language definition, a
domain only needs to provide:

- topological dimension,
- geometric dimension,
- and a parametrization per cell type.

A parametrization is a map out of a cell's coordinate domain, which UFLx can evaluate
and differentiate but does not otherwise interpret. A finite element basis summed
against a mesh's coordinate dofs is one way to describe such a map; a closed form
expression is another, and neither is privileged. Two parametrizations compose when the
first lands where the second starts, so a mesh can be carried onto a surface given in
closed form rather than interpolated through the nodes of a higher degree element.

A domain may contain cells of more than one type, and then it has one parametrization
per type. Integrating over such a domain is integrating over each of its cell types and
adding the results, since each type has its own coordinate domain and its own map out of
it. Pulling such an integral back to the cells' coordinates therefore gives a sum of
integrals, one per cell type.

A domain may also be presented as a union of parametrized pieces, each parametrized in its own
right. That is what the facets of a mesh are: a piece per cell type and local facet,
each the inclusion of the reference facet into the cell followed by the cell's own map.
Integrating over a union is integrating over each piece and adding the results, as for
cell types. A union is not an atlas: its pieces carry no claim of a relation to each
other, where an atlas's parametrizations cover one set and are related by transition maps.

The actual mesh is still external to UFLx. So is which facets a union of them has: that
a facet lies on the boundary is a question about a cell complex, so the domain is named
for the sub-entities it is parametrized by and claims no boundary. A domain carries an opaque
label for a consumer to tell one such selection from another.

### Function spaces

A function space defines a finite-dimensional space of functions. There are three types
of function spaces:

1. standard finite-element spaces, defined by a domain and an element,
1. constant function spaces, defined by a shape and a scalar type,
1. non-finite-element function spaces, defined by domain, shape and scalar type, but no
   element.

Examples:

```python
V = FunctionSpace(domain=domain, element=element)  # standard finite-element space
C = FunctionSpace(shape=(2,), scalar_type="real")  # constant function space
N = FunctionSpace(
    domain=domain, shape=(2,), scalar_type="real"
)  # non-finite-element function space
```

In older UFL examples, `Argument` and `Coefficient` are sometimes created directly from
an element. This is not allowed.

### Terminals

Terminals are the leaves of UFLx expression trees.

The essential terminals are:

- `Argument`, representing an arbitrary basis function of a function space,

  It has a `number` attribute to represent the position of the `Argument` in the form.
  For example, consider a bilinear form `a(u, v)`. The `Argument` representing `u` has
  `number=1`, and the `Argument` representing `v` has `number=0`.

  > Experimental (API breaking) note: Existence of an `Argument` and a `Coefficient` is
  > rooted deep into the legacy UFL. However, more useful concept would be to only allow
  > for `Function` terminals, which can represent both known and unknown functions. It
  > would be a responsibility of the user, when creating a form, to define which of the
  > `Function` terminals are to be considered as "arguments" and which not. This would
  > be much closer to a mathematical notation and would solve many user error of having
  > to treat `Argument` and `Coefficient` differently.

- `Coefficient`, representing a known function in a function space,

- literals such as numbers,

- and geometric quantities such as coordinates, normals, Jacobians, cell diameter, facet
  area, etc.

Example:

```python
v = Argument(function_space=V, name="v", number=0)
u = Argument(function_space=V, name="u", number=1)
f = Coefficient(function_space=V, name="f")
```

> Legacy note: There is no `Constant` in UFLx as a means to represent literal value
> which is not known at compile time.

UFLx supports named symbolic terminals. This allows for better error messages and better
postprocessing of expressions. In addition, it improves the stability of the expression
tree.

> Legacy note: In legacy UFL, symbolic objects do not have names. It requires expensive
> and error prone renumbering of terminals in expressions before lowering and code
> generation.

### Expressions

An expression is either a terminal or an operator applied to other expressions.

Examples of operators are addition, multiplication, indexing, `grad`, `div`, `dot`,
`inner`, `jump`, and `avg`.

Expressions are symbolic. They are not evaluated when they are created.

A minimal expression carries these static properties:

- tensor shape,
- free indices,
- associated domain or domains,
- scalar type (real or complex),
- continuity information,
- (possibly) physical dimension.

These properties are enough to reject many meaningless expressions before lowering or
code generation.

### Measures

A measure is the domain integrated over paired with the density used on it. Neither half
is optional. A function is not a differential form, so integrating one needs a density,
and a density is a density of something.

The density is not stated by the caller. `dx(domain)` is the measure the domain's own
parametrization induces, `sqrt(det g)`, so that lengths and volumes measured in a
parametrization's coordinates agree with the ambient ones. Where the parametrization is
the identity, as it is on R^d or on a region of it, that is one and the measure is the
Lebesgue measure.

There is no way to weigh a measure by a function, and this is deliberate. A factor that
is geometric is the density of some parametrization, so it arrives by presenting the
domain through that parametrization rather than by weighting: the `r` of an axisymmetric
problem is `sqrt(det g)` of a cylindrical parametrization, not a weight on a
two-dimensional one. A factor that is not geometric is part of the physics and belongs
in the integrand. What does lie outside an induced density is a measure of a different
kind rather than a weighted one, such as the Dirac measure a point evaluation integrates
against.

```python
dx(mesh)  # the measure of a mesh
dx(mesh).domain  # mesh
dx(mesh).density  # VolumeElement(mesh)
```

A measure always names its domain: there is no measure of no domain in particular, and
an integral never guesses which domain was meant. Pulling an integral back onto a cell's
coordinates moves its measure onto that cell, which is what change of variables does to
a measure.

`ds` and `dS` are not kinds of measure. In UFL they are, because the measure does not
name its domain and so has to encode it; here an exterior-facet integral is this same
measure over a domain of facets, and an interior-facet integral is it over an interface.

The first of those is constructible. `facets_of(mesh)` is the domain a mesh's facets
make up, and an exterior-facet integral is the one measure over it, with the fields in
the integrand read there through its parametrization:

```python
facets = facets_of(mesh)
inner(restrict(grad(u), facets), UnitNormal(facets)) * dx(facets)
```

Pulling that back gives one integral per facet over the reference facet, each carrying
that facet's volume element as its density. No part of it is a facet measure.

An interface is not constructible yet: it needs a parametrization per side. The `+`/`-`
restriction below belongs to the integrand either way, being a field pulled back along
one side's inclusion, and not to the measure.

### Integrals

An integral is made by multiplying a scalar expression by a measure.

Examples:

```python
v * f * dx(mesh)
inner(grad(v), grad(u)) * dx(mesh)
restrict(v, facets_of(mesh)) * g * dx(facets_of(mesh))  # an exterior facet integral
jump(v) * dS  # once an interface domain exists
```

A valid integral has a scalar-valued integrand with no unresolved free indices.

### Forms

A form is a multi-linear functional $a(u, v, w, \\ldots): U \\times V \\times W \\times
\\ldots \\to {\\mathbb{R}, \\mathbb C}$.

Arguments determine the arity of a form.

> Legacy note: In contrast to the legacy UFL, in UFLx a form does not need to be defined
> as a sum of integrals. For example, the following is a valid UFLx form:
> `a = (u * dx(mesh)) * (v * dx(mesh))`, or
> `a = PointEvaluation(u, (0.5, 0.5)) * PointEvaluation(v, (0.5, 0.5))`.

Examples:

```python
a = inner(grad(v), grad(u)) * dx(mesh)
L = v * f * dx(mesh)
F = a - L
```

> Experimental (API breaking) note: Similar to the `Argument` and `Coefficient`
> terminals, the `Form` in the legacy UFL is defined implicitly, by the presence of
> `Argument` terminals. In UFLx, a `Form` could be an explicit object that can be
> created by the user, e.g.

```python
a = Form(inner(grad(v), grad(u)) * dx(mesh), arguments=[u, v])
```

## 3. Minimal validity rules

The language is small, but it is typed.

- Addition and subtraction require compatible shapes and compatible free indices.
- Multiplication follows scalar, tensor, matrix-vector, matrix-matrix, and
  indexed-expression rules.
- Tensor contractions such as `dot` and `inner` require compatible tensor shapes.
- Spatial derivatives require a spatial domain.
- Facet quantities such as a normal are read off the domain they are a normal to, so a
  facet's normal is `UnitNormal(facets_of(mesh))`, and a normal is only defined in
  codimension one.
- A field on one domain is read on another through a restriction, `restrict(u, facets)`,
  which is the pullback of the field along that domain's parametrization.
- Interior-facet expressions may use side restrictions:
  ```python
  u("+")
  u("-")
  ```
  which are that same restriction along one side's inclusion, once a domain has a parametrization
  per side.
- Operators such as `jump` and `avg` are shorthand for combinations of such
  restrictions.
- An integral is valid only when the final integrand is scalar and has no free indices.
- A form is valid when all of its integrals are valid and their domains, arguments, and
  coefficients are mutually consistent.

## 4. What transformation procedures does the UFLx core provide

There are two key transformation procedures that the UFLx core provides:

1. Automatic/symbolic differentiation.

   Similar to UFL, UFLx provides a symbolic differentiation via a directional (Gateaux)
   derivative. The derivative of an Expression, Integral or a Form is lazy/eagerly
   evaluated. That means the derivative is stored as an expression tree operator, and
   the actual differentiation is performed when the derivative is used in a context
   where it needs to be evaluated, such as lowering or code generation.

   Example:

   ```python
    a = inner(grad(v), grad(u)) * dx(mesh)
    da = derivative(a, u, du)
   ```

1. Expression and integral transformation between different configurations.

   > Legacy note: In legacy UFL, there is a very central perspective on the
   > configuration of where objects are defined. Usually, all expressions are defined in
   > the physical configuration, and UFL applies pullbacks to the entire expression
   > tree.

   In UFLx, configuration where expressions are defined are a consequence of the domain.
   Since the domain is part of the Expression's static attributes, it is propagated
   through the expression tree based on the language rules.

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

This is enough to describe ordinary weak forms such as mass, stiffness, elasticity,
Stokes, mixed Poisson, and DG forms.

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

These are consumers or transformations of the symbolic language, not the core language
itself.

## 7. Extensibility of the attributes

> Legacy note: In legacy UFL, the set of static attributes is fixed and hardcoded. This
> led to many language features being hacked in (e.g. complex numbers support, dual
> spaces, etc.).

In UFLx, the set of static attributes is extensible. This allows for cleaner language
extensions in the future, such as support for physical dimensions, dual spaces, error
estimates, custom quadrature rules, custom integrability/continuity rules, etc.

## 8. Examples

### Poisson equation

```python
parametrization = FiniteElement("Lagrange", triangle, 1, shape=(2,))
element = FiniteElement("Lagrange", triangle, 1)

mesh = Mesh(parametrization)
V = FunctionSpace(mesh, element)

u = TrialFunction(V, name="u")
v = TestFunction(V, name="v")
f = Coefficient(V, name="f")

a = inner(grad(v), grad(u)) * dx(mesh)
L = v * f * dx(mesh)
```

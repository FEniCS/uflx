# Design choices

Documentation of design decisions that have been made

## Naming

File names are plural to allow (eg) `uflx.domains` to be the module and `uflx.domain` to be
the function defined in domains.py that is imported into __init__.py.

Whenever a function to initialise a class is defined, it's name is the same as the class name
but lowercase, with `_x` wherever the class name has (eg) a capital `X` mid name.

## Default values

Class `__init__` functions have no default values. If default values are required, a function
to initialise the class with defaults should be defined.

## Mixins and protocols

A class gets something it did not write in two ways, and they are for different things.

A **mixin** is inherited and supplies a rule that holds of every value. There is one,
`Immutable` in `uflx/graphs/graphs.py`, inherited by every abstract root
(`AbstractExpression`, `AbstractEntity`, `AbstractFiniteElement`, `AbstractFunctionSpace`,
`AbstractDomain`, `AbstractReferenceMap`, `AbstractMeasure`, `AbstractIntegral`), so a
consumer's cell or element inherits it too. It says that a value is built once and never
changed, that its hash is computed once and kept, and that its identity is its class and
what it was built from: a class that says `init_args` gets `__eq__`, `__hash__` and
`__repr__` from it and writes none of them, and one whose identity is deliberately
something else writes its own and says why. A new mixin is for a rule of that kind, one
that holds of everything; code shared between a few classes is a function.

A **protocol** names a question a node can be asked, and is never inherited. It is a
`runtime_checkable` `Protocol`, so a node implements it by having the attribute and an
algorithm asks with `isinstance`; a node that cannot answer is not an instance, and the
asker needs no default. That keeps the layering one-way: a protocol lives in the lowest
module that asks it, which may be below the module that implements it. It is named as
the thing it says about the node, a verb phrase (`ComposesADomain`, `ReadAtAVariable`) or
an adjective (`ComplexValued`, `Commutative`). The protocols are:

- `uflx/graphs/graphs.py`: `GraphNode`, what every node is (`successors`, `init_args`);
  `RepresentedByGraph`, an object whose construction is a graph.
- `uflx/algorithms/`: `PullBackToReference`, a node that pulls itself back;
  `ExtraReplacement`, a node that says what replaces it; `Simplifiable`,
  `SimplifiableInProduct`, `SimplifiableInSum`, `SimplifiableInMatrixProduct` and
  `RightSimplifiableInMatrixProduct`, a node that simplifies alone or against a
  neighbour; `Commutative`, one whose operands may be reordered.
- `uflx/complex.py`: `ComplexValued`. `uflx/geometry.py`: `ExpandableGeometry`.
  `uflx/maps.py`: `IsPushedForward`, `IsPulledBack`.

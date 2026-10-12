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
changed, and that its hash is computed once and kept. It also says that a value's identity
is its class and its `init_args`, and it supplies `__eq__`, `__hash__`, `__repr__` and
`init_args` itself. A new mixin is for a rule of that kind, one
that holds of everything. Code shared between a few classes is a function.

Avoid defining `__eq__` or `__hash__` on a class at all costs. An identity other than the
arguments is one that rebuilding a node can silently change, and equality that reaches
outside the class breaks hashing. A class that seems to need one is usually missing an
argument or keeping one it should not.

`init_args` is read back from the constructor. Each argument `name` is stored as
`self._name`, and the mixin reads those attributes in the order of the constructor's
signature. A read-only property exposes one that should be public. A class with no
constructor of its own has no arguments. Three rules follow.

- Store an argument in the form it should be compared in, such as a tuple rather than a
  list. Then `type(v)(*v.init_args)` rebuilds `v`, which `reconstruct_node` relies on.
- Store every argument. Two values with equal `init_args` are the same value, so
  anything left out is lost.
- Store only what is hashable and compared exactly, all the way down. A tolerance is not
  transitive, so it cannot decide which values are equal. An array is kept as nested
  tuples of its entries.

A derived value, such as a lookup table built from an argument, goes under a name of its
own. A constructor that takes `*args`, `**kwargs` or keyword-only arguments cannot be
called with a tuple, so its class defines `init_args` itself. So does a class that stores
an argument in a way the rules above do not cover, and it says why in a comment.

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

# Language rules

## Conventions

### Normative wording

| Wording | Meaning |
| ------- | ------- |
| **Must** | An implementation or library has to satisfy the rule. |
| **Should** | An implementation or library satisfies the rule unless the note states a reason not to. |
| **May** | A choice left open to the implementation or library. |
| **Present tense**, such as “the value shape is a `tuple`” | A rule of the language. |
| **Not yet covered** | A decided gap that the authors intend to specify later. |
| **Open** | A question the authors have not decided. |
| **Legacy UFL** dropdown | A comparison with legacy UFL, not a rule. |

### Mathematical admissibility

Static type checking establishes the documented Python argument types. Construction
checks the remaining structural and value constraints. Values and rewrites assume the required
mathematical hypotheses; the core neither verifies them nor extends formulas
through singularities.

### Glossary

| Term | Meaning |
| ---- | ------- |
| [Leaf](#leaves-and-subclassing) | A primitive object rather than an operation; it may have graph dependencies. |
| [Terminal](#leaves-and-subclassing) | A leaf that is an expression, such as a literal, `Coordinate`, `Function` or `Normal`. |
| [Construction free function](expressions.md#node-class-or-free-function) | A pure free function that may validate inputs and builds fixed expression compositions or canonical space descriptions. |
| [Pass](#passes) | A pure free function that walks a graph after construction to compute a result or rewrite the graph. |
| [Evaluation traversal](evaluation.md#evaluation-traversal) | A walk through expression operands, excluding dependencies stored in geometry, spaces and labels. |
| [Free occurrence](evaluation.md#free-occurrences) | An occurrence of a name outside its binding scope, on the walk chosen by the operation; construction dependencies govern capture and geometry-dependence checks. |
| [Independence](evaluation.md#independence) | An operand's value is unchanged when the operation's target varies with other inputs fixed. |
| [Domain-free expression](expressions.md#domain-free-expressions) | An expression with no spatial argument, interpreted as a constant field where a domain is required. |
| [Grammar pattern](evaluation.md#grammars) | A structural condition on a node and objects within a bounded number of edges, independent of pass results. |
| [Rule pattern](extensions.md#rule-tables) | A condition selecting a local rewrite from node structure, static types or explicit analyses. |
| [Strategy](evaluation.md#transformations-recipes-and-contracts) | Selection and scheduling of rewrite rules, including traversal and stopping conditions. |
| [Recipe](evaluation.md#transformations-recipes-and-contracts) | A composition of strategies with stated input, output and preservation contracts. |
| [Array rank](expressions.md#abstractexpression), or **rank** | The number of axes in an expression's value shape. |
| [Arity](forms.md#argument-spaces) | The number of argument spaces of a form, counting a product-space argument once. |
| [Structurally equal](#the-base-object-base) | Objects with the same class and equal construction arguments, as tested by `==`. |
| Alpha-equivalent | Expressions differing only in consistent lexical renaming of bound names, with free names and independent parameter signatures preserved; distinct from structural equality. |
| [Mathematically equal](#the-base-object-base) | Expressions denoting the same value, whether or not their graphs are structurally equal. |

## Scope

### What UFLx is and isn't

UFLx is a typed symbolic language for writing variational forms.

The UFLx core language, henceforth the core, is independent of any particular
discretisation, solver or finite-element library. This draft requires vector calculus on
domains in $R^n$ and on manifolds embedded in $R^n$, in the coordinates of
$R^n$. One gradient serves every dimension, so the statement of each concept
does not depend on the dimension.

UFLx is *not* a symbolic language for reasoning about problems in functional
analysis or differential geometry. The core should therefore contain only the
concepts and relations needed to write variational forms as code for later
algorithmic transformation.

A UFLx program is ordinary Python code that builds a graph of symbolic objects,
usually a bilinear form, linear form or functional.

UFLx does not handle concrete data or numerical algorithms, such as meshes,
boundary markers, values of coefficients or assembly procedures. Those belong
to the surrounding problem-solving environment, such as DOLFINx, Firedrake,
FIAT and Basix. UFLx specifies a set of leaf nodes, and the external library
provides the concrete data behind each leaf node.

### Basic principles

- Each core object definition should follow the mathematical definition of
  the object closely.
- The core must not add definitions, properties or terminology outside
  differential geometry and functional analysis.
- The core values strictness and expressivity over brevity. Higher-level APIs
  can add brevity later.
- Use dimension-independent definitions where the mathematics permits.
  Operations such as `Curl` state their dimension-specific cases explicitly.
- The core contains no algorithms solely for a form compiler's convenience.
  Each compiler states the forms it accepts. A pass outside the core must
  lower to those forms, using only core node classes and
  [language extensions](extensions.md#language-extensions) accepted by that compiler.
  For example, a basic form compiler might only accept integrals over a
  reference cell, such as a unit simplex, as in each cell term of a legacy UFL
  `dx` integral.

```{dropdown} Legacy UFL
TSFC lowers a UFL form to GEM, which repeats much of the algebra of UFL, partly because UFL simplifies during construction. The UFLx core does not simplify during construction, so a form compiler may be able to use a lowered UFLx form directly, adding only leaves for tables and kernel inputs. This has not been confirmed with the TSFC or FFCx developers.
```

## Rules of the language

### The base object `Base`

Every object in the directed acyclic graph (DAG) must extend `Base`. A graph
edge means "is defined from"; geometric edges are explicitly called edges of
a triangle or cell. Everything defining the variational problem must be
recoverable from the graph.

| `Base`     |                                                                    |
| ---------- | ------------------------------------------------------------------ |
| Properties | `construction_arguments`, the construction arguments, as in rule 4 |

Objects must satisfy the following rules:

1. A construction argument must be one of a small curated set of types: a base
   object, or exactly an `int`, a finite `float`, a `str`, or a `tuple` whose
   items are of these types. A subclass of these Python types is not accepted,
   so a `bool` is not an `int` and a named tuple is not a `tuple`. Construction
   arguments are deeply immutable. Their hashes must be deterministic across
   processes on the same machine, and they are compared by type and value.
   `1`, `1.0` and `True` are distinct, `0.0` and `-0.0` are distinct, and every
   value is equal to itself. The permitted argument types are checked statically
   under rule 8. Construction raises `TypeError` for a subclass where an exact
   built-in type is required, and `ValueError` for a `float` that is not finite.
2. Construction and all results depend only on the class and construction
   arguments, never on global state.
3. Objects are immutable after construction.
4. Reconstructing an object from its class and construction arguments must give
   an equal object. `construction_arguments` records those arguments in order
   as a `tuple`. A constructor may also store derived data.
5. The edges of an object must be exactly the base objects among its
   construction arguments, including the base objects nested in `tuple`s, in
   order.
6. Any order that the core or a pass relies on must be deterministic across
   processes. A pass, and any other walk of the graph such as hashing,
   serialisation or printing, follows the edges it visits in order (rule 5).
   A pass states which edges it visits. Where unrelated objects must be sorted,
   the sort uses a total order defined by the class and the construction
   arguments alone.
7. Equality (`==`) returns a `bool`: objects are equal exactly when their
   classes and construction arguments are equal. Whether `is` is structural
   is implementation dependent. An implementation in which `is` is
   structural must intern objects, so that equal objects are the same object.
   Structural equality does not decide mathematical equality: $a + b$ and
   $b + a$ are not structurally equal. Mathematical equalities are applied
   in later graph algorithms.
8. Requirements expressible to a static type checker must be checked statically
   and must not be duplicated by runtime validation. This includes argument
   classes, container member types, fixed tuple structure and permitted literal
   alternatives. Callers must satisfy the documented types; runtime behaviour
   for a statically invalid call is not specified. Constructors and construction
   free functions check the remaining structural and value requirements, such
   as compatible shapes, matching extents and distinct binders. They never
   convert inputs. Unless stated otherwise, failure of a runtime value
   requirement raises `ValueError`. Class tables list additional details or
   exceptions. Exact built-in types, including the exclusion of `bool` and
   named tuples in rule 1, and finite floats still require runtime checks.
9. A constructor has no default values, so every construction argument is
   written at construction.
10. A concrete class states an `identifier`, a `str` with the public name of
    the class qualified by the package name, such as `"uflx.Sum"` for the core
    or `"library.UnitSimplex"` for a library. The hash of rule 1, the total
    order of rule 6 and serialisation use the identifier and the construction
    arguments, and never the location of the class in the code. Construction
    raises `TypeError` for a class with no identifier. The package name keeps
    the identifiers of different libraries apart, and a library keeps the
    identifiers of the classes of the library distinct.

Construction validates immutable objects once; every existing object is structurally valid.

### Leaves and subclassing

Only leaves are subclassed. A leaf is a primitive object rather than an
operation, such as the coordinate domain `RD(n)` or a function whose values an
external library supplies. Leaves may have graph dependencies: a `Function`
depends on its space, and a `Coordinate` on its domain. A terminal is a leaf
that is an expression, such as a literal, `Coordinate`, `Function` or `Normal`
of the core, or a terminal of a downstream library.

An operation is represented by a single core class or a construction free
function returning core objects, never a subclass. For example, a construction
free function builds the map of an affinely mapped triangle as a
`DomainMap`, rather than a subclass of `DomainMap`. A language extension is the
one exception (see [Language extensions](extensions.md#language-extensions)).

Concrete core classes are final; libraries subclass only abstract classes.
A parent class never constructs a subclass, though its properties may return
existing subclass instances. A subclass may narrow a property's return type.

### Construction policy

| Kind | Role and permitted inspection | Validation and exceptions |
| ---- | ----------------------------- | ------------------------- |
| Property | Returns stored or derived data; creates no objects | Derived objects may be built and stored during construction |
| Constructor | Creates a node and computes derived data; does not simplify expressions | Checks constraints not established statically, as in rule 8 |
| Construction free function | Builds a fixed composition from expression static types and documented interfaces of explicit non-expression arguments; does not inspect expression structure to choose a composition | Checks constraints not established statically, itself or through constructors; may inspect structure for validation only. Space factories apply canonical forms; `indexed`, `component`, `component_tensor` and `index_sum` apply their empty-index identities; `apply(f, ())` returns an arity-zero `f` |
| Pass | After construction, inspects graph structure and returns information or a rewritten graph | Enforces its own preconditions; constructors validate nodes it builds. Space substitution preserves canonical forms |

`replace`, `partial` and `rebind` are explicit scope-aware transformation passes,
not construction factories. They may inspect declarations and traverse bodies
under their documented substitution rules; none evaluates unrelated lazy nodes.

Free functions and passes are pure. A class method immediately calls a free
function. Constructor arguments are positional only; the type checker rejects
named arguments. Free functions may take named arguments, such as
`domain_map(P=T, psi=psi, Q=RD(2))`.

A stored property returns a construction argument; a derived property depends
on the class and those arguments. Concrete classes inherit the derived
properties of their abstract classes.

### No simplification during construction

The construction policy preserves the user's written operations. For example,
in `(a + b) + c`, the user may need `(a + b)` evaluated before `c` is added
in floating-point arithmetic. Building a fixed composition defines an
operation; it does not simplify its operands. Later passes, chosen by the
user, may rewrite it.

### Passes

A pass is a pure free function that walks a graph after construction and returns
an analysis result or a rewritten graph. A rewrite reconstructs dependent
objects along the visited edges up to the root. Immutable objects allow
results to be cached per object.

Each pass states its traversal and purpose. Passes compose; the caller chooses
which run and in what order. The core substitution pass is
[`replace(e, mapping)`](evaluation.md#substitution-and-binding).

Unmatched nodes and failures follow the [strategy and contract rules](evaluation.md#transformations-recipes-and-contracts)
and the [extension policy](extensions.md#rule-tables). Analyses may report
operation-specific errors, such as inconsistent units.

### Operator overloads

Only `AbstractExpression` defines arithmetic and indexing overloads, listed in
the [table of free functions](reference.md#appendix-free-functions).

To allow external libraries to handle operations such as `expr * measure`,
binary arithmetic overloads return `NotImplemented` for non-expression operands.
Otherwise, overloads call the corresponding factory with operands unchanged,
except for indexing below. The core defines no reflected overloads.

An overload applies one free function to the operands, so the graph follows
the grouping and the precedence of Python, and construction never regroups
(see [No simplification during
construction](#no-simplification-during-construction)). For example,
`a + b + c` is `Sum(Sum(a, b), c)` and `a + (b + c)` is `Sum(a, Sum(b, c))`,
and the two are not equal.

Python passes `f[i]` to `__getitem__` as `i`, but `f[i, j]` as `(i, j)`.
`__getitem__` wraps a single entry in a `tuple` and calls `component` (see
[Partial indexing](expressions.md#partial-indexing)).

`AbstractExpression` defines `__iter__` to raise `TypeError`, because
expressions are not iterable and Python otherwise iterates through
`__getitem__`.

### Numbers in expressions

Every number in an expression is written with literals, such as `Integer(2)`,
`Float(0.5)`, `Integer(1) / Integer(3)` or `Complex(Integer(0), Integer(1))`. A Python `int`, `float`,
`complex`, `Fraction` or NumPy scalar is never converted and is rejected by
static type checking.

With `X = coordinate(RD(2))`, `X[0]**Integer(2)` is accepted, but `X[0]**2`
is a static type error. The same rule applies inside containers:
`tensor((Integer(1), Integer(2)))` is accepted, but `tensor((1, 2))` is not.

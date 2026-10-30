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

## Implementation principles

- Explicit is better than implicit.
- Boilerplate and error-prone patterns should be minimised.
- UFLx should be deterministic on the same machine using the same byte-code.
  Serialisation and pickling should be considered early.
- Performance matters, so the data structures and algorithms of the core
  should be chosen with care and tested for scalability.
- The minimal language should need only a small number of object kinds. When
  designing the language, consider whether objects can be composed rather than
  extended, and what could be layered on top of the core.

## Enforcement of the rules of the language

The rules of the base object in language.md must be enforced to help library
developers. Candidate tools are `dataclass`, deterministic hashing,
hash-consing, `WeakValueDictionary`, `ruff` and
[`semgrep`](https://github.com/semgrep/semgrep).

These rules should be put in place at the start of the project.

How best to enforce the rules is not yet clear. Custom static analysis, such
as `semgrep`, may be more robust than making Python 'well behaved' at runtime.
For example, the purity of a function cannot be checked in Python.

Python randomises the hash of a `str` for each process, so rule 1 of the base
object, a hash that is deterministic across processes, needs a custom hash for
`str` in an implementation in Python.

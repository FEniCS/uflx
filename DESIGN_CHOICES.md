# Design choices

Documentation of design decisions that have been made

## Naming

File names are plural to allow (eg) `uflx.domains` to be the module and `uflx.domain` to
be the function defined in domains.py that is imported into __init__.py.

Whenever a function to initialise a class is defined, it's name is the same as the class
name but lowercase, with `_x` wherever the class name has (eg) a capital `X` mid name.

## Default values

A class `__init__` does not choose a value on the caller's behalf. Where a default is
wanted, define a function or a named instance that initialises the class with it, as
`dx` does for a measure.

`None` is allowed as a default for a field that is not yet known, as against one being
chosen: a label that has not been minted, a variable an expression is not yet bound to,
a point a geometric quantity has not been told. There the field is absent and its type
says so, and requiring every caller to write `None` would say nothing extra.

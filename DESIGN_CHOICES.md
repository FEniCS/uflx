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
`dx(domain)` does for a measure: it supplies the density, which `Measure` does not
choose.

`None` is allowed as a default for a field that is not yet known, as against one being
chosen: a label that has not been minted, a variable an expression is not yet bound to,
a point a geometric quantity has not been told. There the field is absent and its type
says so, and requiring every caller to write `None` would say nothing extra.

## Derived domains

A domain derived from another is constructed by one function per construction and not
one per case of it. `sub_entities_of(domain, codim)` is the sub-entities of a domain at
any codimension, and there is no `facets_of`, `ridges_of` or `peaks_of`: the codimension
names which sub-entities are meant, so a function per codimension would be a second
spelling rather than a second construction. `facets_of` existed and was removed for that
reason.

The pull of such a name is that an exterior facet integral has one in the literature and
in every other library, `ds`. That belongs where the rest of the inherited vocabulary
is, in `uflx.legacy`, which spells it `ds` and builds `sub_entities_of(mesh, 1)`
underneath.

## Protocols of the domain layer

A protocol names a question a node can be asked, is a `runtime_checkable` `Protocol`
asked with `isinstance`, is never inherited, and lives in the lowest module that asks
it, which may be below the one that implements it: `domains.py` holds what `geometry.py`
and `integrals.py` ask of an interface, which `interfaces.py` implements above them. It
is named as the thing it says about the node, a verb phrase or an adjective. Those of
the domain layer are:

- `uflx/domains.py`: `ReadOnADomain`, an expression that says which domain it is read
  on; `ComposesADomain`, a domain presented as another with a map applied;
  `IncludesASubEntity`, a map that says which sub-entity of a cell it includes;
  `MadeOfPieces`, a domain integrated over as the sum over its pieces;
  `EmbedsAnInterface`, a domain that is one embedding of an interface; `CarriesADomain`,
  a map that starts from a whole domain; `StatesAnOutwardNormal`, a domain that says
  which way is out of what it bounds.
- `uflx/variables.py`: `ReadAtAVariable`, a node told where it is read after being
  built. `uflx/finite_elements.py`: `HasAGeometricDimension`.
- `uflx/algorithms/maps.py`: `PullBackToEntity`, `SplitByPiece`, `SplitByCellType`, what
  the pull back asks of an integral.

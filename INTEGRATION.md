# Measures and integrals

What integration means in UFLx, written against differential geometry rather than finite
element usage. A design note: it argues for the model, where `LANGUAGE.md` states it.
Arguments against the alternatives are in the notes at the end, and the vocabulary is
Lee's, cited there.

## 1. What an integral is

There is only one coordinate-free, metric-free integral in differential geometry:

> An `n`-form integrates over an oriented `n`-manifold.

For `ω` an `n`-form on an oriented `n`-manifold `M`, `∫_M ω` is defined, and for a
diffeomorphism `φ: U → M`,

```
∫_M ω = ∫_U φ*ω
```

*exactly*, with nothing left over. Change of variables is not a correction applied to
the integral; it is the statement that `∫` is defined on forms and `φ*` is how forms
move. In coordinates `φ*(dx¹ ∧ … ∧ dxⁿ) = det(J) dX¹ ∧ … ∧ dXⁿ`, so the Jacobian
determinant is not a factor someone multiplies in, it is what the pullback of the
coordinate volume form *is*.

Functions are not forms, so `∫_M f` means nothing until `M` is given a volume form or a
density. Given a metric `g`, in coordinates `X` on `U`:

```
∫_M f dV_g = ∫_U f(φ(X)) √det g(X) dX
```

For an immersion `φ: U ⊂ Rⁿ → R^d` with the Euclidean metric pulled back, `g = JᵀJ` and
the factor is `√det(JᵀJ)`, which for `n = d` is `|det J|`.

Three distinct objects are involved, and they are routinely conflated:

1. the **domain**, an oriented manifold, possibly with a boundary,
1. the **measure**, a volume form or a density on it,
1. the **pullback**, which moves both an integrand and a measure along a map.

A **density** is the absolute value of an `n`-form. It integrates over an unoriented
manifold and is always non-negative. UFLx integrates functions against densities rather
than forms against an orientation, which is the right choice for a mesh, whose cells are
not consistently oriented and whose volumes must be positive. The sign a density
discards is the orientation the parametrization gives, and §4 says what that sign is and
is not good for.

## 2. What a measure is

A measure is the domain integrated over paired with the density used on it:

```
measure = (the domain integrated over, the density on it)
```

Neither half is optional. A function is not a differential form, so integrating one
needs a density, and a density is a density of something.

The density is the one the domain's own parametrization induces, `√det g`, read off the
measure:

```python
dx(omega)  # Measure(omega), whose density is VolumeElement(omega)
```

It is a property of the domain rather than a second argument, and there is no way to
weigh a measure by a function [1]. Where the parametrization is the identity, as on
`R^d` or on a region of it, the density is one and the measure is the Lebesgue measure.
Seen another way the density is `dμ/dλ`, the derivative of the measure with respect to
the Lebesgue measure on its parametrization's parameter region, which is why the
parametrization fixes it.

The domain is stated and never inferred, so `dx` is a function of a domain rather than a
thing that acquires one [2]. `Integral` reads the domain off its measure and rejects a
function living on another.

The domain is any parametrized domain: one that offers a map out of a parameter region.
A reference cell is one kind of parameter region, a region of `R^d` is another, and
`R^d` itself is parametrized by itself. Being integrable over does not imply being made
of cells.

Change of variables acts on the measure as well as the integrand. Pulling an integral
back onto a cell's coordinates multiplies the measure's density into the integrand and
retargets the measure onto that cell, which is §1's identity rather than a factor
injected by an algorithm.

What a measure does not carry is topological selection. A codimension is relative to a
domain, and whether an entity lies on a boundary is a question about a cell complex UFLx
does not own; both are ways of *selecting a derived domain* from a mesh. That selection
is real and necessary, and it belongs to domain construction — see §3.

## 3. One measure, three domains

`dx`, `ds` and `dS` are one measure over three different domains [3].

- `dx` is the induced density on `Ω`.
- `ds` is the induced density on a codimension-one domain of facets, in its own right.
  Its parametrization per facet is the reference facet's inclusion into the cell
  followed by the cell's own map, which is a `ComposedParametrization`.
- `dS` is the induced density on an interface, a codimension-one domain again [4].

The `+`/`-` of UFL is not part of the measure either. It is a *restriction of a field*,
which is a pullback along one side's inclusion map, and it lives in the integrand.

So a codimension and a boundary flag become operators on domains —
`sub_entities_of(omega, codim)` and `interface(gamma, embeddings)` — returning domains
that then need only the one measure. That is the home for a mesh-relative query, and it
makes the derived domain an object carrying a parametrization, a metric, a normal and a
measure like any other.

The derived domain is named for the sub-entities it is parametrized by and not for a
boundary. Which facets it has — those with one adjacent cell, those with two, those a
consumer has marked — is a question about a cell complex UFLx does not own, so
membership is external, as a mesh's cells are. `sub_entities_of` says only how such a
domain is parametrized.

A marked selection of them is a domain and not an id on a measure, which is what UFL
encodes as a subdomain id. `LabelledDomain` is a domain and an opaque label a consumer
selects it by, changing no cell, no parametrization and no coordinate, so changing no
density and no kernel either — only identity, which is the point. It labels each *piece*
rather than the union of them, because splitting an integral retargets it onto one piece
at a time: a label on the union alone would be gone by the time there were kernels to
assign to it, where a label on the pieces rides into the kernel the way the facet's own
geometry does, through the quantities that name the piece. The same mechanism serves
marked cells, so `dx` of a subdomain needs no special case. A label says which entities
a consumer means and nothing about the fields on them, so an integral takes a label off
before asking whether a function belongs to the domain it is over.

The local index of a sub-entity is not a label of this kind. It is part of a
parametrization's identity, picking which inclusion and so which kernel, where a tag
picks which entities a kernel runs on. The inclusion an entity supplies says which
sub-entity it is (`IncludesASubEntity`), so the index is read off a pulled back integral
rather than recovered by matching the inclusion's matrix.

Several facets are parametrized from one reference facet, so the parameter region no
longer identifies the parametrization, and the domain is held as a union of parametrized
pieces. A union is not an atlas: its pieces carry no claim of a relation to each other.
Integrating over one is the sum over its pieces, which generalises the fan out over cell
types, and a piece names one parametrization and so one cell type too.

A region of `R^d` shows this with no cells involved: a box carries a measure, and so
does each of its faces, which is an exterior facet measure rather than a second kind of
measure. Its faces are written out by hand in `test/test_regions.py`, since
`sub_entities_of` reaches a cell's sub-entities through its entity and a region has no
cells; §6 says what asking a domain for its own boundary pieces would take.

## 4. Orientation

A density integrates on an unoriented manifold, which is what makes `|det J|` right for
a mesh. Stokes' theorem is about forms:

```
∫_M dω = ∫_∂M ι* ω
```

with `∂M` carrying the orientation induced from `M`. The divergence theorem is this, and
the outward normal in it is not a free sign — it is the induced boundary orientation
written as a vector.

Three separate things are at stake, and the difference decides how much `ds` costs.

- **A facet's outward conormal needs no orientation.** Outwardness is about the
  interior, not the orientation: `φ` carries interior to interior, so if `n̂·v > 0` for
  outward `v` then `(J⁻ᵀn̂)·(Jv) = n̂·v > 0`. Under a reflection, `J⁻ᵀn̂` gives the
  outward normal where the cross product of the mapped facet's own Jacobian does not. So
  this costs no new information, only a reference normal and a transport rule.
- **A manifold's surface normal needs a genuine external choice** of which side is up,
  and `det J` does not exist there to supply it.
- **`sign(det J)`** says whether a parametrization inverts. It feeds contravariant
  Piola, not normals.

What a normal needs is the domain it is a facet *of*, which is structural rather than
numerical. `test/test_regions.py` makes it exact: a box face's parametrization is an
inclusion whose offset holds one coordinate fixed, a Jacobian does not see an offset, so
the lower and upper faces of one axis have identical Jacobians at every point while
their outward normals are opposite. No function of the Jacobian can tell them apart. The
information belongs to the domain, which is where a region puts it and a cell complex
does not.

The design consequence is narrow: the sign is a property of a map into the domain and
not of a normal. On an interface, `ι₊` and `ι₋` induce opposite orientations on `Γ`, so
once a restriction names a side the normal's sign follows from that name with no
separate orientation terminal. An outward normal asks "outward from which side", a
restriction answers it, and integration by parts becomes statable in the language rather
than in a consumer's conventions.

## 5. The interface, without finite element terminology

Two manifolds `Ω₊` and `Ω₋` meeting along `Γ` of codimension one. Two embeddings

```
ι₊ : Γ → closure(Ω₊)        ι₋ : Γ → closure(Ω₋)
```

A field `u₊` on `Ω₊` restricts to `Γ` as the pullback `ι₊* u₊`, and the jump

```
[u] = ι₊* u₊ − ι₋* u₋
```

typechecks precisely because both terms are fields on the same `Γ`. There is one
integration domain, one measure, and two maps named in the integrand.

Two cases, and they are not the same difficulty.

**Matching: `Γ` is a facet of both sides.** Then `Γ` has two parametrizations out of the
same reference facet, `ψ± = φ± ∘ ι̂±`, differing by a reparametrization of the reference
facet, `ψ₊ = ψ₋ ∘ τ`, fixed by the two cells' vertex orderings. `τ` is an isometry of
the reference facet, so `|det Dτ| = 1`: the two sides induce *the same* density on `Γ`,
which is why one measure serves both. This is a two-map atlas with a genuine transition
map, and it is the one place in UFLx where an atlas is earned rather than decorative.

`Interface` is that atlas: a domain `Γ` and, per side, an *image* of it, which is `Γ`
parametrized by that side's own map. An image is a piece as `sub_entities_of` already
builds them, or a face of a box, so an interface of two cells and an interface of two
regions are one construction. The transition map is never an object. The point of the
`-` side's chart that `X` corresponds to is `ψ₋⁻¹(ψ₊(X))`, written literally as the
preimage through `ψ₋` of the image through `ψ₊`, both of which exist already; the
preimage is the terminal the surrounding library evaluates, as a facet permutation here
and by point location below. No derivative of `τ` is ever wanted, a restriction being
the pullback of a function and the measure being `Γ`'s. The two sides of one mesh share
one `φ`, so `ψ₋` is named by its side (`LabelledParametrization`), which is to a map
what a labelled domain is to a domain: identity and nothing geometric.

**Non-matching: `Γ` is its own domain.** Mortar coupling, a Lagrange multiplier on a
contact surface, two independently meshed subdomains. `Γ` is a facet of neither side and
`ι±` is not combinatorial, so restricting `u₊` to `Γ` means evaluating it at
`φ₊⁻¹(ψ_Γ(X))` — pulling a point of `Γ` into `Ω₊`'s coordinates through a map with no
symbolic inverse.

UFLx has that primitive. `PreimagePoint` is `φ⁻¹` on points, a terminal by design,
standing for a mapping the surrounding library provides, and it takes the
parametrization rather than the domain precisely because which cell the point lands in
is part of the question. So the general interface needs no new geometry: a point of `Γ`
is pushed forward by `Γ`'s own parametrization and pulled back through each side's.

Read in that order, UFL's `+`/`-` is the special case in which `φ₊⁻¹ ∘ ψ_Γ` is known
combinatorially and so costs nothing. In general the map of a side written in charts is
`ι̂_s = φ_s⁻¹ ∘ T_s ∘ ψ_Γ`, from `Γ`'s region into a cell's coordinates, and
`CoordinateRepresentation` is that map: its value is the preimage, the terminal again,
and its Jacobian is the chain rule in closed form. `T_s` is a map the user states, the
identity of sets unless said, which is what makes a periodic boundary an interface whose
other side is a translation away. An image built from it,
`composed_domain(composed_domain(region, CoordinateRepresentation(Γ, φ_s, T_s)), φ_s)`,
has the shape a facet of a mesh has and is read the same way, and being parametrized
from `Γ`'s own region it needs no chart transition. The matching case rested on three
assumptions, that `ι̂_s` lands in a facet, that `ι_s` is the inclusion of sets, and that
each piece of `Γ` meets one cell per side; this drops the first two. The third is where
a supermesh lives, as the pieces of `Γ`, a union whose membership is external as a
mesh's cells are, and the language says nothing about it. An interface immersed in one
mesh, both sides carried by the same map, has sides that are one, and the jump of a
field on that mesh vanishes node for node, which is right: without an enriched space
there is no jump across a line through a cell.

## 6. What is built, and what is not

Built, as `ds` and `dS`:

1. **The derived domain of facets.** `AbstractEntity.sub_entity_inclusion(dim, index)`
   is the reference geometry UFLx cannot supply, so the entity's definer supplies it and
   the default raises; `sub_entity_domain(cell, dim, index)` is a sub-entity of a cell
   as a domain, being that inclusion composed onto the sub-entity's own coordinate
   domain; `UnionDomain` holds several parametrized pieces; and
   `sub_entities_of(omega, codim)` is the union of the pieces, each a `ComposedDomain`
   of the inclusion and the cell's own map. A parametrization says which coordinate
   domain it lands in, which is how a facet's points are known to be named in the cell's
   coordinates. Any codimension, not a facet's: `sub_entities_of(omega, codim)` is the
   general form, and what fails to generalise above codimension one is a claim rather
   than a construction, a normal being a single direction only there.
1. **Restriction as a pullback along a named parametrization.** `restrict(u, domain)` is
   `ι*u`, the field read at the point the domain's parametrization gives. It lives in
   the integrand, and a gradient shows why: inside a restriction the entity gradient is
   still carried by the *cell's* inverse Jacobian, read at the included point, where the
   facet's own would be the tangential gradient.
1. **The fan out over pieces.** `Integral.split_by_piece` gives one integral per piece,
   carrying the measure, the restrictions and the geometric quantities onto it, and runs
   before the cell type split. Pulling back then multiplies the piece's volume element
   into the integrand, which is the factor an exterior facet integral picks up.
1. **An image per side.** `Interface(Γ, embeddings)` is §5's atlas, and the facets of a
   mesh presented twice, `interface(facets, (("+", facets), ("-", facets)))`, is `dS`.
   There is no second spelling of it. A restriction names a side,
   `restrict(u, Γ.embedding("-"))`, and reads the field at that side's chart point
   through that side's inclusion; `jump` and `avg` are the two restrictions' difference
   and mean. An interface of unions fans out over pairs of pieces, which is what
   identifies a kernel, and the pull back is onto the parameter region of the domain's
   one parametrization, so an interface needs no special case there. A mesh and a
   submesh of its entities meet the same way, with no spelling of their own: the
   entities presented by the mesh as `sub_entities_of(mesh, codim)` and by the submesh
   as `sub_entities_of(submesh, 0)`, at codimension zero as much as one, with the map
   between the two meshes' entities as the chart transition. A form mixing fields on
   both is each field restricted to its side.
1. **The outward normal.** `OutwardNormal(side)` is §4's transport rule: the reference
   facet's outward normal, which the entity supplies as `facet_outward_normal`, carried
   by the cell map's inverse transpose at the side's cell point. It needs no orientation
   of `Γ`, is outward under a reflection, and `ν₊ = -ν₋` is checked rather than defined.
   `UnitNormal` stays as the direction a domain of codimension one has on its own.

Not built:

1. **Pull back off a cell.** Change of variables onto a parametrization's parameter
   region is meaningful for any parametrized domain — the identity for `R^d` and for a
   region, the inclusion for a face — and `pull_back_to_entity` now asks for that
   region, but the variable it hands the integrand is still a cell's, so a region
   refuses. An interface of two boxes builds and has its normals, and stops there.
1. **A measure over a domain of topological dimension zero**, which is the counting
   measure a sub-entity domain of full codimension would give. The empty Gram
   determinant is one, but the empty matrix its Jacobian needs cannot be built, as
   `test_the_measure_of_a_point_is_the_counting_measure` records.
1. **A domain asked for its own boundary pieces.** `sub_entities_of` reaches a cell's
   sub-entities through its entity's combinatorics, so it serves a mesh and refuses a
   region, while §3 presents a box's faces as the same idea -- and they are, by hand. A
   domain that answered "the pieces of my boundary" however it liked, a cell type from
   its entity and a box from its bounds, would make that one construction rather than
   two. What this is *not* is a domain knowing it is a subset of another: the extent of
   a reference cell is the element library's, as its vertices are, and nothing here
   reads it. It is also not what the composition check wants, which is the declared
   intent of a map and not containment -- the reference triangle sits inside the
   reference square, so containment would accept a map written for the square being
   applied to the triangle.
1. **Atomic measures**, for a point evaluation and for a quadrature rule as a sum of
   Diracs. A sibling of `Measure` rather than a weighting of one.

## 7. What the design has to pass

Acceptance criteria, as tests rather than prose to agree with.

- **Change of variables holds by construction.** Pulling an integral back moves the
  measure onto the reference cell and multiplies its density into the integrand.
  `test_pulling_back_retargets_the_measure`.
- **A composed domain's measure is `√det g`.**
  `test_the_volume_element_is_the_metrics_gram_determinant`, and on a surface
  `test_the_volume_element_of_a_surface_is_its_area_element`.
- **Mixed cell types split, and the measure goes with them.**
  `test_an_integral_restricts_to_one_cell_type`.
- **A measure needs a parametrization and not a cell.**
  `test_a_measure_on_r_d_is_the_lebesgue_one` and
  `test_a_measure_on_a_regions_boundary_is_an_exterior_facet_measure`.
- **An axisymmetric measure is the volume element of a cylindrical parametrization.**
  Needs trigonometric functions in the expression language, of which `Sqrt` is currently
  the only neighbour.
- **A marked part of a domain is a domain, as far as a kernel.** `ds(3)` and `ds(4)`
  pull back to kernels that still say which selection they are for, and which
  sub-entity. `test_a_marked_tag_is_still_there_in_the_kernels` and
  `test_a_kernel_says_which_sub_entity_it_is_for`.
- **An exterior facet integral is the one measure over a domain of facets.** Pulling
  `inner(restrict(grad(u), facets), UnitNormal(facets)) * dx(facets)` back gives one
  integral per facet over the reference facet, carrying that facet's volume element and
  normal and the cell's inverse Jacobian at the included point.
  `test_pulling_back_an_integral_over_the_facets_gives_one_kernel_per_facet`.
- **A facet measure on a curved domain is arc length.** A triangle laid on the
  paraboloid `z = u^2 + v^2` has facets that are curves in space, and the density UFLx
  gives on one, integrated, is the length of the lifted curve rather than of the flat
  one. `test_the_measure_of_a_curved_facet_integrates_to_its_arc_length`.
- **The facet geometry is the cell's geometry, asked for.** The inclusion an entity
  supplies gives the same facet lengths as the one `test_boundaries.py` writes by hand.
  `test_the_measure_of_a_facet_is_its_length` in `test_facet_domains.py`.
- **One integral, two domains.** `inner(jump(u, Γ), v) * dx(Γ)` builds, the two
  restrictions are distinguishable in the graph, and pulling it back gives one kernel
  per pair of facets, the `-` side read at the chart transition.
  `test_the_jump_is_a_difference_of_restrictions_on_one_domain` and
  `test_pulling_an_interface_integral_back_gives_one_kernel_per_pair`.
- **The jump is side-symmetric in the measure.** Pulling the interface integral back
  through either side's parametrization gives the same density factor.
  `test_the_jump_is_side_symmetric_in_the_measure`.
- **The divergence theorem is statable.** `∫ div(u) dx(omega)` and
  `∫ inner(restrict(u, facets), OutwardNormal(facets)) dx(facets)` are both expressible,
  with the normal's sign fixed by the cell it points out of rather than by a consumer's
  convention, and the two sides of an interface have opposite ones.
  `test_the_outward_normals_of_the_two_sides_are_opposite`.
- **Two regions meet as two cells do.** Two boxes sharing a face make an `Interface`
  with the one measure and two outward normals, and no cell anywhere in it.
  `test_two_boxes_meet_along_a_face_with_no_cells_anywhere`.
- **A side need not be a facet.** A periodic boundary's other side is a translation
  away, checked to the number; two meshes couple through a multiplier on an interface
  mesh of its own; a line of intervals runs through tetrahedra.
  `test_a_periodic_boundary_is_an_interface_whose_other_side_is_a_translation_away`,
  `test_a_mortar_coupling_of_two_meshes_through_a_multiplier` and
  `test_a_well_in_three_dimensions_is_an_interface_of_codimension_two`.

## Notes

1. **Why a measure cannot be weighed by a function.** For any factor one might want to
   put in a measure: if it is geometric, it is the density of some parametrization, so
   it comes from composing, and the thing to change is the domain; if it is not
   geometric, it is physics, and belongs in the integrand. The axisymmetric `r dr dz` is
   the case worth working through, looking most like a weight. The cylindrical
   parametrization `psi(r, z, theta) = (r cos theta, r sin theta, z)` has
   `g = J^T J = diag(1, 1, r^2)`, so `sqrt(det g) = r`: the `r` is the volume element of
   a domain presented through that parametrization, ordinary induced geometry. It looks
   like a weight only once the domain is pretended to be two-dimensional. Stated
   honestly, the axisymmetric measure is the pushforward of the three-dimensional one
   along the projection that forgets `theta`, and the `2 pi` is `∫ dtheta`. A
   parametrization, then a map out of it. Completeness points the same way: every
   absolutely continuous measure is the pushforward of Lebesgue measure under some
   diffeomorphism — in one dimension the inverse of the distribution function, in higher
   dimensions Moser's argument — so a weight adds nothing composition cannot give. For
   `exp(-|x|^2) dx` that map is not elementary, which is a reason not to pretend it is a
   composition one can write down, not a reason to add a weight. What lies outside an
   induced density is a different kind of measure rather than a weighted one: a Dirac
   measure has no density with respect to anything, so it cannot be `w . mu`, and a
   quadrature rule is a sum of Diracs.
1. **Why the domain is stated rather than inferred.** Inference reads the domain off the
   functions in the integrand and requires them to agree, so it cannot express an
   integrand spanning two domains at all — which is exactly what an interface needs. It
   also cannot be contradicted: an integral that guesses its domain is an integral that
   cannot be told it guessed wrong, and where domain equality is structural the guess
   can be wrong silently, two distinct meshes built from the same coordinate element
   comparing equal.
1. **Why UFL has three measures rather than one.** Its measure does not name its domain,
   so what it integrates over has to be encoded in the measure instead, as a codimension
   and a boundary flag. Subdomain ids on the measure are the same conflation from the
   other side: a subdomain is a domain, not an indicator weight on a measure.
1. **Why `dS` as a kind of measure gets it backwards.** The matching case of §5 — where
   the two sides' parametrizations differ by a combinatorial reparametrization — is the
   one that costs nothing, and making `dS` a measure kind hardcodes it into the
   vocabulary. That leaves no room for the non-matching case, which is the one that
   needs the language's help and which UFL cannot express at all.

## References

The vocabulary is taken from John M. Lee's two textbooks, so that a word here means what
it means there.

- J. M. Lee, *Introduction to Smooth Manifolds*, 2nd ed., Graduate Texts in Mathematics
  218, Springer, 2012. Charts, atlases and transition maps (ch. 1), which an `Interface`
  is two of with a transition between them; the *coordinate representation* of a smooth
  map, `F̂ = ψ ∘ F ∘ φ⁻¹` (ch. 2), which `CoordinateRepresentation` is for the map of an
  embedding; immersions and embeddings (ch. 4) and embedded submanifolds (ch. 5), which
  an `Embedding` is taken to be on the user's word; restricting maps to submanifolds
  (ch. 5), which is what `Restriction` is named for; outward-pointing vectors along a
  boundary (ch. 5) and the boundary orientation they induce (ch. 15), which
  `OutwardNormal` writes as a unit vector; pullbacks of covector fields and of forms
  (ch. 11 and 14), which a value map is; Riemannian metrics and the metric a submanifold
  inherits (ch. 13) and the Riemannian volume form (ch. 15), which `MetricTensor` and
  `VolumeElement` are; and integration of forms, change of variables as the pullback of
  a form, Stokes's theorem with the boundary carrying the induced orientation,
  integration on Riemannian manifolds and densities (ch. 16), which §1 and §4 restate.
- J. M. Lee, *Introduction to Riemannian Manifolds*, 2nd ed., Graduate Texts in
  Mathematics 176, Springer, 2018. The unit normal and the second fundamental form of a
  hypersurface (ch. 8, Riemannian submanifolds), of which `UnitNormal` is the first and
  nothing yet is the second.

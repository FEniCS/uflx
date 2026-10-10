# Measures and integrals

What integration means in UFLx, written against differential geometry rather than finite
element usage. A design note: it argues for the model, where `LANGUAGE.md` states it.

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
discards is the orientation the chart gives, and §5 says what that sign is and is not
good for.

## 2. What a measure is

A measure is the domain integrated over paired with the density used on it:

```
measure = (the domain integrated over, the density on it)
```

Neither half is optional. A function is not a differential form, so integrating one
needs a density, and a density is a density of something.

The density is the one the domain's own chart induces, `√det g`, read off the measure:

```python
dx(omega)  # Measure(omega), whose density is VolumeElement(omega)
```

It is a property of the domain rather than a second argument, for the reason in §3.
Where the chart is the identity, as on `R^d` or on a region of it, the density is one
and the measure is the Lebesgue measure. Seen another way the density is `dμ/dλ`, the
derivative of the measure with respect to the Lebesgue measure on its chart's parameter
region, which is why the chart fixes it.

The domain is stated and never inferred, so `dx` is a function of a domain rather than a
thing that acquires one. An integral that guesses its domain from its integrand is an
integral that cannot be told it guessed wrong, and guessing cannot express an integrand
spanning two domains at all, which is what an interface needs. `Integral` therefore
reads the domain off its measure and rejects a function living on another.

The domain is any charted domain: one that offers a map out of a parameter region. A
reference cell is one kind of parameter region, a region of `R^d` is another, and `R^d`
itself is charted by itself. Being integrable over does not imply being made of cells.

Change of variables acts on the measure as well as the integrand. Pulling an integral
back onto a cell's coordinates multiplies the measure's density into the integrand and
retargets the measure onto that cell, which is §1's identity rather than a factor
injected by an algorithm.

What a measure does not carry is topological selection. A codimension is relative to a
domain, and whether an entity lies on a boundary is a question about a cell complex UFLx
does not own; both are ways of *selecting a derived domain* from a mesh. That selection
is real and necessary, and it belongs to domain construction — see §4.

## 3. No weighted measures

There is no way to weigh a measure by a function. For any factor one might want to put
in a measure:

- if it is geometric, it is the density of some chart, so it comes from composing — the
  thing to change is the domain and not the measure;
- if it is not geometric, it is physics, and belongs in the integrand.

The axisymmetric `r dr dz` is the case worth working through, looking most like a
weight. The cylindrical chart `psi(r, z, theta) = (r cos theta, r sin theta, z)` has

```
g = J^T J = diag(1, 1, r^2)        sqrt(det g) = r        |det J| = r
```

so `r` is the volume element of a domain presented through that chart: ordinary induced
geometry. It looks like a weight only once the domain is pretended to be
two-dimensional. Stated honestly, the axisymmetric measure is the pushforward of the
three-dimensional one along the projection that forgets `theta`, and the `2 pi` is
`∫ dtheta`. A chart, then a map out of it. Nothing weighs anything.

Completeness points the same way. Every absolutely continuous measure is the pushforward
of Lebesgue measure under some diffeomorphism — in one dimension the inverse of the
distribution function, in higher dimensions Moser's argument — so a weight adds nothing
composition cannot give. For `exp(-|x|^2) dx` that map is not elementary, which is a
reason not to pretend it is a composition one can write down, not a reason to add a
weight.

What lies outside an induced density is a different kind of measure rather than a
weighted one. A Dirac measure has no density with respect to anything, so it cannot be
`w . mu`, and a quadrature rule is a sum of Diracs. Those are atomic, and would be a
sibling of `Measure`. Restriction to a subdomain is likewise a domain and not an
indicator weight, which is why UFL's subdomain ids on the measure are the same
conflation from the other side.

## 4. `ds` and `dS` are not kinds of measure

In UFL, `dx`, `ds` and `dS` are three measures. Here they are one measure over three
different domains, and they look like three only when the domain is unstated and has to
be encoded in the measure instead.

- `dx` is the induced density on `Ω`.
- `ds` is the induced density on `∂Ω`, a codimension-one domain in its own right. Its
  chart per facet type is the reference facet's inclusion into the cell followed by the
  cell's own map, which is a `ComposedParametrization`.
- `dS` is the induced density on an interface, a codimension-one domain again.

The `+`/`-` of UFL is not part of the measure either. It is a *restriction of a field*,
which is a pullback along one side's inclusion map, and it lives in the integrand.

So a codimension and a boundary flag become operators on domains — `boundary_of(omega)`,
`interface_of(omega_plus, omega_minus)` — returning domains that then need only the one
measure. That is the home for a mesh-relative query, and it makes the derived domain an
object carrying a chart, a metric, a normal and a measure like any other.

A region of `R^d` shows this with no cells involved: a box carries a measure, and so
does each of its faces, which is an exterior facet measure rather than a second kind of
measure.

## 5. Orientation

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
- **`sign(det J)`** says whether a chart inverts. It feeds contravariant Piola, not
  normals.

What a normal needs is the domain it is a facet *of*, which is structural rather than
numerical. `test/test_regions.py` makes it exact: a box face's chart is an inclusion
whose offset holds one coordinate fixed, a Jacobian does not see an offset, so the lower
and upper faces of one axis have identical Jacobians at every point while their outward
normals are opposite. No function of the Jacobian can tell them apart. The information
belongs to the domain, which is where a region puts it and a cell complex does not.

The design consequence is narrow: the sign is a property of a map into the domain and
not of a normal. On an interface, `ι₊` and `ι₋` induce opposite orientations on `Γ`, so
once a restriction names a side the normal's sign follows from that name with no
separate orientation terminal. An outward normal asks "outward from which side", a
restriction answers it, and integration by parts becomes statable in the language rather
than in a consumer's conventions.

## 6. The interface, without finite element terminology

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

**Matching: `Γ` is a facet of both sides.** Then `Γ` has two charts out of the same
reference facet, `ψ± = φ± ∘ ι̂±`, differing by a reparametrization of the reference
facet, `ψ₊ = ψ₋ ∘ τ`, fixed by the two cells' vertex orderings. `τ` is an isometry of
the reference facet, so `|det Dτ| = 1`: the two sides induce *the same* density on `Γ`,
which is why one measure serves both. This is a two-chart atlas with a genuine
transition map, and it is the one place in UFLx where an atlas is earned rather than
decorative. It also wants a chart per `(cell type, side)`, where
`AbstractParametrizedDomain.parametrization` is keyed on cell type alone.

**Non-matching: `Γ` is its own domain.** Mortar coupling, a Lagrange multiplier on a
contact surface, two independently meshed subdomains. `Γ` is a facet of neither side and
`ι±` is not combinatorial, so restricting `u₊` to `Γ` means evaluating it at
`φ₊⁻¹(ψ_Γ(X))` — pulling a point of `Γ` into `Ω₊`'s coordinates through a map with no
symbolic inverse.

UFLx has that primitive. `PulledBackPoint` is `φ⁻¹` on points, a terminal by design,
standing for a Newton solve a consumer performs, and it takes the chart rather than the
domain precisely because which cell the point lands in is part of the question. So the
general interface needs no new geometry: a point of `Γ` is pushed forward by `Γ`'s own
chart and pulled back through each side's.

Read in that order, UFL's `+`/`-` is the special case in which `φ₊⁻¹ ∘ ψ_Γ` is known
combinatorially and so costs nothing. Designing `dS` as a measure kind gets this
backwards: it hardcodes the easy case into the vocabulary and leaves no room for the
other.

## 7. What is not built

1. **The derived domains**, `boundary_of` and `interface_of`. This needs reference
   geometry UFLx cannot supply: `AbstractEntity` knows a facet's vertices but not where
   they are, so the facet inclusion `s ↦ (1 − s, s)` comes from outside, as an abstract
   hook, a consumer-supplied chart, or an unexpandable terminal as in UFL. That fork is
   the real decision and it is not a measure question.
1. **A chart per side**, generalising `parametrization(cell)`, plus restriction as a
   pullback along a named chart. This is `dS` in both cases of §6, and it subsumes
   `+`/`-`.
1. **Orientation carried on a restriction**, giving a signed normal and a statable
   divergence theorem.
1. **Pull back off a cell.** Change of variables onto a chart's parameter region is
   meaningful for any charted domain — the identity for `R^d` and for a region, the
   inclusion for a face — but the variable it hands the integrand is a cell's, so
   `restricted_to` and `pull_back_to_entity` require a cellular domain and say so.
1. **Atomic measures**, for a point evaluation and for a quadrature rule as a sum of
   Diracs. A sibling of `Measure` rather than a weighting of one.

## 8. What the design has to pass

Acceptance criteria, as tests rather than prose to agree with.

- **Change of variables holds by construction.** Pulling an integral back moves the
  measure onto the reference cell and multiplies its density into the integrand.
  `test_pulling_back_retargets_the_measure`.
- **A composed domain's measure is `√det g`.**
  `test_the_volume_element_is_the_metrics_gram_determinant`, and on a surface
  `test_the_volume_element_of_a_surface_is_its_area_element`.
- **Mixed cell types split, and the measure goes with them.**
  `test_an_integral_restricts_to_one_cell_type`.
- **A measure needs a chart and not a cell.**
  `test_a_measure_on_r_d_is_the_lebesgue_one` and
  `test_a_measure_on_a_regions_boundary_is_an_exterior_facet_measure`.
- **An axisymmetric measure is the volume element of a cylindrical chart.** Needs
  trigonometric functions in the expression language, of which `Sqrt` is currently the
  only neighbour.
- **One integral, two domains.** `inner(ι₊* u₊ − ι₋* u₋, v) * dx(interface)` builds, and
  the two restrictions are distinguishable in the graph.
- **The jump is side-symmetric in the measure.** Pulling the interface integral back
  through either side's chart gives the same density factor.
- **The divergence theorem is statable.** `∫ div(u) dx(omega)` and
  `∫ inner(u, n) dx(boundary_of(omega))` are both expressible, with `n`'s sign fixed by
  the boundary's inclusion rather than by a consumer's convention.

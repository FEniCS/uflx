# Measures and integrals

What integration means, what UFLx currently does, and where the two disagree.
`LANGUAGE.md` leaves the measure section as a TODO and `OPEN_ISSUES.md` records that
the volume form is not an object; this note is the argument behind both. It is a design
note, not a specification: nothing here is implemented yet.

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

For an immersion `φ: U ⊂ Rⁿ → R^d` with the Euclidean metric pulled back,
`g = JᵀJ` and the factor is `√det(JᵀJ)`, which for `n = d` is `|det J|`.

Three distinct objects are involved, and they are routinely conflated:

1. the **domain**, an oriented manifold, possibly with a boundary,
1. the **measure**, a volume form or a density on it,
1. the **pullback**, which moves both an integrand and a measure along a map.

A **density** is the absolute value of an `n`-form. It integrates over an unoriented
manifold and is always non-negative. Taking `|det J|` rather than `det J` is the choice
of densities over forms, and it is the right choice for finite elements — a mesh's cells
are not consistently oriented and a cell volume must be positive — but it is a choice,
and discarding the sign discards exactly the information an outward normal needs. The
unexpressible orientation in `OPEN_ISSUES.md` and the unsigned `UnitNormal` are not two
gaps. They are one.

## 2. What UFLx does now

```python
dx = Measure(dim=None, codim=0, boundary_only=False)
```

A measure holds three integers and flags and nothing else. Nothing reads them: they go
into `init_args` and `__repr__`, and `codegeneration` consults `_codim` and
`_boundary_only` only to refuse what it does not handle. The geometry lives in
`Integral.pull_back_to_entity`, which multiplies it into the integrand:

```python
det = abs(JacobianDeterminant(domain))
return Integral(det * integrand, self._measure, ...)
```

Four consequences, each checked against the code on this branch.

**The measure does not say what it integrates over.** `Integral.domain` reads the domain
off the integrand, through `extract_domain`, and asserts that every function agrees. So
the domain is inferred from the data rather than stated, and an integrand mentioning
functions on two domains is rejected:

```python
V1, V2 = function_space(m1, P1), function_space(m2, P1)
inner(Coefficient(V1), TestFunction(V2)) * dx  # AssertionError
```

That is the motivating example of this branch — two subdomains meeting on an interface —
and `dx` cannot hold it. The assert also passes silently when it should not, since
domain equality is structural: two distinct meshes built from the same coordinate
element compare equal (`OPEN_ISSUES.md`, multi-mesh forms).

**Its fields are a mesh query, not geometry.** `codim` is a codimension relative to a
domain the measure does not name. `boundary_only` asks a topological question about a
cell complex UFLx does not own. Neither is a property of a measure; together they are a
way of *selecting a derived domain* from a mesh. That selection is real and necessary,
but it belongs to domain construction.

**Pulling back does not change the measure.** After `pull_back_to_entity` the integral
is over a reference cell, with the volume factor in the integrand, and still carries the
same `Measure` object:

```
measure before: Measure(codim=0)   after: Measure(codim=0)   same object: True
```

A measure that says nothing cannot be made wrong by a transformation, which is the only
reason this is currently harmless. Section 1 says a pullback acts on the measure, and
here it provably does not.

**The volume factor is named for a thing it is not.** `JacobianDeterminant` already
returns `abs(det J)` for a square map and `Sqrt(det(JᵀJ))` for a tall one — it is
`√det g`, a density factor, never a determinant. `pull_back_to_entity` then takes `abs`
of it again, and `Abs(Abs(·))` survives to the expanded, simplified graph:

```
after simplify, Abs count: 2    # Abs(Abs(Sum(...)))
```

Numerically harmless, conceptually not: three different objects — the signed `det J` of
a square map, the density factor `√det g`, and the orientation sign — are wearing one
name, and the one that is missing is the one nothing can currently express.

## 3. Is `dx` correct

No, in three separable ways.

- It does not name its domain, so it relies on inference, cannot express an integrand
  spanning two domains, and cannot be retargeted by a pullback.
- It carries topological selection (`dim`, `codim`, `boundary_only`) that is not a
  property of a measure and that only means anything relative to an external mesh.
- It carries no density, so the one thing a measure is for is supplied by an algorithm
  instead, hardcoded to the induced Riemannian one.

Removing all three leaves a measure that is the pair a measure actually is:

```
measure = (the domain integrated over, the density on it)
```

with the density defaulting to the one the domain's own parametrization induces,
`√det g`. The suggested spelling keeps `f * dx` working for the single-domain case, using
the `None`-means-not-yet-known allowance in `DESIGN_CHOICES.md` rather than a default
chosen in `__init__`:

```python
dx  # Measure(domain=None), the domain to be inferred as today
dx(omega)  # Measure(omega), the domain stated
```

Inference then stays available as sugar for the common case and is documented as sugar,
rather than being the model.

A measure that carries its density also buys something UFL cannot express at all. An
axisymmetric problem integrates against `r dr dz`, and in UFL the `r` is multiplied into
every integrand by hand, where it is indistinguishable from part of the physics. As a
weighted measure it is where it belongs:

```python
dx_axi = Measure(omega, density=x[0] * InducedDensity(omega))
```

The same slot holds a surface measure that is not the induced one, and — stretching
further — a Dirac measure on a point domain, which is what `LANGUAGE.md`'s
`PointEvaluation` is.

## 4. `ds` and `dS` are not kinds of measure

In UFL, `dx`, `ds` and `dS` are three measures. Under §3 they are one measure over three
different domains, and the only reason they look like three is that the domain is not
stated and so has to be encoded in the measure instead.

- `dx` is the induced density on `Ω`.
- `ds` is the induced density on `∂Ω`, a codimension-one domain in its own right. Its
  parametrization per facet type is the reference facet's inclusion into the cell
  followed by the cell's own map — a `ComposedParametrization`, which already exists.
- `dS` is the induced density on an interface, a codimension-one domain again.

The `+`/`-` of UFL is then not part of the measure either. It is a *restriction of a
field*, which is a pullback along one side's inclusion map, and it lives in the
integrand.

So `boundary_only` and `codim` become operators on domains — `boundary_of(omega)`,
`interface_of(omega_plus, omega_minus)` — returning domains that then need only the one
measure. That is the honest home for a mesh-relative query, and it makes the derived
domain an object that can carry a parametrization, a metric, a normal and a measure like
any other.

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
integration domain, one measure, and two maps named in the integrand. The current
language cannot say this, because the only thing it can say about where a function lives
is the single domain its function space carries, and `Integral` demands that they all
agree.

Two cases, and they are not the same difficulty.

**Matching: `Γ` is a facet of both sides.** Then `Γ` has two parametrizations out of the
same reference facet, `ψ± = φ± ∘ ι̂±`, and they differ by a reparametrization of the
reference facet, `ψ₊ = ψ₋ ∘ τ`, fixed by the two cells' vertex orderings. `τ` is an
isometry of the reference facet, so `|det Dτ| = 1`: the two sides induce *the same*
density on `Γ`, which is why one measure legitimately serves both. This is a two-chart
atlas with a genuine transition map, and it is the first place in UFLx where an atlas is
earned rather than decorative. It also breaks an assumption:
`AbstractParametrizedDomain.parametrization(cell)` is a function of cell type alone, and
an interface needs a map per `(cell type, side)`. The signature has to generalise before
`dS` can be anything.

**Non-matching: `Γ` is its own domain.** Mortar coupling, a Lagrange multiplier on a
contact surface, two independently meshed subdomains. Now `Γ` is not a facet of either
side, and `ι±` is not combinatorial. Restricting `u₊` to `Γ` means evaluating it at
`φ₊⁻¹(ψ_Γ(X))` — pulling a point of `Γ` into `Ω₊`'s coordinates through a map that has
no symbolic inverse.

UFLx already has this primitive. `PulledBackPoint` is exactly `φ⁻¹` on points, a
terminal by design, standing for a Newton solve a consumer performs, and it takes the
parametrization rather than the domain precisely because which cell the point lands in
is part of the question. So the general interface needs no new geometry: a point of `Γ`
is pushed forward by `Γ`'s own map and pulled back through each side's.

Read in that order, UFL's `+`/`-` is the special case in which `φ₊⁻¹ ∘ ψ_Γ` is known
combinatorially and so costs nothing, and the general case is the one that was never
expressible. Designing `dS` as a measure kind gets this backwards: it hardcodes the easy
case into the vocabulary and leaves no room for the other.

## 6. Orientation, and why it cannot be deferred much longer

A density integrates on an unoriented manifold, which is what makes `|det J|` right for
a mesh. But Stokes' theorem is about forms:

```
∫_M dω = ∫_∂M ι* ω
```

with `∂M` carrying the orientation induced from `M`. The divergence theorem is this, and
the outward normal in it is *not* a free sign — it is the induced boundary orientation
written as a vector. This is why `UnitNormal` cannot fix its own sign from the tangent
space alone, and why `CellOrientation` exists in UFL as a thing a consumer supplies.

The design consequence is narrow and worth stating: the sign is not a property of a
normal, it is a property of a map into the domain. On an interface, `ι₊` and `ι₋` induce
opposite orientations on `Γ`, so once a restriction names a side the normal's sign is
determined by that name and needs no separate orientation terminal. An outward normal
asks "outward from which side", the restriction already answers it, and integration by
parts becomes statable in the language rather than in a consumer's conventions.

## 7. A staged proposal

Each stage stands alone and is listed with what it costs.

1. **Name the density factor what it is.** Split `JacobianDeterminant` into the signed
   determinant of a square map and a density factor — `VolumeElement`, say, being
   `√det g` — and drop the duplicated `abs`. Pure clarification, small, and it is what
   makes an orientation sign expressible later rather than being absorbed into a name.
   Breaks `codegeneration`'s pattern match on `JacobianDeterminant`.
1. **Give the measure its domain.** `Measure(domain=None)`, `dx(omega)`, and
   `Integral` reconciles the stated domain with the inferred one and errors when they
   disagree. `pull_back_to_entity` retargets the measure onto the reference cell instead
   of leaving it untouched. No behaviour change for single-domain forms.
1. **Give the measure its density.** `Measure(domain, density)`, defaulting through a
   factory to the induced one; `pull_back_to_entity` reads the factor off the measure
   rather than constructing it. The density becomes a successor of the measure in the
   graph, so `replace` and `expand_geometry` reach it with no new traversal.
   Axisymmetric and other weighted measures follow for free.
1. **Retire `dim`, `codim` and `boundary_only`** in favour of domain constructors,
   `boundary_of` and `interface_of`. This needs the reference geometry UFLx cannot
   supply (`OPEN_ISSUES.md`): the facet inclusion `s ↦ (1 − s, s)` has to come from
   outside, as an abstract hook, a consumer-supplied parametrization, or an unexpandable
   terminal as in UFL. That fork is the real decision and it is not a measure question.
1. **Generalise `parametrization(cell)` to admit a map per side**, and add restriction as
   a pullback along a named map. This is `dS`, in both the matching and non-matching
   cases, and it subsumes `+`/`-`.
1. **Carry orientation on the restriction**, giving a signed normal and a statable
   divergence theorem.

## 8. What the design has to pass

Acceptance criteria, as tests to write rather than prose to agree with.

- **Change of variables holds by construction.** Pulling an integral back changes the
  measure's domain and its density together, and the result's domain is the reference
  cell. Today the first half silently does not happen.
- **A weighted measure is distinguishable from physics.** `∫ f r dx` as an axisymmetric
  measure and `∫ (f r) dx` as a weighted integrand are different objects, and the
  measure survives differentiation of the form with respect to a coefficient untouched.
- **A composed domain's measure is `√det g`.** Already covered by
  `test_the_volume_scaling_is_the_metrics_gram_determinant`, which should stop needing
  to apologise for `abs(det J)` and `√det(JᵀJ)` being equal as numbers but not as
  expressions.
- **One integral, two domains.** `inner(ι₊* u₊ − ι₋* u₋, v) * dx(interface)` builds,
  and the two restrictions are distinguishable in the graph.
- **The jump is side-symmetric in the measure.** Pulling the interface integral back
  through either side's chart gives the same density factor.
- **Mixed cell types still split.** `split_by_cell_type` keeps working when the measure
  carries a domain, giving one measure per cell type rather than one shared measure.
- **The divergence theorem is statable.** `∫ div(u) dx(omega)` and
  `∫ inner(u, n) dx(boundary_of(omega))` are both expressible, with `n`'s sign fixed by
  the boundary's inclusion rather than by a consumer's convention.

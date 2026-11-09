# Patch problems

## A mean-constrained DG patch correction

On a nonempty, facet-connected vertex patch $\omega_p$, solve a real DG Poisson
correction with homogeneous Neumann conditions and zero mean:

$$
\begin{aligned}
a_p^N(\delta u_p,v)+\lambda_p\int_{\omega_p}v\,dx &= r_p(v),\\
\int_{\omega_p}\delta u_p\,dx &= 0.
\end{aligned}
$$

Here $r_p$ is the restricted global residual. The scalar multiplier removes
the constant nullspace. The auxiliary operator $a_p^N$ contains cell and internal
SIPG facet terms, with a sufficiently large penalty; it omits patch-boundary
terms and therefore differs from a restriction of the global DG matrix.

### Inputs

The library supplies these objects and their transformation laws. Code uses
core names without imports, as in the other examples.

| Input | Contract |
| ----- | -------- |
| `T` | A fixed local patch configuration, with one reference-tetrahedron slot per cell and `sides = ()`. |
| `S` | A fixed configuration of internal facets, one reference-triangle slot per facet, with sides naming its two distinct parent slots in `T`. |
| `X` | Real continuous vector P1 coordinate space on `T`, shape `(3,)`. |
| `V` | Real scalar DG P1 space on `T`. |
| `Q` | Real scalar space on `S`, constant independently on each facet. |
| `C` | Real scalar constant space: `domain = None`, `value_shape = ()`, `argument_spaces = ()`. |

The solver supplies geometry, dofmaps, penalties including size scaling, and
[side data](../language/topology.md#side-data) pairing the two traces on each
facet. This construction describes one fixed patch configuration. The package
may specialise it for different configurations, or assemble a runtime-sized
patch from single-cell and two-parent facet kernels. Runtime cell counts and
neighbour lookup are not core inputs. `C` supplies one multiplier for the
whole patch; assembly must associate all local constraint contributions with
that same multiplier.

### Form

```python
x = function(X, "patch_coordinates")
D = mesh(x, RD(3))
E = side_domain(D, S, "+")

# Trial unknowns and test functions.
du = function(V, "correction")
lam = function(C, "multiplier")
v = function(V, "test")
mu = function(C, "constraint_test")

du_D = push_forward(du, D, transformation_law(V))
v_D = push_forward(v, D, transformation_law(V))
gdu = gradient(du_D, D)
gv = gradient(v_D, D)

sigma_ref = function(Q, "patch_penalty")
sigma = push_forward(sigma_ref, E, transformation_law(Q))
# Cell and interior-facet terms; gradients precede restriction.
i = index("component")
half = Integer(1) / Integer(2)


def dot(a, b):
    return index_sum(a[i] * b[i], (i,))


cell_term = integral(dot(gdu, gv), D)

plus = side_map(D, E, "+")
minus = side_map(D, E, "-")
n = normal(plus)

du_p, du_m = pull_back(du_D, plus), pull_back(du_D, minus)
v_p, v_m = pull_back(v_D, plus), pull_back(v_D, minus)
gdu_p, gdu_m = pull_back(gdu, plus), pull_back(gdu, minus)
gv_p, gv_m = pull_back(gv, plus), pull_back(gv, minus)

jump_du = du_p - du_m
jump_v = v_p - v_m
avg_dn_du = half * dot(gdu_p + gdu_m, n)
avg_dn_v = half * dot(gv_p + gv_m, n)

facet_term = integral(
    -avg_dn_du * jump_v - avg_dn_v * jump_du + sigma * jump_du * jump_v,
    E,
)
# One patch-wide mean constraint.
constraint_term = integral(lam * v_D + mu * du_D, D)

a_patch = form(
    cell_term + facet_term + constraint_term,
    ((v, mu), (du, lam)),
)

r = function(dual_space(V), "patch_residual")
r_named = name_arguments(r, (v,))
L_patch = form(action(r_named, v, v), ((v, mu),))
```

The mixed arguments are the test pair `(v, mu)` and trial pair `(du, lam)`.
The load has zero constraint component. Assembly produces

$$
\begin{bmatrix}A_p^N&m_p\\m_p^T&0\end{bmatrix}
\begin{bmatrix}\delta u_p\\\lambda_p\end{bmatrix}
=\begin{bmatrix}r_p\\0\end{bmatrix},
\qquad (m_p)_i=\int_{\omega_p}\phi_i\,dx.
$$

The solver supplies each patch's data, assembles and solves this system, then
scatters only the correction. An additive sweep uses the same global residual
for all patches. The zero-mean corrections need a complementary component,
such as a coarse correction, to address the remaining modes.

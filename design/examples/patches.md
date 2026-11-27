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

The patch size may vary between calls. DOLFINx supplies packed local data and
local-to-patch degree-of-freedom maps. UFLx expresses the patch operator using
ordinary arrays, gathers and sums; no runtime topology or mesh-cell selector
is needed.

### Packed inputs

Use DG P1 tetrahedra, with four local degrees of freedom per cell. An interior
facet contribution has eight: four from each parent. Cell and facet kernels
use fixed quadrature rules. Their local contributions are written as array
inputs here to keep the patch expression focused on indexing.

The library supplies the following real arrays and index tables. `Nc`, `Nf`
and `N` are `AbstractExtent`s for the numbers of patch cells, internal facets
and patch degrees of freedom. DOLFINx supplies their values for each call.
For this DG space, `N` has value four times `Nc`; symbolic extent arithmetic
is not needed.

| Input | Shape | Meaning |
| ----- | ----- | ------- |
| `A` | `(Nc, 4, 4)` | Cell stiffness matrices, with test rows and trial columns. |
| `F` | `(Nf, 8, 8)` | Interior SIPG facet matrices, including consistency, symmetry and penalty terms. |
| `m` | `(Nc, 4)` | Cell integrals of the local basis functions. |
| `r` | `(N,)` | Restricted residual in patch degree-of-freedom order. |
| `cell_dofs` | `(Nc, 4)` | `AbstractIndexTable` with extent `N`, mapping cell-local to patch degrees of freedom. |
| `facet_dofs` | `(Nf, 8)` | `AbstractIndexTable` with extent `N`, mapping facet-local to patch degrees of freedom. |

Each internal facet occurs once. Its first four entries refer to the plus
parent and the last four to the minus parent, in the same order as `F`.
These entries refer to the existing patch degrees of freedom, not additional
facet unknowns. Both maps use the same ordering as the corresponding local
matrices. For example, a facet between the first two cells might map to
`(0, 1, 2, 3, 4, 5, 6, 7)`; other patch orderings work unchanged.

The [facet example](coupling.md#an-interior-facet-integral) gives the local
SIPG form used to compute `F`. The cell form computes `A`, and integrating
each cell basis function computes `m`. The package supplies geometry,
orientations and penalties to those kernels and packs their results. The
expression below combines the contributions without requiring a preassembled
patch matrix.

### Patch form

The library supplies real spaces of arrays with no spatial domain:
`V` has value shape `(N,)`, and `C` has value shape `()`. Both have argument
spaces `()`. Functions in `V` represent patch coefficient vectors; functions
in `C` represent the single patch-wide multiplier and its test value.

```python
du, v = function(V, "correction"), function(V, "test")
lam, mu = function(C, "multiplier"), function(C, "constraint_test")

c, f = index("cell"), index("facet")
a, b = index("cell_row"), index("cell_column")
i, j = index("facet_row"), index("facet_column")
p = index("patch_dof")

ca = gathered_index(cell_dofs, (c, a))
cb = gathered_index(cell_dofs, (c, b))
fi = gathered_index(facet_dofs, (f, i))
fj = gathered_index(facet_dofs, (f, j))

cell_term = index_sum(v[ca] * A[c, a, b] * du[cb], (c, a, b))
facet_term = index_sum(v[fi] * F[f, i, j] * du[fj], (f, i, j))
constraint_term = index_sum(m[c, a] * (lam * v[ca] + mu * du[ca]), (c, a))

a_patch = form(cell_term + facet_term + constraint_term, ((v, mu), (du, lam)))
L_patch = form(index_sum(r[p] * v[p], (p,)), ((v, mu),))
```

For example, `du[ca]` reads the patch coefficient for local degree of freedom
`a` of packed cell `c`. The indices `c` and `f` index input arrays; they do
not select cells from an `AbstractTopology`. Their extents come from the
array and table shapes. All indices are summed, so the forms have no free
indices. The same expressions serve different patch sizes, including a
single-cell patch with `Nf = 0`, whose facet sum is zero.

The mixed arguments are the test pair `(v, mu)` and trial pair `(du, lam)`.
The load has zero constraint component. The forms describe

$$
\begin{bmatrix}A_p^N&m_p\\m_p^T&0\end{bmatrix}
\begin{bmatrix}\delta u_p\\\lambda_p\end{bmatrix}
=\begin{bmatrix}r_p\\0\end{bmatrix},
\qquad (m_p)_i=\int_{\omega_p}\phi_i\,dx.
$$

DOLFINx supplies each patch's packed inputs, evaluates the local system,
solves it and scatters only the correction. An additive sweep uses the same
global residual for all patches. The zero-mean corrections need a complementary
component, such as a coarse correction, to address the remaining modes.

For mixed cell types, use one set of arrays and sums per cell type and per
interface configuration, with their fixed local matrix sizes. All index tables
map into the same patch vector. The runtime number of contributions in each
set may vary without changing the expression.

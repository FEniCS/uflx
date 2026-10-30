# Tensor operations

These examples use the core names from [language.md](../language.md) without
imports. Shared helper definitions are in the
[reference library](library.md).

## Transpose

A transpose is a fixed composition of index operations. For a matrix expression
`A`, the core's `transpose(A)` builds:

```python
k = Index(("uflx.transpose", A, 0))
l = Index(("uflx.transpose", A, 1))
component_tensor(A[k, l], (l, k))
```

The [derived indices](../language.md#derived-indices) keep the labels tied to `A`.
Using user indices instead gives a mathematically equal expression with a
different graph. `Indexed` rejects an operand of array rank other than 2,
including a scalar, with `IndexError`.

## Tensor algebra with `EllipsisIndex`

An ellipsis index lets a library define a tensor operation once for every
array rank. With derived ellipsis indices `r` and `s` and an ordinary derived
index `k`, the definitions are:

| Operation | Fixed composition |
| --------- | ----------------- |
| `outer(u, v)` | `component_tensor(u[r] * v[s], (r, s))` |
| `dot(u, v)` | `component_tensor(index_sum(u[r, k] * v[k, s], (k,)), (r, s))` |
| `inner(A, B)` | `index_sum(A[r] * conjugate(B[r]), (r,))` |

These are library operations. Core tensor arithmetic instead builds one index
per axis. For vectors and matrices, `dot(u, v)` and `u @ v` contract the same
axes but are not structurally equal. Nor does an ellipsis composition reduce
when it meets a binder with separate indices, as in `inner(transpose(A), B)`;
that requires a later rewrite.

## Symmetric tensors

Tensor symmetry can live in the data representation. A 2D stiffness tensor
with minor and major symmetries has 16 entries but only six stored values.

Let `T` be a triangle topology, `x` its coordinate function of value shape
`(2,)`, `V` a real vector space of the same value shape on `T`, and `K` a
real space of constants of value shape `(6,)`. The library supplies an
`AbstractIndexValue` table `packed`, of shape `(2,2,2,2)` and extent 6,
selecting the stored value for each tensor entry:

```python
D = image(x, RD(2))
u, v = function(V, "u"), function(V, "v")
i, j, k, l = index("i"), index("j"), index("k"), index("l")
two = Integer(2)
c = function(K, "c")
stiffness = component_tensor(c[gathered_index(packed, (i, j, k, l))], (i, j, k, l))

grad_u, grad_v = gradient(push_forward(u, D), D), gradient(push_forward(v, D), D)
eps_u = (grad_u + transpose(grad_u)) / two
eps_v = (grad_v + transpose(grad_v)) / two
elastic = form(
    integral(index_sum(stiffness[i, j, k, l] * eps_u[k, l] * eps_v[i, j], (i, j, k, l)), D),
    (v, u),
)
```

The kernel receives six values for `c`; gathering reconstructs the full
stiffness tensor. The core needs no symmetry attribute.

```{dropdown} Symmetric directions and Voigt notation
For an energy $C_{01}^2$, unconstrained differentiation gives $2C_{01}$ at
entry `(0,1)` and zero at `(1,0)`. Differentiating in symmetric unit directions
gives $C_{01}$ at both. The results agree in every symmetric direction, but
only the second is symmetric.

Voigt and Mandel representations scale off-diagonal strain entries by 2 and
$\sqrt{2}$ respectively. A library can read either through a table and scale
factors. Differentiation with respect to stored Voigt values differs from
differentiation with respect to the full tensor by a factor of 2 off the diagonal.
```

# Compiling numerical kernels

## Scope and inputs

This page gives non-normative compiler guidance.

[Core normal form](evaluation.md#the-core-normal-form) leaves reference-domain
integrals, function evaluations and permitted geometric operations—such as
terminal gradients and pull-backs—for the consumer.
A compiler can use core normal form or another suitable representation.

Expression optimisation may be provided by UFLx as a reusable facility, with
compilers choosing when to invoke it, which transformations to explore and how
to assess execution cost. It may run before or after discretisation. Such
optimisation is separate from establishing core normal form and follows the
[optimisation contract](optimisation.md).

The core does not specify an ABI or perform global assembly.

## Library, compiler and runtime

| Component | Responsibility |
| --- | --- |
| Cell and tabulation library, such as [Basix](https://docs.fenicsproject.org/basix/main/python/index.html) | Define reference cells and finite elements; provide basis values and derivatives at requested points, quadrature rules and DOF transformations |
| Compiler | Request the reference data needed by the form; build and optimise expressions; generate kernels and choose data layouts |
| Solver runtime, such as DOLFINx | Supply mesh coordinates, coefficient values, constants and entity information; call kernels and assemble their contributions |

Table generation, code generation and data loading may be separate stages.
Tables may be embedded constants, separately loaded data or runtime inputs.

## Compiler choices

| Concern | Possible treatment |
| --- | --- |
| Integration | Replace an integral by a weighted sum of point evaluations, or use another integration method |
| Function evaluation | Request basis values and derivatives from the tabulation library; embed tables, arrange runtime tabulation or generate evaluation code |
| Inputs and outputs | Represent coefficients, geometry and tables as arrays; choose output axes for a local tensor or an operator action |
| Symbolic work | Simplify indexing and identity contractions, share repeated expressions, or rearrange sums and products |
| Execution | Choose loops, temporary storage, vectorisation and parallel work |
| Numerical behaviour | Choose precision, reduction order and permitted floating-point transformations |

Tabulation may preserve tensor-product factors rather than produce a dense
table. Exploiting those factors is a compiler choice.

## An MLIR example

This example shows how a symbolic expression for a local matrix becomes code
that computes the matrix entries.

Consider the unit-coefficient mass form for a real scalar Lagrange element.
Let $\psi$ map reference cell `K` to a full-dimensional physical cell. The
identity [transformation law](geometry.md#transformationlaw) gives
$\varphi_i\circ\psi=\phi_i$ for the physical basis. After choosing the same
basis for both form arguments, the local matrix is

```python
A = component_tensor(integral(phi[i] * phi[j] * rho, K), (i, j))
```

| Symbol | Meaning | Shape | Domain |
| --- | --- | --- | --- |
| `K` | Reference cell | — | — |
| `phi` | Library terminal containing the real reference basis functions, with no form argument slots | `(ndof,)` | `K` |
| `rho` | Integration density, with no form argument slots or free indices | `()` | `K` or `None` if spatially constant |
| `i`, `j` | Distinct free indices, each of extent `ndof` | — | — |

With $J(\xi)$ the Jacobian of $\psi$, `rho` is the core-normal-form expansion
of [`density(ψ)`](geometry.md#jacobian-and-density):

$$
\rho(\xi) = |\det J(\xi)|
$$

| Expression | Shape | Free indices | Domain |
| --- | --- | --- | --- |
| `phi[i] * phi[j] * rho` | `()` | `i, j` | `K` |
| `integral(phi[i] * phi[j] * rho, K)` | `()` | `i, j` | `None` |
| `A` | `(ndof, ndof)` | None | `None` |

After quadrature and tabulation,

$$
A_{ij} \approx \sum_q w_q\,\rho_q\,B_{qi}B_{qj},
$$

or, in the same index notation,

```python
A_q = component_tensor(
    index_sum(w[q] * rho_q[q] * B[q, i] * B[q, j], (q,)),
    (i, j),
)
```

At reference quadrature points $\xi_q$, the compiler uses:

| Array | Shape | Values | Source |
| --- | --- | --- | --- |
| `B` | `(nq, ndof)` | Basis values $\phi_i(\xi_q)$ | Tabulation library, evaluated at the chosen points |
| `w` | `(nq,)` | Quadrature weights | Quadrature library, for the requested rule |
| `rho_q` | `(nq,)` | Densities $\rho(\xi_q)$ | Generated geometry computation using runtime coordinate DOFs |

`A_q` also satisfies core normal form.

### Geometry preparation

Take a linear scalar element on an affine triangle.
The reference vertices are $(0,0)$, $(1,0)$ and $(0,1)$, with basis
$\phi=(1-\xi-\eta,\xi,\eta)$. Use quadrature points $(1/6,1/6)$,
$(2/3,1/6)$ and $(1/6,2/3)$, each with weight $1/6$.

This example represents the basis and quadrature tables as constant globals.
The solver runtime supplies physical vertices $X_0,X_1,X_2$ through
`coordinate_dofs`. Geometry preparation computes

$$
J = [X_1-X_0,\;X_2-X_0], \qquad \rho = |\det J|.
$$

The cell density can be reused across contractions. A backend may batch
geometry preparation across cells or fuse preparation with a contraction.

```mlir
func.func @cell_density(%coordinate_dofs: memref<3x2xf64>) -> f64 {
  %c0 = arith.constant 0 : index
  %c1 = arith.constant 1 : index
  %c2 = arith.constant 2 : index

  %x0 = memref.load %coordinate_dofs[%c0, %c0] : memref<3x2xf64>
  %y0 = memref.load %coordinate_dofs[%c0, %c1] : memref<3x2xf64>
  %x1 = memref.load %coordinate_dofs[%c1, %c0] : memref<3x2xf64>
  %y1 = memref.load %coordinate_dofs[%c1, %c1] : memref<3x2xf64>
  %x2 = memref.load %coordinate_dofs[%c2, %c0] : memref<3x2xf64>
  %y2 = memref.load %coordinate_dofs[%c2, %c1] : memref<3x2xf64>
  %j00 = arith.subf %x1, %x0 : f64
  %j01 = arith.subf %x2, %x0 : f64
  %j10 = arith.subf %y1, %y0 : f64
  %j11 = arith.subf %y2, %y0 : f64
  %ad = arith.mulf %j00, %j11 : f64
  %bc = arith.mulf %j01, %j10 : f64
  %det = arith.subf %ad, %bc : f64
  %rho = math.absf %det : f64
  return %rho : f64
}
```

### Matrix contraction

The contraction consumes the prepared density and constant reference tables.
The output buffer follows a UFCx-like accumulation convention: the caller zeros
`A` before the first contribution. `%w` denotes quadrature weights.

```mlir
memref.global "private" constant @basis : memref<3x3xf64> = dense<[
  [0.6666666666666666, 0.1666666666666667, 0.1666666666666667],
  [0.1666666666666667, 0.6666666666666666, 0.1666666666666667],
  [0.1666666666666667, 0.1666666666666667, 0.6666666666666666]
]>
memref.global "private" constant @weights : memref<3xf64> =
  dense<0.1666666666666667>

func.func @mass_matrix(%A: memref<3x3xf64>, %rho: f64) {
  %B = memref.get_global @basis : memref<3x3xf64>
  %w = memref.get_global @weights : memref<3xf64>

  linalg.generic {
    indexing_maps = [
      affine_map<(i, j, q) -> (q)>,     // w[q]
      affine_map<(i, j, q) -> ()>,      // Cell-constant rho
      affine_map<(i, j, q) -> (q, i)>,  // B[q, i]
      affine_map<(i, j, q) -> (q, j)>,  // B[q, j]
      affine_map<(i, j, q) -> (i, j)>   // A[i, j]
    ],
    iterator_types = ["parallel", "parallel", "reduction"]
  } ins(%w, %rho, %B, %B : memref<3xf64>, f64,
                          memref<3x3xf64>, memref<3x3xf64>)
    outs(%A : memref<3x3xf64>) {
    ^bb0(%weight: f64, %density: f64, %bi: f64, %bj: f64, %acc: f64):
      %scale = arith.mulf %weight, %density : f64
      %left = arith.mulf %scale, %bi : f64
      %term = arith.mulf %left, %bj : f64
      %next = arith.addf %acc, %term : f64
      linalg.yield %next : f64
  }
  return
}
```


| MLIR construct | Role |
| --- | --- |
| `indexing_maps` | Select each operand's entry at `(i, j, q)`; `%B` appears twice with different maps |
| `iterator_types` | Mark `i, j` as independent output dimensions and `q` as a reduction |
| `outs(%A ...)` | Accumulate into the existing output entries; initialisation belongs to the caller |
| Region arguments | Supply the weight, cell density, two basis values and current output accumulator |
| `linalg.yield` | Return the updated accumulator |

The empty map `()` broadcasts the cell-constant density. These are array indexing
maps, distinct from the cell map $\psi$. Explicit loops and accumulators
using [SCF](https://mlir.llvm.org/docs/Dialects/SCFDialect/) are another option.

This representation exposes the array accesses and reduction to MLIR. Loop
order, temporary storage and vectorisation remain compiler choices.

# Local fields and library bindings

## Local and assembled DOFs

Two P1 triangle slots have six independent local coefficients. A library map
can associate them with four assembled DOFs:

```text
cell slot:           0        1
local coefficient:   0  1  2  3  4  5
assembled DOF:       0  1  2  0  2  3
```

Gathering repeats shared values; scattering sums their contributions. The local
matrix remains 6 by 6, while the assembled matrix is 4 by 4. A different map can
keep all six DOFs independent without changing the expression.

## Air, liquid and solid

For stationary heat conduction, use three volume calculations and two interface
calculations:

```text
volumes:     air cell       liquid cell       solid cell
interfaces:  air cell -- liquid cell     liquid cell -- solid cell
```

Each role is one local triangle. The library selects actual cells independently
for each invocation; the two liquid interface neighbours need not coincide.
No three-cell patch or description of physical region shapes is required.

The following declarations use the existing reference-library constructors:

```python
K, I = ReferenceSimplex(2), ReferenceSimplex(1)

# One local triangle role per material, with library vertex labels.
Ta = StaticTopology((K,), ((0, 1, 2),))
Tl = StaticTopology((K,), ((1, 3, 2),))
Ts = StaticTopology((K,), ((3, 4, 2),))

# Each interface has one local interval and two declared parent selections.
Sal = StaticTopology((I,), ((1, 2),))
Sal_sides = static_sides(Sal, ((("air", Ta, 0), ("liquid", Tl, 0)),))
Sls = StaticTopology((I,), ((3, 2),))
Sls_sides = static_sides(Sls, ((("liquid", Tl, 0), ("solid", Ts, 0)),))
```

The library supplies vector coordinate spaces `Xa`, `Xl`, `Xs` on these
one-cell topologies, each with value shape `(2,)`. Their functions stay symbolic:

```python
xa, xl, xs = function(Xa, "xa"), function(Xl, "xl"), function(Xs, "xs")
Da = mesh(xa, EuclideanSpace(2))
Dl = mesh(xl, EuclideanSpace(2))
Ds = mesh(xs, EuclideanSpace(2))

Eal = mesh(pull_back(xa, side_map(Ta, Sal, *Sal_sides["air"].entries)), EuclideanSpace(2))
Els = mesh(pull_back(xl, side_map(Tl, Sls, *Sls_sides["liquid"].entries)), EuclideanSpace(2))
```

Let `Ua`, `Ul`, `Us` be scalar P1 temperatures pushed forward to `Da`, `Dl`,
`Ds`, with test fields `Wa`, `Wl`, `Ws`. The three volume terms are
$\int_{D_m} k_m\nabla U_m\cdot\nabla W_m$, for material $m$.
Each temperature has three local basis DOFs.

Finite interface conductance allows temperature jumps. For air–liquid exchange:

```python
ma = side_map(Da, Eal, *Sal_sides["air"].entries)
ml = side_map(Dl, Eal, *Sal_sides["liquid"].entries)
jump_U = pull_back(Ua, ma) - pull_back(Ul, ml)
jump_W = pull_back(Wa, ma) - pull_back(Wl, ml)
a_air_liquid = form(integral(hal * jump_U * jump_W, Eal), ((va, vl), (ua, ul)))
```

Here `va`, `vl`, `ua`, `ul` are the underlying test and trial functions.
Liquid–solid exchange has the same form on `Els` with conductance `hls`.
The volume tensors are 3 by 3; each interface tensor contains four 3 by 3 blocks.
Thus five integration contributions update seven assembled material blocks:

```text
                air  liquid solid
      air       [ *     *     0 ]
      liquid    [ *     *     * ]
      solid     [ 0     *     * ]
```

The package binds liquid temperature and geometry in each calculation:

| Calculation | Selected liquid cell |
| --- | --- |
| Liquid volume | A volume cell |
| Air–liquid interface | The liquid neighbour of the selected air interface facet |
| Liquid–solid interface | The liquid neighbour of the selected solid interface facet |

Each invocation has its own input environment and output DOF maps. A local
function may be reused across calculations, or distinct functions may bind to
the same library field. The block system does not prescribe a kernel decomposition: five kernels, one
per integration domain, or eleven, splitting each interface into four material
blocks, can produce the same assembled system. The form compiler and library
choose the kernel decomposition and corresponding assembly routine.

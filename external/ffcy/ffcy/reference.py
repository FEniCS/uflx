# NumPy and Basix implementations of actions, used to check generated code.

from itertools import product

import basix
import numpy as np
import numpy.typing as npt
from basix import CellType


def tp_element(cell: CellType, degree: int) -> basix.finite_element.FiniteElement:
    variant = basix.LagrangeVariant.gll_warped
    return basix.create_tp_element(basix.ElementFamily.P, cell, degree, variant)


def random_cells(
    cell: CellType, num_cells: int, rng: np.random.Generator, scale: float = 0.15
) -> npt.NDArray[np.float64]:
    """Coordinate dofs of perturbed reference cells in tensor product order."""
    vertices = tp_element(cell, 1).points
    return vertices + scale * rng.uniform(-1, 1, (num_cells, *vertices.shape))


def geometry_degree(x: npt.NDArray[np.float64]) -> int:
    """The degree of coordinate dofs with (k + 1)^d nodes per cell."""
    d = x.shape[2]
    return round(x.shape[1] ** (1 / d)) - 1


def raise_geometry(
    cell: CellType, x: npt.NDArray[np.float64], degree: int
) -> npt.NDArray[np.float64]:
    """The coordinate dofs of degree `degree` of the same cells as degree 1 dofs `x`.

    Each cell's multilinear map is evaluated at the nodes of the higher degree
    element, so the cells keep their shape.
    """
    nodes = tp_element(cell, degree).points
    phi = tp_element(cell, 1).tabulate(0, nodes)[0, :, :, 0]
    return np.einsum("av,cvg->cag", phi, x)


def jacobian(cell, x, quadrature_degree):
    """J[c, q, g, r], the derivative of x_g along reference axis r."""
    points, _ = basix.make_quadrature(cell, quadrature_degree)
    d = x.shape[2]
    dphi = tp_element(cell, geometry_degree(x)).tabulate(1, points)[1 : d + 1, :, :, 0]
    return np.einsum("cag,rqa->cqgr", x, dphi)


def jacobian_determinant(cell, x, quadrature_degree):
    return np.linalg.det(jacobian(cell, x, quadrature_degree))


def mass_action(cell, x, u, degree, quadrature_degree):
    points, weights = basix.make_quadrature(cell, quadrature_degree)
    phi = tp_element(cell, degree).tabulate(0, points)[0, :, :, 0]
    scale = weights * np.abs(jacobian_determinant(cell, x, quadrature_degree))
    return (scale * (u @ phi.T)) @ phi


def laplace_action(cell, x, u, degree, quadrature_degree):
    points, weights = basix.make_quadrature(cell, quadrature_degree)
    d = x.shape[2]
    dphi = tp_element(cell, degree).tabulate(1, points)[1 : d + 1, :, :, 0]
    J = jacobian(cell, x, quadrature_degree)
    # Physical gradients are J^-T times reference gradients.
    grad_phi = np.einsum("cqrg,rqi->cqgi", np.linalg.inv(J), dphi)
    grad_u = np.einsum("cqgi,ci->cqg", grad_phi, u)
    scale = weights * np.abs(np.linalg.det(J))
    return np.einsum("cq,cqgi,cqg->ci", scale, grad_phi, grad_u)


def helmholtz_action(cell, x, u, degree, quadrature_degree):
    return laplace_action(cell, x, u, degree, quadrature_degree) + mass_action(
        cell, x, u, degree, quadrature_degree
    )


def box_mesh(
    cell: CellType, degree: int, n: int, rng: np.random.Generator, scale: float = 0.15
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int32], int]:
    """A box of n^d cells with perturbed interior vertices.

    Returns the coordinate dofs, the dofmap in tensor product order and the number of
    dofs. Global dofs are numbered row-major on the lattice of all cells' dof points.
    """
    d = len(basix.cell.topology(cell)) - 1
    grid = np.stack(np.meshgrid(*[np.arange(n + 1)] * d, indexing="ij"), axis=-1) / n
    interior = (slice(1, n),) * d
    grid[interior] += scale / n * rng.uniform(-1, 1, grid[interior].shape)

    cells = np.array(list(product(range(n), repeat=d)))
    corners = tp_element(cell, 1).points.astype(int)
    x = np.array([[grid[tuple(c + v)] for v in corners] for c in cells])

    # The position of each dof point on its cell's 1D lattice.
    points = tp_element(cell, degree).points
    lattice = np.stack(
        [np.unique(points[:, a].round(12), return_inverse=True)[1] for a in range(d)], axis=1
    )
    size = degree * n + 1
    strides = size ** np.arange(d - 1, -1, -1)
    dofmap = (cells[:, None, :] * degree + lattice[None, :, :]) @ strides
    return x, dofmap.astype(np.int32), size**d


def transpose_dofmap(
    dofmap: npt.NDArray[np.int32], num_dofs: int, width: int
) -> npt.NDArray[np.int32]:
    """For each dof, the flat E-vector indices c * ndofs + i that map to it, padded with -1.

    `width` is the most cells that can share a dof, 2^d on a box mesh.
    """
    transpose = np.full((num_dofs, width), -1, dtype=np.int32)
    filled = np.zeros(num_dofs, dtype=int)
    for e, dof in enumerate(dofmap.ravel()):
        transpose[dof, filled[dof]] = e
        filled[dof] += 1
    return transpose


def assembled_action(action, cell, x, dofmap, num_dofs, v, degree, quadrature_degree):
    """The action on an assembled vector: gather v onto cells, apply, and sum back."""
    y = np.zeros(num_dofs)
    np.add.at(y, dofmap, action(cell, x, v[dofmap], degree, quadrature_degree))
    return y

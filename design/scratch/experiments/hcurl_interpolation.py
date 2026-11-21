"""Interpolate vector P1 into Nedelec on an affine physical triangle."""

import basix
import numpy as np

from uflx_compiler import compile_form, evaluate
from uflx_spec_interpreter import (
    RD,
    cell_map,
    directional_derivative,
    domain_map,
    function,
    image,
    integral,
    num,
    pull_back,
    push_forward,
    tensor,
)
from uflx_toy.basix_cells import basix_cell
from uflx_toy.basix_spaces import basix_space
from uflx_toy.interpolation import interpolate
from uflx_toy.nodes import VertexCoordinates
from uflx_toy.topology import static_topology


def dot(a, b):
    """Contract two real vectors."""
    return a[0] * b[0] + a[1] * b[1]


def point_value(f, point, cell):
    """Construct a domain-free symbolic point evaluation."""
    position = tensor(tuple(num(float(x)) for x in point))
    return integral(pull_back(f, domain_map(RD(0), position, cell)), RD(0))


def mapped_values(source, target, coefficients, points, J):
    """Apply Basix's interpolation operator with the target's Piola mappings."""
    operator = basix.compute_interpolation_operator(source, target)
    dofs = operator @ (coefficients @ J).ravel()
    reference = np.einsum("qac,a->qc", target.tabulate(0, points)[0], dofs)
    return reference @ np.linalg.inv(J)


def main():
    """Compile interpolation, its energy and directional derivatives."""
    cell = basix_cell(basix.CellType.triangle)
    topology = static_topology((cell,), ((0, 1, 2),))
    space = basix_space("P1-vector", topology, basix.ElementFamily.P, 1, shape=(2,))
    u, h = function(space, "u"), function(space, "h")
    source = space.element(0)
    vertices = ((0.7, -0.4), (2.4, -0.6), (1.1, 0.9))
    D = image(VertexCoordinates(topology, vertices), RD(2))
    C = cell_map(D, 0)
    J = (np.asarray(vertices)[1:] - vertices[0]).T
    target = basix_space("Nedelec", topology, basix.ElementFamily.N1E, 1)
    field = push_forward(u, D, space.law)
    reduced = interpolate(field, target)
    rng = np.random.default_rng(42)
    coefficients, direction = rng.normal(size=(2, source.dim, 2))
    env = {"coefficients:u": coefficients.ravel(), "coefficients:h": direction.ravel()}
    points = np.array([[0.13, 0.21], [0.6, 0.1], [0.0, 0.0], [0.2, 0.8]])
    compiled = [compile_form(point_value(pull_back(reduced, C), p, cell), 4) for p in points]
    actual = np.array([evaluate(e, env) for e in compiled])
    expected = mapped_values(source, target.element(0), coefficients, points, J)
    np.testing.assert_allclose(actual, expected, atol=1e-12)
    print(f"Physical P1-vector -> Nedelec: max error {np.max(np.abs(actual - expected)):.3e}")

    # Every lowest-order Nedelec field must be reproduced exactly.
    p = source.points
    exact_coefficients = np.column_stack((2 - 3 * p[:, 1], -1 + 3 * p[:, 0])) @ np.linalg.inv(J)
    actual = np.array(
        [evaluate(e, {"coefficients:u": exact_coefficients.ravel()}) for e in compiled]
    )
    expected = np.column_stack((2 - 3 * points[:, 1], -1 + 3 * points[:, 0])) @ np.linalg.inv(J)
    np.testing.assert_allclose(actual, expected, atol=1e-12)
    print("Nedelec field reproduction: passed")

    q, weights = basix.make_quadrature(basix.CellType.triangle, 4)
    energy = integral(dot(reduced, reduced) / num(2), D)
    derivative = directional_derivative(energy, u, u, h)
    second = directional_derivative(derivative, u, u, h)
    values = mapped_values(source, target.element(0), coefficients, q, J)
    variations = mapped_values(source, target.element(0), direction, q, J)
    expected = abs(np.linalg.det(J)) * np.array(
        (
            np.einsum("q,qc,qc->", weights, values, values) / 2,
            np.einsum("q,qc,qc->", weights, values, variations),
            np.einsum("q,qc,qc->", weights, variations, variations),
        )
    )
    for label, expression, wanted in zip(
        ("energy", "first derivative", "second derivative"),
        (energy, derivative, second),
        expected,
        strict=True,
    ):
        actual = evaluate(compile_form(expression, 4), env)
        np.testing.assert_allclose(actual, wanted, atol=1e-12)
        print(f"Compiled {label}: error {abs(actual - wanted):.3e}")


if __name__ == "__main__":
    main()

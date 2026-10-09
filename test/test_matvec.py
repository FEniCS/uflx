"""Test matrix-vector products."""

from uflx.expressions import MatVec, RealScalar
from uflx.tensors import Matrix, Vector


def test_non_square_matvec():
    """A (3, 2) matrix times a vector of length 2 is a vector of length 3."""
    a = Matrix([[RealScalar(float(3 * i + j)) for j in range(2)] for i in range(3)])
    x = Vector([RealScalar(1.0), RealScalar(2.0)])

    product = MatVec(a, x)
    assert product.value_shape == (3,)
    for i in range(3):
        assert product.component(i).as_float() == (3 * i) * 1.0 + (3 * i + 1) * 2.0

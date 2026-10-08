"""Test uflx.tensors.Matrix's compute_inverse/compute_determinant."""

from itertools import product

import numpy as np
import pytest

from uflx.expressions import RealScalar
from uflx.tensors import FlattenedTensorMap, Matrix


def _to_matrix(values: np.ndarray) -> Matrix:
    return Matrix([[RealScalar(float(v)) for v in row] for row in values])


def _to_numpy(matrix: Matrix) -> np.ndarray:
    rows, cols = matrix.value_shape
    return np.array([[matrix.component(i, j).as_float() for j in range(cols)] for i in range(rows)])


@pytest.mark.parametrize("n", [1, 2, 3])
def test_compute_inverse(n):
    """compute_inverse() should match numpy's inverse for 1x1, 2x2 and 3x3."""
    rng = np.random.default_rng(0)
    values = rng.uniform(0.5, 2.0, (n, n))

    inverse = _to_numpy(_to_matrix(values).compute_inverse())

    np.testing.assert_allclose(inverse, np.linalg.inv(values), rtol=1e-12)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_compute_determinant(n):
    """compute_determinant() should match numpy's determinant for 1x1, 2x2 and 3x3."""
    rng = np.random.default_rng(1)
    values = rng.uniform(0.5, 2.0, (n, n))

    det = _to_matrix(values).compute_determinant().as_float()

    np.testing.assert_allclose(det, np.linalg.det(values), rtol=1e-12)


def test_compute_inverse_not_implemented_for_4x4():
    """compute_inverse() should still raise a clear error for unsupported sizes."""
    matrix = _to_matrix(np.eye(4))
    with pytest.raises(NotImplementedError):
        matrix.compute_inverse()


@pytest.mark.parametrize("shape", [(2, 1), (3, 1), (3, 2), (1, 2), (1, 3), (2, 3)])
def test_compute_determinant_non_square(shape):
    """For a non-square matrix, compute_determinant() is the pseudo-determinant."""
    rng = np.random.default_rng(2)
    values = rng.uniform(0.5, 2.0, shape)
    gram = values.T @ values if shape[0] > shape[1] else values @ values.T

    det = _to_matrix(values).compute_determinant().as_float()

    np.testing.assert_allclose(det, np.sqrt(np.linalg.det(gram)), rtol=1e-12)


def test_compute_determinant_interval_in_2d():
    """The Jacobian (3, 4)^T of an interval in 2D scales lengths by 5."""
    det = _to_matrix(np.array([[3.0], [4.0]])).compute_determinant().as_float()

    assert det == pytest.approx(5.0, rel=1e-14)


@pytest.mark.parametrize("shape", [(4,), (5, 3), (2, 3, 4), (2, 3, 1, 2)])
def test_flattened_tensor_map_matches_numpy(shape):
    """The flat index of an entry is that of a C-ordered numpy array."""
    for indices in product(*(range(n) for n in shape)):
        expected = int(np.ravel_multi_index(indices, shape))
        assert FlattenedTensorMap(indices, shape[1:]).flat_index == expected


@pytest.mark.parametrize("index, trailing_shape", [((7, 2), (3,)), ((1, 2, 3), (3, 4)), ((5,), ())])
def test_flattened_tensor_map_any_number_of_rows(index, trailing_shape):
    """Get the flat index of an entry `index` in a (Z, trailing) shaped tensor."""
    exact_flat_index = 0
    for i in range(len(index)):
        exact_flat_index += index[i] * np.prod(trailing_shape[i:])
    assert FlattenedTensorMap(index, trailing_shape).flat_index == int(exact_flat_index)


def test_flattened_tensor_map_is_a_scalar():
    """An entry is a scalar, without components or successors."""
    entry = FlattenedTensorMap((1, 2), (3,))
    assert entry.value_shape == ()
    assert entry.successors == set()
    with pytest.raises(ValueError, match="scalar"):
        entry.component(0)


def test_flattened_tensor_map_invalid():
    """Invalid indices and extents raise."""
    with pytest.raises(ValueError, match="one more"):
        FlattenedTensorMap((0,), (3,))
    with pytest.raises(ValueError, match="one more"):
        FlattenedTensorMap((), ())
    with pytest.raises(ValueError, match="positive"):
        FlattenedTensorMap((0, 0), (0,))
    with pytest.raises(IndexError):
        FlattenedTensorMap((0, 3), (3,))
    with pytest.raises(IndexError):
        FlattenedTensorMap((-1, 0), (3,))


def test_flattened_tensor_map_equality():
    """Entries are equal if their indices and extents are."""
    entry = FlattenedTensorMap((1, 2), (3,))
    assert entry == FlattenedTensorMap((1, 2), (3,))
    assert hash(entry) == hash(FlattenedTensorMap((1, 2), (3,)))
    assert entry != FlattenedTensorMap((2, 1), (3,))
    assert entry != FlattenedTensorMap((1, 2), (4,))

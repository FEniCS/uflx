"""Tensors."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, TypeAlias, cast

from uflx.expressions import AbstractExpression, RealScalar, Sqrt, expression_sum
from uflx.graphs import GraphNode

NestedSequence: TypeAlias = AbstractExpression | Sequence["NestedSequence"]
NestedTuple: TypeAlias = AbstractExpression | tuple["NestedTuple", ...]


class FlattenedTensorMap(AbstractExpression):
    """A scalar entry of a row-major tensor, given by its position in the flattened 1D array.

    The extent of the first dimension does not affect the position of an entry, so it is not
    an input: the entry is the same in a tensor with any number of rows.

    Args:
        indices: The N indices of the entry.
        trailing_shape: The N - 1 extents of the tensor from the second dimension onwards.
    """

    def __init__(self, indices: tuple[int, ...], trailing_shape: tuple[int, ...]):
        """Initialise."""
        if len(indices) != len(trailing_shape) + 1:
            raise ValueError("Number of indices must be one more than the number of extents.")
        if any(n <= 0 for n in trailing_shape):
            raise ValueError(f"Extents {trailing_shape} must be positive.")
        if indices[0] < 0 or any(i < 0 or i >= n for i, n in zip(indices[1:], trailing_shape)):
            raise IndexError(f"Indices {indices} out of range for extents {trailing_shape}.")

        self._indices = indices
        self._trailing_shape = trailing_shape

    @property
    def flat_index(self) -> int:
        """The position of the entry in the flattened array."""
        index = self._indices[0]
        for i, n in zip(self._indices[1:], self._trailing_shape):
            index = index * n + i
        return index

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._indices, self._trailing_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class Tensor(AbstractExpression):
    """A general tensor."""

    def __init__(self, entries: NestedSequence, shape: tuple[int, ...] | None = None):
        """Initalise.

        Args:
            entries: The entries, nested as deeply as the rank
            shape: The shape, where the entries cannot give one: an axis
                of length zero says nothing about the axes below it, so a
                matrix with no rows needs telling how many columns it
                has. It must agree with the entries as far as they go.

        Raises:
            ValueError: If the shape disagrees with the entries
        """

        def to_shape_and_tuple(items) -> tuple[tuple[int, ...], NestedTuple]:
            if isinstance(items, AbstractExpression):
                # A non-scalar entry is stacked: its shape is appended.
                return items.value_shape, items
            s: tuple[int, ...] | None = None
            t = []
            for i in items:
                sub_s, sub_t = to_shape_and_tuple(i)
                if s is None:
                    s = sub_s
                if s != sub_s:
                    raise ValueError(f"Tensor entries of shapes {s} and {sub_s} cannot be stacked.")
                t.append(sub_t)
            if s is None:
                # Nothing to stack, so this axis is of length zero and there
                # are no axes below it. A matrix with no rows cannot be read
                # this way, since nothing says how many columns it has, so
                # the shape is given instead.
                return (0,), ()
            return (len(t), *s), tuple(t)

        self._shape, self._entries = to_shape_and_tuple(entries)
        if shape is not None:
            # An axis of length zero says nothing about the axes below it:
            # a matrix with no rows has some number of columns, and the
            # entries cannot say which. So a given shape may extend the
            # one read off the entries past its first zero, and must
            # agree with it up to there.
            zeros = [i for i, length in enumerate(self._shape) if length == 0]
            known = zeros[0] + 1 if zeros else len(self._shape)
            if shape[:known] != self._shape[:known]:
                raise ValueError(f"These entries are of shape {self._shape}, not {shape}.")
            self._shape = shape

    def __repr__(self):
        """Representation."""
        return f"Tensor({self._entries})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

        def extract_successors(items: NestedTuple) -> set[GraphNode]:
            if isinstance(items, GraphNode):
                return {items}
            return set().union(*[extract_successors(i) for i in items])

        return extract_successors(self._entries)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object.

        The shape is one of them: two tensors of the same (empty) entries
        and different shapes are different tensors, and a tensor rebuilt
        from its arguments has to come back the same.
        """
        return self._entries, self._shape

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""

        def extract_component(items: NestedTuple, indices: tuple[int, ...]) -> AbstractExpression:
            if isinstance(items, AbstractExpression):
                return items.component(*indices) if len(indices) > 0 else items
            assert len(indices) > 0
            return extract_component(items[indices[0]], indices[1:])

        return extract_component(self._entries, indices)

    def _as_a_matrix(self) -> tuple[int, int]:
        """This tensor's shape, as a matrix's.

        The algebra below is a matrix's, so it asks for the rank it needs
        rather than being defined on a class of that rank: a tensor is
        told apart by its shape, and a rank two tensor is a matrix.

        Raises:
            ValueError: If this tensor is not of rank two
        """
        if len(self._shape) != 2:
            raise ValueError(f"A tensor of shape {self._shape} is not a matrix.")
        rows, cols = self._shape
        return rows, cols

    def transpose(self) -> Tensor:
        """Get the transpose of this tensor, read as a matrix."""
        rows, cols = self._as_a_matrix()
        entries = cast(tuple[tuple[AbstractExpression, ...], ...], self._entries)
        return Tensor(
            [[entries[i][j] for i in range(rows)] for j in range(cols)], shape=(cols, rows)
        )

    def matmat(self, other: Tensor) -> Tensor:
        """Compute a matrix-matrix product.

        Raises:
            ValueError: If either tensor is not a matrix, or the shapes do
                not meet
        """
        rows, cols = self._as_a_matrix()
        other_rows, other_cols = other._as_a_matrix()
        if cols != other_rows:
            raise ValueError(
                f"A {rows}x{cols} matrix does not multiply a {other_rows}x{other_cols} one."
            )
        return Tensor(
            [
                [
                    expression_sum(
                        self.component(i, k) * other.component(k, j) for k in range(cols)
                    )
                    for j in range(other_cols)
                ]
                for i in range(rows)
            ],
            shape=(rows, other_cols),
        )

    def compute_inverse(self) -> Tensor:
        """Compute the inverse of this tensor, read as a matrix."""
        rows, cols = self._as_a_matrix()
        if rows > cols:
            return self.transpose().matmat(self).compute_inverse().matmat(self.transpose())
        elif rows < cols:
            return self.transpose().matmat(self.matmat(self.transpose()).compute_inverse())
        else:
            entries = cast(tuple[tuple[AbstractExpression, ...], ...], self._entries)
            match rows:
                case 0:
                    return Tensor([], shape=(0, 0))
                case 1:
                    [[a]] = entries
                    assert isinstance(a, AbstractExpression)
                    return Tensor([[RealScalar(1.0) / a]])
                case 2:
                    [[a, b], [c, d]] = entries
                    assert isinstance(a, AbstractExpression)
                    assert isinstance(b, AbstractExpression)
                    assert isinstance(c, AbstractExpression)
                    assert isinstance(d, AbstractExpression)
                    det = a * d - b * c
                    return Tensor([[d / det, -b / det], [-c / det, a / det]])
                case 3:
                    [[a, b, c], [d, e, f], [g, h, i]] = entries
                    assert isinstance(a, AbstractExpression)
                    assert isinstance(b, AbstractExpression)
                    assert isinstance(c, AbstractExpression)
                    assert isinstance(d, AbstractExpression)
                    assert isinstance(e, AbstractExpression)
                    assert isinstance(f, AbstractExpression)
                    assert isinstance(g, AbstractExpression)
                    assert isinstance(h, AbstractExpression)
                    assert isinstance(i, AbstractExpression)
                    det = a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g)
                    # Inverse = adjugate / det, adjugate = cofactor matrix transposed.
                    return Tensor(
                        [
                            [
                                (e * i - f * h) / det,
                                (c * h - b * i) / det,
                                (b * f - c * e) / det,
                            ],
                            [
                                (f * g - d * i) / det,
                                (a * i - c * g) / det,
                                (c * d - a * f) / det,
                            ],
                            [
                                (d * h - e * g) / det,
                                (b * g - a * h) / det,
                                (a * e - b * d) / det,
                            ],
                        ]
                    )
                case _:
                    raise NotImplementedError(f"Inverting of {rows}x{rows} not implemented.")

    def compute_determinant(self) -> AbstractExpression:
        """Compute the determinant of this tensor, read as a matrix.

        For a non-square matrix A this is the pseudo-determinant sqrt(det(A^T A)) (tall A) or
        sqrt(det(A A^T)) (wide A), i.e. the volume scaling of the reference-to-physical map of a
        manifold cell.

        An empty matrix gives one. A block diagonal determinant is the
        product of its blocks', so adding a 0 by 0 block must change
        nothing: det(diag(A, I0)) = det(A) det(I0) = det(A). numpy and
        MATLAB answer one too.

        Raises:
            ValueError: If this tensor is not a matrix
        """
        rows, cols = self._as_a_matrix()
        if rows == 0 or cols == 0:
            # The empty Gram determinant is one: a map out of no coordinates
            # scales no volume, so a domain of topological dimension zero has
            # a density of one, which is the counting measure.
            return RealScalar(1.0)
        if rows > cols:
            return Sqrt(self.transpose().matmat(self).compute_determinant())
        elif rows < cols:
            return Sqrt(self.matmat(self.transpose()).compute_determinant())
        else:
            entries = cast(tuple[tuple[AbstractExpression, ...], ...], self._entries)
            match rows:
                case 0:
                    return RealScalar(1.0)
                case 1:
                    [[a]] = entries
                    assert isinstance(a, AbstractExpression)
                    return a
                case 2:
                    [[a, b], [c, d]] = entries
                    assert isinstance(a, AbstractExpression)
                    assert isinstance(b, AbstractExpression)
                    assert isinstance(c, AbstractExpression)
                    assert isinstance(d, AbstractExpression)
                    return a * d - b * c
                case _:
                    [[a, b, c], [d, e, f], [g, h, i]] = entries
                    assert isinstance(a, AbstractExpression)
                    assert isinstance(b, AbstractExpression)
                    assert isinstance(c, AbstractExpression)
                    assert isinstance(d, AbstractExpression)
                    assert isinstance(e, AbstractExpression)
                    assert isinstance(f, AbstractExpression)
                    assert isinstance(g, AbstractExpression)
                    assert isinstance(h, AbstractExpression)
                    assert isinstance(i, AbstractExpression)
                    return a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g)


def tensor(entries: NestedSequence, shape: tuple[int, ...] | None = None) -> Tensor:
    """Create a tensor of the given entries.

    Args:
        entries: The entries, nested as deeply as the rank
        shape: The shape, where the entries cannot give one

    Returns:
        That tensor
    """
    return Tensor(entries, shape)


def _of_rank(
    entries: NestedSequence, rank: int, shape: tuple[int, ...] | None, name: str
) -> Tensor:
    """Create a tensor, insisting on a rank.

    There is one tensor class, told apart by its shape, so a vector and a
    matrix are ranks rather than kinds. Asking for one is still worth it
    at the point of writing a literal, where a mis-nested entry is a
    mistake to catch early.

    Raises:
        ValueError: If the entries are not of that rank
    """
    built = Tensor(entries, shape)
    if len(built.value_shape) != rank:
        raise ValueError(
            f"A {name} is of rank {rank}, and these entries are of shape {built.value_shape}."
        )
    return built


def vector(entries: NestedSequence, shape: tuple[int, ...] | None = None) -> Tensor:
    """Create a tensor of rank one.

    Args:
        entries: The entries, which must be scalars
        shape: The shape, where the entries cannot give one

    Returns:
        That tensor
    """
    return _of_rank(entries, 1, shape, "vector")


def matrix(entries: NestedSequence, shape: tuple[int, ...] | None = None) -> Tensor:
    """Create a tensor of rank two.

    Args:
        entries: The rows, whose entries must be scalars
        shape: The shape, where the entries cannot give one, as a matrix
            with no rows cannot

    Returns:
        That tensor
    """
    return _of_rank(entries, 2, shape, "matrix")


def zero(shape: tuple[int, ...]) -> RealScalar | Tensor:
    """Create a tensor full of zeros.

    Args:
        shape: The value shape of the zero tensor to build -- () for a
            scalar, or any tuple of positive extents for a
            correctly-shaped all-zero tensor.

    Returns:
        RealScalar(0.0) if shape == (), else a tensor of that shape built
        from nested RealScalar(0.0) leaves.

    Raises:
        ValueError: if any extent in `shape` is not a positive integer --
            eg a zero or negative extent, which the underlying
            tensor representation cannot express.
    """
    if any(not isinstance(n, int) or n <= 0 for n in shape):
        raise ValueError(
            f"Cannot create a zero tensor of shape {shape}: every extent must be "
            "a positive integer."
        )
    if shape == ():
        return RealScalar(0.0)
    if len(shape) == 1:
        return vector([RealScalar(0.0) for _ in range(shape[0])])
    if len(shape) == 2:
        return matrix(
            [[RealScalar(0.0) for _ in range(shape[1])] for _ in range(shape[0])],
            shape=shape,
        )

    def build_entries(shape: tuple[int, ...]) -> NestedSequence:
        """Recursively build a zero tensor."""
        if len(shape) == 1:
            return [RealScalar(0.0) for _ in range(shape[0])]
        return [build_entries(shape[1:]) for _ in range(shape[0])]

    return Tensor(build_entries(shape))


class Identity(AbstractExpression):
    """Identity matrix."""

    def __init__(self, size: int):
        """Initialise."""
        self.size = size

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.size, self.size)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.size,)

    def __repr__(self) -> str:
        """Representation."""
        return f"Identity({self.size})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        assert len(indices) == 2 and all(i < self.size for i in indices)
        return RealScalar(1.0 if indices[0] == indices[1] else 0.0)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        return other

    def simplified_matrix_product_right(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        return other

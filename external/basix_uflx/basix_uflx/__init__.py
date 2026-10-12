"""Basix UFLx interface."""

from abc import abstractmethod
from collections.abc import Sequence
from typing import Any
from warnings import warn

import basix
import numpy as np
import numpy.typing as npt
import uflx

from basix_uflx.finite_element import AbstractFiniteElement, BlockedElement, MixedElement
from basix_uflx.utils import number_of_derivatives

__all__ = [
    "blocked_element",
    "custom_element",
    "element",
    "mixed_element",
    "quadrature_element",
    "real_element",
    "wrap_element",
]


def convert_map(basix_map: basix.MapType) -> uflx.maps.AbstractReferenceMap:
    """Convert a basix map to a UFLx map."""
    match basix_map:
        case basix.MapType.identity:
            return uflx.maps.IdentityReferenceMap()
        case _:
            raise NotImplementedError()


class BasixCell(uflx.entities.AbstractEntity):
    """A cell defined by Basix."""

    def __init__(self, basix_cell: basix.CellType):
        self._basix_cell = basix_cell

    @property
    def topological_dimension(self) -> int:
        """Topological dimension of the entity."""
        return len(basix.cell.topology(self._basix_cell)) - 1

    @property
    def name(self) -> str:
        """Name of the entity type."""
        return self._basix_cell.name

    def sub_entities(self, dim: int) -> list[uflx.entities.AbstractEntity]:
        """Get a list of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            A list of sub-entities of the given dimension.
        """
        return [BasixCell(et) for et in basix.cell.subentity_types(self._basix_cell)[dim]]

    def sub_entity_vertices(self, dim: int) -> list[list[int]]:
        """Get lists of the vertices of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            A list of lists of vertices of sub-entities of the given dimension.
        """
        return [c[0] for c in basix.cell.sub_entity_connectivity(self._basix_cell)[dim]]


def _as_tuples(data: Any) -> Any:
    """Convert nested lists to nested tuples so they can be hashed."""
    if isinstance(data, list):
        return tuple(_as_tuples(i) for i in data)
    return data


ArrayData = tuple[tuple[int, ...], tuple[float, ...]]


def _array_data(array: npt.ArrayLike) -> ArrayData:
    """Store an array as its shape and its entries, so that an empty array keeps its shape."""
    a = np.asarray(array)
    return tuple(a.shape), tuple(a.ravel().tolist())


def _array(data: ArrayData, dtype: str) -> npt.NDArray:
    shape, entries = data
    return np.array(entries, dtype=dtype).reshape(shape)


class _WrappedBasixElement(AbstractFiniteElement):
    """An element that Basix tabulates, built from its arguments."""

    _basix_element: basix.finite_element.FiniteElement

    @property
    def cell(self) -> uflx.entities.AbstractEntity:
        return BasixCell(self._basix_element.cell_type)

    @property
    def real_valued(self) -> bool:
        return True

    @property
    def dim(self) -> int:
        return self._basix_element.dim

    @property
    def lagrange_superdegree(self) -> int:
        d = self._basix_element.embedded_superdegree
        if d is None:
            raise ValueError("Element does not have a well-defined Lagrange superdegree")
        return d

    @property
    def reference_map(self) -> uflx.maps.AbstractReferenceMap:
        return convert_map(self._basix_element.map_type)

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        return tuple(self._basix_element.value_shape)

    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        return np.asarray(self._basix_element.tabulate(derivatives, np.asarray(points)))


class BasixElement(_WrappedBasixElement):
    """An element from one of Basix's families.

    Users should not call this class's initialiser directly, but should
    use the `element` function instead.
    """

    def __init__(
        self,
        family: basix.ElementFamily,
        cell: basix.CellType,
        degree: int,
        lagrange_variant: basix.LagrangeVariant,
        dpc_variant: basix.DPCVariant,
        discontinuous: bool,
        dof_ordering: tuple[int, ...],
        dtype: str,
    ):
        """Initialise."""
        self._family = family
        self._cell = cell
        self._degree = degree
        self._lagrange_variant = lagrange_variant
        self._dpc_variant = dpc_variant
        self._discontinuous = discontinuous
        self._dof_ordering = dof_ordering
        self._dtype = dtype
        self._basix_element = basix.create_element(
            family,
            cell,
            degree,
            lagrange_variant,
            dpc_variant,
            discontinuous,
            dof_ordering=list(dof_ordering) or None,
            dtype=np.dtype(dtype),
        )


class CustomBasixElement(_WrappedBasixElement):
    """A custom Basix element, given by its coefficients and interpolation data.

    Users should not call this class's initialiser directly, but should
    use the `custom_element` function instead.
    """

    def __init__(
        self,
        cell: basix.CellType,
        value_shape: tuple[int, ...],
        wcoeffs: ArrayData,
        x: tuple[tuple[ArrayData, ...], ...],
        M: tuple[tuple[ArrayData, ...], ...],
        interpolation_nderivs: int,
        map_type: basix.MapType,
        sobolev_space: basix.SobolevSpace,
        discontinuous: bool,
        embedded_subdegree: int,
        embedded_superdegree: int,
        polyset_type: basix.PolysetType,
        dtype: str,
    ):
        """Initialise."""
        self._cell = cell
        self._value_shape = value_shape
        self._wcoeffs = wcoeffs
        self._x = x
        self._M = M
        self._interpolation_nderivs = interpolation_nderivs
        self._map_type = map_type
        self._sobolev_space = sobolev_space
        self._discontinuous = discontinuous
        self._embedded_subdegree = embedded_subdegree
        self._embedded_superdegree = embedded_superdegree
        self._polyset_type = polyset_type
        self._dtype = dtype
        self._basix_element = basix.create_custom_element(
            cell,
            value_shape,
            _array(wcoeffs, dtype),
            [[_array(j, dtype) for j in i] for i in x],
            [[_array(j, dtype) for j in i] for i in M],
            interpolation_nderivs,
            map_type,
            sobolev_space,
            discontinuous,
            embedded_subdegree,
            embedded_superdegree,
            polyset_type,
            dtype=np.dtype(dtype),
        )


class QuadratureElement(AbstractFiniteElement):
    """A quadrature element."""

    def __init__(
        self,
        cell: basix.CellType,
        points: npt.NDArray[np.floating],
        weights: npt.NDArray[np.floating],
        reference_map: uflx.maps.AbstractReferenceMap,
        degree: int | None = None,
        dtype: npt.DTypeLike = np.float64,
    ):
        """Initialise the element."""
        self._points = _as_tuples(np.asarray(points, dtype=dtype).tolist())
        self._weights = _as_tuples(np.asarray(weights, dtype=dtype).tolist())
        self._cell = cell
        self._reference_map = reference_map
        self._degree = len(self._points) if degree is None else degree
        self._dtype = np.dtype(dtype).name

    @property
    def cell(self) -> uflx.entities.AbstractEntity:
        return BasixCell(self._cell)

    @property
    def real_valued(self) -> bool:
        return True

    @property
    def dim(self) -> int:
        return len(self._points)

    @property
    def lagrange_superdegree(self) -> int:
        return self._degree  # TODO: this is not right

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        return ()

    @property
    def reference_map(self) -> uflx.maps.AbstractReferenceMap:
        return self._reference_map

    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        if derivatives > 0:
            raise ValueError("Cannot take derivatives of Quadrature element.")

        points = np.asarray(points)
        own = np.asarray(self._points, dtype=self._dtype)
        if points.shape != own.shape or not np.allclose(points, own):
            raise ValueError("Mismatch of tabulation points and element points.")
        npts = points.shape[0]
        table = np.eye(npts, npts).reshape([1, npts, npts, 1])
        return table


class RealElement(AbstractFiniteElement):
    """A real element."""

    def __init__(self, cell: basix.CellType, value_shape: tuple[int, ...]):
        """Initialise the element."""
        self._cell = cell
        self._value_shape = tuple(value_shape)

    @property
    def cell(self) -> uflx.entities.AbstractEntity:
        return BasixCell(self._cell)

    @property
    def real_valued(self) -> bool:
        return True

    @property
    def dim(self) -> int:
        return self.reference_value_size

    @property
    def lagrange_superdegree(self) -> int:
        return 0

    @property
    def reference_map(self) -> uflx.maps.AbstractReferenceMap:
        if self._value_shape == ():
            return uflx.maps.IdentityReferenceMap()
        else:
            return uflx.maps.BlockedReferenceMap(
                uflx.maps.IdentityReferenceMap(), self._value_shape
            )

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        return self._value_shape

    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        points = np.asarray(points)
        table = np.zeros(
            [
                number_of_derivatives(derivatives, self.cell),
                points.shape[0],
                self.reference_value_size,
                self.reference_value_size,
            ]
        )
        table[0, :, :, :] = 1.0
        return table


def element(
    family: basix.ElementFamily | str,
    cell: basix.CellType | str,
    degree: int,
    lagrange_variant: basix.LagrangeVariant | str = basix.LagrangeVariant.unset,
    dpc_variant: basix.DPCVariant | str = basix.DPCVariant.unset,
    discontinuous: bool = False,
    shape: tuple[int, ...] | None = None,
    symmetry: bool | None = None,
    dof_ordering: list[int] | None = None,
    dtype: npt.DTypeLike | None = None,
) -> AbstractFiniteElement:
    """Create a UFLx compatible element using Basix.

    Args:
        family: Element family/type.
        cell: Element cell type.
        degree: Degree of the finite element.
        lagrange_variant: Variant of Lagrange to be used.
        dpc_variant: Variant of DPC to be used.
        discontinuous: If ``True``, the discontinuous version of the
            element is created.
        shape: Value shape of the element. For scalar-valued families,
            this can be used to create vector and tensor elements.
        symmetry: Set to ``True`` if the tensor is symmetric. Valid for
            rank 2 elements only.
        dof_ordering: Ordering of dofs for ``ElementDofLayout``.
        dtype: Floating point data type.

    Returns:
        A finite element.

    """
    # Conversion of string arguments to types
    if isinstance(cell, str):
        cell = basix.CellType[cell]
    if isinstance(family, str):
        if family.startswith("Discontinuous "):
            family = family[14:]
            discontinuous = True
        if family in ["DP", "DG", "DQ"]:
            family = "P"
            discontinuous = True
        if family == "DPC":
            discontinuous = True
        family = basix.finite_element.string_to_family(family, cell.name)
    if isinstance(lagrange_variant, str):
        lagrange_variant = basix.LagrangeVariant[lagrange_variant]
    if isinstance(dpc_variant, str):
        dpc_variant = basix.DPCVariant[dpc_variant]

    # Default variant choices
    EF = basix.ElementFamily
    if lagrange_variant == basix.LagrangeVariant.unset:
        if family == EF.P:
            lagrange_variant = basix.LagrangeVariant.gll_warped
        elif family in [EF.RT, EF.N1E]:
            lagrange_variant = basix.LagrangeVariant.legendre
        elif family in [EF.serendipity, EF.BDM, EF.N2E]:
            lagrange_variant = basix.LagrangeVariant.legendre

    if dpc_variant == basix.DPCVariant.unset:
        if family in [EF.serendipity, EF.BDM, EF.N2E]:
            dpc_variant = basix.DPCVariant.legendre
        elif family == EF.DPC:
            dpc_variant = basix.DPCVariant.diagonal_gll

    ufl_e = BasixElement(
        family,
        cell,
        degree,
        lagrange_variant,
        dpc_variant,
        discontinuous,
        tuple(dof_ordering or ()),
        np.dtype(dtype or np.float64).name,
    )

    if shape is None or shape == ufl_e.reference_value_shape:
        if symmetry is not None:
            raise ValueError("Cannot pass a symmetry argument to this element.")
        return ufl_e
    else:
        return blocked_element(ufl_e, shape=shape, symmetry=symmetry)


def custom_element(
    cell_type: basix.CellType,
    reference_value_shape: Sequence[int],
    wcoeffs: npt.ArrayLike,
    x: Sequence[Sequence[npt.ArrayLike]],
    M: Sequence[Sequence[npt.ArrayLike]],
    interpolation_nderivs: int,
    map_type: basix.MapType,
    sobolev_space: basix.SobolevSpace,
    discontinuous: bool,
    embedded_subdegree: int,
    embedded_superdegree: int,
    polyset_type: basix.PolysetType = basix.PolysetType.standard,
    dtype: npt.DTypeLike | None = None,
) -> AbstractFiniteElement:
    """Create a UFLx compatible custom Basix element.

    Args:
        cell_type: The cell type
        reference_value_shape: The reference value shape of the element
        wcoeffs: Matrices for the kth value index containing the
            expansion coefficients defining a polynomial basis spanning
            the polynomial space for this element. Shape is
            ``(dim(finite element polyset), dim(Legenre polynomials))``.
        x: Interpolation points. Indices are ``(tdim, entity index,
            point index, dim)``.
        M: The interpolation matrices. Indices are ``(tdim, entity
            index, dof, vs, point_index, derivative)``.
        interpolation_nderivs: The number of derivatives that need to be
            used during interpolation.
        map_type: The type of map to be used to map values from the
            reference to a cell.
        sobolev_space: Underlying Sobolev space for the element.
        discontinuous: Indicates whether or not this is the
            discontinuous version of the element.
        embedded_subdegree: The highest degree ``n`` such that a
            Lagrange (or vector Lagrange) element of degree ``n`` is a
            subspace of this element.
        embedded_superdegree: The highest degree of a polynomial in this
            element's polyset.
        polyset_type: Polyset type for the element.
        dtype: Floating point data type.

    Returns:
        A custom finite element.
    """
    return CustomBasixElement(
        cell_type,
        tuple(reference_value_shape),
        _array_data(wcoeffs),
        tuple(tuple(_array_data(j) for j in i) for i in x),
        tuple(tuple(_array_data(j) for j in i) for i in M),
        interpolation_nderivs,
        map_type,
        sobolev_space,
        discontinuous,
        embedded_subdegree,
        embedded_superdegree,
        polyset_type,
        np.dtype(dtype or np.float64).name,
    )


def mixed_element(elements: Sequence[AbstractFiniteElement]) -> AbstractFiniteElement:
    """Create a UFLx compatible mixed element from a list of elements.

    Args:
        elements: The list of elements

    Returns:
        A mixed finite element.
    """
    return MixedElement(list(elements))


def quadrature_element(
    cell: str | basix.CellType,
    value_shape: Sequence[int] = (),
    scheme: str | None = None,
    degree: int | None = None,
    points: npt.ArrayLike | None = None,
    weights: npt.ArrayLike | None = None,
    reference_map: uflx.maps.AbstractReferenceMap | None = None,
    symmetry: bool | None = None,
    dtype: npt.DTypeLike = np.float64,
) -> AbstractFiniteElement:
    """Create a quadrature element.

    When creating this element, either the quadrature scheme and degree
    must be input or the quadrature points and weights must be.

    Args:
        cell: Cell to create the element on.
        value_shape: Value shape of the element.
        scheme: Quadrature scheme.
        degree: Quadrature degree.
        points: Quadrature points.
        weights: Quadrature weights.
        reference_map: Reference map
        symmetry: Set to ``True`` if the tensor is symmetric. Valid for
            rank 2 elements only.
        dtype: Data type of quadrature points and weights

    Returns:
        A 'quadrature' finite element.
    """
    if isinstance(cell, str):
        cell = basix.CellType[cell]

    if reference_map is None:
        if tuple(value_shape) == ():
            reference_map = uflx.maps.IdentityReferenceMap()
        else:
            reference_map = uflx.maps.BlockedReferenceMap(
                uflx.maps.IdentityReferenceMap(), tuple(value_shape)
            )

    if points is None:
        assert weights is None
        assert degree is not None
        if scheme is None:
            points, weights = basix.make_quadrature(cell, degree)
        else:
            points, weights = basix.make_quadrature(
                cell, degree, rule=basix.quadrature.string_to_type(scheme)
            )

    assert points is not None
    assert weights is not None

    e = QuadratureElement(
        cell, np.asarray(points), np.asarray(weights), reference_map, degree, dtype=dtype
    )
    if tuple(value_shape) == ():
        if symmetry is not None:
            raise ValueError("Cannot pass a symmetry argument to this element.")
        return e
    else:
        return blocked_element(e, shape=tuple(value_shape), symmetry=symmetry)


def real_element(
    cell: basix.CellType | str,
    value_shape: Sequence[int],
) -> AbstractFiniteElement:
    """Create a real element.

    Args:
        cell: Cell to create the element on.
        value_shape: Value shape of the element.

    Returns:
        A 'real' finite element.

    """
    if isinstance(cell, str):
        cell = basix.CellType[cell]

    return RealElement(cell, tuple(value_shape))


def blocked_element(
    sub_element: AbstractFiniteElement,
    shape: Sequence[int],
    symmetry: bool | None = None,
) -> AbstractFiniteElement:
    """Create a UFLx compatible blocked element.

    Args:
        sub_element: Element used for each block.
        shape: Value shape of the element. For scalar-valued families,
            this can be used to create vector and tensor elements.
        symmetry: Set to ``True`` if the tensor is symmetric. Valid for
            rank 2 elements only.

    Returns:
        A blocked finite element.
    """
    if len(sub_element.reference_value_shape) != 0:
        raise ValueError("Cannot create a blocked element containing a non-scalar element.")

    return BlockedElement(sub_element, shape=tuple(shape), symmetry=symmetry)


def wrap_element(element: basix.finite_element.FiniteElement) -> AbstractFiniteElement:
    """Wrap a Basix element as a Basix UFLx element.

    The element is described by what Basix reports about it, so wrapping
    an element gives the same value as building it with `element` or
    `custom_element`.
    """
    dtype = np.dtype(element.dtype).name
    if element.family == basix.ElementFamily.custom:
        return custom_element(
            element.cell_type,
            element.value_shape,
            element.wcoeffs,
            element.x,
            element.M,
            element.interpolation_nderivs,
            element.map_type,
            element.sobolev_space,
            element.discontinuous,
            element.embedded_subdegree,
            element.embedded_superdegree,
            element.polyset_type,
            dtype,
        )
    return BasixElement(
        element.family,
        element.cell_type,
        element.degree,
        element.lagrange_variant,
        element.dpc_variant,
        element.discontinuous,
        tuple(element.dof_ordering),
        dtype,
    )

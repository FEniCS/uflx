# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test domains."""

import pytest
from conftest import (
    BlockedIdentityMappedElement,
    Interval,
    NonIdentityMappedElement,
    Tetrahedron,
    Triangle,
)

from uflx import composed_domain, parametrized_domain
from uflx.domains import (
    RD,
    AbstractCellularDomain,
    AbstractCoordinateDomain,
    EntityDomain,
    entity_domain,
    labelled_domain,
    union_domain,
)
from uflx.expressions import RealScalar
from uflx.geometry import VolumeElement
from uflx.parametrizations import AffineParametrization
from uflx.points import Point


def test_has_affine_parametrization_true_for_degree_one_simplex(lagrange_element):
    """A degree 1 Lagrange coordinate element on a simplex cell is an affine map."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    assert domain.has_affine_parametrization

    domain = parametrized_domain(lagrange_element("tetrahedron", 1, (3,)))
    assert domain.has_affine_parametrization


def test_has_affine_parametrization_false_for_higher_degree(lagrange_element):
    """A higher degree Lagrange coordinate element is not an affine map."""
    domain = parametrized_domain(lagrange_element("triangle", 2, (2,)))
    assert not domain.has_affine_parametrization


def test_has_affine_parametrization_false_for_tensor_product_cell(lagrange_element):
    """A degree 1 Lagrange coordinate element on a non-simplex cell is not affine.

    Even at degree 1, a quadrilateral or hexahedron's coordinate map is
    multilinear, not affine, since these cells aren't simplices.
    """
    domain = parametrized_domain(lagrange_element("quadrilateral", 1, (2,)))
    assert not domain.has_affine_parametrization

    domain = parametrized_domain(lagrange_element("hexahedron", 1, (3,)))
    assert not domain.has_affine_parametrization


@pytest.mark.parametrize("cell", ["interval", "triangle", "quadrilateral", "tetrahedron"])
def test_entity_domain_dimensions_agree(cell, lagrange_element):
    """An entity domain's geometry is the identity, so its dimensions agree."""
    (entity,) = parametrized_domain(lagrange_element(cell, 1, (3,))).cell_types
    domain = entity_domain(entity)

    assert isinstance(domain, AbstractCoordinateDomain)
    assert domain.geometric_dimension == domain.topological_dimension
    assert domain.geometric_dimension == entity.topological_dimension
    assert domain.cell_types == (entity,)


def test_an_entity_domains_parametrization_is_the_identity(lagrange_element):
    """An entity domain's points are already coordinates, so nothing maps them."""
    (entity,) = parametrized_domain(lagrange_element("triangle", 1, (2,))).cell_types
    domain = entity_domain(entity)

    assert isinstance(domain, AbstractCellularDomain)
    parametrization = domain.parametrization(entity)
    assert parametrization.is_identity
    assert parametrization.is_affine
    assert parametrization.source == domain
    assert parametrization.target_dimension == domain.geometric_dimension
    assert domain.has_affine_parametrization


def test_an_entity_domain_only_parametrizes_its_own_entity(lagrange_element):
    """Asking for another cell's map is a caller mistake."""
    (triangle,) = parametrized_domain(lagrange_element("triangle", 1, (2,))).cell_types
    (interval,) = parametrized_domain(lagrange_element("interval", 1, (1,))).cell_types

    with pytest.raises(ValueError, match="not the entity"):
        entity_domain(triangle).parametrization(interval)


@pytest.mark.parametrize("dim", range(1, 4))
def test_rd_equality(dim):
    """Separately built ambient domains of the same dimension are equal."""
    assert RD(dim) == RD(dim)
    assert hash(RD(dim)) == hash(RD(dim))
    assert RD(dim) != RD(dim + 1)


def test_entity_domain_equality(lagrange_element):
    """Separately built entity domains on the same entity are equal."""
    (triangle,) = parametrized_domain(lagrange_element("triangle", 1, (2,))).cell_types
    (interval,) = parametrized_domain(lagrange_element("interval", 1, (1,))).cell_types

    assert EntityDomain(triangle) == EntityDomain(triangle)
    assert hash(EntityDomain(triangle)) == hash(EntityDomain(triangle))
    assert EntityDomain(triangle) != EntityDomain(interval)


def test_parametrized_domain_equality(lagrange_element):
    """Domains with the same geometric description are equal."""
    a = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    b = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    c = parametrized_domain(lagrange_element("triangle", 2, (2,)))

    assert a == b
    assert hash(a) == hash(b)
    assert a != c


def test_an_entity_domain_is_not_a_parametrized_domain(lagrange_element):
    """The two kinds of domain never compare equal, whatever their dimensions."""
    domain = parametrized_domain(lagrange_element("interval", 1, (1,)))
    (entity,) = domain.cell_types

    assert EntityDomain(entity) != domain
    assert domain != EntityDomain(entity)


def test_parametrization_must_be_identity_mapped(lagrange_element):
    """A parametrization's values are the ambient coordinates, so no mapping is allowed."""
    (triangle,) = parametrized_domain(lagrange_element("triangle", 1, (2,))).cell_types

    with pytest.raises(ValueError, match="identity mapped"):
        parametrized_domain(NonIdentityMappedElement(triangle, 1, (2,)))


def test_parametrization_may_block_an_identity_map(lagrange_element):
    """A vector element blocks a scalar identity map, which still leaves values untouched.

    This is what an element library gives for vector Lagrange, so rejecting it
    would reject every real mesh.
    """
    (triangle,) = parametrized_domain(lagrange_element("triangle", 1, (2,))).cell_types
    domain = parametrized_domain(BlockedIdentityMappedElement(triangle, 1, (2,)))

    assert domain.geometric_dimension == 2
    assert domain.topological_dimension == 2


def in_the_plane(matrix, offset):
    """A domain of the plane, being the reference interval mapped into it."""
    source = entity_domain(Interval())
    return composed_domain(source, AffineParametrization(source, RD(2), matrix, offset))


def test_a_union_takes_the_dimensions_of_its_pieces():
    """Named in the pieces' coordinates, of the pieces' topological dimension."""
    pieces = [
        in_the_plane(((1.0,), (0.0,)), (0.0, 0.0)),
        in_the_plane(((0.0,), (1.0,)), (0.0, 0.0)),
    ]

    union = union_domain(pieces)

    assert union.geometric_dimension == 2
    assert union.topological_dimension == 1
    assert union.pieces == tuple(pieces)


def test_a_union_of_pieces_of_different_dimension_has_no_topological_dimension():
    """As a domain of mixed cell types has none: there is no one answer."""
    union = union_domain([in_the_plane(((1.0,), (0.0,)), (0.0, 0.0)), entity_domain(Triangle())])

    assert union.topological_dimension is None
    assert union.geometric_dimension == 2


def test_a_union_needs_pieces_in_the_same_coordinates():
    """A piece of the plane and a piece of space are not one domain."""
    with pytest.raises(ValueError, match="same coordinates"):
        union_domain([entity_domain(Triangle()), entity_domain(Tetrahedron())])

    with pytest.raises(ValueError, match="union of no pieces"):
        union_domain([])


def test_a_union_is_told_apart_by_its_pieces():
    """It carries no label of its own: a selection is marked on the pieces.

    A label on the union alone would be dropped the moment an integral
    over it were split onto one piece, which is where a consumer needs
    it, so the pieces carry it instead.
    """
    pieces = [in_the_plane(((1.0,), (0.0,)), (0.0, 0.0))]
    marked = [labelled_domain(piece, 3) for piece in pieces]

    assert union_domain(pieces) == union_domain(pieces)
    assert hash(union_domain(pieces)) == hash(union_domain(pieces))
    assert union_domain(marked) != union_domain(pieces)
    assert union_domain(marked) == union_domain([labelled_domain(p, 3) for p in pieces])
    assert len({union_domain(pieces), union_domain(marked), union_domain(marked)}) == 2


def test_a_labelled_domain_differs_only_in_identity():
    """Same cells, same parametrizations, same coordinates, told apart by the label."""
    domain = in_the_plane(((1.0,), (0.0,)), (0.0, 0.0))

    marked = labelled_domain(domain, "inflow")

    assert marked.domain == domain
    assert marked.label == "inflow"
    assert marked.cell_types == domain.cell_types
    assert marked.parametrizations == domain.parametrizations
    assert marked.coordinate_domain == domain.coordinate_domain
    assert marked.geometric_dimension == domain.geometric_dimension
    assert marked.topological_dimension == domain.topological_dimension
    assert marked != domain
    assert marked != labelled_domain(domain, "outflow")
    assert marked.on_cell_type(Interval()) == labelled_domain(
        domain.on_cell_type(Interval()), "inflow"
    )


def test_a_union_holds_every_pieces_parametrization():
    """Several parametrizations out of one region is why a union exists.

    It is not asked to choose between them. Needing exactly one parametrization is
    a question about evaluating a quantity, which geometry asks of the
    point it is read at.
    """
    pieces = [
        in_the_plane(((1.0,), (0.0,)), (0.0, 0.0)),
        in_the_plane(((0.0,), (1.0,)), (0.0, 0.0)),
    ]
    interval = entity_domain(Interval())

    parametrizations = union_domain(pieces).parametrizations

    assert parametrizations == pieces[0].parametrizations + pieces[1].parametrizations
    assert [parametrization.source for parametrization in parametrizations] == [interval, interval]

    with pytest.raises(ValueError, match="2 parametrizations out of"):
        VolumeElement(union_domain(pieces), Point([RealScalar(0.5)], interval)).expand_geometry()

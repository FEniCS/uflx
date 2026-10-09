# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test domains."""

import pytest

from uflx import parametrized_domain
from uflx.domains import RD, AbstractCoordinateDomain, EntityDomain, entity_domain


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


def test_entity_domain_has_no_parametrization(lagrange_element):
    """An entity domain carries no parametrization, unlike a parametrized domain."""
    (entity,) = parametrized_domain(lagrange_element("triangle", 1, (2,))).cell_types

    assert not hasattr(entity_domain(entity), "parametrization")


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

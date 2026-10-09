"""Test the C generator."""

import pytest
from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, inner
from uflx.expressions import RealScalar, Sum

import uflx_codegeneration
from uflx_codegeneration import symbols
from uflx_codegeneration.c import CGenerator
from uflx_codegeneration.nodes import ArrayEntry


def test_core_classes_are_not_modified():
    """Generating C code does not add methods to classes in the UFLx core."""
    assert not hasattr(Sum, "generate_c")
    assert not hasattr(RealScalar, "generate_c")


def test_custom_generator(lagrange_element):
    """A custom generator can be passed to generate."""

    class CountingGenerator(CGenerator):
        def __init__(self):
            self.count = 0

        def generate(self, node):
            self.count += 1
            return super().generate(node)

    space = function_space(
        coordinate_element(lagrange_element("triangle", 1, (2,))), lagrange_element("triangle", 1)
    )
    form = inner(TrialFunction(space), TestFunction(space)) * dx

    symbols.global_variable_namer.reset()
    code, _ = uflx_codegeneration.generate(form)
    generator = CountingGenerator()
    symbols.global_variable_namer.reset()
    custom_code, _ = uflx_codegeneration.generate(form, generator=generator)

    assert generator.count > 0
    assert custom_code == code


def test_object_outside_the_core():
    """An object implementing GenerateC generates its own code."""

    class Twice:
        def __init__(self, operand):
            self.operand = operand

        def generate_c(self, generator: CGenerator) -> str:
            return f"(2 * {generator.generate(self.operand)})"

    assert CGenerator().generate(Twice(ArrayEntry("w", (0,)))) == "(2 * w[0])"


def test_unsupported_node():
    """Generating code for an unsupported node raises an error."""
    with pytest.raises(NotImplementedError):
        CGenerator().generate(object())

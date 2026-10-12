"""Checks that apply to every test."""

import pytest
from uflx.graphs.graphs import _Interning, describe


@pytest.fixture(autouse=True)
def every_value_rebuilds_to_itself(monkeypatch):
    """Check that every value a test builds is rebuilt from its init_args as the same object."""
    build = _Interning.__call__

    def checked(cls, *args, **kwargs):
        value = build(cls, *args, **kwargs)
        assert build(cls, *value.init_args) is value, describe(value)
        return value

    monkeypatch.setattr(_Interning, "__call__", checked)

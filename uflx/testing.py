"""Checks for a test suite that builds UFLx values.

A consumer turns them on for its whole suite with one line in its
top-level ``conftest.py``::

    pytest_plugins = ["uflx.testing"]
"""

import pytest

from uflx.graphs.graphs import _Interning


@pytest.fixture(autouse=True)
def every_value_rebuilds_to_itself(monkeypatch):
    """Check that every value a test builds rebuilds from its init_args to itself.

    Outside a test only the first value of each class is checked. A
    value that does not rebuild to itself raises ``TypeError`` where it
    is built.
    """
    monkeypatch.setattr(_Interning, "check_every_value", True)

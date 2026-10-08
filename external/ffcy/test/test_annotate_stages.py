from collections import Counter

import pytest
from basix import CellType
from xdsl.dialects import tensor
from xdsl.dialects.linalg import ops as linalg

from ffcy.pipeline import apply, compile_form
from ffcy.transforms.annotate_stages import BUFFER, SINK, STAGE, AnnotateStages
from ffcy.transforms.linalg_builders import PENCIL
from ffcy.transforms.split_cells import SplitCells
from forms.laplace import laplace_form
from forms.mass import mass_form


@pytest.mark.parametrize(
    ("forms", "sinks", "fused"),
    [
        # Evaluation along x and y, then evaluation along z, the flux and integration
        # along z in one stage, then integration along y and x. Laplace also sums its
        # terms in the last stage.
        (mass_form, {0: {2: 1}, 1: {3: 1}, 2: {4: 1}, 3: {3: 1}, 4: {2: 1}}, {2}),
        (laplace_form, {0: {2: 3}, 1: {3: 3}, 2: {4: 3}, 3: {3: 3}, 4: {2: 1}}, {2, 4}),
    ],
)
def test_annotate_stages(forms, sinks, fused):
    """Thread-local chains share a stage, and every op writes into its own buffer."""
    module = apply(
        compile_form(forms(2, CellType.hexahedron), 5, True, expand_results=True),
        SplitCells(4),
        AnnotateStages(),
    )
    (action,) = [f for f in module.ops if f.sym_name.data == "action"]

    found: dict[int, Counter[int]] = {}
    members = Counter()
    for op in action.walk():
        if isinstance(op, linalg.GenericOp):
            if SINK in op.attributes:
                stage = op.attributes[SINK].value.data
                found.setdefault(stage, Counter())[op.attributes[PENCIL].value.data] += 1
            else:
                members[op.attributes[STAGE].value.data] += 1
    assert found == {stage: Counter(axes) for stage, axes in sinks.items()}
    assert set(members) == fused

    empties = [op for op in action.walk() if isinstance(op, tensor.EmptyOp)]
    assert all(len(list(op.tensor.uses)) == 1 for op in empties)
    assert len({op.attributes[BUFFER].value.data for op in empties}) == len(empties)

"""Generate the reference-cell figure as SVG and PDF.

Run with Python and Matplotlib; outputs are written beside this script.
"""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def main():
    """Draw the two local reference cells and their selected edge maps."""
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "svg.fonttype": "path",
            "svg.hashsalt": "uflx-reference-slots",
            "pdf.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    fig.patch.set_facecolor("white")
    for slot, (ax, labels, colour) in enumerate(
        zip(axes, ((0, 1, 2), (1, 3, 2)), ("#2563a6", "#b45309"), strict=True)
    ):
        ax.set_aspect("equal")
        ax.set_xlim(-0.28, 1.35)
        ax.set_ylim(-0.3, 1.35)
        ax.axis("off")
        ax.add_patch(
            Polygon(
                [(0, 0), (1, 0), (0, 1)], facecolor="#f3f6fa", edgecolor="#64748b", linewidth=1.5
            )
        )
        start = (1, 0) if slot == 0 else (0, 0)
        ax.plot([start[0], 0], [start[1], 1], color=colour, linewidth=3)
        # Arrows show increasing interval coordinate, not normal orientation.
        arrow_start = (0.7, 0.3) if slot == 0 else (0, 0.3)
        arrow_end = (0.35, 0.65) if slot == 0 else (0, 0.65)
        ax.annotate(
            "",
            xy=arrow_end,
            xytext=arrow_start,
            arrowprops={"arrowstyle": "-|>", "color": colour, "lw": 2.5, "mutation_scale": 17},
        )
        for (x, y), label, offset, alignment in zip(
            ((0, 0), (1, 0), (0, 1)),
            labels,
            ((0, -13), (0, -13), (0, 12)),
            ("top", "top", "bottom"),
            strict=True,
        ):
            ax.plot(x, y, "o", color="#1e293b", markersize=4.5)
            ax.annotate(
                str(label),
                (x, y),
                xytext=offset,
                textcoords="offset points",
                ha="center",
                va=alignment,
                weight="bold",
                color="#1e293b",
            )
        ax.set_title(
            "Cell A" if slot == 0 else "Cell B", loc="left", weight="bold", fontsize=13, pad=8
        )
        label = "inclusion_a" if slot == 0 else "inclusion_b"
        ax.text(0.42, 0.78, label, color=colour, fontsize=11, fontfamily="DejaVu Sans Mono")
        ax.text(
            0.5,
            -0.27,
            r"Reference coordinates $(\xi,\eta)$",
            ha="center",
            color="#475569",
            fontsize=10,
        )
    fig.subplots_adjust(left=0.04, right=0.97, bottom=0.06, top=0.86, wspace=0.18)
    output = Path(__file__).parent
    for extension in ("svg", "pdf"):
        metadata = {"Date": None} if extension == "svg" else {"CreationDate": None, "ModDate": None}
        fig.savefig(output / f"reference-slots.{extension}", metadata=metadata)
    plt.close(fig)


if __name__ == "__main__":
    main()

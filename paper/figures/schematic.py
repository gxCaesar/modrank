"""Drawing primitives for schematic panels, so a flowchart is built rather than drawn by hand.

Three reasons this is a module and not inline code in one figure script.

WHY PRIMITIVES. A flowchart assembled from raw `ax.add_patch` calls acquires its geometry from
whatever numbers were typed at each call site, and the first time a box is widened every arrow
attached to it has to be found and moved by hand. Here a box knows its own anchors and an arrow
asks for them, so moving a box moves its arrows.

WHY NO TYPED NUMERALS. Every quantity printed on a schematic in this manuscript is loaded from the
results JSON it belongs to. A number typed into a figure script is a number no gate reads, and this
project has already paid twice for a value that lived in exactly one place.

WHY THE COLLISION CHECK. A schematic degrades silently: text that overflows its box still renders,
and looks deliberate. `Box.fit_fontsize` shrinks to fit and RAISES below a floor rather than
emitting something unreadable, because the failure this guards against is a 5 pt label nobody
noticed until a referee did.
"""

from __future__ import annotations

import matplotlib.patches as mpatches
import numpy as np

from style import INK

MIN_PT = 6.0  # nothing in a figure of this manuscript renders below this


class Box:
    """A rounded rectangle in AXES coordinates that knows its own edge anchors."""

    def __init__(self, ax, x, y, w, h, facecolor="none", edgecolor=INK, lw=0.8,
                 radius=0.012, ls="-", zorder=2, alpha=1.0):
        self.ax, self.x, self.y, self.w, self.h = ax, x, y, w, h
        self.patch = mpatches.FancyBboxPatch(
            (x, y), w, h, boxstyle=mpatches.BoxStyle("Round", pad=0, rounding_size=radius),
            facecolor=facecolor, edgecolor=edgecolor, linewidth=lw, linestyle=ls,
            zorder=zorder, alpha=alpha, transform=ax.transAxes, clip_on=False)
        ax.add_patch(self.patch)

    # -- anchors, so arrows never carry their own coordinates
    @property
    def left(self):
        return (self.x, self.y + self.h / 2)

    @property
    def right(self):
        return (self.x + self.w, self.y + self.h / 2)

    @property
    def top(self):
        return (self.x + self.w / 2, self.y + self.h)

    @property
    def bottom(self):
        return (self.x + self.w / 2, self.y)

    @property
    def centre(self):
        return (self.x + self.w / 2, self.y + self.h / 2)

    def text(self, s, size=7, weight="normal", dy=0.0, colour=INK, ha="center", va="center",
             x=None, linespacing=1.25):
        cx, cy = self.centre
        self.ax.text(x if x is not None else cx, cy + dy, s, size=size, weight=weight, color=colour,
                     ha=ha, va=va, transform=self.ax.transAxes, zorder=4,
                     linespacing=linespacing)

    def fit_fontsize(self, s, start=7.0, chars_per_inch=None, ax_width_in=None):
        """Largest size at or below `start` whose longest line fits, never below MIN_PT.

        Deliberately crude -- a character-width estimate, not a render-and-measure loop -- because
        the point is to CATCH the overflow, and a rule that raises is worth more than one that is
        exact. If it cannot fit at MIN_PT it raises, so the figure fails loudly instead of
        shipping something unreadable.
        """
        longest = max(len(line) for line in s.split("\n"))
        if chars_per_inch is None:
            chars_per_inch = 1.0 / 0.0062          # ~ Arial average advance at 1 pt, in inches
        width_in = self.w * (ax_width_in if ax_width_in is not None else 6.30)
        size = min(start, width_in * chars_per_inch / max(longest, 1))
        if size < MIN_PT:
            raise ValueError(
                "text does not fit at the %.1f pt floor: %d chars in %.2f in -> %.2f pt\n  %r"
                % (MIN_PT, longest, width_in, size, s.split("\n")[0]))
        return size


def arrow(ax, p0, p1, colour=INK, lw=0.9, ls="-", rad=0.0, head=2.2, zorder=3, alpha=1.0):
    ax.annotate("", xy=p1, xytext=p0, xycoords=ax.transAxes, textcoords=ax.transAxes,
                zorder=zorder,
                arrowprops=dict(arrowstyle="-|>", color=colour, linewidth=lw, linestyle=ls,
                                alpha=alpha, mutation_scale=head * 3.2,
                                connectionstyle="arc3,rad=%.3f" % rad,
                                shrinkA=1.5, shrinkB=1.5))


def label(ax, x, y, s, size=7, weight="normal", colour=INK, ha="center", va="center",
          style="normal", linespacing=1.25, zorder=5):
    return ax.text(x, y, s, size=size, weight=weight, color=colour, ha=ha, va=va, style=style,
                   transform=ax.transAxes, zorder=zorder, linespacing=linespacing)


def panel_letter(ax, s, x=-0.005, y=1.0):
    ax.text(x, y, s, size=9, weight="bold", color=INK, ha="left", va="top",
            transform=ax.transAxes, zorder=6)


def blank(ax):
    """A drawing surface: no axes, no ticks, unit coordinates, nothing clipped."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    return ax


def bracket(ax, x0, x1, y, s, size=6.5, colour=INK, drop=0.035, lw=0.7, above=True):
    """A horizontal brace with a label, for `these three things together are X`."""
    d = drop if above else -drop
    ax.plot([x0, x0, x1, x1], [y - d * 0.45, y, y, y - d * 0.45], color=colour, lw=lw,
            transform=ax.transAxes, clip_on=False, zorder=3)
    ax.text((x0 + x1) / 2, y + (0.012 if above else -0.028), s, size=size, color=colour,
            ha="center", va="bottom" if above else "top", transform=ax.transAxes, zorder=5)


def assert_no_overlap(boxes, tol=0.004):
    """Boxes on one panel must not overlap. A schematic whose boxes collide still renders."""
    for i, a in enumerate(boxes):
        for b in boxes[i + 1:]:
            dx = min(a.x + a.w, b.x + b.w) - max(a.x, b.x)
            dy = min(a.y + a.h, b.y + b.h) - max(a.y, b.y)
            if dx > tol and dy > tol:
                raise ValueError("schematic boxes overlap by %.3f x %.3f in axes units" % (dx, dy))

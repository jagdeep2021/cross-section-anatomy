"""Review and correct the auto pass, image by image, in napari.

  python review.py <anatomy_dir>          (the folder run_batch.py wrote)

Clicks (in the default pan/zoom mode):
  Shift-click        add a cell (fills the dark lumen out to its walls)
  Ctrl-click         delete the cell under the cursor
  Alt-click          add a metaxylem vessel here; Alt-click a vessel to remove it
Paint tools (napari's own, top-left once a layer is selected):
  region layer       paint 2 = stele, 1 = root, 0 = outside; fill bucket for big areas
  cells layer        brush / eraser to fix a boundary, picker to see an id
Keys:  Shift-R recompute layers   Shift-S save + next   . next   , previous
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import napari
import numpy as np
from magicgui.widgets import ComboBox, Container, Label, PushButton
from napari.utils.colormaps import DirectLabelColormap

import pipeline as P

# C1..C12, inner to outer; repeats past 12
_PAL = ["#3b4cc0", "#1f9e89", "#6ece58", "#fde725", "#f89540", "#e16462",
        "#b12a90", "#7e03a8", "#2c7fb8", "#7fcdbb", "#d95f0e", "#636363"]
LAYER_COLORS = {k: _PAL[(k - 1) % len(_PAL)] for k in range(1, 41)}

HELP = ("Shift-click: add cell   Ctrl-click: delete cell   Alt-click: add/remove metaxylem\n"
        "region layer: paint 2 = stele, 1 = root   Shift-R: recompute   Shift-S: save + next")


class Reviewer:
    def __init__(self, anatomy_dir: Path, show: bool = True):
        self.dir = Path(anatomy_dir)
        self.settings = json.loads((self.dir / "settings.json").read_text())
        self.image_dir = Path(self.settings["image_dir"])
        self.angle = float(self.settings.get("quadrant_angle_deg", 0))
        self.files = sorted((self.dir / "masks").glob("*.npz"))
        if not self.files:
            sys.exit(f"No masks in {self.dir / 'masks'}; run run_batch.py first.")
        self.viewer = napari.Viewer(title="Root cross-section review", show=show)
        self.i = None
        self._build_dock()
        first = next((k for k, f in enumerate(self.files) if P.load_state(f)[3] != "reviewed"), 0)
        self.open(first)

    # ---------------------------------------------------------- UI
    def _build_dock(self):
        self.pick = ComboBox(choices=self._choices(), label="image")
        self.pick.changed.connect(lambda v: self.open(self._choices().index(v)))
        self.info = Label(value="")
        b = {n: PushButton(text=n) for n in ("< Prev", "Next >", "Next unreviewed", "Recompute layers", "Save + next", "Reset to auto")}
        b["< Prev"].changed.connect(lambda: self.open(max(self.i - 1, 0)))
        b["Next >"].changed.connect(lambda: self.open(min(self.i + 1, len(self.files) - 1)))
        b["Next unreviewed"].changed.connect(self.next_unreviewed)
        b["Recompute layers"].changed.connect(self.recompute)
        b["Save + next"].changed.connect(self.save_next)
        b["Reset to auto"].changed.connect(self.reset)
        self.dock = Container(widgets=[self.pick, *b.values(), self.info, Label(value=HELP)])
        self.viewer.window.add_dock_widget(self.dock, area="right", name="Review")
        v = self.viewer
        v.bind_key("Shift-R", lambda _: self.recompute(), overwrite=True)
        v.bind_key("Shift-S", lambda _: self.save_next(), overwrite=True)
        v.bind_key(".", lambda _: self.open(min(self.i + 1, len(self.files) - 1)), overwrite=True)
        v.bind_key(",", lambda _: self.open(max(self.i - 1, 0)), overwrite=True)
        v.mouse_drag_callbacks.append(self._click)

    def _choices(self):
        return [f.stem + ("  ✓" if P.load_state(f)[3] == "reviewed" else "") for f in self.files]

    # ---------------------------------------------------------- image I/O
    def open(self, i: int):
        if i == self.i:
            return
        self.i = i
        f = self.files[i]
        z = np.load(f)
        self.rel = str(z["image_file"])
        self.bar = int(z["scale_bar_px"]) if "scale_bar_px" in z else 0
        self.gray = P.load_gray(self.image_dir / self.rel)
        cells, region, mx, self.status = P.load_state(f)
        v = self.viewer
        v.layers.clear()
        v.add_image(self.gray, name="image")
        self.region = v.add_labels(region.astype(np.uint8), name="region (2 stele, 1 root)", opacity=1,
                                   colormap=DirectLabelColormap(color_dict={None: "transparent", 0: "transparent",
                                                                            1: (1, 1, 0, 0.06), 2: (0, 0.9, 1, 0.3)}))
        self.layers = v.add_labels(np.zeros_like(cells), name="cortical layers", opacity=0.5,
                                   colormap=DirectLabelColormap(color_dict={None: "gray", 0: "transparent", **LAYER_COLORS}))
        self.mx = v.add_labels(mx, name="metaxylem", opacity=0.6,
                               colormap=DirectLabelColormap(color_dict={None: "red", 0: "transparent"}))
        self.cells = v.add_labels(cells, name="cells", opacity=0.9,
                                  colormap=DirectLabelColormap(color_dict={None: "yellow", 0: "transparent"}))
        self.cells.contour = 1
        self.air = v.add_labels(np.zeros_like(cells), name="air spaces", opacity=0.8,
                                colormap=DirectLabelColormap(color_dict={None: "magenta", 0: "transparent"}))
        self.layers.editable = False
        self.air.editable = False
        self.quads = v.add_shapes(name="quadrants", edge_color="white", edge_width=1.5, face_color="transparent")
        v.layers.selection.active = self.cells
        self.recompute()
        with self.pick.changed.blocked():
            self.pick.choices = self._choices()
            self.pick.value = self._choices()[i]

    def save(self):
        P.save_state(self.files[self.i], self.cells.data.astype(np.uint16), self.region.data.astype(np.uint8),
                     self.mx.data.astype(np.uint16), status="reviewed", image_file=self.rel, scale_bar_px=self.bar)
        self.status = "reviewed"

    def save_next(self):
        self.save()
        self.next_unreviewed()

    def next_unreviewed(self):
        n = len(self.files)
        for k in range(1, n + 1):
            j = (self.i + k) % n
            if P.load_state(self.files[j])[3] != "reviewed":
                self.open(j)
                return
        self.info.value = "All images reviewed. Run summarise.py"

    def reset(self):
        gray, cells, region, mx = P.propose(self.image_dir / self.rel, self.settings.get("diameter"))
        self.cells.data, self.region.data, self.mx.data = cells, region, mx
        self.recompute()

    # ---------------------------------------------------------- edits
    def recompute(self):
        props, _, tr = P.derive(self.cells.data, self.region.data, self.mx.data, self.angle)
        self.layers.data = P.layer_image(self.cells.data, props)
        cy, cx = P.centre_of(self.region.data)
        h, w = self.gray.shape
        a = np.radians(self.angle)
        L = 0.5 * max(h, w)
        ends = [np.array([[cy - L * np.sin(a + t), cx + L * np.cos(a + t)], [cy, cx]]) for t in (0, np.pi / 2, np.pi, 3 * np.pi / 2)]
        self.quads.data = []
        self.quads.add_lines(ends, edge_color="white", edge_width=1.5)
        qt = P.quadrant_traits(self.cells.data, self.region.data, self.mx.data, props, self.angle)
        self.air.data = np.isin(self.cells.data, props.label[props.tissue == "air_space"].values).astype(np.uint8)
        self.info.value = (f"{self.i + 1}/{len(self.files)}  {self.rel}  [{self.status}]\n"
                           f"layers {tr['n_cortical_layers']}   metaxylem {tr['mx_n']}   cortex cells {tr['n_cortex_cells']}   air spaces {tr['n_air_spaces']}\n"
                           f"layers by quadrant Q1-Q4: {'/'.join(str(n) for n in qt.n_cortical_layers)}")

    def _click(self, viewer, event):
        mods = set(event.modifiers)
        if not mods & {"Shift", "Control", "Alt"} or event.button != 1:
            return
        y, x = (int(round(c)) for c in self.cells.world_to_data(event.position)[-2:])
        if not (0 <= y < self.gray.shape[0] and 0 <= x < self.gray.shape[1]):
            return
        self.click(y, x, mods)

    def click(self, y, x, mods):
        cells, mx = self.cells.data, self.mx.data
        if "Shift" in mods:
            P.add_cell(self.gray, cells, y, x)
        elif "Control" in mods:
            if cells[y, x]:
                cells[cells == cells[y, x]] = 0
        elif "Alt" in mods:
            if mx[y, x]:
                mx[mx == mx[y, x]] = 0
            else:
                P.add_mx(self.gray, mx, cells, y, x)
                self.mx.data = mx
        self.cells.data = cells
        self.mx.refresh()
        self.recompute()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    Reviewer(Path(sys.argv[1]))
    napari.run()

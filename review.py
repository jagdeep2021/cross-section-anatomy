"""Review and correct the auto pass, image by image, in napari.

  python review.py <anatomy_dir>          (the folder run_batch.py wrote)

Two click modes (dropdown in the Review panel, or Shift-M to switch):

  Fix cells
    Shift-click      add a cell (fills the dark lumen out to its walls)
    Ctrl-click       delete the cell under the cursor
    Alt-click        add a metaxylem vessel here; Alt-click a vessel to remove it
    region layer     paint 2 = stele, 1 = root, 0 = outside (napari brush / fill bucket)

  Draw cell (button or Shift-D)   trace a missing cell freehand; it takes any overlapping
                   pixels (each pixel belongs to one cell, so nothing is counted twice).
                   In File paths mode the drawn cell also joins the path.
  Draw metaxylem (button or Shift-X)   same, for a vessel, but a red vessel it touches is merged
                   into it (one vessel, union area); it also clears cell pixels under it.

  File paths  (one spoke per quadrant; its cells, stele -> out, are files 1..n)
    click            add / remove the cell from its quadrant's path (zig-zag is fine)
    Shift-click      move this quadrant's spoke through the click (path is re-proposed)

Keys:  Shift-M switch mode   Shift-R recompute   Shift-S save + next   . next   , previous
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import napari
import numpy as np
from magicgui.widgets import ComboBox, Container, Label, PushButton
from napari.utils.colormaps import DirectLabelColormap
from scipy import ndimage as ndi

import pipeline as P

_PAL = ["#3b4cc0", "#1f9e89", "#6ece58", "#fde725", "#f89540", "#e16462",
        "#b12a90", "#7e03a8", "#2c7fb8", "#7fcdbb", "#d95f0e", "#636363"]
RING_COLORS = {k: _PAL[(k - 1) % len(_PAL)] for k in range(1, 41)}
MODES = ("Fix cells", "File paths")

HELP = {"Fix cells": "Shift-click: add cell   Ctrl-click: delete cell   Alt-click: add/remove metaxylem\n"
                     "region layer: paint 2 = stele, 1 = root   Shift-D: draw a cell   Shift-X: draw a metaxylem vessel",
        "File paths": "click a cell: add/remove it from the path (files numbered from the stele out)\n"
                      "Shift-click: move this quadrant's spoke here   Shift-D: draw a missing cell (joins the path)"}


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
        self.mode = MODES[0]
        self._build_dock()
        first = next((k for k, f in enumerate(self.files) if P.load_state(f)[3] != "reviewed"), 0)
        self.open(first)

    # ---------------------------------------------------------- UI
    def _build_dock(self):
        self.pick = ComboBox(choices=self._choices(), label="image")
        self.pick.changed.connect(lambda v: self.open(self._choices().index(v)))
        self.mode_box = ComboBox(choices=MODES, value=self.mode, label="click mode")
        self.mode_box.changed.connect(self.set_mode)
        self.info = Label(value="")
        self.help = Label(value=HELP[self.mode])
        b = {n: PushButton(text=n) for n in ("Draw cell", "Draw metaxylem", "< Prev", "Next >", "Next unreviewed", "Recompute", "Save + next",
                                             "Re-propose file paths", "Reset image to auto")}
        b["Draw cell"].changed.connect(lambda: self.start_drawing("cell"))
        b["Draw metaxylem"].changed.connect(lambda: self.start_drawing("mx"))
        b["< Prev"].changed.connect(lambda: self.open(max(self.i - 1, 0)))
        b["Next >"].changed.connect(lambda: self.open(min(self.i + 1, len(self.files) - 1)))
        b["Next unreviewed"].changed.connect(self.next_unreviewed)
        b["Recompute"].changed.connect(self.recompute)
        b["Save + next"].changed.connect(self.save_next)
        b["Re-propose file paths"].changed.connect(self.repropose_paths)
        b["Reset image to auto"].changed.connect(self.reset)
        self.dock = Container(widgets=[self.pick, self.mode_box, *b.values(), self.info, self.help])
        self.viewer.window.add_dock_widget(self.dock, area="right", name="Review")
        v = self.viewer
        v.bind_key("Shift-M", lambda _: self.set_mode(MODES[1 - MODES.index(self.mode)]), overwrite=True)
        v.bind_key("Shift-R", lambda _: self.recompute(), overwrite=True)
        v.bind_key("Shift-D", lambda _: self.start_drawing("cell"), overwrite=True)
        v.bind_key("Shift-X", lambda _: self.start_drawing("mx"), overwrite=True)
        v.bind_key("Shift-S", lambda _: self.save_next(), overwrite=True)
        v.bind_key(".", lambda _: self.open(min(self.i + 1, len(self.files) - 1)), overwrite=True)
        v.bind_key(",", lambda _: self.open(max(self.i - 1, 0)), overwrite=True)
        v.mouse_drag_callbacks.append(self._click)

    def set_mode(self, mode):
        self.mode = mode
        self.help.value = HELP[mode] + "\n\nKeys:  Shift-M switch mode   Shift-D draw cell   Shift-X draw metaxylem   Shift-R redraw   Shift-S save + next   . next   , previous"
        with self.mode_box.changed.blocked():
            self.mode_box.value = mode
        self.rings.visible = mode == "Fix cells"

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
        self.rings = v.add_labels(np.zeros_like(cells), name="rings (gap check)", opacity=0.4,
                                  colormap=DirectLabelColormap(color_dict={None: "gray", 0: "transparent", **RING_COLORS}))
        self.mx = v.add_labels(mx, name="metaxylem", opacity=0.6,
                               colormap=DirectLabelColormap(color_dict={None: "red", 0: "transparent"}))
        self.path = v.add_labels(np.zeros_like(cells), name="file path", opacity=0.55,
                                 colormap=DirectLabelColormap(color_dict={None: "orangered", 0: "transparent"}))
        self.cells = v.add_labels(cells, name="cells", opacity=0.9,
                                  colormap=DirectLabelColormap(color_dict={None: "yellow", 0: "transparent"}))
        self.cells.contour = 1
        self.air = v.add_labels(np.zeros_like(cells), name="air spaces", opacity=0.8,
                                colormap=DirectLabelColormap(color_dict={None: "magenta", 0: "transparent"}))
        self.quads = v.add_shapes(name="quadrants + spokes", face_color="transparent")
        self.draw = v.add_shapes(name="draw a cell", face_color=[1, 1, 1, 0.2], edge_color="white", edge_width=1)
        self.draw.events.data.connect(self._on_draw)
        self.nums = v.add_points(np.empty((0, 2)), name="file numbers", size=0,
                                 text={"string": "{file}", "color": "white", "size": 9})
        for lyr in (self.rings, self.air, self.path):
            lyr.editable = False
        v.layers.selection.active = self.cells
        props, _, _ = P.derive(cells, region, mx, self.angle)
        saved = P.load_paths(f)
        if saved is None:
            self.spokes, self.path_q, self.path_cells = P.auto_paths(cells, region, props, self.angle)
        else:
            self.spokes, self.path_q, self.path_cells = (np.array(a) for a in saved)
        self.set_mode(self.mode)
        self.recompute()
        with self.pick.changed.blocked():
            self.pick.choices = self._choices()
            self.pick.value = self._choices()[i]

    def save(self):
        P.save_state(self.files[self.i], self.cells.data.astype(np.uint16), self.region.data.astype(np.uint8),
                     self.mx.data.astype(np.uint16), status="reviewed", image_file=self.rel, scale_bar_px=self.bar,
                     spokes=self.spokes, path_q=np.asarray(self.path_q, np.int16),
                     path_cells=np.asarray(self.path_cells, np.int32))
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
        self.repropose_paths()

    def repropose_paths(self):
        props, _, _ = P.derive(self.cells.data, self.region.data, self.mx.data, self.angle)
        self.spokes, self.path_q, self.path_cells = P.auto_paths(self.cells.data, self.region.data, props, self.angle)
        self.recompute()

    # ---------------------------------------------------------- display
    def recompute(self):
        cells, region = self.cells.data, self.region.data
        props, _, tr = P.derive(cells, region, self.mx.data, self.angle)
        self.props = props
        self.rings.data = P.layer_image(cells, props)
        self.air.data = np.isin(cells, props.label[props.tissue == "air_space"].values).astype(np.uint8)
        ft = P.file_table(cells, region, props, self.path_q, self.path_cells)
        self.path.data = np.isin(cells, ft.label.values).astype(np.uint8) if len(ft) else np.zeros_like(cells, np.uint8)
        pr = props.set_index("label")
        if len(ft):
            self.nums.data = np.c_[pr.loc[ft.label, "y"].values, pr.loc[ft.label, "x"].values]
            self.nums.features = {"file": ft.file.values.astype(int)}
            self.nums.text = {"string": "{file}", "color": "white", "size": 9}
        else:
            self.nums.data = np.empty((0, 2))

        cy, cx = P.centre_of(region)
        L = 0.5 * max(self.gray.shape)
        line = lambda t, c: np.array([[cy, cx], [cy - L * np.sin(np.radians(t)), cx + L * np.cos(np.radians(t))]])
        self.quads.data = []
        self.quads.add_lines([line(self.angle + 90 * k, 0) for k in range(4)], edge_color="white", edge_width=1.2)
        self.quads.add_lines([line(t, 0) for t in self.spokes], edge_color="orange", edge_width=1.5)

        per_q = ft.groupby("quadrant").file.max().reindex([1, 2, 3, 4]).fillna(0).astype(int) if len(ft) else [0] * 4
        self.info.value = (f"{self.i + 1}/{len(self.files)}  {self.rel}  [{self.status}]\n"
                           f"files along the path, Q1-Q4:  {'/'.join(str(n) for n in per_q)}\n"
                           f"metaxylem {tr['mx_n']}   cortex cells {tr['n_cortex_cells']}   air spaces {tr['n_air_spaces']}")

    # ---------------------------------------------------------- clicks
    def _pos(self, event):
        y, x = (int(round(c)) for c in self.cells.world_to_data(event.position)[-2:])
        return (y, x) if 0 <= y < self.gray.shape[0] and 0 <= x < self.gray.shape[1] else None

    def _click(self, viewer, event):
        if event.button != 1:
            return
        if self.viewer.layers.selection.active is getattr(self, "draw", None):
            return   # clicks are drawing a cell
        mods = set(event.modifiers)
        pos = self._pos(event)
        if pos is None:
            return
        if self.mode == "Fix cells":
            if mods & {"Shift", "Control", "Alt"}:
                self.click(*pos, mods)
            return
        if mods & {"Control", "Alt"}:
            return
        dragged = False   # a plain click in path mode must not fire when the user pans
        yield
        while event.type == "mouse_move":
            dragged = True
            yield
        if not dragged:
            self.path_click(*pos, mods)

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

    def start_drawing(self, target="cell"):
        """Freehand outline -> new cell or metaxylem vessel. Hold the mouse and trace; release
        closes the shape. Pressing the same button/key again stops."""
        drawing = self.viewer.layers.selection.active is self.draw and self.draw.mode != "pan_zoom"
        if drawing and getattr(self, "draw_target", None) == target:
            return self.stop_drawing()
        self.draw_target = target
        self.viewer.layers.selection.active = self.draw
        self.draw.mode = "add_polygon_lasso"
        what = "a cell" if target == "cell" else "a METAXYLEM vessel"
        key = "Shift-D" if target == "cell" else "Shift-X"
        self.help.value = (f"DRAWING {what}: hold the mouse and trace; release to finish. Draw as many as needed.\n"
                           f"{key} again to stop.")

    def stop_drawing(self):
        self.draw.mode = "pan_zoom"
        self.viewer.layers.selection.active = self.cells
        self.set_mode(self.mode)

    def _on_draw(self, event):
        if getattr(event, "action", None) != "added" or not len(self.draw.data):
            return
        masks = self.draw.to_masks(mask_shape=self.cells.data.shape)
        cells, mx = self.cells.data, self.mx.data
        target = cells if getattr(self, "draw_target", "cell") == "cell" else mx
        centre = P.centre_of(self.region.data)
        for m in masks:
            if m.sum() < 10:
                continue
            new = int(target.max()) + 1
            hit = {int(l): int((target == l).sum()) for l in np.unique(target[m]) if l}
            target[m] = new   # one owner per pixel: the drawn object takes any overlap
            for l, before in hit.items():
                if target is mx:   # vessels never touch, so an overlapped vessel is the same vessel: merge
                    target[target == l] = new
                elif (target == l).sum() < 0.25 * before:   # a cell that lost most of itself is gone, not a sliver
                    target[target == l] = 0
            if target is mx:
                cells[m] = 0   # a vessel is not a cell
            elif self.mode == "File paths":
                cy, cx = np.argwhere(m).mean(0)
                self.path_q = np.r_[self.path_q, int(P.quadrant_of(cy, cx, centre, self.angle))].astype(np.int16)
                self.path_cells = np.r_[self.path_cells, new].astype(np.int32)
        self.mx.data = mx
        self.cells.data = cells
        with self.draw.events.data.blocker():
            self.draw.data = []
        self.recompute()

    def path_click(self, y, x, mods):
        centre = P.centre_of(self.region.data)
        if "Shift" in mods:  # move this quadrant's spoke through the click
            q = int(P.quadrant_of(y, x, centre, self.angle))
            t = float(np.degrees(np.arctan2(-(y - centre[0]), x - centre[1])))
            self.spokes[q - 1] = t
            keep = self.path_q != q
            new = P.path_along(self.cells.data, self.region.data, self.props, t)
            self.path_q = np.r_[self.path_q[keep], np.full(len(new), q)].astype(np.int16)
            self.path_cells = np.r_[self.path_cells[keep], new].astype(np.int32)
        else:
            lbl = int(self.cells.data[y, x])
            if not lbl:
                return
            if lbl in self.path_cells:
                keep = self.path_cells != lbl
                self.path_q, self.path_cells = self.path_q[keep], self.path_cells[keep]
            else:
                cyx = ndi.center_of_mass(self.cells.data == lbl)
                q = int(P.quadrant_of(cyx[0], cyx[1], centre, self.angle))
                self.path_q = np.r_[self.path_q, q].astype(np.int16)
                self.path_cells = np.r_[self.path_cells, lbl].astype(np.int32)
        self.recompute()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    Reviewer(Path(sys.argv[1]))
    napari.run()

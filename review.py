"""Review and correct the auto pass, image by image, in napari.

  python review.py <anatomy_dir>          (the folder run_batch.py wrote)

Pick what you are editing (keys 1-4; its name follows the cursor):
  1 Cell   2 Metaxylem   3 Stele   4 File path

The same four gestures work on every target:
  click              add here      cell / vessel: fills the lumen out to its walls
                                   stele: the cell under the cursor becomes stele
                                   file path: the cell joins its quarter's path
  Cmd-click          remove what is under the cursor        (Ctrl-click on Windows/Linux)
  Shift-drag         draw an outline -> add it               (cell: new cell; vessel: merges with
                                   vessels it touches; stele: area becomes stele; path: cells inside join)
  Cmd-Shift-drag     draw an outline -> remove everything inside it
  drag               pan   ·   scroll  zoom
  Alt-click          (file path) move this quarter's spoke through the click

One owner per pixel: a drawn cell takes the pixels it overlaps (a neighbour left with < 25% of
its area is removed), so nothing is counted twice. A vessel clears the cell pixels under it.

Keys:  1-4 target   Shift-H hide overlays   Shift-S save + next   , previous   . next
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import napari
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from magicgui.widgets import CheckBox, ComboBox, Container, FloatSlider, Label, PushButton
from napari.utils.colormaps import DirectLabelColormap
from scipy import ndimage as ndi
from skimage.draw import polygon as fill_polygon

import pattern_plots as PP
import pipeline as P

_PAL = ["#3b4cc0", "#1f9e89", "#6ece58", "#fde725", "#f89540", "#e16462",
        "#b12a90", "#7e03a8", "#2c7fb8", "#7fcdbb", "#d95f0e", "#636363"]
RING_COLORS = {k: _PAL[(k - 1) % len(_PAL)] for k in range(1, 41)}
TARGETS = {"Cell": "yellow", "Metaxylem": "red", "Stele": "cyan", "File path": "orange"}
HELP = ("click  add here\n⌘-click  remove\nShift-drag  draw outline → add\n⌘-Shift-drag  draw outline → remove\n"
        "drag  pan · scroll  zoom\nFile path: Alt-click moves the spoke\n\n"
        "1-4 target · Shift-H hide overlays\nShift-S save + next · , / . prev / next")


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
        self.target = "Cell"
        self._build_dock()
        first = next((k for k, f in enumerate(self.files) if P.load_state(f)[3] != "reviewed"), 0)
        self.open(first)

    # ---------------------------------------------------------- UI
    def _build_dock(self):
        self.pick = ComboBox(choices=self._choices(), label="image")
        self.pick.changed.connect(lambda v: self.open(self._choices().index(v)))
        self.target_box = ComboBox(choices=list(TARGETS), value=self.target, label="editing")
        self.target_box.changed.connect(self.set_target)
        self.img_op = FloatSlider(value=1.0, min=0.1, max=1.0, step=0.05, label="image")
        self.img_op.changed.connect(lambda v: setattr(self.viewer.layers["image"], "opacity", v))
        self.show_rings = CheckBox(value=False, text="show rings (gap check)")
        self.show_rings.changed.connect(lambda v: setattr(self.rings, "visible", v))
        prev, nxt, save = PushButton(text="< Prev"), PushButton(text="Next >"), PushButton(text="Save + next")
        prev.changed.connect(lambda: self.open(max(self.i - 1, 0)))
        nxt.changed.connect(lambda: self.open(min(self.i + 1, len(self.files) - 1)))
        save.changed.connect(self.save_next)
        self.info = Label(value="")
        self.dock = Container(widgets=[self.pick, self.target_box,
                                       Container(widgets=[prev, nxt], layout="horizontal", labels=False),
                                       save, self.img_op, self.show_rings, self.info, Label(value=HELP)])
        for w in self.dock:
            if isinstance(w, Label):
                w.native.setWordWrap(True)
        self.dock.native.setMaximumWidth(320)
        # live cortex pattern: cell size per quarter, by file and by distance from the stele
        self.fig = Figure(figsize=(3.2, 3.9), dpi=90, layout="constrained")
        self.prof_axes = self.fig.subplots(2, 1)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setMinimumHeight(330)
        lay = self.dock.native.layout()
        lay.insertWidget(lay.count() - 1, self.canvas)
        self.viewer.window.add_dock_widget(self.dock, area="right", name="Review")
        qv = self.viewer.window._qt_viewer   # napari's own layer panels are not needed
        for d in (qv.dockLayerList, qv.dockLayerControls):
            d.setVisible(False)

        v = self.viewer
        for k, t in enumerate(TARGETS, 1):
            v.bind_key(str(k), lambda _, t=t: self.set_target(t), overwrite=True)
        v.bind_key("Shift-H", lambda _: self.toggle_overlays(), overwrite=True)
        v.bind_key("Shift-S", lambda _: self.save_next(), overwrite=True)
        v.bind_key(".", lambda _: self.open(min(self.i + 1, len(self.files) - 1)), overwrite=True)
        v.bind_key(",", lambda _: self.open(max(self.i - 1, 0)), overwrite=True)
        v.mouse_drag_callbacks.append(self._mouse)
        v.mouse_move_callbacks.append(self._follow)
        v.cursor.style = "crosshair"

    def set_target(self, t):
        self.target = t
        with self.target_box.changed.blocked():
            self.target_box.value = t
        if hasattr(self, "tag"):
            self.tag.text = {"string": t, "color": TARGETS[t], "size": 11, "translation": [-14, 0]}

    def toggle_overlays(self):
        over = [l for l in self.viewer.layers if l.name not in ("image", "cursor") and l is not self.rings]
        show = not any(l.visible for l in over)
        for l in over:
            l.visible = show

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
        v.add_image(self.gray, name="image", opacity=self.img_op.value)
        lab = lambda data, name, cmap, op: v.add_labels(data, name=name, opacity=op,
                                                       colormap=DirectLabelColormap(color_dict=cmap))
        self.region = lab(region.astype(np.uint8), "region", {None: "transparent", 0: "transparent",
                                                               1: (1, 1, 0, 0.06), 2: (0, 0.9, 1, 0.3)}, 1)
        self.rings = lab(np.zeros_like(cells), "rings", {None: "gray", 0: "transparent", **RING_COLORS}, 0.4)
        self.mx = lab(mx, "metaxylem", {None: "red", 0: "transparent"}, 0.6)
        self.path = lab(np.zeros_like(cells), "file path", {None: "orangered", 0: "transparent"}, 0.55)
        self.cells = lab(cells, "cells", {None: "yellow", 0: "transparent"}, 0.9)
        self.cells.contour = 1
        self.air = lab(np.zeros_like(cells), "air spaces", {None: "magenta", 0: "transparent"}, 0.8)
        self.quads = v.add_shapes(name="quadrants + spokes", face_color="transparent")
        self.outline = v.add_shapes(name="outline", edge_color="white", edge_width=1.5, face_color="transparent")
        self.nums = v.add_points(np.empty((0, 2)), name="file numbers", size=0)
        self.tag = v.add_points(np.empty((0, 2)), name="cursor", size=0)
        for lyr in v.layers[1:]:
            lyr.editable = False
        self.rings.visible = self.show_rings.value
        v.layers.selection.active = self.cells
        props, _, _ = P.derive(cells, region, mx, self.angle)
        saved = P.load_paths(f)
        if saved is None:
            self.spokes, self.path_q, self.path_cells = P.auto_paths(cells, region, props, self.angle)
        else:
            self.spokes, self.path_q, self.path_cells = (np.array(a) for a in saved)
        self.set_target(self.target)
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
        n = len(self.files)
        for k in range(1, n + 1):
            j = (self.i + k) % n
            if P.load_state(self.files[j])[3] != "reviewed":
                return self.open(j)
        self.info.value = "All images reviewed. Run summarise.py"

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
        line = lambda t: np.array([[cy, cx], [cy - L * np.sin(np.radians(t)), cx + L * np.cos(np.radians(t))]])
        self.quads.data = []
        self.quads.add_lines([line(self.angle + 90 * k) for k in range(4)], edge_color="white", edge_width=1.2)
        self.quads.add_lines([line(t) for t in self.spokes], edge_color="orange", edge_width=1.5)
        self.draw_profile(ft)
        per_q = ft.groupby("quadrant").file.max().reindex([1, 2, 3, 4]).fillna(0).astype(int) if len(ft) else [0] * 4
        self.info.value = (f"{self.i + 1}/{len(self.files)}  ·  {self.status}\n"
                           f"Files Q1-Q4:  {'/'.join(str(n) for n in per_q)}\n"
                           f"Metaxylem {tr['mx_n']}  ·  air spaces {tr['n_air_spaces']}")

    def draw_profile(self, ft):
        """Live cortex pattern: per quarter, along the spoke and all cortical cells (same as the result figures)."""
        s = P.um_per_px(self.dir, self.rel, self.bar, self.settings)   # samples.xlsx override, else the bar
        um = not np.isnan(s)
        s = s if um else 1
        for ax in self.prof_axes:
            ax.clear()
        PP.quadrant_profiles(self.prof_axes, ft.assign(area=ft.area_px * s * s),
                             self.props.assign(area=self.props.area_px * s * s), "area",
                             "Cell area (µm²)" if um else "Cell area (px)", fontsize=7, compact=True)
        for ax in self.prof_axes:
            ax.tick_params(labelsize=7)
            ax.xaxis.label.set_size(8), ax.yaxis.label.set_size(8), ax.title.set_size(8)
        self.canvas.draw_idle()

    # ---------------------------------------------------------- mouse
    def _pos(self, event):
        y, x = (int(round(c)) for c in self.cells.world_to_data(event.position)[-2:])
        return (y, x) if 0 <= y < self.gray.shape[0] and 0 <= x < self.gray.shape[1] else None

    def _follow(self, viewer, event):
        """The target's name rides next to the cursor."""
        self.tag.data = np.array([event.position[-2:]])
        self.tag.text = {"string": self.target, "color": TARGETS[self.target], "size": 11, "translation": [-14, 0]}

    def _mouse(self, viewer, event):
        if event.button != 1:
            return
        mods = set(event.modifiers)
        remove = "Control" in mods          # Qt maps the Mac Cmd key to Control
        start = self._pos(event)
        if start is None:
            return
        if "Shift" in mods:                 # draw an outline
            self.viewer.camera.mouse_pan = False
            pts = [start]
            yield
            while event.type == "mouse_move":
                p = self._pos(event)
                if p and p != pts[-1]:
                    pts.append(p)
                    if len(pts) % 3 == 0:
                        self.outline.data = [np.array(pts)]
                        self.outline.shape_type = ["path"]
                yield
            self.viewer.camera.mouse_pan = True
            self.outline.data = []
            if len(pts) >= 3:
                self.outline_edit(np.array(pts), remove)
            return
        dragged = False                     # plain drag pans; a click without movement edits
        yield
        while event.type == "mouse_move":
            dragged = True
            yield
        if dragged:
            return
        if "Alt" in mods and self.target == "File path":
            self.move_spoke(*start)
        else:
            self.click_edit(*start, remove)

    # ---------------------------------------------------------- edits (same rules for every target)
    def _quadrant(self, mask):
        cy, cx = np.argwhere(mask).mean(0)
        return int(P.quadrant_of(cy, cx, P.centre_of(self.region.data), self.angle))

    def _path_add(self, labels):
        for l in labels:
            if l and l not in self.path_cells:
                self.path_q = np.r_[self.path_q, self._quadrant(self.cells.data == l)].astype(np.int16)
                self.path_cells = np.r_[self.path_cells, l].astype(np.int32)

    def _path_remove(self, labels):
        keep = ~np.isin(self.path_cells, list(labels))
        self.path_q, self.path_cells = self.path_q[keep], self.path_cells[keep]

    def _take(self, arr, m, merge):
        """Give pixels m to a new object in arr; one owner per pixel."""
        new = int(arr.max()) + 1
        hit = {int(l): int((arr == l).sum()) for l in np.unique(arr[m]) if l}
        arr[m] = new
        for l, before in hit.items():
            if merge:
                arr[arr == l] = new                 # vessels never touch: overlap = same vessel
            elif (arr == l).sum() < 0.25 * before:
                arr[arr == l] = 0                   # no slivers
        return new

    def click_edit(self, y, x, remove):
        cells, mx, region, t = self.cells.data, self.mx.data, self.region.data, self.target
        lbl = int(cells[y, x])
        if t == "Cell":
            if remove and lbl:
                cells[cells == lbl] = 0
            elif not remove:
                P.add_cell(self.gray, cells, y, x)
        elif t == "Metaxylem":
            if remove and mx[y, x]:
                mx[mx == mx[y, x]] = 0
            elif not remove and not mx[y, x]:
                P.add_mx(self.gray, mx, cells, y, x)
        elif t == "Stele":
            m = cells == lbl if lbl else P.lumen_at(self.gray, y, x, 8 * P.typical_cell_area(cells))
            region[m & (region >= 1) if remove else m] = 1 if remove else 2
        elif t == "File path" and lbl:
            self._path_remove([lbl]) if remove or lbl in self.path_cells else self._path_add([lbl])
        self._refresh()

    def outline_edit(self, pts, remove):
        rr, cc = fill_polygon(pts[:, 0], pts[:, 1], self.gray.shape)
        m = np.zeros(self.gray.shape, bool)
        m[rr, cc] = True
        if m.sum() < 10:
            return
        cells, mx, region, t = self.cells.data, self.mx.data, self.region.data, self.target
        inside = [int(l) for l in np.unique(cells[m]) if l and (cells[m] == l).sum() >= 0.5 * (cells == l).sum()]
        if t == "Cell":
            if remove:
                cells[m] = 0
            else:
                self._take(cells, m, merge=False)
        elif t == "Metaxylem":
            if remove:
                mx[m] = 0
            else:
                self._take(mx, m, merge=True)
                cells[m] = 0                        # a vessel is not a cell
        elif t == "Stele":
            region[m & (region >= 1)] = 1 if remove else 2     # stele stays inside the root
        elif t == "File path":
            if remove:
                self._path_remove(inside)
            else:
                empty = m & (cells == 0)            # an outline over a gap draws the missing cell
                if empty.sum() > 0.3 * m.sum():
                    inside.append(self._take(cells, m, merge=False))
                self._path_add(inside)
        self._refresh()

    def move_spoke(self, y, x):
        centre = P.centre_of(self.region.data)
        q = int(P.quadrant_of(y, x, centre, self.angle))
        t = float(np.degrees(np.arctan2(-(y - centre[0]), x - centre[1])))
        self.spokes[q - 1] = t
        keep = self.path_q != q
        new = P.path_along(self.cells.data, self.region.data, self.props, t)
        self.path_q = np.r_[self.path_q[keep], np.full(len(new), q)].astype(np.int16)
        self.path_cells = np.r_[self.path_cells[keep], new].astype(np.int32)
        self.recompute()

    def _refresh(self):
        for lyr in (self.cells, self.mx, self.region):
            lyr.refresh()
        self.recompute()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    Reviewer(Path(sys.argv[1]))
    napari.run()

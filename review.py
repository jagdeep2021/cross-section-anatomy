"""The app: measure an image set, set the scale, check outlines, make the tables. In napari.

  python review.py [<anatomy_dir>]

Buttons:  Open folder...     the experiment folder, its raw/ (images) or anatomy/. Images not measured yet
                             are outlined first (anatomy/ next to the images; progress in the panel;
                             stop/resume any time), then they open.
          um / pixel         the scale: always typed by you (scale bars are not read)
          Make results       tables + figures -> results/ next to anatomy/

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
import os
import subprocess
import sys
from pathlib import Path

import napari
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from magicgui.widgets import ComboBox, Container, FloatSlider, FloatSpinBox, Label, PushButton
from napari.qt.threading import thread_worker
from napari.utils.colormaps import DirectLabelColormap
from qtpy.QtCore import QEvent, QSize, Qt
from qtpy.QtWidgets import QFileDialog, QFrame, QHBoxLayout, QLabel, QMessageBox, QScrollArea, QWidget
from scipy import ndimage as ndi
from skimage.draw import polygon as fill_polygon

import pandas as pd
import pattern_plots as PP
import pipeline as P
import run_batch as RB
import summarise as SUM

TARGETS = {"Cell": "yellow", "Metaxylem": "red", "Stele": "cyan", "File path": "orange"}
CMD = "⌘" if sys.platform == "darwin" else "Ctrl"   # napari's "Control" is Cmd on a Mac
_KEYS = [("Click", "add"), (f"{CMD}-click", "remove (stele: the whole stele)"), ("Shift-drag", "draw outline, add"),
         (f"{CMD}-Shift-drag", "draw outline, remove"), ("Drag / scroll", "move / zoom"),
         ("Alt-click", "move spoke (file path)"), ("1 2 3 4", "cell, metaxylem, stele, file path"),
         ("Shift-H", "hide / show outlines"), ("Shift-S", "save + next"), (", .", "previous / next image")]
HELP = ("<table cellspacing='0' cellpadding='2' style='font-size:12px'>" +
        "".join(f"<tr><td style='color:#ffd479; padding-right:10px'><b>{k}</b></td><td>{v}</td></tr>"
                for k, v in _KEYS) + "</table>")


class WrapBar(QWidget):
    """The top bar: its boxes sit in one row, and move to the next row when the window is too narrow."""
    GAP = 6

    def __init__(self, boxes):
        super().__init__()
        self.boxes, self.dock = boxes, None
        for b in boxes:
            b.setParent(self)

    def _place(self, width):
        x = y = row_h = 0
        for b in self.boxes:
            s = b.sizeHint()
            if x and x + s.width() > width - 2 * self.GAP:
                x, y, row_h = 0, y + row_h + self.GAP, 0
            b.setGeometry(self.GAP + x, self.GAP // 2 + y, s.width(), s.height())
            x += s.width() + self.GAP
            row_h = max(row_h, s.height())
        h = y + row_h + self.GAP
        if h != self.height():
            self.setFixedHeight(h)
            if self.dock is not None:            # napari's dock would otherwise stretch far past the bar
                self.dock.setFixedHeight(h + 4)

    def resizeEvent(self, e):
        self._place(e.size().width())
        super().resizeEvent(e)

    def event(self, e):                          # a box changed size (e.g. "All Maize" -> "All Setaria")
        if e.type() == QEvent.Type.LayoutRequest:
            self._place(self.width())
        return super().event(e)

    def minimumSizeHint(self):
        return QSize(max(b.sizeHint().width() for b in self.boxes) + 2 * self.GAP, self.height())

    def sizeHint(self):
        return QSize(sum(b.sizeHint().width() + self.GAP for b in self.boxes) + self.GAP, self.height())


class Reviewer:
    def __init__(self, anatomy_dir: Path | None = None, show: bool = True):
        self.dir, self.files, self.i, self.busy = None, [], None, False
        self.viewer = napari.Viewer(title="Root cross-section anatomy", show=show)
        self.target = "Cell"
        self._build_dock()
        if anatomy_dir:
            self.load_experiment(anatomy_dir)
        else:
            self.info.value = ("Start: Open folder... and pick the experiment folder (or its raw folder). "
                               "Images not measured yet are measured first.")

    def load_experiment(self, d):
        d = Path(d)
        if not (d / "settings.json").exists() and (d / "anatomy" / "settings.json").exists():
            d = d / "anatomy"                          # the experiment folder was picked
        if not (d / "settings.json").exists():
            self.info.value = f"No measurements in {d}. Use Open folder... on the experiment or raw folder."
            return
        self.dir = d
        self.settings = json.loads((d / "settings.json").read_text())
        self.image_dir = P.image_dir(d, self.settings)
        self.angle = float(self.settings.get("quadrant_angle_deg", 0))
        self.files = sorted((d / "masks").glob("*.npz"))
        self.i = None
        if not self.files:
            self.info.value = "No images measured yet."
            return
        with self.pick.changed.blocked():
            self.pick.choices = self._choices()
        first = next((k for k, f in enumerate(self.files) if P.load_state(f)[3] != "reviewed"), 0)
        self.open(first)

    # ---------------------------------------------------------- measure / scale / summarise
    def _folder(self, title):
        start = str(self.dir.parent) if self.dir is not None else str(Path.home())
        return QFileDialog.getExistingDirectory(self.dock.native, title, start) or None

    def _run(self, work, done, progress=None):
        """Run a slow job in the background; the panel shows its progress."""
        if self.busy:
            self.info.value = "Busy: wait until the current job is done."
            return
        self.busy = True
        w = work()
        w.yielded.connect(progress or (lambda m: setattr(self.info, "value", m)))
        w.returned.connect(lambda res: (setattr(self, "busy", False), done(res)))
        w.errored.connect(lambda e: (setattr(self, "busy", False), self._error(e)))
        w.start()

    def open_folder(self):
        """One button for everything: pick the experiment folder, its raw/ (images) or its anatomy/.
        Images not measured yet are measured first (resumes where it stopped), then the images open."""
        d = self._folder("Choose the experiment folder (or its raw or anatomy folder)")
        if not d:
            return
        self.info.value = f"Looking in {d} ..."
        anatomy, images = RB.locate(d)
        if images is None and anatomy is None:
            self.info.value = (f"No images found in\n{d}\nPick the folder that holds the images "
                               f"({', '.join(P.IMAGE_EXT)}) or the experiment folder above it.")
            return
        todo = RB.n_to_measure(images, anatomy) if images is not None else 0
        if todo == 0:
            return self.load_experiment(anatomy)
        ok = QMessageBox.question(self.dock.native, "Measure images",
                                  f"Measure {todo} image(s) in\n{images}\n\nResults go to\n{anatomy}\n\n"
                                  "Each image takes about 10 s to 2 min. Images appear as they are done.")
        if ok != QMessageBox.StandardButton.Yes:
            self.info.value = "Cancelled."
            return
        if (anatomy / "settings.json").exists():
            self.load_experiment(anatomy)            # show what is measured already while the rest runs

        @thread_worker
        def work():
            yield (f"{todo} image(s) to measure in\n{images}\nLoading the outlining model (the first time "
                   f"this downloads ~1 GB: a few minutes)...")
            gen = RB.batch(images, anatomy)
            try:
                while True:
                    i, n, line = next(gen)
                    yield (i, n, f"Measuring image {i} of {n}\n{line}\n(you can close the app; Open folder "
                                 f"continues where it stopped)")
            except StopIteration as stop:
                return stop.value

        def progress(msg):
            if isinstance(msg, str):
                self.info.value = msg
                return
            i, n, text = msg
            self.viewer.title = f"Measuring {i} of {n} - Root cross-section anatomy"
            self._add_new(anatomy)
            self.info.value = text

        def done(out):
            self.viewer.title = "Root cross-section anatomy"
            self._add_new(out)
            self.info.value = f"All {len(self.files)} images measured. Now type the um / pixel scale, then Make results."

        self._run(work, done, progress)

    def _add_new(self, anatomy):
        """Images measured so far appear in the list; the first one opens straight away."""
        if self.dir is None or Path(anatomy).resolve() != self.dir.resolve():
            if (Path(anatomy) / "settings.json").exists() and any((Path(anatomy) / "masks").glob("*.npz")):
                self.load_experiment(anatomy)
            return
        files = sorted((self.dir / "masks").glob("*.npz"))
        if len(files) != len(self.files):
            cur = self.files[self.i] if self.i is not None else None
            self.files = files
            with self.pick.changed.blocked():
                self.pick.choices = self._choices()
                if cur is not None:
                    self.i = files.index(cur)
                    self.pick.value = self._choices()[self.i]
            if cur is None:
                self.open(0)

    def _error(self, e):
        import traceback
        txt = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        log = Path(__file__).with_name("app_log.txt")
        log.write_text(txt)
        print(txt, flush=True)
        self.info.value = f"Stopped: {e}\n(details in {log.name} in the tool folder; send it to Jagdeep)"

    def set_scale(self, all_in_group: bool):
        if self.dir is None or self.i is None:
            return
        f = SUM.write_samples(self.dir)
        s = pd.read_excel(f, dtype={"image_file": str, "plant": str})
        rows = s.image_file.str.split("/").str[0] == self.rel.split("/")[0] if all_in_group else s.image_file == self.rel
        s.loc[rows, "um_per_px"] = self.scale.value if self.scale.value > 0 else np.nan
        try:
            s.to_excel(f, index=False)
        except PermissionError:
            self.info.value = "samples.xlsx is open in Excel: close it, then set the scale again."
            return
        self.recompute()
        self.info.value += f"\nScale {self.scale.value:g} um/px set for {int(rows.sum())} image(s)."

    def summarise(self):
        if self.dir is None:
            return

        @thread_worker
        def work():
            yield "Making the tables and figures... (1-3 min)"
            SUM.main(self.dir.resolve())
            return self.dir.parent / "results"

        def done(res):
            self.info.value = f"Done. Tables and figures are in\n{res}"
            if sys.platform == "win32":
                os.startfile(res)
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(res)])

        self._run(work, done)

    # ---------------------------------------------------------- UI
    def _build_dock(self):
        opn, summ = PushButton(text="Open folder..."), PushButton(text="Make results")
        opn.changed.connect(self.open_folder)
        summ.changed.connect(self.summarise)
        self.scale = FloatSpinBox(value=0, min=0, max=1000, step=0.0001, label="um / pixel")
        self.scale.native.setDecimals(4)
        self.scale_grp, scale_img = PushButton(text="All"), PushButton(text="This image")
        self.scale_grp.changed.connect(lambda: self.set_scale(True))
        scale_img.changed.connect(lambda: self.set_scale(False))
        self.pick = ComboBox(choices=[], label="image")
        self.pick.changed.connect(lambda v: v and self.open(self._choices().index(v)))
        self.target_box = ComboBox(choices=list(TARGETS), value=self.target, label="editing")
        self.target_box.changed.connect(self.set_target)
        # two sliders instead of a box per layer: 0 hides, 1 is full strength
        self.img_op = FloatSlider(value=1.0, min=0, max=1.0, step=0.05, readout=False)
        self.img_op.changed.connect(lambda v: "image" in self.viewer.layers and
                                    setattr(self.viewer.layers["image"], "opacity", v))
        self.ov_op = FloatSlider(value=1.0, min=0, max=1.0, step=0.05, readout=False)
        self.ov_op.changed.connect(lambda v: self._apply_overlays())
        for s_ in (self.img_op, self.ov_op):
            s_.native.setFixedWidth(70)
        prev, nxt, save = PushButton(text="<"), PushButton(text=">"), PushButton(text="Save + next")
        prev.changed.connect(lambda: self.step(-1))
        nxt.changed.connect(lambda: self.step(1))
        save.changed.connect(self.save_next)
        self.info = Label(value="")
        self.info.native.setMinimumHeight(54)
        # Top bar: every control, left to right in the order of the work (the run and its scale | the image |
        # editing it). Right panel: what you look at (the numbers, the cortex pattern, what is shown, the keys).
        self.pick.native.setMinimumContentsLength(10)
        self.scale.native.setMinimumWidth(120)
        self.scale.native.lineEdit().setMinimumWidth(70)    # room for all 4 decimals (0.7194)

        def group(*items):                       # controls that belong together share one outlined box
            box = QFrame()
            box.setObjectName("grp")
            box.setStyleSheet("#grp {border: 1px solid #5c6370; border-radius: 6px;}")
            gl = QHBoxLayout(box)
            gl.setContentsMargins(6, 3, 6, 3)
            gl.setSpacing(4)
            for w in items:
                gl.addWidget(QLabel(w) if isinstance(w, str) else w.native)
            return box
        bar = WrapBar([group(opn, summ),
                       group("µm/px", self.scale, self.scale_grp, scale_img),
                       group(self.pick, prev, nxt, save),
                       group("editing", self.target_box),
                       group("image", self.img_op, "outlines", self.ov_op)])
        bar.dock = self.viewer.window.add_dock_widget(bar, area="top", name="Controls")
        bar.dock.setTitleBarWidget(QWidget())    # no close button: the bar can't be lost by accident
        bar.dock.setFixedHeight(bar.height() + 4)
        help_lbl = QLabel(HELP)
        help_lbl.setTextFormat(Qt.TextFormat.RichText)
        self.dock = Container(widgets=[self.info])
        self.info.native.setWordWrap(True)
        self.info.native.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.dock.native.setMinimumWidth(240)   # no maximum: drag the panel edge to make it wider or narrower
        # live cortex pattern: cell size per quarter, by file and by distance from the stele
        self.fig = Figure(figsize=(3.2, 3.9), dpi=90, layout="constrained")
        self.prof_axes = self.fig.subplots(2, 1)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setMinimumHeight(300)
        lay = self.dock.native.layout()
        lay.setSpacing(4)
        lay.addWidget(self.canvas, stretch=1)            # the graphs take all the spare height
        lay.addWidget(help_lbl)
        scroll = QScrollArea()                   # a short screen scrolls the panel instead of cutting it off
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidget(self.dock.native)
        self.viewer.window.add_dock_widget(scroll, area="right", name="Cortex pattern")
        qv = self.viewer.window._qt_viewer   # napari's own layer panels are not needed
        for d in (qv.dockLayerList, qv.dockLayerControls):
            d.setVisible(False)
        self.viewer.window._qt_window.statusBar().hide()   # napari's own key hints (1 = eraser ...) contradict ours

        v = self.viewer
        for k, t in enumerate(TARGETS, 1):
            v.bind_key(str(k), lambda _, t=t: self.set_target(t), overwrite=True)
        v.bind_key("Shift-H", lambda _: self.toggle_overlays(), overwrite=True)
        v.bind_key("Shift-S", lambda _: self.save_next(), overwrite=True)
        v.bind_key(".", lambda _: self.step(1), overwrite=True)
        v.bind_key(",", lambda _: self.step(-1), overwrite=True)
        v.mouse_drag_callbacks.append(self._mouse)
        v.mouse_move_callbacks.append(self._follow)
        v.cursor.style = "crosshair"

    def step(self, k):
        if self.files and self.i is not None:
            self.open(min(max(self.i + k, 0), len(self.files) - 1))

    def set_target(self, t):
        self.target = t
        with self.target_box.changed.blocked():
            self.target_box.value = t
        if hasattr(self, "tag"):
            self.tag.text = {"string": t, "color": TARGETS[t], "size": 11, "translation": [-14, 0]}

    def toggle_overlays(self):
        if self.i is None:
            return
        if self.ov_op.value > 0:
            self._ov_last, self.ov_op.value = self.ov_op.value, 0
        else:
            self.ov_op.value = getattr(self, "_ov_last", 1.0)

    def _apply_overlays(self):
        if self.i is None:
            return
        a = self.ov_op.value
        for lyr, base in self._base.items():
            lyr.opacity, lyr.visible = base * a, a > 0

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
        self.gray = P.load_gray(self.image_dir / self.rel)
        cells, region, mx, self.status = P.load_state(f)
        v = self.viewer
        v.layers.clear()
        v.add_image(self.gray, name="image", opacity=self.img_op.value)
        lab = lambda data, name, cmap, op: v.add_labels(data, name=name, opacity=op,
                                                       colormap=DirectLabelColormap(color_dict=cmap))
        self.region = lab(region.astype(np.uint8), "region", {None: "transparent", 0: "transparent",
                                                               1: (1, 1, 0, 0.06), 2: (0, 0.9, 1, 0.3)}, 1)
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
        self._base = {lyr: lyr.opacity for lyr in v.layers[1:] if lyr is not self.tag}   # each layer's full strength
        self._apply_overlays()
        v.layers.selection.active = self.cells
        props, _, _ = P.derive(cells, region, mx, self.angle)
        saved = P.load_paths(f)
        if saved is None:
            self.spokes, self.path_q, self.path_cells = P.auto_paths(cells, region, props, self.angle)
        else:
            self.spokes, self.path_q, self.path_cells = (np.array(a) for a in saved)
        self.set_target(self.target)
        self.recompute()
        um = P.um_per_px(self.dir, self.rel)
        with self.scale.changed.blocked():
            self.scale.value = 0 if np.isnan(um) else um
        self.scale_grp.text = f"All {self.rel.split('/')[0]}" if "/" in self.rel else "All"
        with self.pick.changed.blocked():
            self.pick.choices = self._choices()
            self.pick.value = self._choices()[i]

    def save(self):
        P.save_state(self.files[self.i], self.cells.data.astype(np.uint16), self.region.data.astype(np.uint8),
                     self.mx.data.astype(np.uint16), status="reviewed", image_file=self.rel,
                     spokes=self.spokes, path_q=np.asarray(self.path_q, np.int16),
                     path_cells=np.asarray(self.path_cells, np.int32))
        self.status = "reviewed"

    def save_next(self):
        if self.i is None:
            return
        self.save()
        n = len(self.files)
        for k in range(1, n + 1):
            j = (self.i + k) % n
            if P.load_state(self.files[j])[3] != "reviewed":
                return self.open(j)
        self.info.value = "All images checked. Click Make results."

    # ---------------------------------------------------------- display
    def recompute(self):
        cells, region = self.cells.data, self.region.data
        props, _, tr = P.derive(cells, region, self.mx.data, self.angle)
        self.props = props
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
        s = P.um_per_px(self.dir, self.rel)   # samples.xlsx, given by the user
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
        if self.i is None:
            return
        self.tag.data = np.array([event.position[-2:]])
        self.tag.text = {"string": self.target, "color": TARGETS[self.target], "size": 11, "translation": [-14, 0]}

    def _mouse(self, viewer, event):
        if event.button != 1 or self.i is None:
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
            if remove:
                region[region == 2] = 1             # the whole stele goes; redraw it with Shift-drag
                self._refresh()
                self.info.value += "\nStele removed. Shift-drag around the stele to draw it again."
                return
            else:
                m = cells == lbl if lbl else P.lumen_at(self.gray, y, x, 8 * P.typical_cell_area(cells))
                region[m] = 2
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
    r = Reviewer(Path(sys.argv[1]) if len(sys.argv) > 1 else None)
    # napari reopens where it last closed, which can be off a screen that is no longer attached:
    # always open on the main screen, full size
    win = r.viewer.window._qt_window
    win.setGeometry(napari.qt.get_qapp().primaryScreen().availableGeometry())
    win.showMaximized()
    win.raise_()
    napari.run()

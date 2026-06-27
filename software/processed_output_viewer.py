"""
processed_output_viewer.py
==========================
可视化 processed_output.csv：勾选/取消各列曲线显示。
支持全部三个参数：hbo、hbr、cyt
用法:
  python processed_output_viewer.py
  python processed_output_viewer.py --input result_table/2026-06-24_14-19-17/processed_output.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyqtgraph as pg
from PyQt5 import QtCore, QtWidgets

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "result_table" / "2026-06-27_14-05-21" / "processed_output.csv"

TIME_CANDIDATES = ("Time", "Time (s)")


def _resolve_time_column(columns: list[str]) -> str:
    for name in TIME_CANDIDATES:
        if name in columns:
            return name
    raise ValueError(f"未找到时间列（期望 {TIME_CANDIDATES}），现有列: {columns}")


def load_processed_csv(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError("CSV 为空。")

    time_col = _resolve_time_column(list(df.columns))
    times = df[time_col].to_numpy(dtype=float)
    channels: dict[str, np.ndarray] = {}
    for col in df.columns:
        if col == time_col:
            continue
        channels[str(col)] = df[col].to_numpy(dtype=float)
    if not channels:
        raise ValueError("除时间列外没有其他数据列。")
    return times, channels


class ProcessedOutputViewer(QtWidgets.QMainWindow):
    def __init__(self, initial_path: Path | None = None) -> None:
        super().__init__()
        self.setWindowTitle("processed_output 可视化")
        self.resize(1200, 720)

        self._csv_path: Path | None = None
        self._times: np.ndarray | None = None
        self._channels: dict[str, np.ndarray] = {}
        self._checkboxes: dict[str, QtWidgets.QCheckBox] = {}
        self._curve_items: dict[str, pg.PlotDataItem] = {}

        self._build_ui()
        if initial_path is not None and initial_path.is_file():
            self._load_file(initial_path)

    def _build_ui(self) -> None:
        central = QtWidgets.QWidget()
        root = QtWidgets.QHBoxLayout(central)

        controls = QtWidgets.QWidget()
        controls.setMaximumWidth(360)
        layout = QtWidgets.QVBoxLayout(controls)
        layout.setContentsMargins(10, 10, 10, 10)

        self.path_label = QtWidgets.QLabel("未加载文件")
        self.path_label.setWordWrap(True)
        open_btn = QtWidgets.QPushButton("打开 CSV…")
        open_btn.clicked.connect(self._open_file_dialog)

        self.info_label = QtWidgets.QLabel("采样点: —")
        self.status_label = QtWidgets.QLabel("")
        self.status_label.setStyleSheet("color: #c0392b;")
        self.status_label.setWordWrap(True)

        preset_row = QtWidgets.QHBoxLayout()
        show_all_btn = QtWidgets.QPushButton("全部显示")
        hide_all_btn = QtWidgets.QPushButton("全部隐藏")
        ssr_only_btn = QtWidgets.QPushButton("仅 SSR")
        show_all_btn.clicked.connect(lambda: self._set_all_visible(True))
        hide_all_btn.clicked.connect(lambda: self._set_all_visible(False))
        ssr_only_btn.clicked.connect(self._show_ssr_only)
        preset_row.addWidget(show_all_btn)
        preset_row.addWidget(hide_all_btn)
        preset_row.addWidget(ssr_only_btn)

        self.channel_box = QtWidgets.QGroupBox("曲线显示")
        self.channel_layout = QtWidgets.QVBoxLayout(self.channel_box)
        self.channel_scroll = QtWidgets.QScrollArea()
        self.channel_scroll.setWidgetResizable(True)
        self.channel_inner = QtWidgets.QWidget()
        self.channel_inner_layout = QtWidgets.QVBoxLayout(self.channel_inner)
        self.channel_inner_layout.addStretch()
        self.channel_scroll.setWidget(self.channel_inner)

        layout.addWidget(self.path_label)
        layout.addWidget(open_btn)
        layout.addWidget(self.info_label)
        layout.addLayout(preset_row)
        layout.addWidget(self.channel_scroll, stretch=1)
        layout.addWidget(self.status_label)

        self.plot_widget = pg.PlotWidget()
        self.plot_widget.showGrid(x=True, y=True, alpha=0.25)
        self.plot_widget.setLabel("bottom", "Time (s)")
        self.plot_widget.setLabel("left", "Concentration change (mM)")
        self.plot_widget.addLegend(offset=(10, 10))

        root.addWidget(controls, stretch=0)
        root.addWidget(self.plot_widget, stretch=1)
        self.setCentralWidget(central)

    def _open_file_dialog(self) -> None:
        start_dir = str(self._csv_path.parent if self._csv_path else SCRIPT_DIR)
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "选择 processed_output.csv",
            start_dir,
            "CSV Files (*.csv);;All Files (*)",
        )
        if path:
            self._load_file(Path(path))

    def _clear_channel_controls(self) -> None:
        while self.channel_inner_layout.count() > 1:
            item = self.channel_inner_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self._checkboxes.clear()

    def _load_file(self, path: Path) -> None:
        try:
            times, channels = load_processed_csv(path)
        except Exception as exc:
            self.status_label.setText(f"加载失败：{exc}")
            return

        self._csv_path = path
        self._times = times
        self._channels = channels
        self.path_label.setText(str(path))
        self.info_label.setText(f"采样点: {len(times)} | 通道: {len(channels)}")
        self.status_label.setText("")

        self._clear_channel_controls()
        colors = [
            "#c0392b",
            "#2980b9",
            "#27ae60",
            "#8e44ad",
            "#d35400",
            "#16a085",
            "#2c3e50",
            "#e67e22",
            "#1abc9c",
        ]
        for idx, name in enumerate(channels):
            cb = QtWidgets.QCheckBox(name)
            color = colors[idx % len(colors)]
            cb.setStyleSheet(f"QCheckBox {{ color: {color}; font-weight: 500; }}")
            cb.setChecked(name.endswith("_ssr_hbo"))
            cb.toggled.connect(self._update_plot)
            self._checkboxes[name] = cb
            self.channel_inner_layout.insertWidget(idx, cb)

        self._rebuild_curves()
        self._update_plot()

    def _set_all_visible(self, visible: bool) -> None:
        for cb in self._checkboxes.values():
            cb.blockSignals(True)
            cb.setChecked(visible)
            cb.blockSignals(False)
        self._update_plot()

    def _show_ssr_only(self) -> None:
        for name, cb in self._checkboxes.items():
            cb.blockSignals(True)
            cb.setChecked("_ssr_" in name)
            cb.blockSignals(False)
        self._update_plot()

    def _rebuild_curves(self) -> None:
        self.plot_widget.clear()
        self.plot_widget.addLegend(offset=(10, 10))
        self._curve_items.clear()
        if self._times is None:
            return

        colors = [
            "#c0392b",
            "#2980b9",
            "#27ae60",
            "#8e44ad",
            "#d35400",
            "#16a085",
            "#2c3e50",
            "#e67e22",
            "#1abc9c",
        ]
        for idx, name in enumerate(self._channels):
            color = colors[idx % len(colors)]
            pen = pg.mkPen(
                color=color,
                width=2,
                style=QtCore.Qt.DashLine if "_ssr_" in name else QtCore.Qt.SolidLine,
            )
            item = self.plot_widget.plot(
                self._times,
                self._channels[name],
                pen=pen,
                name=name,
            )
            item.setVisible(self._checkboxes[name].isChecked())
            self._curve_items[name] = item

    def _update_plot(self) -> None:
        if not self._curve_items:
            return
        for name, item in self._curve_items.items():
            cb = self._checkboxes.get(name)
            if cb is not None:
                item.setVisible(cb.isChecked())


def main() -> None:
    parser = argparse.ArgumentParser(description="可视化 processed_output.csv")
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="processed_output.csv 路径",
    )
    args = parser.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    initial = args.input if args.input.is_file() else None
    viewer = ProcessedOutputViewer(initial_path=initial)
    viewer.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

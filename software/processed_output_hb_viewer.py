"""
processed_output_hb_viewer.py
=============================
可视化 HbO/HbR 版 processed_output_hb.csv：勾选/取消各列曲线显示。

默认打开:
    software/result_table/2026-06-24_14-25-15/processed_output_hb.csv

用法:
    python software/processed_output_hb_viewer.py
    python software/processed_output_hb_viewer.py --input software/result_table/某次/processed_output_hb.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from PyQt5 import QtWidgets

from processed_output_viewer import ProcessedOutputViewer

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "result_table" / "2026-06-24_15-50-55" / "processed_output_hb.csv"


class ProcessedOutputHbViewer(ProcessedOutputViewer):
    """HbO/HbR 专用默认入口，复用 processed_output_viewer 的列选择绘图逻辑。"""

    def __init__(self, initial_path: Path | None = None) -> None:
        super().__init__(initial_path=initial_path)
        self.setWindowTitle("processed_output_hb 可视化")
        self.plot_widget.setLabel("left", "HbO / HbR concentration change")


def main() -> None:
    parser = argparse.ArgumentParser(description="可视化 processed_output_hb.csv")
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="processed_output_hb.csv 路径",
    )
    args = parser.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    initial = args.input if args.input.is_file() else None
    viewer = ProcessedOutputHbViewer(initial_path=initial)
    viewer.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

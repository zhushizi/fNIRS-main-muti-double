"""
双接收源五波长实时 HbO / HbR / Cyt 曲线查看器。

从串口持续读取五波长光强，用前几秒为 S1-D1 / S1-D2 分别建立基线，
再用广义 MBLL 输出 HbO、HbR 和探索版 Cyt/oxCCO 相对变化。
"""

from __future__ import annotations

import sys
import time
from collections import deque

import numpy as np
import nirsimple.preprocessing as nsp
import pyqtgraph as pg
import serial
from PyQt5 import QtCore, QtWidgets
from scipy.signal import butter, resample_poly, sosfiltfilt

from config import (
    ACK_TIMEOUT_SECONDS,
    BAUD_RATE,
    CHROMOPHORE_TYPES,
    DETECTOR_CHANNELS,
    DEFAULT_INTENSITY_MA,
    MAX_RETRIES,
    MBLL_DEFAULT_AGE,
    SERIAL_PORT,
    TIMEOUT,
    WAVELENGTH_CHANNELS,
)
from fNIRS_processing import generalized_mbll
from protocol import FrameReader, build_command_frame, parse_data_frame, send_frame_with_ack

# 用前若干秒的平均光强作为 OD 基线
BASELINE_SECONDS = 5.0
# MBLL 更新节流：至少间隔多少秒再算下一对（避免 GUI 卡顿）
UPDATE_INTERVAL_S = 0.1
# 绘图时间窗（秒）
WINDOW_SECONDS = 60.0
BP_LOW_HZ = 0.05
BP_HIGH_HZ = 0.5
BP_ORDER = 4
BP_TARGET_FS = 20.0


def butter_bandpass_sos(lowcut: float, highcut: float, fs: float, order: int = 4):
    """构造带通滤波器的 SOS 形式。"""
    nyq = 0.5 * fs
    if nyq <= 0 or lowcut >= highcut or highcut >= nyq:
        return None
    return butter(order, [lowcut / nyq, highcut / nyq], btype="band", output="sos")


def smart_bandpass(
    data: np.ndarray,
    fs: float,
    lowcut: float = BP_LOW_HZ,
    highcut: float = BP_HIGH_HZ,
    order: int = BP_ORDER,
    target_fs: float = BP_TARGET_FS,
) -> np.ndarray:
    """
    对 OD 数据做稳健带通。

    当采样率过高时，先降采样再滤波，最后升采样回来，
    可以减少数值不稳定和不必要的计算量。
    """
    if data.shape[1] < max(16, 3 * (order + 1) + 1):
        return data

    if fs > target_fs + 1:
        decim = int(round(fs / target_fs))
        fs_ds = fs / decim
        data_ds = resample_poly(data, up=1, down=decim, axis=1)
    else:
        decim, fs_ds, data_ds = 1, fs, data

    sos = butter_bandpass_sos(lowcut, highcut, fs_ds, order)
    if sos is None or data_ds.shape[1] < max(16, 3 * (order + 1) + 1):
        return data

    padlen = min(data_ds.shape[1] - 1, 3 * (order + 1))
    if padlen <= 0:
        return data

    data_bp = sosfiltfilt(sos, data_ds, axis=1, padtype="odd", padlen=padlen)
    if decim > 1:
        data_bp = resample_poly(data_bp, up=decim, down=1, axis=1)
    return data_bp


def open_serial() -> serial.Serial:
    return serial.Serial(SERIAL_PORT, baudrate=BAUD_RATE, timeout=TIMEOUT)


class SerialReaderThread(QtCore.QThread):
    """后台串口读线程，与 adc_live 相同逻辑：发启动命令、持续读 0x02、回 ACK、发信号。"""
    newData = QtCore.pyqtSignal(float, int, int, int, str)

    def __init__(self, ser: serial.Serial, parent=None):
        super().__init__(parent)
        self.ser = ser
        self.reader = FrameReader(ser)
        self.running = True
        self.start_time = time.time()

    def _emit(self, sample):
        t = time.time() - self.start_time
        self.newData.emit(
            t,
            sample.detector_code,
            sample.wavelength_code,
            sample.value,
            sample.channel_name or "",
        )

    def _start_stream(self) -> bool:
        cmd = build_command_frame(True, DEFAULT_INTENSITY_MA)
        deadline = time.time() + max(ACK_TIMEOUT_SECONDS * 10, 0.2)
        for attempt in range(MAX_RETRIES + 1):
            self.ser.write(cmd)
            while time.time() < deadline:
                frame = self.reader.read_frame(timeout_seconds=min(TIMEOUT, max(0.01, deadline - time.time())))
                if frame is None:
                    continue
                if frame.frame_type == 0x02:
                    self._emit(parse_data_frame(frame))
                    return True
            deadline = time.time() + max(ACK_TIMEOUT_SECONDS * 10, 0.2)
        return False

    def run(self):
        started = self._start_stream()
        if not started:
            print("[hbo_hbr_live] Start command sent (max attempts reached). Keep waiting for data frames...")
        while self.running:
            frame = self.reader.read_frame(timeout_seconds=TIMEOUT)
            if frame is None or frame.frame_type != 0x02:
                continue
            self._emit(parse_data_frame(frame))
        try:
            send_frame_with_ack(self.ser, self.reader, build_command_frame(False, DEFAULT_INTENSITY_MA))
        except Exception:
            pass

    def stop(self):
        self.running = False


class MainWindow(QtWidgets.QWidget):
    """实时 HbO / HbR / Cyt 坐标图。"""

    def __init__(self):
        super().__init__()
        pg.setConfigOption("background", "w")
        pg.setConfigOption("foreground", "k")
        self.setWindowTitle("S1-D1 / S1-D2 HbO / HbR / Cyt 实时")

        self.required_keys = {(det.code, wl.code) for det in DETECTOR_CHANNELS for wl in WAVELENGTH_CHANNELS}
        self.latest_values: dict[tuple[int, int], float] = {}
        self.baseline_values: dict[tuple[int, int], list[float]] = {key: [] for key in self.required_keys}
        self.references: dict[tuple[int, int], float] | None = None
        self.last_mbll_time = 0.0
        self.time_od = deque(maxlen=3000)
        self.od_history = {
            det.code: {wl.code: deque(maxlen=3000) for wl in WAVELENGTH_CHANNELS}
            for det in DETECTOR_CHANNELS
        }

        self.curve_data = {}
        colors = ["red", "blue", "darkGreen", "magenta", "cyan", "black"]

        layout = QtWidgets.QVBoxLayout(self)
        self.plot = pg.PlotWidget(title="HbO / HbR / Cyt 实时")
        self.plot.showGrid(x=True, y=True)
        self.plot.setLabel("bottom", "Time (s)")
        self.plot.setLabel("left", "Concentration (Δ, Cyt UCL extinction)")
        self.plot.addLegend()
        self.plot.setXRange(0, WINDOW_SECONDS, padding=0)
        idx = 0
        for det in DETECTOR_CHANNELS:
            for chrom in CHROMOPHORE_TYPES:
                key = (det.code, chrom)
                self.curve_data[key] = {
                    "time": deque(maxlen=3000),
                    "value": deque(maxlen=3000),
                    "curve": self.plot.plot(
                        pen=pg.mkPen(colors[idx % len(colors)], width=2),
                        name=f"{det.name}_{chrom}",
                    ),
                }
                idx += 1
        layout.addWidget(self.plot)

        self.status = QtWidgets.QLabel("建立基线…")
        layout.addWidget(self.status)

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self._update_plot)
        self.timer.start(50)

    @QtCore.pyqtSlot(float, int, int, int, str)
    def on_new_data(self, elapsed: float, detector_code: int, wavelength_code: int, value: int, channel_name: str):
        if value <= 0:
            return
        key = (detector_code, wavelength_code)
        if key not in self.required_keys:
            return

        value_f = float(value)
        self.latest_values[key] = value_f
        if self.references is None and elapsed <= BASELINE_SECONDS:
            self.baseline_values[key].append(value_f)

        if self.references is None and elapsed >= BASELINE_SECONDS:
            if all(self.baseline_values[k] for k in self.required_keys):
                self.references = {
                    k: max(float(np.mean(v)), 1.0) for k, v in self.baseline_values.items()
                }
                self.status.setText("基线就绪，开始计算 HbO/HbR/Cyt")
            else:
                self.status.setText("基线数据不足，继续等待完整 2x5 波长周期")
                return

        if self.references is None or set(self.latest_values.keys()) != self.required_keys:
            return
        if (elapsed - self.last_mbll_time) < UPDATE_INTERVAL_S:
            return

        self._compute_concentrations(elapsed)

    def _compute_concentrations(self, elapsed: float):
        self.time_od.append(elapsed)
        wavelengths = [wl.mbll_nm for wl in WAVELENGTH_CHANNELS]
        dpfs = [nsp.get_dpf(wl, MBLL_DEFAULT_AGE) for wl in wavelengths]

        for det in DETECTOR_CHANNELS:
            for wl in WAVELENGTH_CHANNELS:
                key = (det.code, wl.code)
                intensity = max(self.latest_values[key], 1.0)
                reference = self.references[key] if self.references is not None else intensity
                self.od_history[det.code][wl.code].append(-np.log10(intensity / reference))

            od_matrix = np.vstack(
                [np.asarray(self.od_history[det.code][wl.code], dtype=float) for wl in WAVELENGTH_CHANNELS]
            )
            od_times = np.asarray(self.time_od, dtype=float)
            dt = np.median(np.diff(od_times)) if od_times.size >= 2 else np.nan
            fs = 1.0 / dt if np.isfinite(dt) and dt > 0 else 1.0
            od_filt = smart_bandpass(od_matrix, fs, lowcut=BP_LOW_HZ, highcut=BP_HIGH_HZ, order=BP_ORDER)
            try:
                delta_c = generalized_mbll(od_filt[:, -1:], wavelengths, dpfs, det.distance_cm)
            except Exception as exc:
                self.status.setText(f"MBLL 计算失败: {exc}")
                return
            for idx, chrom in enumerate(CHROMOPHORE_TYPES):
                curve = self.curve_data[(det.code, chrom)]
                curve["time"].append(elapsed)
                curve["value"].append(float(delta_c[idx, -1]))

        self.last_mbll_time = elapsed
        self.status.setText(f"HbO/HbR/Cyt updated t={elapsed:.2f}s (Cyt UCL extinction)")

    def _update_plot(self):
        t_max = 0.0
        for series in self.curve_data.values():
            if not series["time"]:
                continue
            series["curve"].setData(list(series["time"]), list(series["value"]))
            t_max = max(t_max, series["time"][-1])
        if t_max > 0:
            x_min = max(0.0, t_max - WINDOW_SECONDS)
            self.plot.setXRange(x_min, max(t_max, WINDOW_SECONDS), padding=0)


def main():
    print(f"[hbo_hbr_live] Opening {SERIAL_PORT} @ {BAUD_RATE} ...")
    ser = open_serial()
    ser.reset_input_buffer()
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    thread = SerialReaderThread(ser)
    thread.newData.connect(window.on_new_data)
    thread.start()
    app.aboutToQuit.connect(lambda: (thread.stop(), thread.wait(2000), ser.close()))
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

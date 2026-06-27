"""
双接收源五波长 fNIRS 采集与处理主脚本。

完整链路如下：
1. 通过 26B 协议采集 S1-D1 / S1-D2 五波长原始数据
2. 写入 all_groups.csv
3. 对每个接收源/波长组合做阈值/低通/RMS 预处理
4. 聚合成 2 接收源 x 5 波长光强宽表
5. 计算 OD -> 短距回归校正 -> 广义 MBLL -> HbO/HbR/Cyt
6. 输出 processed_output.csv
"""

from __future__ import annotations

import csv
import os
import threading
import time
from datetime import datetime

import nirsimple.preprocessing as nsp
import numpy as np
import pandas as pd
import serial
from scipy.signal import butter, filtfilt, resample_poly, sosfiltfilt
from tabulate import tabulate

from config import (
    ACK_TIMEOUT_SECONDS,
    BAUD_RATE,
    BP_HIGH_HZ,
    BP_LOW_HZ,
    BP_ORDER,
    BP_TARGET_FS_HZ,
    CHROMOPHORE_TYPES,
    CYT_DIFFERENCE_EXTINCTION,
    DETECTOR_CHANNELS,
    DEFAULT_INTENSITY_MA,
    DEFAULT_STREAM_ENABLED,
    INTENSITY_COLUMNS,
    MAX_RETRIES,
    MBLL_DEFAULT_AGE,
    RAW_OUTPUT_CSV,
    SERIAL_PORT,
    TIMEOUT,
    WAVELENGTH_CHANNELS,
    WAVELENGTH_OFF_CODE,
    detector_channel_by_code,
    intensity_column_for,
    intensity_columns_for_detector,
    wavelength_channel_by_code,
)
from protocol import (
    FrameReader,
    build_command_frame,
    parse_data_frame,
    send_frame_with_ack,
)


def open_serial() -> serial.Serial:
    """按配置打开串口。"""
    return serial.Serial(SERIAL_PORT, baudrate=BAUD_RATE, timeout=TIMEOUT)


def frame_to_hex(frame_bytes: bytes) -> str:
    """把协议帧转成十六进制字符串，便于现场联调。"""
    return " ".join(f"{b:02X}" for b in frame_bytes)


def _start_enter_listener(stop_event: threading.Event) -> threading.Thread:
    """起一个后台线程监听回车，用于结束采集。"""
    def _wait_for_enter() -> None:
        try:
            input("Press Enter to stop capture and continue processing.\n")
            stop_event.set()
        except EOFError:
            return

    thread = threading.Thread(target=_wait_for_enter, daemon=True)
    thread.start()
    return thread


def start_stream_until_data(
    ser: serial.Serial,
    reader: FrameReader,
    intensity_ma: int,
) -> object | None:
    """
    发送启动命令并等待第一帧数据。

    有些固件版本对单次启动命令不敏感；这里复用实时 ADC 脚本的策略：
    最多发送 MAX_RETRIES + 1 次，收到第一帧 0x02 即认为启动成功。
    """
    command = build_command_frame(DEFAULT_STREAM_ENABLED, intensity_ma)
    startup_timeout = max(ACK_TIMEOUT_SECONDS * 10, 0.2)

    for attempt in range(MAX_RETRIES + 1):
        print(f"TX start command attempt {attempt + 1}: {frame_to_hex(command)}")
        ser.write(command)
        deadline = time.time() + startup_timeout
        while time.time() < deadline:
            frame = reader.read_frame(timeout_seconds=min(TIMEOUT, max(0.01, deadline - time.time())))
            if frame is None:
                continue
            if frame.frame_type == 0x02:
                print("Stream became active; first data frame received.")
                return frame
            print(f"Ignored startup frame type: 0x{frame.frame_type:02X}")

    print("No data frame received during startup; keep waiting in capture loop.")
    return None


def threshold_filter(
    df: pd.DataFrame,
    signal_col: str = "Value",
    lower_threshold: int = 50000,
    upper_threshold: int = 300000,
    zero_level: int = 170000,
) -> pd.DataFrame:
    """对采样值做阈值抑制，超限值直接替换为 zero_level。"""
    filtered = df.copy()
    filtered[signal_col] = np.where(
        (df[signal_col] < lower_threshold) | (df[signal_col] > upper_threshold),
        zero_level,
        df[signal_col],
    )
    return filtered


def butter_lowpass_filter(
    df: pd.DataFrame,
    cutoff_hz: float,
    fs: float,
    order: int = 4,
    signal_col: str = "Value",
) -> pd.DataFrame:
    """按接收源/波长组合分别低通，避免不同 LED/PD 的电平互相污染。"""
    if len(df) < max(12, order * 3):
        return df.copy()

    filtered = df.copy()
    nyquist = 0.5 * fs
    if nyquist <= 0 or cutoff_hz >= nyquist:
        return filtered

    b, a = butter(order, cutoff_hz / nyquist, btype="low", analog=False)
    group_cols = [col for col in ("DetectorId", "Wavelength") if col in filtered.columns]
    grouped = filtered.groupby(group_cols, sort=False) if group_cols else [(None, filtered)]
    for _, group in grouped:
        if len(group) < max(12, order * 3):
            continue
        padlen = min(len(group) - 1, 3 * (max(len(a), len(b)) - 1))
        if padlen <= 0:
            continue
        filtered.loc[group.index, signal_col] = filtfilt(
            b,
            a,
            group[signal_col].astype(float),
            padlen=padlen,
        )
    return filtered


def sliding_window_rms(
    df: pd.DataFrame,
    signal_col: str = "Value",
    wavelength_col: str = "Wavelength",
    detector_col: str = "DetectorId",
    remove_dc: bool = False,
) -> pd.DataFrame:
    """
    按 (DetectorId, Wavelength) 连续段切段，并把每一段压成一行 RMS 代表值。

    这样后续在做双波长配对时，每个波长块只保留一个物理有效样本，
    避免把同一块里的重复 RMS 再展开成多行。
    """
    if df.empty:
        return df.copy()

    keys = df[[detector_col, wavelength_col]].astype(int).to_numpy()
    change_points = np.where(np.any(np.diff(keys, axis=0) != 0, axis=1))[0] + 1
    segments = np.split(np.arange(len(df)), change_points)
    rows = [] # 存储每一段数据的 RMS 代表值

    for segment in segments:
        if len(segment) == 0:
            continue
        segment_idx = segment.tolist()
        segment_data = df.loc[segment_idx, signal_col].astype(float)
        if remove_dc:
            segment_data = segment_data - segment_data.mean()
        rms_val = float(np.sqrt(np.mean(np.square(segment_data))))
        rep_row = df.loc[segment_idx[0]].copy()
        rep_row[signal_col] = rms_val
        if "Time (s)" in df.columns:
            rep_row["Time (s)"] = float(df.loc[segment_idx, "Time (s)"].mean())
        rows.append(rep_row)

    return pd.DataFrame(rows).reset_index(drop=True)


def aggregate_wavelength_cycles(
    df: pd.DataFrame,
    signal_col: str = "Value",
    mode_col: str = "Wavelength",
) -> pd.DataFrame:
    """
    将 RMS 后的 2 接收源 x 5 波长段合并为一行。

    在时序上每凑齐 DETECTOR_CHANNELS x WAVELENGTH_CHANNELS 即输出一行；
    与具体交错顺序弱相关，仍以 payload 中的 DetectorId/Wavelength 为准。
    """
    required_keys = {(det.code, wl.code) for det in DETECTOR_CHANNELS for wl in WAVELENGTH_CHANNELS}
    pending: dict[tuple[int, int], pd.Series] = {}
    rows: list[dict[str, float]] = []

    def flush_pending() -> None:
        if set(pending.keys()) != required_keys:
            return
        times = [float(row["Time (s)"]) for row in pending.values()]
        out: dict[str, float] = {"Time (s)": float(np.mean(times))}
        for det in DETECTOR_CHANNELS:
            for wl in WAVELENGTH_CHANNELS:
                out[intensity_column_for(det.name, wl.emitter_nm)] = float(
                    pending[(det.code, wl.code)][signal_col]
                )
        rows.append(out)
        pending.clear()

    for _, row in df.reset_index(drop=True).iterrows():
        code = int(row[mode_col])
        detector_code = int(row["DetectorId"])
        if wavelength_channel_by_code(code) is None or detector_channel_by_code(detector_code) is None:
            continue
        key = (detector_code, code)
        if key in pending:
            # 周期未凑齐却遇到同一 detector/wavelength，说明中间缺帧；重启当前周期。
            pending.clear()
        pending[key] = row
        if len(pending) == len(required_keys):
            flush_pending()

    return pd.DataFrame(rows)


def stack_intensities_for_detector(df: pd.DataFrame, detector) -> np.ndarray:
    """按 WAVELENGTH_CHANNELS 顺序堆叠为 (n_wavelengths, n_timepoints)。"""
    return np.vstack(
        [df[col].to_numpy(dtype=float) for col in intensity_columns_for_detector(detector)]
    )


def short_separation_regression(
    long_od: np.ndarray,
    short_od: np.ndarray,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """
    用短距 OD 的波动成分作为浅层参考，逐波长回归校正长距 OD。

    返回:
    - corrected_long_od: 校正后的长距 OD，形状与 long_od 相同
    - beta: 每个波长的浅层耦合系数，形状为 (n_wavelengths,)
    """
    if long_od.shape != short_od.shape:
        raise ValueError("long_od and short_od must have the same shape.")
    if long_od.ndim != 2:
        raise ValueError("OD matrices must be shaped as (n_wavelengths, n_timepoints).")

    short_centered = short_od - np.mean(short_od, axis=1, keepdims=True)
    long_centered = long_od - np.mean(long_od, axis=1, keepdims=True)
    denominator = np.sum(short_centered * short_centered, axis=1)
    numerator = np.sum(long_centered * short_centered, axis=1)
    beta = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=float),
        where=denominator > eps,
    )
    corrected = long_od - beta[:, np.newaxis] * short_centered
    return corrected, beta


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
    target_fs: float = BP_TARGET_FS_HZ,
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


def _hemoglobin_extinctions(wavelengths: list[float], table: str = "wray") -> np.ndarray:
    """读取 nirsimple 内置 HbO/HbR 消光系数表，并插值到五个波长。"""
    ex_path = os.path.join(os.path.dirname(nsp.__file__), "tables", f"{table}.csv")
    ex_df = pd.read_csv(ex_path)
    table_wls = ex_df["lambda"].to_numpy(dtype=float)
    hbo = ex_df["hbo"].to_numpy(dtype=float)
    hbr = ex_df["hbr"].to_numpy(dtype=float)
    return np.column_stack(
        [
            np.interp(wavelengths, table_wls, hbo),
            np.interp(wavelengths, table_wls, hbr),
        ]
    )


def _extinction_matrix(wavelengths: list[float], table: str = "wray") -> np.ndarray:
    """
    构造 HbO/HbR/Cyt 三色团消光系数矩阵。

    Cyt 使用 UCL-NIR-Spectra 的 cytochrome oxidase difference extinction
    spectrum，原始单位 OD / cm / mM；这里转为 OD / cm / M，与 HbO/HbR 一致。
    """
    hb_ex = _hemoglobin_extinctions(wavelengths, table)
    cyt = np.asarray([CYT_DIFFERENCE_EXTINCTION[float(wl)] for wl in wavelengths], dtype=float) * 1000.0
    return np.column_stack([hb_ex, cyt])


def _cyt_extinction_vector(wavelengths: list[float]) -> np.ndarray:
    """返回 Cyt/oxCCO 差分消光系数向量，单位与 HbO/HbR 保持一致。"""
    return np.asarray([CYT_DIFFERENCE_EXTINCTION[float(wl)] for wl in wavelengths], dtype=float) * 1000.0


def generalized_mbll_hb(
    delta_od: np.ndarray,
    wavelengths: list[float],
    dpfs: list[float],
    distance_cm: float,
    table: str = "wray",
) -> np.ndarray:
    """用五波长最小二乘只反演 HbO/HbR，返回 (2, n_timepoints)。"""
    hb_ex = _hemoglobin_extinctions(wavelengths, table)
    pathlength = np.asarray(dpfs, dtype=float) * float(distance_cm)
    a_hb = hb_ex * pathlength[:, np.newaxis]
    return np.linalg.pinv(a_hb) @ delta_od


def generalized_mbll(
    delta_od: np.ndarray,
    wavelengths: list[float],
    dpfs: list[float],
    distance_cm: float,
    table: str = "wray",
) -> np.ndarray:
    """
    分步广义 MBLL：先稳定反演 HbO/HbR，再从 Hb 残差中估计 Cyt。

    返回行顺序仍为 (HbO, HbR, Cyt)，以保持 processed_output.csv 列名兼容。
    """
    hb_ex = _hemoglobin_extinctions(wavelengths, table)
    cyt_ex = _cyt_extinction_vector(wavelengths)
    pathlength = np.asarray(dpfs, dtype=float) * float(distance_cm)
    a_hb = hb_ex * pathlength[:, np.newaxis]
    a_cyt = (cyt_ex * pathlength)[:, np.newaxis]

    hb_delta_c = np.linalg.pinv(a_hb) @ delta_od
    hb_fit_od = a_hb @ hb_delta_c
    residual_od = delta_od - hb_fit_od

    # 只使用 Cyt 光谱中不能被 HbO/HbR 解释的部分，降低 cyt 的不确定性对 HbO/HbR 的污染。
    hb_projection = a_hb @ np.linalg.pinv(a_hb)
    cyt_residual_basis = a_cyt - hb_projection @ a_cyt
    cyt_delta_c = np.linalg.pinv(cyt_residual_basis) @ residual_od
    return np.vstack([hb_delta_c, cyt_delta_c])


def process_csv_dataset(
    input_csv: str,
    output_csv: str,
    age: int = MBLL_DEFAULT_AGE,
    molar_ext_coeff_table: str = "wray",
    bp_low: float = BP_LOW_HZ,
    bp_high: float = BP_HIGH_HZ,
    bp_order: int = BP_ORDER,
) -> None:
    """从聚合后的 CSV 计算最终的 HbO/HbR/Cyt。"""
    df = pd.read_csv(input_csv)
    required_cols = ["Time (s)", *INTENSITY_COLUMNS]
    if df.empty or any(col not in df.columns for col in required_cols):
        print("Insufficient or invalid interleaved data for processing.")
        return

    times = df["Time (s)"].to_numpy(dtype=float)
    if len(times) < 2:
        print("Need at least two time samples to run MBLL.")
        return

    dt = np.mean(np.diff(times))
    fs = 1.0 / dt if dt > 0 else 1.0

    wavelengths = [ch.mbll_nm for ch in WAVELENGTH_CHANNELS]
    dpfs = [nsp.get_dpf(wl, age) for wl in wavelengths]
    outputs: list[np.ndarray] = []
    processed_headers: list[str] = []
    table_data = []
    detector_delta_od: dict[str, np.ndarray] = {}

    for detector in DETECTOR_CHANNELS:
        samples = stack_intensities_for_detector(df, detector)
        delta_od = nsp.intensities_to_od_changes(samples)
        detector_delta_od[detector.name] = delta_od
        delta_od_filt = smart_bandpass(delta_od, fs, lowcut=bp_low, highcut=bp_high, order=bp_order)
        delta_c = generalized_mbll(
            delta_od_filt,
            wavelengths,
            dpfs,
            detector.distance_cm,
            table=molar_ext_coeff_table,
        )
        outputs.append(delta_c)
        for chrom_idx, chrom_type in enumerate(CHROMOPHORE_TYPES):
            processed_headers.append(f"{detector.name}_{chrom_type}")
            table_data.append([detector.name, chrom_type, f"{delta_c[chrom_idx, -1]:.4e}"])

    if len(DETECTOR_CHANNELS) >= 2:
        short_detector = min(DETECTOR_CHANNELS, key=lambda det: det.distance_cm)
        long_detector = max(DETECTOR_CHANNELS, key=lambda det: det.distance_cm)
        if short_detector.name != long_detector.name:
            corrected_od, beta = short_separation_regression(
                detector_delta_od[long_detector.name],
                detector_delta_od[short_detector.name],
            )
            # 对校正后的 OD 数据做稳健带通。
            corrected_od_filt = smart_bandpass(
                corrected_od,
                fs,
                lowcut=bp_low,
                highcut=bp_high,
                order=bp_order,
            )
            corrected_c = generalized_mbll(
                corrected_od_filt,
                wavelengths,
                dpfs,
                long_detector.distance_cm,
                table=molar_ext_coeff_table,
            )
            outputs.append(corrected_c)
            corrected_prefix = f"{long_detector.name}_ssr"
            for chrom_idx, chrom_type in enumerate(CHROMOPHORE_TYPES):
                processed_headers.append(f"{corrected_prefix}_{chrom_type}")
                table_data.append([corrected_prefix, chrom_type, f"{corrected_c[chrom_idx, -1]:.4e}"])
            beta_text = ", ".join(
                f"{wl.emitter_nm:g}nm={wl_beta:.4g}" for wl, wl_beta in zip(WAVELENGTH_CHANNELS, beta)
            )
            print(
                f"Short-separation regression: {short_detector.name} ({short_detector.distance_cm:g} cm) "
                f"used as shallow reference for {long_detector.name} ({long_detector.distance_cm:g} cm)."
            )
            print(f"SSR beta by wavelength: {beta_text}")

    delta_c_all = np.vstack(outputs)
    header_out = ["Time"] + processed_headers

    with open(output_csv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(header_out)
        n_cols = min(len(times), delta_c_all.shape[1])
        for col_idx in range(n_cols):
            writer.writerow([times[col_idx]] + list(delta_c_all[:, col_idx]))

    print(f"Post-processing complete. Output saved to '{output_csv}'.")
    print("\nExample: Processed concentrations at the final time sample (Cyt uses UCL extinction):")
    print(tabulate(table_data, headers=["Channel", "Type", "Concentration"]))


def capture_data(
    csv_filename: str = RAW_OUTPUT_CSV,
    stop_on_enter: bool = True,
    duration_seconds: float | None = None,
    intensity_ma: int = DEFAULT_INTENSITY_MA,
) -> None:
    """
    采集双接收源五波长原始数据并写 CSV。

    每收到一帧 0x02 数据帧就：
    1. 解析波长/传感器/采样值
    2. 记录到 all_groups.csv
    """
    ser = open_serial() # 打开串口
    reader = FrameReader(ser) # 创建帧读取器
    stop_event = threading.Event() 
    listener = _start_enter_listener(stop_event) if stop_on_enter else None # 创建监听器

    try:
        ser.reset_input_buffer() # 清空串口接收缓冲区里的残留数据
        first_frame = start_stream_until_data(ser, reader, intensity_ma)

        with open(csv_filename, mode="w", newline="", encoding="utf-8") as csvfile:
            writer = csv.writer(csvfile)
            # 原始采样窄表：时间 + 接收源号 + 通道名 + 波长编号 + 数值
            writer.writerow(["Time (s)", "DetectorId", "Channel", "Wavelength", "Value"])

            print("Starting dual-detector five-wavelength raw ADC logging (seconds elapsed)...")
            start_time = time.time()
            last_wait_log = start_time
            rows_written = 0

            while True:
                if stop_event.is_set():
                    print("Stop requested by user.")
                    break
                if duration_seconds is not None and (time.time() - start_time) >= duration_seconds:
                    print("Capture duration reached.")
                    break

                if first_frame is not None:
                    frame = first_frame
                    first_frame = None
                else:
                    frame = reader.read_frame(timeout_seconds=TIMEOUT)
                if frame is None:
                    now = time.time()
                    if now - last_wait_log >= 1.0:
                        print(
                            f"Waiting for data frames... elapsed={now - start_time:.1f}s "
                            f"rows={rows_written}"
                        )
                        last_wait_log = now
                    continue
                if frame.frame_type != 0x02:
                    print(f"Ignored frame type: 0x{frame.frame_type:02X}")
                    continue

                sample = parse_data_frame(frame)
                elapsed_time = round(time.time() - start_time, 6)
                writer.writerow(
                    [
                        elapsed_time,
                        sample.detector_code,
                        sample.channel_name or "",
                        sample.wavelength_code,
                        sample.value,
                    ]
                )
                rows_written += 1
                csvfile.flush()
                print(
                    f"{elapsed_time:.3f}s - value={sample.value} "
                    f"wl={int(sample.wavelength_code)} detector={int(sample.detector_code)} "
                    f"channel={sample.channel_name or 'unknown'}"
                )
    finally:
        try:
            stop_command = build_command_frame(False, intensity_ma)
            print(f"TX stop command: {frame_to_hex(stop_command)}")
            send_frame_with_ack(ser, reader, stop_command)
        except Exception:
            pass
        ser.close()
        if listener is not None and listener.is_alive():
            stop_event.set()


# 三个结果表统一放在此目录下，每次运行占一个以时间命名的子文件夹
RESULT_TABLE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result_table")


def run_pipeline() -> None:
    """一键跑完整链路：采集 -> 预处理 -> 配对 -> MBLL -> CSV 输出。三个表写入 result_table/<时间>/。"""
    os.makedirs(RESULT_TABLE_DIR, exist_ok=True)
    run_id = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_dir = os.path.join(RESULT_TABLE_DIR, run_id)
    os.makedirs(output_dir, exist_ok=True)
    raw_path = os.path.join(output_dir, "all_groups.csv")
    interleaved_path = os.path.join(output_dir, "interleaved_output.csv")
    processed_path = os.path.join(output_dir, "processed_output.csv")
    print(f"本次结果将保存到: {output_dir}")

    # 采集数据
    capture_data(csv_filename=raw_path, stop_on_enter=True)

    # 读取数据
    df = pd.read_csv(raw_path)
    if df.empty or len(df) < 2:
        print("No enough raw rows captured; skipping processing.")
        return

    # 过滤数据
    n_raw = len(df)
    df = df[df["Wavelength"] != WAVELENGTH_OFF_CODE].reset_index(drop=True) # 去掉未点亮的数据
    dropped = n_raw - len(df) # 计算去掉的数据量
    if dropped:
        print(f"Dropped {dropped} raw row(s) with Wavelength=OFF (0x00); not used for dual-wavelength pairing.")

    if df.empty or len(df) < 2:
        print("No enough non-OFF samples after filtering; skipping processing.")
        return

    dt = df["Time (s)"].diff().mean() # 计算时间戳的差值的平均值
    fs = 1.0 / dt if pd.notna(dt) and dt > 0 else 1.0 # 计算采样率
    print(f"采样率: {fs} Hz")

    # # 阈值截断
    # df = threshold_filter(df)

    # 低通滤波
    df = butter_lowpass_filter(df=df, cutoff_hz=1.0, fs=fs, order=4)

    # 分段 RMS
    df = sliding_window_rms(df=df)

    # 多波长周期聚合（列名由 config.WAVELENGTH_CHANNELS 决定）
    final_df = aggregate_wavelength_cycles(df, mode_col="Wavelength")
    if final_df.empty:
        print("No complete wavelength cycles were formed; skipping MBLL.")
        return

    # 按配对后的真实平均间隔重建等间隔时间戳，避免采集抖动影响后续基于 fs 的滤波。
    pair_times = final_df["Time (s)"].to_numpy(dtype=float)
    pair_dt = np.diff(pair_times)
    valid_pair_dt = pair_dt[np.isfinite(pair_dt) & (pair_dt > 0)]
    increment = float(np.mean(valid_pair_dt)) if valid_pair_dt.size else 0.001
    if "Time (s)" in final_df.columns:
        final_df = final_df.drop(columns=["Time (s)"])
    final_df.insert(0, "Time (s)", [i * increment for i in range(len(final_df))])

    # 统一保留6位小数，减少浮点表示伪差。
    final_df["Time (s)"] = final_df["Time (s)"].round(6)

    # 写入文件
    final_df.to_csv(interleaved_path, index=False)
    print(final_df.head(20))

    # 计算最终的 HbO/HbR/Cyt
    process_csv_dataset(interleaved_path, processed_path)


if __name__ == "__main__":
    run_pipeline()

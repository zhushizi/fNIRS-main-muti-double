"""
双接收源五波长 fNIRS 采集与处理主脚本（HbO/HbR 版）。

完整链路如下：
1. 通过 26B 协议采集 S1-D1 / S1-D2 五波长原始数据
2. 写入 all_groups.csv
3. 对每个接收源/波长组合做低通/RMS 预处理
4. 聚合成 2 接收源 x 5 波长光强宽表，写入 interleaved_output_hb.csv
5. 计算 OD -> 短距回归校正 -> 五波长 MBLL -> HbO/HbR
6. 输出 processed_output_hb.csv

本版本不计算 Cyt/oxCCO，其余处理链路与 fNIRS_processing.py 保持一致。
"""

from __future__ import annotations

import csv
import os
from datetime import datetime
from pathlib import Path

import nirsimple.preprocessing as nsp
import numpy as np
import pandas as pd
from tabulate import tabulate

from config import (
    BP_HIGH_HZ,
    BP_LOW_HZ,
    BP_ORDER,
    DETECTOR_CHANNELS,
    INTENSITY_COLUMNS,
    MBLL_DEFAULT_AGE,
    WAVELENGTH_CHANNELS,
    WAVELENGTH_OFF_CODE,
)
from fNIRS_processing import (
    aggregate_wavelength_cycles,
    butter_lowpass_filter,
    capture_data,
    short_separation_regression,
    sliding_window_rms,
    smart_bandpass,
    stack_intensities_for_detector,
)

HEMOGLOBIN_TYPES: tuple[str, str] = ("hbo", "hbr")
RESULT_TABLE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result_table")


def _hemoglobin_extinctions(wavelengths: list[float], table: str = "wray") -> np.ndarray:
    """读取 nirsimple 内置 HbO/HbR 消光系数表，并插值到五个波长。"""
    ex_path = Path(nsp.__file__).resolve().parent / "tables" / f"{table}.csv"
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


def generalized_mbll_hb(
    delta_od: np.ndarray,
    wavelengths: list[float],
    dpfs: list[float],
    distance_cm: float,
    table: str = "wray",
) -> np.ndarray:
    """用五波长最小二乘只反演 HbO/HbR，返回 (2, n_timepoints)。"""
    ex = _hemoglobin_extinctions(wavelengths, table)
    pathlength = np.asarray(dpfs, dtype=float) * float(distance_cm)
    a_matrix = ex * pathlength[:, np.newaxis]
    return np.linalg.pinv(a_matrix) @ delta_od


def process_csv_dataset_hb(
    input_csv: str,
    output_csv: str,
    age: int = MBLL_DEFAULT_AGE,
    molar_ext_coeff_table: str = "wray",
    bp_low: float = BP_LOW_HZ,
    bp_high: float = BP_HIGH_HZ,
    bp_order: int = BP_ORDER,
) -> None:
    """从聚合后的 CSV 计算最终 HbO/HbR，不输出 Cyt。"""
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
        delta_c = generalized_mbll_hb(
            delta_od_filt,
            wavelengths,
            dpfs,
            detector.distance_cm,
            table=molar_ext_coeff_table,
        )
        outputs.append(delta_c)
        for chrom_idx, chrom_type in enumerate(HEMOGLOBIN_TYPES):
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
            corrected_od_filt = smart_bandpass(
                corrected_od,
                fs,
                lowcut=bp_low,
                highcut=bp_high,
                order=bp_order,
            )
            corrected_c = generalized_mbll_hb(
                corrected_od_filt,
                wavelengths,
                dpfs,
                long_detector.distance_cm,
                table=molar_ext_coeff_table,
            )
            outputs.append(corrected_c)
            corrected_prefix = f"{long_detector.name}_ssr"
            for chrom_idx, chrom_type in enumerate(HEMOGLOBIN_TYPES):
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
    print("\nExample: Processed HbO/HbR concentrations at the final time sample:")
    print(tabulate(table_data, headers=["Channel", "Type", "Concentration"]))


def run_pipeline() -> None:
    """一键跑完整链路：采集 -> 预处理 -> 配对 -> HbO/HbR CSV 输出。"""
    os.makedirs(RESULT_TABLE_DIR, exist_ok=True)
    run_id = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_dir = os.path.join(RESULT_TABLE_DIR, run_id)
    os.makedirs(output_dir, exist_ok=True)
    raw_path = os.path.join(output_dir, "all_groups.csv")
    interleaved_path = os.path.join(output_dir, "interleaved_output_hb.csv")
    processed_path = os.path.join(output_dir, "processed_output_hb.csv")
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
    df = df[df["Wavelength"] != WAVELENGTH_OFF_CODE].reset_index(drop=True)
    dropped = n_raw - len(df)
    if dropped:
        print(f"Dropped {dropped} raw row(s) with Wavelength=OFF (0x00); not used for wavelength pairing.")

    if df.empty or len(df) < 2:
        print("No enough non-OFF samples after filtering; skipping processing.")
        return

    dt = df["Time (s)"].diff().mean()
    fs = 1.0 / dt if pd.notna(dt) and dt > 0 else 1.0
    print(f"采样率: {fs} Hz")

    # 低通滤波，保持与 fNIRS_processing.py 一致。
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

    # 写入 Hb 版中间表
    final_df.to_csv(interleaved_path, index=False)
    print(final_df.head(20))

    # 计算最终 HbO/HbR
    process_csv_dataset_hb(interleaved_path, processed_path)


if __name__ == "__main__":
    run_pipeline()

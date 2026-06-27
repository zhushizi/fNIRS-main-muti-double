"""
从 all_groups.csv 离线计算五波长 HbO/HbR（不计算 Cyt）。

默认输入:
    software/result_table/2026-06-24_14-44-47/all_groups.csv

输出默认写在同一目录:
    interleaved_output_hb.csv
    processed_output_hb.csv

处理流程与 fNIRS_processing.py 保持一致:
    all_groups.csv -> 去 OFF -> raw 低通 -> 分段 RMS -> 五波长周期聚合
    -> OD -> 带通 -> 五波长 HbO/HbR MBLL -> SSR HbO/HbR
用于确认cyt的引入是否影响hbo和hbr的准确性
"""

from __future__ import annotations

import argparse
import csv
import os
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
    short_separation_regression,
    sliding_window_rms,
    smart_bandpass,
    stack_intensities_for_detector,
)

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "result_table" / "2026-06-27_10-30-58" / "all_groups.csv"
HEMOGLOBIN_TYPES: tuple[str, str] = ("hbo", "hbr")


def _hemoglobin_extinctions(wavelengths: list[float], table: str = "wray") -> np.ndarray:
    """读取 nirsimple 内置 HbO/HbR 消光系数表，并插值到当前五个波长。"""
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


def process_interleaved_hb(
    input_csv: Path,
    output_csv: Path,
    age: int = MBLL_DEFAULT_AGE,
    molar_ext_coeff_table: str = "wray",
    bp_low: float = BP_LOW_HZ,
    bp_high: float = BP_HIGH_HZ,
    bp_order: int = BP_ORDER,
) -> None:
    """从聚合后的 interleaved CSV 计算 HbO/HbR，不输出 Cyt。"""
    df = pd.read_csv(input_csv)
    required_cols = ["Time (s)", *INTENSITY_COLUMNS]
    if df.empty or any(col not in df.columns for col in required_cols):
        raise ValueError("Insufficient or invalid interleaved data for processing.")

    times = df["Time (s)"].to_numpy(dtype=float)
    if len(times) < 2:
        raise ValueError("Need at least two time samples to run MBLL.")

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
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["Time", *processed_headers])
        n_cols = min(len(times), delta_c_all.shape[1])
        for col_idx in range(n_cols):
            writer.writerow([times[col_idx], *list(delta_c_all[:, col_idx])])

    print(f"Post-processing complete. Output saved to '{output_csv}'.")
    print("\nExample: Processed HbO/HbR concentrations at the final time sample:")
    print(tabulate(table_data, headers=["Channel", "Type", "Concentration"]))


def process_all_groups_hb(
    input_csv: Path,
    interleaved_csv: Path | None = None,
    output_csv: Path | None = None,
) -> tuple[Path, Path]:
    """直接从 all_groups.csv 离线生成 HbO/HbR 结果表。"""
    if interleaved_csv is None:
        interleaved_csv = input_csv.with_name("interleaved_output_hb.csv")
    if output_csv is None:
        output_csv = input_csv.with_name("processed_output_hb.csv")

    df = pd.read_csv(input_csv)
    if df.empty or len(df) < 2:
        raise ValueError("No enough raw rows for processing.")

    n_raw = len(df)
    df = df[df["Wavelength"] != WAVELENGTH_OFF_CODE].reset_index(drop=True)
    dropped = n_raw - len(df)
    if dropped:
        print(f"Dropped {dropped} raw row(s) with Wavelength=OFF (0x00); not used for wavelength pairing.")

    if df.empty or len(df) < 2:
        raise ValueError("No enough non-OFF samples after filtering.")

    dt = df["Time (s)"].diff().mean()
    fs = 1.0 / dt if pd.notna(dt) and dt > 0 else 1.0
    print(f"采样率: {fs} Hz")

    df = butter_lowpass_filter(df=df, cutoff_hz=1.0, fs=fs, order=4)
    df = sliding_window_rms(df=df)

    final_df = aggregate_wavelength_cycles(df, mode_col="Wavelength")
    if final_df.empty:
        raise ValueError("No complete wavelength cycles were formed; skipping MBLL.")

    pair_times = final_df["Time (s)"].to_numpy(dtype=float)
    pair_dt = np.diff(pair_times)
    valid_pair_dt = pair_dt[np.isfinite(pair_dt) & (pair_dt > 0)]
    increment = float(np.mean(valid_pair_dt)) if valid_pair_dt.size else 0.001
    if "Time (s)" in final_df.columns:
        final_df = final_df.drop(columns=["Time (s)"])
    final_df.insert(0, "Time (s)", [i * increment for i in range(len(final_df))])
    final_df["Time (s)"] = final_df["Time (s)"].round(6)

    interleaved_csv.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(interleaved_csv, index=False)
    print(f"Interleaved Hb input saved to '{interleaved_csv}'.")
    print(final_df.head(20))

    process_interleaved_hb(interleaved_csv, output_csv)
    return interleaved_csv, output_csv


def main() -> None:
    parser = argparse.ArgumentParser(description="从 all_groups.csv 离线计算五波长 HbO/HbR，不计算 Cyt。")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="输入 all_groups.csv 路径")
    parser.add_argument(
        "--interleaved",
        type=Path,
        default=None,
        help="输出中间 interleaved CSV，默认写到输入同目录 interleaved_output_hb.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="输出 HbO/HbR CSV，默认写到输入同目录 processed_output_hb.csv",
    )
    args = parser.parse_args()

    process_all_groups_hb(args.input, interleaved_csv=args.interleaved, output_csv=args.output)


if __name__ == "__main__":
    main()

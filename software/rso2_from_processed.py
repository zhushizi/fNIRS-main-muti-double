from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path(__file__).resolve().parent / "result_table" / "2026-06-29_14-43-08" / "processed_output.csv"

# (显示名, HbO 列名, HbR 列名)
RSO2_CHANNEL_SPECS: tuple[tuple[str, str, str], ...] = (
    ("S1_D1", "S1_D1_hbo", "S1_D1_hbr"),
    ("S1_D2", "S1_D2_hbo", "S1_D2_hbr"),
    ("S1_D1_ssr", "S1_D1_ssr_hbo", "S1_D1_ssr_hbr"),
)


@dataclass(frozen=True)
class Rso2Series:
    label: str
    times: np.ndarray
    rso2_pct: np.ndarray
    hbo_abs_uM: np.ndarray
    hbr_abs_uM: np.ndarray
    baseline_mean_rso2_pct: float


def _resolve_time_column(columns: list[str]) -> str:
    for name in ("Time", "Time (s)"):
        if name in columns:
            return name
    raise ValueError(f"未找到时间列（期望 Time 或 Time (s)），现有列: {columns}")


def _require_columns(df: pd.DataFrame, *col_names: str) -> None:
    missing = [col for col in col_names if col not in df.columns]
    if missing:
        raise ValueError(f"输入 CSV 缺少列: {missing}。现有列: {list(df.columns)}")


def compute_rso2_series(
    times: np.ndarray,
    delta_hbo: np.ndarray,
    delta_hbr: np.ndarray,
    *,
    label: str,
    baseline_start_s: float,
    baseline_end_s: float,
    baseline_hbt_uM: float,
    baseline_rso2_pct: float,
) -> Rso2Series:
    """由 ΔHbO/ΔHbR 与固定基线假设估算 rSO2 时间序列。"""
    times = np.asarray(times, dtype=float)
    delta_hbo = np.asarray(delta_hbo, dtype=float)
    delta_hbr = np.asarray(delta_hbr, dtype=float)

    baseline_mask = (times >= baseline_start_s) & (times <= baseline_end_s)
    if not np.any(baseline_mask):
        raise ValueError(
            f"[{label}] 在 {baseline_start_s:.1f}s ~ {baseline_end_s:.1f}s 内没有基线样本，"
            "请调整基线时间窗口。"
        )

    baseline_hbt_M = baseline_hbt_uM * 1e-6
    baseline_hbo_abs_M = baseline_hbt_M * (baseline_rso2_pct / 100.0)
    baseline_hbr_abs_M = baseline_hbt_M - baseline_hbo_abs_M

    baseline_delta_hbo_mean = float(np.mean(delta_hbo[baseline_mask]))
    baseline_delta_hbr_mean = float(np.mean(delta_hbr[baseline_mask]))

    hbo_abs_M = baseline_hbo_abs_M + (delta_hbo - baseline_delta_hbo_mean)
    hbr_abs_M = baseline_hbr_abs_M + (delta_hbr - baseline_delta_hbr_mean)
    hbo_abs_M = np.maximum(hbo_abs_M, 0.0)
    hbr_abs_M = np.maximum(hbr_abs_M, 0.0)
    hbt_abs_M = hbo_abs_M + hbr_abs_M

    with np.errstate(divide="ignore", invalid="ignore"):
        rso2_pct = np.where(hbt_abs_M > 0, 100.0 * hbo_abs_M / hbt_abs_M, np.nan)

    return Rso2Series(
        label=label,
        times=times,
        rso2_pct=rso2_pct,
        hbo_abs_uM=hbo_abs_M * 1e6,
        hbr_abs_uM=hbr_abs_M * 1e6,
        baseline_mean_rso2_pct=float(np.nanmean(rso2_pct[baseline_mask])),
    )


def compute_rso2_table(
    input_csv: Path,
    output_csv: Path | None = None,
    baseline_start_s: float = 30.0,
    baseline_end_s: float = 60.0,
    baseline_hbt_uM: float = 80.0,
    baseline_rso2_pct: float = 65.0,
    channels: tuple[tuple[str, str, str], ...] = RSO2_CHANNEL_SPECS,
) -> tuple[Path, list[Rso2Series]]:
    df = pd.read_csv(input_csv)
    if df.empty:
        raise ValueError("输入 CSV 为空。")

    time_col = _resolve_time_column(list(df.columns))
    times = df[time_col].to_numpy(dtype=float)

    series_list: list[Rso2Series] = []
    for label, hbo_col, hbr_col in channels:
        _require_columns(df, hbo_col, hbr_col)
        series_list.append(
            compute_rso2_series(
                times,
                df[hbo_col].to_numpy(dtype=float),
                df[hbr_col].to_numpy(dtype=float),
                label=label,
                baseline_start_s=baseline_start_s,
                baseline_end_s=baseline_end_s,
                baseline_hbt_uM=baseline_hbt_uM,
                baseline_rso2_pct=baseline_rso2_pct,
            )
        )

    result_mask = times > baseline_end_s
    if not np.any(result_mask):
        raise ValueError(f"在 {baseline_end_s:.1f}s 之后没有可输出的 rSO2 数据。")

    out_data: dict[str, np.ndarray] = {"Time": times[result_mask]}
    for item in series_list:
        out_data[f"{item.label}_rso2_pct"] = item.rso2_pct[result_mask]

    out_df = pd.DataFrame(out_data)

    if output_csv is None:
        output_csv = input_csv.with_name("rso2_output.csv")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_csv, index=False)

    print(f"输入文件: {input_csv}")
    print(f"输出文件: {output_csv}")
    print(
        "基线设定: "
        f"HbT={baseline_hbt_uM:.2f} uM, "
        f"rSO2={baseline_rso2_pct:.2f}%, "
        f"窗口={baseline_start_s:.1f}s~{baseline_end_s:.1f}s"
    )
    for item in series_list:
        post = item.rso2_pct[result_mask]
        baseline_mask = (item.times >= baseline_start_s) & (item.times <= baseline_end_s)
        print(f"\n[{item.label}]")
        print(
            "  基线窗内: "
            f"HbO={np.mean(item.hbo_abs_uM[baseline_mask]):.3f} uM, "
            f"HbR={np.mean(item.hbr_abs_uM[baseline_mask]):.3f} uM, "
            f"rSO2={item.baseline_mean_rso2_pct:.3f}%"
        )
        print(
            f"  {baseline_end_s:.1f}s 后: "
            f"min={np.nanmin(post):.3f}%, "
            f"max={np.nanmax(post):.3f}%, "
            f"mean={np.nanmean(post):.3f}%"
        )

    return output_csv, series_list


def plot_rso2(
    series_list: list[Rso2Series],
    *,
    baseline_start_s: float,
    baseline_end_s: float,
    output_png: Path | None = None,
    show: bool = True,
) -> None:
    """临时绘图：同图对比 S1_D1 / S1_D2 / S1_D1_ssr 的 rSO2。"""
    if not series_list:
        return

    colors = {
        "S1_D1": "#c0392b",
        "S1_D2": "#2980b9",
        "S1_D1_ssr": "#27ae60",
    }

    fig, ax = plt.subplots(figsize=(12, 5))
    for item in series_list:
        ax.plot(
            item.times,
            item.rso2_pct,
            label=item.label,
            linewidth=1.4,
            color=colors.get(item.label),
        )

    ax.axvspan(baseline_start_s, baseline_end_s, color="#f1c40f", alpha=0.25, label="baseline window")
    ax.axvline(baseline_end_s, color="#7f8c8d", linestyle="--", linewidth=1.0, label="output starts")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("rSO2 (%)")
    ax.set_title("Regional rSO2: S1_D1 / S1_D2 / S1_D1_ssr")
    ax.grid(alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()

    if output_png is not None:
        output_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_png, dpi=150)
        print(f"图像已保存: {output_png}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="根据 processed_output.csv 计算 S1_D1 / S1_D2 / S1_D1_ssr 的 rSO2，并可临时绘图。"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="输入 processed_output.csv 路径",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="输出 rSO2 CSV 路径，默认与输入同目录下的 rso2_output.csv",
    )
    parser.add_argument(
        "--plot-png",
        type=Path,
        default=None,
        help="可选：保存 rSO2 对比图 PNG 路径，默认与 CSV 同目录 rso2_plot.png",
    )
    parser.add_argument("--baseline-start", type=float, default=30.0, help="基线起始时间（秒）")
    parser.add_argument("--baseline-end", type=float, default=60.0, help="基线结束时间（秒）")
    parser.add_argument("--baseline-hbt", type=float, default=80.0, help="固定基线 HbT（uM）")
    parser.add_argument("--baseline-rso2", type=float, default=65.0, help="固定基线 rSO2（%%）")
    parser.add_argument("--no-plot", action="store_true", help="不弹出 matplotlib 窗口")
    args = parser.parse_args()

    output_csv, series_list = compute_rso2_table(
        input_csv=args.input,
        output_csv=args.output,
        baseline_start_s=args.baseline_start,
        baseline_end_s=args.baseline_end,
        baseline_hbt_uM=args.baseline_hbt,
        baseline_rso2_pct=args.baseline_rso2,
    )

    plot_png = args.plot_png
    if plot_png is None and not args.no_plot:
        plot_png = output_csv.with_name("rso2_plot.png")

    if not args.no_plot or plot_png is not None:
        plot_rso2(
            series_list,
            baseline_start_s=args.baseline_start,
            baseline_end_s=args.baseline_end,
            output_png=plot_png,
            show=not args.no_plot,
        )


if __name__ == "__main__":
    main()

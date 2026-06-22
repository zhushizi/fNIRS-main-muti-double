# 双接收源五波长 fNIRS 处理链路说明

本文档说明当前 `software/fNIRS_processing.py` 的单通道处理链路。  
当前版本已不再使用旧的 64B 多组协议，而是采用 **26B 带帧头/校验/ACK** 的双接收源五波长协议。

---

## 1. 总览

```text
26B 串口数据帧
  -> 帧头同步 / 长度校验 / checksum 校验
  -> 回 0x03 ACK
  -> capture_data -> all_groups.csv
  -> threshold_filter
  -> butter_lowpass_filter
  -> sliding_window_rms
  -> aggregate_wavelength_cycles
  -> interleaved_output.csv
  -> process_csv_dataset
      -> 2 接收源 x 5 波长矩阵
      -> intensities_to_od_changes
      -> smart_bandpass
      -> generalized_mbll
  -> processed_output.csv
```

---

## 2. 串口协议层

### 2.1 串口参数

- 波特率：`115200`
- 帧头：`0x55 0xAA`
- 固定帧长：`26B`
- 类型：
  - `0x01` 指令帧
  - `0x02` 数据帧
  - `0x03` 应答帧

### 2.2 校验规则

- 校验范围：从 **长度字节** 到 **payload 最后一个字节**
- 校验方式：**累加和低字节**

### 2.3 ACK / 重发

- ACK 超时：`10ms`
- 最多重发：`2` 次
- 上位机收到 `0x02` 数据帧后会立即回 `0x03` ACK

### 2.4 指令帧 `0x01`

payload 共 21 字节，目前使用：

| 字节 | 含义 |
|------|------|
| 0 | 启停：`0x00` 停止，`0x01` 启动 |
| 1 | 光强单字节：`0x00~0xFF` |
| 2 | 保留，当前填 0 |
| 3-20 | 保留，当前填 0 |

默认光强：`0xC8`（200）

### 2.5 数据帧 `0x02`

payload 共 21 字节，目前使用：

| 字节 | 含义 |
|------|------|
| 0 | 波长编号：`0x00=未点亮`，`0x01=850nm`，`0x02=810nm`，`0x03=770nm`，`0x04=730nm`，`0x05=700nm` |
| 1 | 接收源编号：`0x01=PD1/S1_D1`，`0x02=PD2/S1_D2` |
| 2-5 | 采集值，4 字节无符号大端 |
| 6-19 | 保留 |
| 20 | 保留位，当前固定 `0x00` |

单个数据帧只携带 **一个接收源、一个波长下的一次采样值**。

说明：

- 当前上位机里，`payload byte0` 表示 LED 波长，`payload byte1` 表示 PD 接收源。
- **`0x00` 行在 `run_pipeline` 读入 `all_groups.csv` 后即丢弃**，不参与 RMS 切段与周期聚合。
- 预期上报周期为每个波长依次 PD1 / PD2，共 `5 x 2 = 10` 帧。

---

## 3. 原始采集层

### 3.1 `capture_data`

`capture_data()` 的职责：

1. 打开串口
2. 发送启动命令帧
3. 按帧头同步读取 26B 数据帧
4. 校验通过后回 ACK
5. 将数据写入 `all_groups.csv`
6. 采集结束后发送停止命令帧

### 3.2 `all_groups.csv` 结构

```text
Time (s),DetectorId,Channel,Wavelength,Value
```

字段说明：

- `Time (s)`：接收该帧时的相对时间
- `DetectorId`：接收源编号，`1=PD1/S1_D1`，`2=PD2/S1_D2`
- `Channel`：通道名
- `Wavelength`：波长编号（流水线中会去掉 `0`）
- `Value`：当前采样值

当前几何包含 `S1_D1=3.0cm` 与 `S1_D2=1.0cm`。

---

## 4. 预处理层

输入：`all_groups.csv`  
输出：`interleaved_output.csv`

### 4.1 阈值滤波 `threshold_filter`

- 针对 `S1_D1` 做阈值过滤
- 越界样本先置空，再插值补齐

### 4.2 低通滤波 `butter_lowpass_filter`

- 对 `S1_D1` 沿时间轴做低通
- 数据太短时自动跳过，避免 `filtfilt` 报错

### 4.3 分段 RMS `sliding_window_rms`

- 依据 `Wavelength` 变化切段
- 每段把 `S1_D1` 用 RMS 代表值回填整段

目的仍然是：

- 用一个稳定值代表同一波长时隙
- 方便后续 660 / 940 配对

### 4.4 周期聚合 `aggregate_wavelength_cycles`

- 将 RMS 后的 `(DetectorId, Wavelength)` 段按周期聚合
- 凑齐 `S1_D1/S1_D2 x 850/810/770/730/700nm` 后输出一行
- 若周期内遇到重复组合，认为缺帧并重启当前周期

### 4.5 `interleaved_output.csv` 结构

```text
Time (s),S1_D1_850,S1_D1_810,S1_D1_770,S1_D1_730,S1_D1_700,S1_D2_850,...
```

每一行表示一组已经对齐好的 2 接收源 x 5 波长样本。

---

## 5. 后处理层

入口函数：`process_csv_dataset`

输入：`interleaved_output.csv`  
输出：`processed_output.csv`

### 5.1 五波长矩阵

每个接收源独立从 `interleaved_output.csv` 中取 5 个波长列，组合为形状 `(5, N)` 的样本矩阵。

### 5.2 `build_channel_info`

当前几何固定为：

- 通道名：`S1_D1`
- 距离：`3.0 cm`（来自 `config.py`）
- 波长：`660 / 940`

### 5.3 OD 转换

```python
delta_od = nsp.intensities_to_od_changes(samples)
```

### 5.4 带通滤波 `smart_bandpass`

- 默认带通：`0.05 ~ 0.1 Hz`
- 样本太短时自动跳过滤波

### 5.5 广义 MBLL

使用五波长最小二乘 / 伪逆求解：

- `S1_D1` 距离 `3.0cm`
- `S1_D2` 距离 `1.0cm`
- 色团输出 `hbo / hbr / cyt`

最终输出：

- `S1_D1_hbo`
- `S1_D1_hbr`
- `S1_D1_cyt`
- `S1_D2_hbo`
- `S1_D2_hbr`
- `S1_D2_cyt`

### 5.6 `processed_output.csv` 结构

```text
Time,S1_D1_hbo,S1_D1_hbr,S1_D1_cyt,S1_D2_hbo,S1_D2_hbr,S1_D2_cyt
```

这仍然保持旧版的输出风格：
第一列时间，后面采用 `{ChannelName}_{Type}` 形式。

注意：当前 `cyt` 使用 `config.CYT_DIFFERENCE_EXTINCTION`，来源为 UCL-NIR-Spectra 的 cytochrome oxidase 差分消光谱（OD / cm / mM），矩阵内会转换为 OD / cm / M 以匹配 HbO/HbR；结果仍建议结合实验标定验证。

---

## 6. 关键配置

位置：`software/config.py`

重点参数：

- `SERIAL_PORT`
- `BAUD_RATE = 115200`
- `ACK_TIMEOUT_SECONDS = 0.01`
- `MAX_RETRIES = 2`
- `DEFAULT_INTENSITY_MA = 200`
- `SOURCE_DETECTOR_DISTANCE_CM = 3.0`

---

## 7. 与旧版的核心差异

| 项目 | 旧版 | 当前版 |
|------|------|--------|
| 协议 | 64B 裸帧 | 26B 带帧头/校验/ACK |
| 通道数 | 8组 x 3 路 | 1 路 `S1_D1` |
| 波长组织 | 由多组 emitter 状态切段 | 单通道按 `Wavelength` 的 660/940 配对 |
| 原始 CSV | `G0...G7` 宽表 | `Time,S1_D1,Wavelength` |
| 处理中间表 | 33 列宽表 | 单通道双波长配对表 |
| 几何 | 短距/长距混合 | 单距离 `3.0 cm` |

---

## 8. 运行建议

1. 先确认 `config.py` 中串口端口正确
2. 用 `python visualizer.py` 验证控制命令与 ACK
3. 用 `python adc_live.py` 查看实时波形
4. 用 `python fNIRS_processing.py` 完成一次采集和处理

如果 `processed_output.csv` 行数太少，优先检查：

- 是否确实收到了 `0x02` 数据帧
- `Wavelength` 是否在 660 / 940 之间切换
- ACK 是否正常返回
- 原始采样是否足够形成多个配对样本

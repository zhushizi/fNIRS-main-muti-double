# fNIRS 双接收源五波长上位机软件

本仓库现在只保留 **上位机软件部分**，并已切换到新的 **26B 双接收源五波长串口协议**。  
软件通过串口与设备通信，接收 `S1_D1` / `S1_D2` 两个接收源的五波长数据（850 / 810 / 770 / 730 / 700nm），完成采集、周期聚合、广义 MBLL 处理，并输出 `hbo / hbr / cyt`。

## 当前协议摘要

- **波特率**：`115200`
- **帧头**：`0x55 0xAA`
- **固定帧长**：`0x1A`（26 字节）
- **帧类型**：
  - `0x01`：指令帧（上位机 -> 下位机）
  - `0x02`：数据帧（下位机 -> 上位机）
  - `0x03`：应答帧（双向 ACK）
- **payload**：`21B`
- **校验**：从长度字节到 payload 末尾的累加和低字节
- **ACK 规则**：10ms 超时，最多重发 2 次

## 当前几何与算法假设

- **物理通道**：`S1_D1` 与 `S1_D2`
- **源探距离**：`S1_D1 = 3.0 cm`，`S1_D2 = 1.0 cm`
- **波长编码**（数据帧 payload byte0）：
  - `0x00`：未点亮
  - `0x01`：850nm
  - `0x02`：810nm
  - `0x03`：770nm
  - `0x04`：730nm
  - `0x05`：700nm
- **接收源编码**（数据帧 payload byte1）：
  - `0x01`：PD1 / `S1_D1`
  - `0x02`：PD2 / `S1_D2`
- **输出风格**：`processed_output.csv` 保持 `Time + {ChannelName}_{Type}` 样式，例如 `S1_D1_hbo`、`S1_D1_hbr`、`S1_D1_cyt`
- **Cyt 说明**：当前 Cyt/oxCCO 使用 UCL-NIR-Spectra 的 cytochrome oxidase 差分消光谱，矩阵内统一到 `OD/cm/M`；结果仍建议结合实验标定验证。

## 快速开始

1. 安装依赖

```bash
cd software
pip install -r requirements.txt
```

2. 配置串口

编辑 `software/config.py`：

- `SERIAL_PORT`
- `BAUD_RATE`
- `DETECTOR_CHANNELS`
- `WAVELENGTH_CHANNELS`
- `CYT_DIFFERENCE_EXTINCTION`

3. 运行方式

- **采集 + 处理**：`python fNIRS_processing.py`
- **实时 ADC 曲线**：`python adc_live.py`
- **控制面板**：`python visualizer.py`
- **演示模式控制面板**：`python visualizer.py demo`

## 输出文件

- `all_groups.csv`
  - 原始窄表采集数据
  - 列：`Time (s), DetectorId, Channel, Wavelength, Value`
- `interleaved_output.csv`
  - 2 接收源 x 5 波长聚合后的中间结果
  - 列：`Time (s), S1_D1_850 ... S1_D2_700`
- `processed_output.csv`
  - 广义 MBLL 输出
  - 列风格：`Time, S1_D1_hbo, S1_D1_hbr, S1_D1_cyt, S1_D2_hbo, ...`

## 主要文件

- `software/config.py`：串口和协议配置
- `software/protocol.py`：26B 协议收发、校验、ACK、重发
- `software/fNIRS_processing.py`：采集、周期聚合、HbO/HbR/Cyt 反演
- `software/adc_live.py`：双接收源五波长实时 ADC 显示
- `software/hbo_hbr_live.py`：实时 HbO/HbR/Cyt 显示
- `software/visualizer.py`：轻量控制面板后端
- `software/index.html`：轻量控制面板前端

更详细的数据链路说明见 `software/fNIRS_processing_pipeline.md`。

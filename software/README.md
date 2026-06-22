# fNIRS Dual-Detector Five-Wavelength Software

This directory contains the host-side software for the **dual-detector five-wavelength** fNIRS workflow.  
It targets the **26-byte framed serial protocol** and assumes one source with two detector channels, `S1_D1` and `S1_D2`, sampled under five wavelengths (`850 / 810 / 770 / 730 / 700nm`).

## Core Features

- **Framed serial protocol**
  - `0x55 0xAA` frame header
  - fixed `26B` frame length
  - checksum validation
  - `0x03` ACK handling with retry support

- **Dual-detector acquisition**
  - one intensity value per data frame
  - wavelength code and detector id included in the payload
  - raw capture stored in `all_groups.csv`

- **Dual-detector five-wavelength processing**
  - threshold filtering
  - low-pass filtering
  - segment RMS by detector/wavelength code
  - 2 detector x 5 wavelength cycle aggregation
  - generalized MBLL for HbO / HbR / Cyt

- **Visualization**
  - `adc_live.py`: live ADC plot for all detector/wavelength combinations
  - `hbo_hbr_live.py`: live HbO / HbR / Cyt plot
  - `adc_animation.py`: replay raw CSV
  - `mBLL_animation.py`: replay processed CSV
  - `visualizer.py`: lightweight control dashboard

## Protocol Summary

### Common Frame Layout

| Field | Size |
|------|------|
| Header | 2 bytes (`0x55 0xAA`) |
| Length | 1 byte (`0x1A`) |
| Type | 1 byte |
| Payload | 21 bytes |
| Checksum | 1 byte |

### Frame Types

- `0x01`: command frame
- `0x02`: data frame
- `0x03`: ACK frame

### Command Payload (`0x01`)

| Byte | Meaning |
|------|---------|
| 0 | start/stop (`0x00` stop, `0x01` start) |
| 1 | intensity (`0x00`~`0xFF`) |
| 2 | reserved, currently `0x00` |
| 3-20 | reserved, currently `0x00` |

### Data Payload (`0x02`)

| Byte | Meaning |
|------|---------|
| 0 | wavelength code (`0x00=off`, `0x01=850nm`, `0x02=810nm`, `0x03=770nm`, `0x04=730nm`, `0x05=700nm`) |
| 1 | detector id (`0x01=PD1/S1_D1`, `0x02=PD2/S1_D2`) |
| 2-5 | sampled value, unsigned 32-bit big-endian |
| 6-19 | reserved |
| 20 | reserved, currently fixed `0x00` |

`fNIRS_processing.py` drops `0x00` wavelength rows before cycle aggregation. The expected cycle is five wavelengths, each followed by PD1 and PD2 samples.

### ACK Handling

- timeout: `10ms`
- retries: `2`

## Configuration

Edit `config.py` before running:

- `SERIAL_PORT`
- `BAUD_RATE`
- `TIMEOUT`
- `DEFAULT_INTENSITY_MA`
- `DETECTOR_CHANNELS`
- `WAVELENGTH_CHANNELS`
- `CYT_DIFFERENCE_EXTINCTION`

## Main Scripts

### `fNIRS_processing.py`

Captures raw frames, writes `all_groups.csv`, creates the 2x5 wavelength `interleaved_output.csv`, and computes HbO/HbR/Cyt in `processed_output.csv`.

### `adc_live.py`

Starts the stream and displays all detector/wavelength values live on one plot.

### `visualizer.py`

Runs a minimal Flask dashboard at `http://127.0.0.1:8050` for:

- start/stop commands
- intensity updates
- latest packet inspection
- ACK status feedback

Run demo mode with:

```bash
python visualizer.py demo
```

## CSV Outputs

### `all_groups.csv`

```text
Time (s),DetectorId,Channel,Wavelength,Value
```

### `interleaved_output.csv`

```text
Time (s),S1_D1_850,S1_D1_810,S1_D1_770,S1_D1_730,S1_D1_700,S1_D2_850,...
```

### `processed_output.csv`

```text
Time,S1_D1_hbo,S1_D1_hbr,S1_D1_cyt,S1_D2_hbo,S1_D2_hbr,S1_D2_cyt
```

`cyt` uses the UCL-NIR-Spectra cytochrome oxidase difference extinction spectrum (`OD / cm / mM`), converted internally to `OD / cm / M` to match HbO/HbR. Experimental calibration is still recommended before interpreting it as a validated absolute concentration.

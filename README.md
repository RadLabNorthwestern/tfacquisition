# TF Acquisition & Plotting

Two Python scripts for measuring and reviewing the RF **Transfer Function (TF)** of
an implanted-lead test fixture using a **Keysight E5063A ENA** vector network
analyzer (VNA). They implement the *TF Measurement Procedure* SOP:

| Script | Purpose |
|--------|---------|
| [`recall_and_acquire_s21.py`](recall_and_acquire_s21.py) | Interactively acquire one averaged S-parameter sweep per z-axis probe position and save each as a native Touchstone `.s2p` file. |
| [`plot_tf_zaxis.py`](plot_tf_zaxis.py) | Read the saved `.s2p` files from a run, plot `\|S21\|` (dB) and phase (deg) versus probe position at one or more target frequencies, and write per-frequency CSVs and the plot into the run's SOP output folders. |

The acquisition script talks to the instrument over VISA/LAN and drives the sweep;
the plot script is offline — it reads the mirrored `.s2p` files and writes the
processed CSVs and overlay plot back into the run's SOP folders.

---

## 1. Requirements

- **Python 3.8+**
- Packages: see [`requirements.txt`](requirements.txt) — `pyvisa`, `numpy`, `matplotlib`
- A **VISA runtime** (Keysight IO Libraries Suite or NI-VISA) for `pyvisa` to reach
  the instrument. *Only needed for acquisition, not for plotting.*
- Network line-of-sight to the E5063A ENA.

```sh
python -m venv .venv
# Windows:
.venv\Scripts\activate
pip install -r requirements.txt
```

---

## 2. Acquisition — `recall_and_acquire_s21.py`

### What it does
1. Connects to the VNA and prints `*IDN?`.
2. (Optional) Recalls a saved instrument state `.sta` file. **By default this is
   OFF** (`RECALL_STATE = False`) so a live 2-port calibration is preserved —
   recalling a state invalidates calibration.
3. Reads back and prints the loaded setup (sweep, points, averaging, IF BW, source
   power, cal state) and warns if it does not match the SOP expectations.
4. Prompts once for device/run identification and builds a nested folder path (the
   "SOP data tree") **on both the VNA disk and a local mirror**.
5. Prompts for the **coordinate map** relating the fixture ruler reading to the lead
   position, and writes machine-readable geometry to `acquisition_info.json`.
6. Enters an **interactive loop**: for each z-axis point it shows the target
   coordinate and filename, waits for you to position the probe and press Enter,
   then triggers one averaged sweep, saves a `.s2p` on the VNA, mirrors that exact
   file into the run's `02_raw/` folder, and appends a line to `run_log.txt`. The
   per-point readout reports `\|S21\|` at `REPORT_FREQ` (not the sweep peak).
7. Optionally runs a **return-to-home drift check**: re-measure the first point and
   confirm `\|S21\|` still matches the baseline within 0.5 dB.
8. If any points were saved, **auto-launches `plot_tf_zaxis.py`** on the run folder
   so the CSVs and overlay plot are generated without a separate step.

### Configure before running
Edit the constants at the top of the script:

| Constant | Meaning |
|----------|---------|
| `inst_address` | VISA resource string, e.g. `TCPIP0::165.124.9.141::inst0::INSTR`. If it changed, run `ipconfig` on the VNA to find its IP. |
| `state_file` | Path to the `.sta` state file (on the VNA disk). |
| `RECALL_STATE` | `False` normally. Set `True` **only** to reload the state *before* calibrating. |
| `VNA_ROOT` | Root folder for saves on the analyzer, e.g. `D:\TF_measurements_workflow`. |
| `LOCAL_ROOT` | Root folder for the local mirror on this PC. Must match the plot script's `LOCAL_ROOT`. |
| `S2P_PORTS` | Ports in the Touchstone file — `(1, 2)`; Port 1 = Receive, Port 2 = Transmit. |
| `REPORT_FREQ` | Frequency (Hz) shown in the per-point `\|S21\|` readout and `run_log.txt`, e.g. `63.6e6`. Display only — the full sweep is always saved. |
| `SETUP_ONLY` | `True` = recall + verify the setup, then stop (use to test the connection). `False` = full acquisition. |

Sweep range, real/imaginary format, −5 dBm power, 500 Hz IF BW and 5 averages come
from the recalled/loaded instrument state — the script verifies but does **not**
re-set them.

### Run it
```sh
python recall_and_acquire_s21.py
```
Then answer the prompts:
- **Device/run identification** → builds the save-path folders.
- **Coordinate mapping** → START ruler & lead positions, step size, whether each
  coordinate increases per step, and an optional END ruler position.

**Loop controls:** `Enter` = acquire & advance · `b` = back one point ·
`f` = skip forward without measuring · `q` = quit.

### Filenames & folder layout
Files are auto-named `\[lead]MM_\[ruler]MM.s2p`, e.g. `40MM_930MM.s2p`. Whole numbers
are written plain; decimals use `p` for the dot (e.g. `152p5MM`). Each run folder
uses the SOP subfolder tree (naming convention from Fuchang's `tf_ops_single_exc.py`):

```
LOCAL_ROOT/
  <AIMD>/<Lead>/<Termination>/<IPG>/<Serial>/[<Excitation>/]<YYYYMMDD_RunN_Initials>/
      run_log.txt                     # human-readable audit log of the run
      01_form/                        # reserved for run paperwork
      02_raw/                         # native Touchstone sweeps + geometry
          0MM_890MM.s2p
          5MM_895MM.s2p
          ...
          acquisition_info.json       # machine-readable geometry (see §4)
      03_raw_processed/               # measured S21 vs distance, per freq (CSV; from the plotter)
      04_interpolated_extrapolated/   # full-length TF curve, per freq (CSV; from the plotter)
      05_plots/                       # overlay PNG (from the plotter)
```

`02_raw/` is filled by the acquisition script; `03`–`05` are filled when
`plot_tf_zaxis.py` runs (automatically at the end of a run, or manually later).

---

## 3. Plotting — `plot_tf_zaxis.py`

### What it does
Reads every `\[lead]MM_\[ruler]MM.s2p` file in a run's `02_raw/` folder, extracts S21
at each `TARGET_FREQS` entry, and plots two stacked panels vs probe position:
**`\|S21\|` (dB)** on top and **phase (deg)** below, with each target frequency
overlaid as its own curve. Use it for the SOP *Data Continuity & Physical Integrity
Check* — look for sharp jumps or discontinuities that indicate the lead touched the
receive ring.

It writes into the run's SOP output folders (one CSV **per target frequency**, plus a
single overlay figure):
- `03_raw_processed/` — `<run>__<freq>__original.csv`: the measured S21-vs-distance
  points actually used (source file, position, real/imag, `\|S21\|` dB, phase deg).
- `04_interpolated_extrapolated/` — `<run>__<freq>__interp_extrap.csv`: the
  full-length curve (measured points plus the extrapolated ends), with a `region`
  column tagging each row `measured` or `extrapolated`. Matches what the plot draws.
- `05_plots/` — `TF_zaxis_<freqs>MHz.png`: the overlay figure, also shown on screen.

You can point it at the **run root**, its **`02_raw/`** folder, or a **flat folder** of
`.s2p` files — it resolves the layout automatically and writes outputs next to `02_raw/`.

If `acquisition_info.json` is present (lead x-axis), the plotter linearly
**extrapolates** dashed segments from the measured range out to the full physical
lead length so the curve — and the `04_...` CSV — spans the whole lead.

### Configure
| Constant | Meaning |
|----------|---------|
| `RUN_DIR` | Run folder to plot. Leave `''` to auto-pick the most recent acquisition under `LOCAL_ROOT`. |
| `LOCAL_ROOT` | Same root as the acquisition script. |
| `TARGET_FREQS` | List of frequencies in **Hz** to read S21 at, e.g. `[63.6e6, 123.2e6]`. Each becomes one overlaid curve. The script snaps to the nearest measured point and reports what it used. |
| `X_AXIS` | `'lead'` or `'ruler'` — which coordinate is the x-axis. `'lead'` requires `acquisition_info.json`. |

### Run it
Usually launched automatically at the end of an acquisition. To run it standalone:
```sh
python plot_tf_zaxis.py                      # auto-pick latest run, or use RUN_DIR
python plot_tf_zaxis.py "C:\path\to\run_dir" # explicit run root or 02_raw folder
```
Output: per-frequency CSVs in `03_raw_processed/` and `04_interpolated_extrapolated/`,
plus `05_plots/TF_zaxis_<freqs>MHz.png`.

---

## 4. Reading the data

### `.s2p` Touchstone files
Plain text, one measured frequency per row. Lines starting with `!` are comments;
the line starting with `#` is the format header:

```
# Hz S RI R 50
```
means frequencies in **Hz**, **S-parameters** in **R**eal/**I**maginary pairs,
reference impedance **50 Ω**. Each data row has 9 columns:

```
freq  S11r S11i  S21r S21i  S12r S12i  S22r S22i
```

**S21** (transmission, Receive←Transmit) is columns 4–5. To reconstruct it:

```python
import numpy as np

def read_s21(path):
    freqs, s21 = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(('!', '#')):
                continue
            v = [float(x) for x in line.split()]
            freqs.append(v[0])                 # Hz
            s21.append(complex(v[3], v[4]))    # real + j*imag
    return np.array(freqs), np.array(s21)

f, s21 = read_s21("40MM_930MM.s2p")
mag_db = 20 * np.log10(np.abs(s21))            # magnitude in dB
phase_deg = np.degrees(np.angle(s21))          # phase in degrees
```

### `acquisition_info.json`
Machine-readable geometry for a run, written into `02_raw/` by the acquisition script
and consumed by the plotter. Key fields:

| Field | Meaning |
|-------|---------|
| `lead_length_mm` | Physical lead/rod length. |
| `plot_start_lead_mm`, `plot_end_lead_mm` | Full lead span the plot should cover (`0` … length). |
| `step_mm` | Step magnitude between points. |
| `ruler_start_mm`, `lead_start_mm` | START coordinates of the run. |
| `ruler_step_mm`, `lead_step_mm` | Signed per-step deltas (direction included). |
| `ruler_end_mm`, `planned_lead_end_mm` | Planned END coordinates (may be `null` if open-ended). |
| `start_unmeasured_mm`, `end_unmeasured_mm` | Distances from each physical lead end that were not measured. |
| `acquired_lead_min_mm`, `acquired_lead_max_mm` | Actual measured span (added after the run finishes). |

### Processed CSV outputs (`03_raw_processed/`, `04_interpolated_extrapolated/`)
Written by `plot_tf_zaxis.py`, one file per target frequency (label like `63p6MHz`):

- **`03_raw_processed/<run>__<freq>__original.csv`** — the measured points used for
  that frequency. Columns: `source_file`, `<lead|ruler>_position_mm`, `S21_real`,
  `S21_imag`, `S21_mag_dB`, `S21_phase_deg`.
- **`04_interpolated_extrapolated/<run>__<freq>__interp_extrap.csv`** — the
  full-length curve. Columns: `<lead|ruler>_position_mm`, `region`
  (`measured`/`extrapolated`), `S21_mag_dB`, `S21_phase_deg`.

### `run_log.txt`
Human-readable audit log at the run-folder root: a header block (operator, device
fields including serial number, stimulus, coordinate map, geometry) followed by one
line per saved point (timestamp, point index, lead & ruler positions, filename,
`\|S21\|` in dB at `REPORT_FREQ` and that frequency, and whether the local mirror
succeeded), then the return-to-home result and a finish line.

---

## 5. Notes & gotchas

- **Calibration:** a valid 2-port cal must be ON, or `.s2p` saves fail (E5063A error
  57). The acquisition script warns if correction is OFF. Recalling a `.sta` state
  invalidates the cal — keep `RECALL_STATE = False` for normal runs.
- **`LOCAL_ROOT` must match** between the two scripts, or the plotter won't find the
  data the acquisition script mirrored.
- **Path length:** the nested SOP folder plus filename can approach the Windows
  260-char limit — and the plotter's `04_interpolated_extrapolated/<run>__<freq>__…csv`
  paths are the longest of all. The acquisition script warns above ~200 chars. Prefer
  short field values, keep `LOCAL_ROOT` shallow, or enable Windows long-path support.
- This repo contains only the two scripts and their docs — the actual measurement
  data (`.s2p`, logs, PNGs) is written under `LOCAL_ROOT` elsewhere and is
  intentionally **not** tracked here (see `.gitignore`).

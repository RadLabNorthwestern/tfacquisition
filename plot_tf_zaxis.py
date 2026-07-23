"""
Post-run TF plot: Magnitude & Phase vs z-axis from the local .s2p mirror.

Reads every  [lead]MM_[ruler]MM.s2p  file in a run's 02_raw folder, pulls S21 at
the TARGET_FREQS, and plots |S21| (dB) and phase (deg) against probe position with
all target frequencies overlaid. Use it for the SOP "Data Continuity & Physical
Integrity Check" - look for sharp jumps or discontinuities that indicate the
lead touched the receive ring.

Outputs follow the SOP data tree created by recall_and_acquire_s21.py:
  03_raw_processed/             <run>__<freq>__original.csv       (measured points)
  04_interpolated_extrapolated/ <run>__<freq>__interp_extrap.csv  (full-length curve)
  05_plots/                     TF_zaxis_<freqs>MHz.png            (overlay plot)

Point this script at either the run root, its 02_raw folder, or a flat folder of
.s2p files - it resolves the layout automatically.
"""

import os
import re
import csv
import glob
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Config

# Local run folder holding the .s2p files. Leave as '' to auto-pick the most
# recent acquisition under LOCAL_ROOT. A command-line arg overrides both.
RUN_DIR = r''

# Root the acquisition script mirrors into (must match recall_and_acquire_s21.py).
LOCAL_ROOT = r'C:\Users\TAZ5297\Documents\TFCapture\TF_data'

# Frequencies (Hz) to read S21 at. Each is plotted as its own overlaid curve.
TARGET_FREQS = [63.6e6, 123.2e6]

# 'lead' or 'ruler' - which coordinate to use as the x-axis.
X_AXIS = 'lead'

# SOP folder tree (must match recall_and_acquire_s21.py).
RAW_FOLDER           = '02_raw'
RAW_PROCESSED_FOLDER = '03_raw_processed'
INTERP_EXTRAP_FOLDER = '04_interpolated_extrapolated'
PLOTS_FOLDER         = '05_plots'

###############################################################################

_UNIT = {'hz': 1.0, 'khz': 1e3, 'mhz': 1e6, 'ghz': 1e9}
# Coordinates are written with 'p' for the decimal point (e.g. 152p5MM), no dots.
_NAME = re.compile(r'(-?\d+(?:p\d+)?)MM_(-?\d+(?:p\d+)?)MM', re.IGNORECASE)


def parse_touchstone(path):
    """Return (freq_hz, s21_complex) from a 2-port real/imaginary .s2p file."""
    fscale = 1.0
    freqs, s21 = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            if line.startswith('#'):
                toks = line.lower().split()
                fscale = _UNIT[toks[1]]
                if 'ri' not in toks:
                    print(f'  WARNING: {os.path.basename(path)} is not RI format')
                continue
            v = [float(x) for x in line.split()]
            # freq  S11r S11i  S21r S21i  S12r S12i  S22r S22i
            freqs.append(v[0] * fscale)
            s21.append(complex(v[3], v[4]))
    return np.array(freqs), np.array(s21)


def resolve_layout(given):
    """Given any of {run root, its 02_raw folder, a flat .s2p folder}, return
    (run_root, raw_dir). Outputs (03/04/05) are written under run_root."""
    given = os.path.normpath(given)
    if os.path.basename(given).lower() == RAW_FOLDER.lower():
        return os.path.dirname(given), given            # given IS the 02_raw folder
    sub_raw = os.path.join(given, RAW_FOLDER)
    if os.path.isdir(sub_raw):
        return given, sub_raw                           # given is the run root
    return given, given                                 # flat folder of .s2p files


def freq_label(hz):
    """Folder-safe frequency label, e.g. 63.6 MHz -> '63p6MHz'."""
    return f'{hz / 1e6:g}'.replace('.', 'p') + 'MHz'


def latest_run_dir(root):
    """Folder containing the most recently written .s2p under root (or None)."""
    s2ps = glob.glob(os.path.join(root, '**', '*.s2p'), recursive=True)
    if not s2ps:
        return None
    return os.path.dirname(max(s2ps, key=os.path.getmtime))


def choose_run_dir_dialog(initialdir):
    """Open a native folder-picker and return the chosen path ('' if cancelled)."""
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception as ex:  # noqa: BLE001 - fall back to a typed path below
        print(f'  (file browser unavailable: {ex})')
        return ''
    root = tk.Tk()
    root.withdraw()                    # hide the empty root window
    root.attributes('-topmost', True)  # bring the dialog to the front
    path = filedialog.askdirectory(
        title='Select the run folder (or its 02_raw folder) with the .s2p files',
        initialdir=initialdir or os.getcwd())
    root.destroy()
    return path or ''


def parse_name(path):
    """Return (lead_mm, ruler_mm) parsed from the filename."""
    m = _NAME.search(os.path.basename(path))
    if not m:
        return None
    return (float(m.group(1).lower().replace('p', '.')),
            float(m.group(2).lower().replace('p', '.')))


def load_acquisition_info(*dirs):
    """Load the machine-readable geometry written by the acquisition script,
    searching each candidate directory in order."""
    for d in dirs:
        path = os.path.join(d, 'acquisition_info.json')
        if os.path.isfile(path):
            with open(path, encoding='utf-8') as fh:
                return json.load(fh)
    return None


def extension_positions(measured_x, target_start, target_end, step):
    """Return regularly spaced positions needed outside the measured interval."""
    if step <= 0:
        raise ValueError('Acquisition step must be positive')
    tol = step * 1e-6
    left = np.arange(measured_x[0] - step, target_start - tol, -step)[::-1]
    right = np.arange(measured_x[-1] + step, target_end + tol, step)
    # Include a boundary that is not an integer number of steps away.
    if target_start < measured_x[0] - tol and (not len(left) or left[0] > target_start + tol):
        left = np.insert(left, 0, target_start)
    if target_end > measured_x[-1] + tol and (not len(right) or right[-1] < target_end - tol):
        right = np.append(right, target_end)
    return left, right


def numpy_extrapolate(x_new, x, y):
    """NumPy linear interpolation with two-point linear endpoint extrapolation."""
    if len(x) < 2:
        return np.full(len(x_new), y[0])

    result = np.interp(x_new, x, y)
    left = x_new < x[0]
    right = x_new > x[-1]
    result[left] = y[0] + (x_new[left] - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
    result[right] = y[-1] + (x_new[right] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return result


def main():
    if len(sys.argv) > 1:
        in_dir = sys.argv[1].strip('" ')
    elif RUN_DIR:
        in_dir = RUN_DIR
    else:
        # Open a folder browser, starting at the most recent acquisition so the
        # right run is one click away (but any folder can be chosen).
        initialdir = latest_run_dir(LOCAL_ROOT) or LOCAL_ROOT
        print('Select the run folder in the file browser...')
        in_dir = choose_run_dir_dialog(initialdir)
        if not in_dir:
            in_dir = input('No folder selected - enter run folder path: ').strip('" ')
    if not in_dir or not os.path.isdir(in_dir):
        sys.exit(f'Not a folder: {in_dir}')

    run_root, raw_dir = resolve_layout(in_dir)

    # Dedupe (Windows glob is case-insensitive, so *.s2p can match *.S2P twice)
    files = {os.path.normcase(p): p for p in glob.glob(os.path.join(raw_dir, '*.s2p'))}
    files = sorted(files.values())
    if not files:
        sys.exit(f'No .s2p files found in {raw_dir}')

    records = []  # (x, [s21 at each target freq], source_name)
    used = None   # actual nearest grid frequencies, reported once
    for fp in files:
        name = parse_name(fp)
        if name is None:
            print(f'  skipping (unrecognized name): {os.path.basename(fp)}')
            continue
        lead, ruler = name
        f, s21 = parse_touchstone(fp)
        idxs = [int(np.argmin(np.abs(f - tf))) for tf in TARGET_FREQS]
        if used is None:                        # report what actually got picked
            used = [f[i] for i in idxs]
            for tf, uf in zip(TARGET_FREQS, used):
                note = '' if abs(uf - tf) < 1e3 else f'  (nearest to {tf/1e6:g} MHz requested)'
                print(f'  using {uf/1e6:g} MHz{note}')
        records.append((lead if X_AXIS == 'lead' else ruler,
                        [s21[i] for i in idxs],
                        os.path.basename(fp)))

    records.sort(key=lambda r: r[0])
    if not records:
        sys.exit('No files had recognizable lead/ruler coordinates.')
    x = np.array([r[0] for r in records])

    info = load_acquisition_info(raw_dir, run_root)
    left_x = right_x = np.array([])
    if X_AXIS == 'lead' and info is not None:
        target_start = float(info['plot_start_lead_mm'])
        target_end = float(info['plot_end_lead_mm'])
        measured_steps = np.diff(x)
        inferred_step = float(np.median(measured_steps)) if len(measured_steps) else 1.0
        step = abs(float(info.get('lead_step_mm', info.get('step_mm', inferred_step))))
        if target_start > x[0] or target_end < x[-1]:
            sys.exit(f'Acquisition geometry {target_start:g}..{target_end:g} mm does not '
                     f'contain measured range {x[0]:g}..{x[-1]:g} mm.')
        left_x, right_x = extension_positions(x, target_start, target_end, step)
        print(f'  requested lead range: {target_start:g}..{target_end:g} mm '
              f'(length {info["lead_length_mm"]:g} mm; unmeasured ends '
              f'{x[0] - target_start:g}/{target_end - x[-1]:g} mm)')
        print(f'  extrapolating {len(left_x)} point(s) before and {len(right_x)} after '
              'using the nearest two measured points')
    elif X_AXIS == 'lead':
        sys.exit('Missing acquisition_info.json. This plotter requires geometry from '
                 'the current acquisition workflow.')

    # Output folders (created next to 02_raw, under the run root).
    proc_dir = os.path.join(run_root, RAW_PROCESSED_FOLDER)
    interp_dir = os.path.join(run_root, INTERP_EXTRAP_FOLDER)
    plots_dir = os.path.join(run_root, PLOTS_FOLDER)
    for d in (proc_dir, interp_dir, plots_dir):
        os.makedirs(d, exist_ok=True)
    run_name = os.path.basename(os.path.normpath(run_root)) or 'run'
    axis_col = f'{X_AXIS}_position_mm'

    xlabel = f'{X_AXIS.capitalize()} position (mm)'
    plt.figure(figsize=(10, 6))
    ax_mag = plt.subplot(2, 1, 1)
    ax_ph = plt.subplot(2, 1, 2)

    for j, uf in enumerate(used):
        s = np.array([r[1][j] for r in records])
        label = f'{uf/1e6:g} MHz'
        flabel = freq_label(uf)
        mag_db = 20 * np.log10(np.abs(s))
        phase_unwrapped = np.degrees(np.unwrap(np.angle(s)))
        mag_line, = ax_mag.plot(x, mag_db, marker='.', markersize=5, label=label)
        phase_line, = ax_ph.plot(x, phase_unwrapped, marker='.', markersize=5, label=label)
        color = mag_line.get_color()

        ext = {'left': (np.array([]), np.array([])), 'right': (np.array([]), np.array([]))}
        for side, x_ext in (('left', left_x), ('right', right_x)):
            if not len(x_ext):
                continue
            mag_ext = numpy_extrapolate(x_ext, x, mag_db)
            phase_ext = numpy_extrapolate(x_ext, x, phase_unwrapped)
            ext[side] = (mag_ext, phase_ext)
            join_x = np.r_[x_ext, x[0]] if side == 'left' else np.r_[x[-1], x_ext]
            join_mag = np.r_[mag_ext, mag_db[0]] if side == 'left' else np.r_[mag_db[-1], mag_ext]
            join_phase = (np.r_[phase_ext, phase_unwrapped[0]] if side == 'left'
                          else np.r_[phase_unwrapped[-1], phase_ext])
            ax_mag.plot(join_x, join_mag, '--', color=color, alpha=.8)
            ax_ph.plot(join_x, join_phase, '--', color=phase_line.get_color(), alpha=.8)

        # 03_raw_processed: the measured S21 vs distance actually used for this freq.
        orig_csv = os.path.join(proc_dir, f'{run_name}__{flabel}__original.csv')
        with open(orig_csv, 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(['source_file', axis_col, 'S21_real', 'S21_imag',
                        'S21_mag_dB', 'S21_phase_deg'])
            for r, sv, md, ph in zip(records, s, mag_db, phase_unwrapped):
                w.writerow([r[2], r[0], sv.real, sv.imag, md, ph])
        print(f'  saved {orig_csv}')

        # 04_interpolated_extrapolated: the full-length curve (measured points plus
        # the extrapolated ends), matching exactly what the plot draws.
        mag_left, phase_left = ext['left']
        mag_right, phase_right = ext['right']
        x_full = np.concatenate([left_x, x, right_x])
        mag_full = np.concatenate([mag_left, mag_db, mag_right])
        phase_full = np.concatenate([phase_left, phase_unwrapped, phase_right])
        region = (['extrapolated'] * len(left_x) + ['measured'] * len(x)
                  + ['extrapolated'] * len(right_x))
        interp_csv = os.path.join(interp_dir, f'{run_name}__{flabel}__interp_extrap.csv')
        with open(interp_csv, 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow([axis_col, 'region', 'S21_mag_dB', 'S21_phase_deg'])
            for xv, rg, md, ph in zip(x_full, region, mag_full, phase_full):
                w.writerow([xv, rg, md, ph])
        print(f'  saved {interp_csv}')

    ax_mag.set(xlabel=xlabel, ylabel='|S21| (dB)',
               title=f'TF Magnitude vs z-axis  ({len(x)} measured points)')
    ax_mag.grid(True)
    ax_mag.legend()
    ax_ph.set(xlabel=xlabel, ylabel='Phase (deg)', title='TF Phase vs z-axis')
    ax_ph.grid(True)
    ax_ph.legend()

    plt.tight_layout()
    tag = '_'.join(f'{uf/1e6:g}'.replace('.', 'p') for uf in used)
    out = os.path.join(plots_dir, f'TF_zaxis_{tag}MHz.png')
    plt.savefig(out, dpi=150)
    print(f'Saved {out}')
    plt.show()


if __name__ == '__main__':
    main()

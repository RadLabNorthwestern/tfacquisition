"""
Post-run TF plot: Magnitude & Phase vs z-axis from the local .s2p mirror.

Reads every  [lead]MM_[ruler]MM.s2p  file in a run folder, pulls S21 at the
TARGET_FREQS, and plots |S21| (dB) and phase (deg) against probe position with
all target frequencies overlaid. Use it for the SOP "Data Continuity & Physical
Integrity Check" - look for sharp jumps or discontinuities that indicate the
lead touched the receive ring.
"""

import os
import re
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


def latest_run_dir(root):
    """Folder containing the most recently written .s2p under root (or None)."""
    s2ps = glob.glob(os.path.join(root, '**', '*.s2p'), recursive=True)
    if not s2ps:
        return None
    return os.path.dirname(max(s2ps, key=os.path.getmtime))


def parse_name(path):
    """Return (lead_mm, ruler_mm) parsed from the filename."""
    m = _NAME.search(os.path.basename(path))
    if not m:
        return None
    return (float(m.group(1).lower().replace('p', '.')),
            float(m.group(2).lower().replace('p', '.')))


def load_acquisition_info(run_dir):
    """Load the machine-readable geometry written by the acquisition script."""
    path = os.path.join(run_dir, 'acquisition_info.json')
    if not os.path.isfile(path):
        return None
    with open(path, encoding='utf-8') as fh:
        return json.load(fh)


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
        run_dir = sys.argv[1].strip('" ')
    elif RUN_DIR:
        run_dir = RUN_DIR
    else:
        run_dir = latest_run_dir(LOCAL_ROOT)
        if run_dir:
            print(f'Latest acquisition: {run_dir}')
        else:
            run_dir = input('No .s2p found under LOCAL_ROOT - enter run folder: ').strip('" ')
    if not os.path.isdir(run_dir):
        sys.exit(f'Not a folder: {run_dir}')

    # Dedupe (Windows glob is case-insensitive, so *.s2p can match *.S2P twice)
    files = {os.path.normcase(p): p for p in glob.glob(os.path.join(run_dir, '*.s2p'))}
    files = sorted(files.values())
    if not files:
        sys.exit(f'No .s2p files found in {run_dir}')

    records = []  # (x, [s21 at each target freq])
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
                        [s21[i] for i in idxs]))

    records.sort(key=lambda r: r[0])
    if not records:
        sys.exit('No files had recognizable lead/ruler coordinates.')
    x = np.array([r[0] for r in records])

    info = load_acquisition_info(run_dir)
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

    xlabel = f'{X_AXIS.capitalize()} position (mm)'
    plt.figure(figsize=(10, 6))
    ax_mag = plt.subplot(2, 1, 1)
    ax_ph = plt.subplot(2, 1, 2)

    for j, uf in enumerate(used):
        s = np.array([r[1][j] for r in records])
        label = f'{uf/1e6:g} MHz'
        mag_db = 20 * np.log10(np.abs(s))
        phase_unwrapped = np.degrees(np.unwrap(np.angle(s)))
        mag_line, = ax_mag.plot(x, mag_db, marker='.', markersize=5, label=label)
        phase_line, = ax_ph.plot(x, phase_unwrapped, marker='.', markersize=5, label=label)
        color = mag_line.get_color()
        for side, x_ext in (('left', left_x), ('right', right_x)):
            if not len(x_ext):
                continue
            mag_ext = numpy_extrapolate(x_ext, x, mag_db)
            phase_ext = numpy_extrapolate(x_ext, x, phase_unwrapped)
            join_x = np.r_[x_ext, x[0]] if side == 'left' else np.r_[x[-1], x_ext]
            join_mag = np.r_[mag_ext, mag_db[0]] if side == 'left' else np.r_[mag_db[-1], mag_ext]
            join_phase = (np.r_[phase_ext, phase_unwrapped[0]] if side == 'left'
                          else np.r_[phase_unwrapped[-1], phase_ext])
            ax_mag.plot(join_x, join_mag, '--', color=color, alpha=.8)
            ax_ph.plot(join_x, join_phase, '--', color=phase_line.get_color(), alpha=.8)

    ax_mag.set(xlabel=xlabel, ylabel='|S21| (dB)',
               title=f'TF Magnitude vs z-axis  ({len(x)} measured points)')
    ax_mag.grid(True)
    ax_mag.legend()
    ax_ph.set(xlabel=xlabel, ylabel='Phase (deg)', title='TF Phase vs z-axis')
    ax_ph.grid(True)
    ax_ph.legend()

    plt.tight_layout()
    tag = '_'.join(f'{uf/1e6:g}'.replace('.', 'p') for uf in used)
    out = os.path.join(run_dir, f'TF_zaxis_{tag}MHz.png')
    plt.savefig(out, dpi=150)
    print(f'Saved {out}')
    plt.show()


if __name__ == '__main__':
    main()

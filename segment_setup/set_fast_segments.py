"""
Rebuild the E5063A segment sweep FROM SCRATCH for a faster TF acquisition.

The previous grid (Frequencies_lock_pm1MHz_0p1MHz.sta) expanded every center into
a +/-1 MHz band at 0.1 MHz spacing, so a lot of time was spent on frequency points
that aren't of interest. This replaces the whole segment table with:

  * BANDS_MHZ         - only the narrow ranges we actually care about, swept at
                        STEP_HZ spacing.
  * SINGLE_POINT_FREQS_HZ - the two top frequencies, one point each (as before).

It reads the live table only to inherit the per-segment IF BW / power *format*
(and reuse the global IF BW / power values), shows a before/after preview, asks
for confirmation, writes the new table, verifies the grid, and saves the grid
back OVER the existing state file (STATE_FILE) that recall_and_acquire_s21.py
already uses -- so there is no new state to recall.

NOTE: changing the stimulus INVALIDATES the calibration. After running this you
must re-run the 2-port cal on this grid before acquiring. recall_and_acquire_s21.py
needs no change (its state_file already points at STATE_FILE, and the new
stimulus is already live on the instrument).
"""

import numpy as np
import pyvisa

inst_address = 'TCPIP0::165.124.9.141::inst0::INSTR'

# The instrument state file to OVERWRITE in place with the new grid. This is the
# same .sta that recall_and_acquire_s21.py recalls (its state_file), so nothing
# there needs to change and there is no new state to recall. Confirm this matches
# the state file you currently have loaded before running.
STATE_FILE = r'D:\TF_measurements_workflow\Frequencies_lock_pm1MHz_0p1MHz.sta'

###############################################################################
# Config -- the new acquisition grid

# Narrow ranges of interest (MHz), inclusive start .. stop. Swept at STEP_HZ.
BANDS_MHZ = [
    (23.0, 24.0),
    (50.0, 51.6),
    (63.0, 64.0),
    (123.0, 124.0),
    (127.0, 128.0),
]

STEP_HZ = 0.1e6            # point spacing inside each band (0.1 MHz = same as before)

# The two top frequencies, kept as a single point each (as they were originally).
SINGLE_POINT_FREQS_HZ = [
    297.0e6,
    447.0e6,
]

###############################################################################

rm = pyvisa.ResourceManager()
inst = rm.open_resource(inst_address)
inst.timeout = 30000
print(inst.query('*IDN?').strip())
inst.write('*CLS')

# --- Read the live segment-table format (flags) so we write a compatible table
if not inst.query(':SENSe1:SWEep:TYPE?').strip().upper().startswith('SEGM'):
    raise SystemExit('Sweep type is not SEGMent - aborting so nothing is clobbered.')

vals = [float(x) for x in inst.query(':SENSe1:SEGMent:DATA?').split(',')]
nflags = int(vals[0])
flags = [int(v) for v in vals[1:1 + nflags]]          # [startstop, ifbw, pow, delay, time]

# Values to stamp on every new segment when the per-segment flag is enabled.
# Reuse the instrument's current global IF BW / power (the SOP 500 Hz / -5 dBm).
ifbw_val = float(inst.query(':SENSe1:BWIDth?'))
power_val = float(inst.query(':SOURce1:POWer?'))
print(f'Reusing IF BW = {ifbw_val:.0f} Hz, power = {power_val:.1f} dBm for all segments.')


def band_points(start, stop):
    return int(round((stop - start) / STEP_HZ)) + 1


# --- Build the planned segment list: (start_hz, stop_hz, npts, kind) ----------
plan = []
for lo_mhz, hi_mhz in BANDS_MHZ:
    lo, hi = lo_mhz * 1e6, hi_mhz * 1e6
    plan.append((lo, hi, band_points(lo, hi), 'band'))
for f in SINGLE_POINT_FREQS_HZ:
    plan.append((f, f, 1, 'single'))

plan.sort(key=lambda p: p[0])

# --- Sanity checks ------------------------------------------------------------
for (a0, a1, *_), (b0, b1, *_) in zip(plan, plan[1:]):
    if b0 <= a1:
        raise SystemExit(f'Segments overlap: {a0/1e6:g}..{a1/1e6:g} and '
                         f'{b0/1e6:g}..{b1/1e6:g} MHz. Fix BANDS_MHZ / '
                         'SINGLE_POINT_FREQS_HZ.')

# --- Build the raw table (force start/stop mode, keep the other flags) --------
out_flags = [0] + flags[1:]
out = [float(nflags)] + [float(f) for f in out_flags] + [float(len(plan))]
for start, stop, npts, _kind in plan:
    seg = [start, stop, float(npts)]
    if flags[1]:
        seg.append(ifbw_val)
    if flags[2]:
        seg.append(power_val)
    if len(flags) > 3 and flags[3]:
        seg.append(0.0)
    if len(flags) > 4 and flags[4]:
        seg.append(0.0)
    out += seg

# --- Preview -----------------------------------------------------------------
print(f'\n{"segment (MHz)":>22}   {"pts":>4}   kind')
print('-' * 44)
total = 0
for start, stop, npts, kind in plan:
    total += npts
    label = f'{start/1e6:g}' if kind == 'single' else f'{start/1e6:g} .. {stop/1e6:g}'
    print(f'{label:>22}   {npts:>4}   {kind}')
print('-' * 44)
print(f'Total points: {total}   (was ~147 on the +/-1 MHz grid)')

if input('\nWrite this segment table to the VNA? [y/N]: ').strip().lower() != 'y':
    raise SystemExit('Aborted - nothing written.')


def fmt(v):
    return str(int(round(v))) if abs(v - round(v)) < 1e-6 else repr(v)


inst.write(':SENSe1:SEGMent:DATA ' + ','.join(fmt(v) for v in out))
inst.query('*OPC?')
err = inst.query(':SYSTem:ERRor?').strip()
if not err.startswith(('+0', '0')):
    raise SystemExit(f'Write failed: {err}')

# --- Verify the grid: every single point + every band edge is hit -------------
freqs = np.array(inst.query_ascii_values(':SENSe1:FREQuency:DATA?'))
print(f'\nWritten. Measured points now: {len(freqs)}')
worst = 0.0
for start, stop, npts, kind in plan:
    for target in ({start} if kind == 'single' else {start, stop}):
        hit = freqs[int(np.argmin(np.abs(freqs - target)))]
        worst = max(worst, abs(hit - target))
print(f'Largest miss (single points + band edges): {worst:.3f} Hz  (should be ~0)')

# --- Overwrite the existing state file in place (state-only) ------------------
inst.write(':MMEMory:STORe:STYPe STATe')      # state only, no cal data (matches the original)
inst.write(f':MMEMory:STORe "{STATE_FILE}"')
inst.query('*OPC?')
serr = inst.query(':SYSTem:ERRor?').strip()
if not serr.startswith(('+0', '0')):
    raise SystemExit(f'State save failed: {serr}')
print(f'\nSaved new grid (state-only) over: {STATE_FILE}')

print('\nStimulus changed -> the cal is now invalid. Next:')
print('  1. Re-run the 2-port cal on this grid (the new stimulus is already live).')
print('  2. recall_and_acquire_s21.py needs no change -- its state_file already')
print(f'     points at {STATE_FILE}.')

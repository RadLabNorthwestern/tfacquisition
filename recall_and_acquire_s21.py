"""
Interactive Transfer Function (TF) acquisition for the Keysight E5063A ENA.

Workflow (matches the TF Measurement Procedure SOP):
  1. Recall Frequencies_lock.sta and verify the loaded stimulus / setup.
  2. Enter device/run info once; the script builds the SOP data-tree folder on
     both the VNA and a local mirror.
  3. For each z-axis position: move the receive probe, press Enter. The script
     triggers one averaged sweep, saves a native Touchstone .s2p on the VNA,
     mirrors that exact file to the local disk, and appends to a run log.
     Files are auto-named  [lead]MM_[ruler]MM.s2p  (e.g. 40MM_150MM.s2p).
  4. Optional return-to-home drift check against the first (baseline) point.

Real/Imaginary SnP format, -5 dBm, 500 Hz IFBW and 5 averages come from the
recalled state (enforced by the Calibration Log), so they are not re-set here.
"""

import os
import re
import sys
import json
import numpy as np
import pyvisa
from datetime import datetime, date

###############################################################################
# Config

inst_address = 'TCPIP0::165.124.9.141::inst0::INSTR'
# If the address changed: on the VNA, open command prompt and run "ipconfig".

state_file = r'D:\TF_measurements_workflow\Frequencies_lock_pm1MHz_0p1MHz.sta'

# Recalling the .sta reloads stimulus and INVALIDATES any active calibration.
# Per the SOP, recall the state manually ONCE before calibrating, then run this
# script with RECALL_STATE = False so the live 2-port cal is preserved.
RECALL_STATE = False

# Roots of the SOP data tree. VNA_ROOT is on the analyzer; LOCAL_ROOT is on this PC.
VNA_ROOT   = r'D:\TF_measurements_workflow'
LOCAL_ROOT = r'C:\Users\TAZ5297\Documents\TFCapture\TF_data'

# Ports that make up the 2-port Touchstone file (Port 1 = Receive, Port 2 = Transmit).
S2P_PORTS = (1, 2)

# Frequency (Hz) reported in the per-point |S21| readout and run log. The full
# sweep is always saved to the .s2p; this only selects which bin to display
# (nearest measured point). Set to the Larmor frequency of interest.
REPORT_FREQ = 63.6e6

# True  -> recall + verify the setup, then stop (use this to test the connection).
# False -> run the full interactive acquisition.
SETUP_ONLY = False

# SOP folder tree created inside each run folder (naming convention from
# Fuchang's tf_ops_single_exc.py). Raw .s2p files land in 02_raw; the remaining
# folders are created ready for downstream processing.
FORM_FOLDER          = '01_form'
RAW_FOLDER           = '02_raw'
RAW_PROCESSED_FOLDER = '03_raw_processed'
INTERP_EXTRAP_FOLDER = '04_interpolated_extrapolated'
PLOTS_FOLDER         = '05_plots'

###############################################################################
# Small console-input helpers


def ask(prompt, default=None, cast=str):
    """Prompt until a valid value is entered. Blank uses default if one exists."""
    suffix = f' [{default}]' if default is not None else ''
    while True:
        raw = input(f'{prompt}{suffix}: ').strip()
        if raw == '':
            if default is not None:
                return default if cast is str else cast(default)
            print('  (required)')
            continue
        try:
            return cast(raw)
        except ValueError:
            print(f'  invalid, expected {cast.__name__}')


def fmt_mm(x):
    """Format a mm coordinate for filenames: int when whole (40), else 'p' decimal (152p5)."""
    if abs(x - round(x)) < 1e-6:
        return str(int(round(x)))
    return f'{x:.1f}'.replace('.', 'p')


_RESERVED = {'CON', 'PRN', 'AUX', 'NUL',
             *(f'COM{i}' for i in range(1, 10)),
             *(f'LPT{i}' for i in range(1, 10))}


def safe(seg):
    """Turn a free-text field into a safe single path segment (folder name)."""
    seg = seg.strip()
    seg = seg.replace('.', 'p')                        # dots -> 'p'  (1.5T -> 1p5T)
    seg = re.sub(r'[<>:"/\\|?*]+', '_', seg)           # illegal Windows chars
    seg = re.sub(r'[\x00-\x1f]+', '_', seg)            # control chars
    seg = re.sub(r'\s+', '_', seg)                     # whitespace -> underscore
    seg = re.sub(r'_+', '_', seg)                      # collapse repeated underscores
    seg = seg.strip('_-')                              # no leading/trailing _ or -
    if seg.upper() in _RESERVED:                       # Windows reserved device names
        seg += '_'
    return seg or 'NA'


###############################################################################
# Connect

rm = pyvisa.ResourceManager()
inst = rm.open_resource(inst_address)
inst.timeout = 60000          # ms; an averaged wide sweep + file save can take a while
print(inst.query('*IDN?').strip())
inst.write('*CLS')            # clear any stale errors from previous runs

###############################################################################
# 1) Optionally recall the state (only do this BEFORE calibrating)

if RECALL_STATE:
    inst.write(f':MMEMory:LOAD "{state_file}"')
    inst.query('*OPC?')       # block until the recall completes
    err = inst.query(':SYSTem:ERRor?').strip()
    if not err.startswith(('+0', '0')):
        raise RuntimeError(f'Instrument reported an error after recall: {err}')
    print('Recalled state (this invalidates any prior calibration).')

###############################################################################
# 2) Verify the setup that was actually loaded (against the SOP expectations)

start   = float(inst.query(':SENSe1:FREQuency:STARt?'))
stop    = float(inst.query(':SENSe1:FREQuency:STOP?'))
swtype  = inst.query(':SENSe1:SWEep:TYPE?').strip()
fdata   = np.array(inst.query_ascii_values(':SENSe1:FREQuency:DATA?'))
points  = len(fdata)                                    # actual measured points
lin_pts = int(float(inst.query(':SENSe1:SWEep:POINts?')))
par1   = inst.query(':CALCulate1:PARameter1:DEFine?').strip()
avgc   = int(float(inst.query(':SENSe1:AVERage:COUNt?')))
avg_on = int(float(inst.query(':SENSe1:AVERage:STATe?')))
ifbw   = float(inst.query(':SENSe1:BWIDth?'))
power  = float(inst.query(':SOURce1:POWer?'))
corr   = int(float(inst.query(':SENSe1:CORRection:STATe?')))

print('\n=== Loaded setup ===')
print(f'  Sweep type : {swtype}')
if swtype.upper().startswith('SEGM'):
    print(f'  Points     : {points} measured (segment sweep; the inactive '
          f'linear-points setting is {lin_pts})')
else:
    print(f'  Sweep      : {start/1e6:.6g} -> {stop/1e6:.6g} MHz, {points} points')
if points <= 24:
    print('  Freqs (MHz): ' + ', '.join(f'{f/1e6:g}' for f in fdata))
print(f'  Tr1 param  : {par1}   (expected S21)')
print(f'  Averaging  : {"ON" if avg_on else "OFF"}, count {avgc}   (SOP: ON, 5)')
print(f'  IF BW      : {ifbw:.0f} Hz   (SOP: 500)')
print(f'  Source pwr : {power:.1f} dBm   (SOP: -5)')
print(f'  2-port cal : {"ON" if corr else "OFF"}   (required for .s2p save)')

if par1.replace('"', '').upper() != 'S21':
    print('  WARNING: Tr1 is not S21 - check the recalled state before proceeding.')
if not corr:
    print('  WARNING: error correction is OFF - .s2p saves will fail (error 57). '
          'Complete the 2-port cal before acquiring.')

if SETUP_ONLY:
    print('\nSETUP_ONLY is set: state recalled and verified. Stopping before acquisition.')
    sys.exit(0)

###############################################################################
# 3) Device / run identification -> builds the SOP data-tree folder

print('\n=== Device / run identification (builds the save path) ===')
aimd     = ask('AIMD type (e.g. CIED, DBS, SCS, VNS, etc.)')
lead     = ask('Lead (vendor_model_length (mm))')
term     = ask('Termination condition (Need to follow the format: Abandoned_Capped, Abandoned_Uncapped, Full_system, etc)')
ipg      = ask('IPG model (Need to follow the format: IPG_manufacturer_model. If none, write "NOIPG")')
serial   = ask('Serial number')
excit    = ask('Excitation / branches (blank if single lead)', default='')
run_date = ask('Date (YYYYMMDD), hit enter for default', default=date.today().strftime('%Y%m%d'))
run_no   = ask('Run number')
initials = ask('Operator initials')
lead_length = ask('Lead / rod length (mm)', cast=float)

if lead_length <= 0:
    sys.exit('Lead length must be positive.')

# Assemble the relative path segments (skip the excitation level for single leads)
parts = [safe(aimd), safe(lead), safe(term), safe(ipg), safe(serial)]
if excit.strip():
    parts.append(safe(excit))
parts.append(f'{safe(run_date)}_Run{safe(run_no)}_{safe(initials)}')

vna_dir   = VNA_ROOT.rstrip('\\') + '\\' + '\\'.join(parts)
local_dir = os.path.join(LOCAL_ROOT, *parts)

print('\nVNA folder   :', vna_dir)
print('Local folder :', local_dir)
if max(len(vna_dir), len(local_dir)) > 200:
    print('WARNING: assembled path is long (>200 chars); with the added filename it '
          'may hit the Windows 260-char limit. Consider shorter field values.')
if input('Proceed with these paths? [Enter]=yes, q=quit: ').strip().lower() == 'q':
    sys.exit(0)

# Create the tree locally and on the VNA (best-effort; "already exists" is ignored)
os.makedirs(local_dir, exist_ok=True)
raw_dir = os.path.join(local_dir, RAW_FOLDER)
os.makedirs(raw_dir, exist_ok=True)
os.makedirs(os.path.join(local_dir, FORM_FOLDER), exist_ok=True)
os.makedirs(os.path.join(local_dir, RAW_PROCESSED_FOLDER), exist_ok=True)
os.makedirs(os.path.join(local_dir, INTERP_EXTRAP_FOLDER), exist_ok=True)
os.makedirs(os.path.join(local_dir, PLOTS_FOLDER), exist_ok=True)

_cur = VNA_ROOT.rstrip('\\')
for _seg in parts:
    _cur = _cur + '\\' + _seg
    inst.write(f':MMEMory:MDIRectory "{_cur}"')
    inst.query(':SYSTem:ERRor?')     # discard "directory already exists"

log_path = os.path.join(local_dir, 'run_log.txt')


def log(line):
    with open(log_path, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


###############################################################################
# 4) Coordinate mapping (fixture ruler <-> lead position)

print('\n=== Coordinate mapping ===')
ruler_start = ask('START ruler position (mm)', cast=float)
lead_start  = ask('START lead position (mm)', cast=float)
step_mag    = ask('Step size (mm)', default='5', cast=float)
ruler_up    = ask('Does the RULER reading increase each step? (y/n)', default='y').lower().startswith('y')
lead_up     = ask('Does the LEAD position increase each step? (y/n)', default='y').lower().startswith('y')
end_raw     = ask('END ruler position (mm), blank = open-ended', default='')

ruler_end  = float(end_raw) if end_raw != '' else None
ruler_step = step_mag if ruler_up else -step_mag
lead_step  = step_mag if lead_up else -step_mag

# Machine-readable acquisition geometry for downstream plotting.  Lead position
# zero and lead_length are the physical ends.  The unmeasured endpoint distances
# are derived from the coordinate map instead of asking the operator for them.
planned_lead_end = None
if ruler_end is not None:
    planned_lead_end = lead_start + ((ruler_end - ruler_start) / ruler_step) * lead_step
metadata_path = os.path.join(raw_dir, 'acquisition_info.json')
metadata = {
    'schema_version': 1,
    'lead_length_mm': lead_length,
    'plot_start_lead_mm': 0.0,
    'plot_end_lead_mm': lead_length,
    'step_mm': step_mag,
    'ruler_start_mm': ruler_start,
    'lead_start_mm': lead_start,
    'start_unmeasured_mm': min(lead_start, lead_length - lead_start),
    'ruler_step_mm': ruler_step,
    'lead_step_mm': lead_step,
    'ruler_end_mm': ruler_end,
    'planned_lead_end_mm': planned_lead_end,
}
if planned_lead_end is not None:
    metadata['end_unmeasured_mm'] = min(planned_lead_end,
                                        lead_length - planned_lead_end)
with open(metadata_path, 'w', encoding='utf-8') as fh:
    json.dump(metadata, fh, indent=2)
    fh.write('\n')

if ruler_end is not None:
    span = ruler_end - ruler_start
    if span == 0:
        print('Planned: 1 point')
    elif (span > 0) == (ruler_step > 0):
        print(f'Planned: {int(abs(span)/step_mag + 1e-6) + 1} points, '
              f'ruler {fmt_mm(ruler_start)} -> {fmt_mm(ruler_end)} mm')
    else:
        print('WARNING: END is on the opposite side of START for the chosen ruler '
              'direction - it will never be reached. Check the direction answers.')

###############################################################################
# 5) Arm the instrument for single, deterministic sweeps

inst.write(':CALCulate1:PARameter1:SELect')     # make S21 (Tr1) the active trace for readback
inst.write(':TRIGger:SEQuence:SOURce BUS')      # only sweep on our bus trigger
inst.write(':INITiate1:CONTinuous ON')          # stay armed / waiting for the bus trigger
if avg_on:
    inst.write(':TRIGger:SEQuence:AVERage ON')  # one trigger completes the whole N-average

inst.write(f':MMEMory:STORe:SNP:TYPE:S2P {S2P_PORTS[0]},{S2P_PORTS[1]}')
inst.query(':SYSTem:ERRor?')                    # clear the error queue before the run

freq = np.array(inst.query_ascii_values(':SENSe1:FREQuency:DATA?'))


def sweep_s21():
    """Trigger one fresh averaged sweep and return complex S21 (numpy array)."""
    if avg_on:
        inst.write(':SENSe1:AVERage:CLEar')     # fresh N-average at THIS position
    inst.write(':TRIGger:SEQuence:SINGle')
    inst.query('*OPC?')                          # block until sweep (+ averaging) is done
    sdata = np.array(inst.query_ascii_values(':CALCulate1:DATA:SDATa?'))
    return sdata[0::2] + 1j * sdata[1::2]


def save_s2p(vna_filepath, local_filepath):
    """Save a native .s2p on the VNA, then mirror that exact file to local disk."""
    inst.write(f':MMEMory:STORe:SNP:DATA "{vna_filepath}"')
    inst.query('*OPC?')
    e = inst.query(':SYSTem:ERRor?').strip()
    if not e.startswith(('+0', '0')):
        raise RuntimeError(f'VNA save failed for {vna_filepath}: {e}')
    # Transfer the just-saved file back byte-for-byte (keeps all 4 S-params + RI format)
    try:
        raw = inst.query_binary_values(f':MMEMory:TRANsfer? "{vna_filepath}"',
                                       datatype='B', container=bytearray)
        with open(local_filepath, 'wb') as fh:
            fh.write(bytes(raw))
        return True
    except Exception as ex:  # noqa: BLE001 - mirror is best-effort; VNA copy is the record
        print(f'  WARNING: local mirror failed ({ex}); the VNA copy is saved.')
        return False


###############################################################################
# 6) Run-log header

log('=' * 72)
log('TF acquisition run log')
log(f'Started      : {datetime.now():%Y-%m-%d %H:%M:%S}')
log(f'Operator     : {initials}')
log(f'AIMD         : {aimd}')
log(f'Lead         : {lead}')
log(f'Termination  : {term}')
log(f'IPG          : {ipg}')
log(f'Serial       : {serial}')
if excit.strip():
    log(f'Excitation   : {excit}')
log(f'VNA folder   : {vna_dir}')
log(f'Local folder : {local_dir}')
log(f'Stimulus     : {start/1e6:.6g}-{stop/1e6:.6g} MHz, {points} pts, '
    f'{avgc} avg, {ifbw:.0f} Hz IFBW, {power:.1f} dBm')
log(f'Coord map    : start ruler={fmt_mm(ruler_start)} lead={fmt_mm(lead_start)}, '
    f'ruler_step={ruler_step:+g}, lead_step={lead_step:+g}, end_ruler={ruler_end}')
log(f'Lead geometry: length={lead_length:g} mm, plot_lead=0..{lead_length:g} mm, '
    f'start_unmeasured={metadata["start_unmeasured_mm"]:g} mm, '
    f'planned_end_unmeasured={metadata.get("end_unmeasured_mm")}')
log(f'Metadata     : {metadata_path}')
log('-' * 72)
log(f'{"timestamp":<19}  {"pt":>3}  {"lead":>7}  {"ruler":>7}  '
    f'{"file":<22}  {"pk|S21|dB":>9}  {"@MHz":>8}  local')

###############################################################################
# 7) Interactive acquisition loop

baseline_db = None
baseline_desc = ''
saved = 0
i = 0

print('\nControls:  [Enter] = acquire & advance    b = back one    '
      'f = skip forward    q = quit\n')

while True:
    ruler = ruler_start + i * ruler_step
    lead_pos = lead_start + i * lead_step
    fname = f'{fmt_mm(lead_pos)}MM_{fmt_mm(ruler)}MM.s2p'
    vna_path = vna_dir + '\\' + fname
    local_path = os.path.join(raw_dir, fname)

    past_end = ruler_end is not None and (
        (ruler_step > 0 and ruler > ruler_end + 1e-6) or
        (ruler_step < 0 and ruler < ruler_end - 1e-6))
    tag = '   <-- past END point' if past_end else ''

    cmd = input(f'Point {i + 1}:  lead={fmt_mm(lead_pos)}mm  ruler={fmt_mm(ruler)}mm  ->  '
                f'{fname}{tag}\n  position probe, [Enter] to acquire (b/f/q): ').strip().lower()

    if cmd == 'q':
        break
    if cmd == 'b':
        i = max(0, i - 1)
        continue
    if cmd == 'f':
        i += 1
        continue
    if cmd != '':
        continue  # unrecognized key -> re-prompt same point

    try:
        s = sweep_s21()
        mirrored = save_s2p(vna_path, local_path)
    except RuntimeError as e:
        print(f'  ERROR: {e}\n  (nothing saved for this point - fix and retry)')
        continue

    mag_db = 20 * np.log10(np.abs(s))
    pk = int(np.argmin(np.abs(freq - REPORT_FREQ)))
    print(f'  saved {fname}   |S21| = {mag_db[pk]:.2f} dB @ {freq[pk]/1e6:.4g} MHz'
          f'{"" if mirrored else "   (VNA only)"}')

    log(f'{datetime.now():%Y-%m-%d %H:%M:%S}  {i + 1:>3}  {fmt_mm(lead_pos):>7}  '
        f'{fmt_mm(ruler):>7}  {fname:<22}  {mag_db[pk]:>9.2f}  {freq[pk]/1e6:>8.4g}  '
        f'{"yes" if mirrored else "NO"}')

    if baseline_db is None:
        baseline_db = mag_db
        baseline_desc = f'lead={fmt_mm(lead_pos)}mm ruler={fmt_mm(ruler)}mm'
        print('  (stored as baseline for the return-to-home drift check)')

    saved += 1
    i += 1

###############################################################################
# 8) Return-to-home drift check (SOP: |S21| must match baseline within 0.5 dB)

print(f'\n{saved} point(s) saved to {raw_dir}')

if baseline_db is not None:
    ans = input('\nReturn-to-home drift check? Move the probe back to the START '
                'coordinate, then [Enter] to measure (s = skip): ').strip().lower()
    if ans != 's':
        s = sweep_s21()
        now_db = 20 * np.log10(np.abs(s))
        d = np.abs(now_db - baseline_db)
        k = int(np.argmax(d))
        verdict = 'PASS' if d[k] <= 0.5 else 'FAIL'
        msg = (f'max |dS21| vs baseline ({baseline_desc}) = {d[k]:.2f} dB '
               f'@ {freq[k]/1e6:.4g} MHz  ->  {verdict}  (tol 0.5 dB)')
        print('  ' + msg)
        log('-' * 72)
        log(f'{datetime.now():%Y-%m-%d %H:%M:%S}  return-to-home  {msg}')

log(f'Finished     : {datetime.now():%Y-%m-%d %H:%M:%S}  ({saved} points)')

# Record the actual acquired span and its unmeasured distance from each physical
# end. This remains accurate if the operator quits before/after the planned END.
if saved:
    acquired_leads = []
    for name in os.listdir(raw_dir):
        m = re.match(r'(-?\d+(?:p\d+)?)MM_', name, re.IGNORECASE)
        if m and name.lower().endswith('.s2p'):
            acquired_leads.append(float(m.group(1).lower().replace('p', '.')))
    if acquired_leads:
        metadata['acquired_lead_min_mm'] = min(acquired_leads)
        metadata['acquired_lead_max_mm'] = max(acquired_leads)
        metadata['start_unmeasured_mm'] = metadata['acquired_lead_min_mm']
        metadata['end_unmeasured_mm'] = lead_length - metadata['acquired_lead_max_mm']
        with open(metadata_path, 'w', encoding='utf-8') as fh:
            json.dump(metadata, fh, indent=2)
            fh.write('\n')

###############################################################################
# 9) Auto-launch the z-axis plot for this run

if saved:
    import subprocess
    plot_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'plot_tf_zaxis.py')
    if os.path.isfile(plot_script):
        print('\nLaunching z-axis plot for this run...')
        try:
            # Non-blocking: the plot window opens while this script exits. The
            # run folder is passed explicitly so no file browser is needed here.
            subprocess.Popen([sys.executable, plot_script, local_dir])
        except Exception as ex:  # noqa: BLE001 - plotting is best-effort
            print(f'  Could not launch plot automatically ({ex}).')
            print(f'  Run it manually: python "{plot_script}" "{local_dir}"')
    else:
        print(f'\nPlot script not found next to this one ({plot_script}); '
              'skipping auto-plot.')

print('\nDone.')

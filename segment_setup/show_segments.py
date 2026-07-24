"""
Read-only dump of the live segment sweep on the Keysight E5063A ENA.

Prints a decoded, human-readable segment table (start/stop, center, points,
per-segment IF BW / power) plus the raw :SENSe1:SEGMent:DATA? string, so the
current acquisition grid can be reviewed before rebuilding it with
set_fast_segments.py in this folder.

Does NOT recall, modify, or save anything on the VNA.
"""

import numpy as np
import pyvisa

inst_address = 'TCPIP0::165.124.9.141::inst0::INSTR'

rm = pyvisa.ResourceManager()
inst = rm.open_resource(inst_address)
inst.timeout = 30000
print(inst.query('*IDN?').strip())

swtype = inst.query(':SENSe1:SWEep:TYPE?').strip()
print(f'\nSweep type: {swtype}')

# --- Global setup (segments may override IF BW / power per-segment) -----------
print(f'Global IF BW : {float(inst.query(":SENSe1:BWIDth?")):.0f} Hz')
print(f'Global power : {float(inst.query(":SOURce1:POWer?")):.1f} dBm')
print(f'Averaging    : {int(float(inst.query(":SENSe1:AVERage:COUNt?")))} '
      f'({"ON" if int(float(inst.query(":SENSe1:AVERage:STATe?"))) else "OFF"})')

freqs = np.array(inst.query_ascii_values(':SENSe1:FREQuency:DATA?'))
print(f'\nTotal measured points: {len(freqs)}')

# --- Decode the segment table -------------------------------------------------
# Layout: nflags, flags[startstop, ifbw, pow, delay, time], nseg, then nseg rows
# of [f0, f1, npts, (ifbw), (pow), (delay), (time)] where the optional columns
# are present only when the matching flag is 1. mode: flags[0] 0=start/stop, 1=center/span.
if swtype.upper().startswith('SEGM'):
    vals = [float(x) for x in inst.query(':SENSe1:SEGMent:DATA?').split(',')]
    nflags = int(vals[0])
    flags = [int(v) for v in vals[1:1 + nflags]]
    nseg = int(vals[1 + nflags])
    per_seg = 3 + sum(flags[1:])          # 3 base cols + optional ifbw/pow/delay/time
    mode = flags[0]

    print(f'\nSegments: {nseg}   (mode: {"center/span" if mode else "start/stop"}, '
          f'per-seg flags ifbw={flags[1]} pow={flags[2]} '
          f'delay={flags[3] if len(flags) > 3 else 0} '
          f'time={flags[4] if len(flags) > 4 else 0})')

    hdr = f'{"#":>3}  {"start(MHz)":>10}  {"stop(MHz)":>10}  {"center(MHz)":>11}  {"pts":>5}'
    if flags[1]:
        hdr += f'  {"IFBW(Hz)":>9}'
    if flags[2]:
        hdr += f'  {"pow(dBm)":>8}'
    print('\n' + hdr)
    print('-' * len(hdr))

    idx = 2 + nflags
    total = 0
    for s in range(nseg):
        seg = vals[idx:idx + per_seg]
        idx += per_seg
        if mode == 1:                     # center/span
            center, span = seg[0], seg[1]
            start, stop = center - span / 2, center + span / 2
        else:                             # start/stop
            start, stop = seg[0], seg[1]
            center = (start + stop) / 2
        npts = int(seg[2])
        total += npts
        col = 3
        row = (f'{s + 1:>3}  {start/1e6:>10.5g}  {stop/1e6:>10.5g}  '
               f'{center/1e6:>11.5g}  {npts:>5}')
        if flags[1]:
            row += f'  {seg[col]:>9.0f}'; col += 1
        if flags[2]:
            row += f'  {seg[col]:>8.1f}'
        print(row)

    print('-' * len(hdr))
    print(f'Sum of segment points: {total}')

    print('\nRaw :SENSe1:SEGMent:DATA? ->')
    print(inst.query(':SENSe1:SEGMent:DATA?').strip())
else:
    print('\nSweep is not in SEGMent mode; no segment table to show.')
    print('Frequencies (MHz):')
    for i, f in enumerate(freqs):
        print(f'  [{i}] {f/1e6:.6g}')

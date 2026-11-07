"""Start queued jobs one at a time, only while enough RAM is free (heavy jobs peak at ~3.5-4 GB each).
    python logs/memqueue.py  (from formula/)"""
import os
import subprocess
import sys
import time

import ctypes


class _MS(ctypes.Structure):
    _fields_ = [('dwLength', ctypes.c_ulong), ('dwMemoryLoad', ctypes.c_ulong), ('ullTotalPhys', ctypes.c_ulonglong),
                ('ullAvailPhys', ctypes.c_ulonglong), ('ullTotalPageFile', ctypes.c_ulonglong),
                ('ullAvailPageFile', ctypes.c_ulonglong), ('ullTotalVirtual', ctypes.c_ulonglong),
                ('ullAvailVirtual', ctypes.c_ulonglong), ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]


def free_gb():
    m = _MS(); m.dwLength = ctypes.sizeof(_MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPhys / 2**30

JOBS = [('wrr80_%d' % i, ['python', 'run_wrr.py', '80', '156', str(i)]) for i in (0, 1, 5, 6, 7)] + \
       [('wref50_%d' % i, ['python', 'run_wref.py', '50', '96', str(i)]) for i in (3, 5, 6, 7, 8, 9)]
NEED_GB, MAX_HEAVY = 7.0, 6


def heavy_running():
    out = subprocess.run(['powershell', '-NoProfile', '-Command',
                          "(Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object { $_.CommandLine "
                          "-match 'run_wrr.py 80|run_wref.py 50' }).Count"], capture_output=True, text=True).stdout
    return int(out.strip() or 0)


procs = []
for name, cmd in JOBS:
    while free_gb() < NEED_GB or heavy_running() >= MAX_HEAVY:
        time.sleep(60)
    log = open(os.path.join('logs', name + '.log'), 'w')
    procs.append(subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT))
    print(time.strftime('%H:%M'), 'started', name, flush=True)
    time.sleep(180)          # let it reach its working size before judging free memory again
for p in procs:
    p.wait()
print(time.strftime('%H:%M'), 'all done; exit codes', [p.returncode for p in procs], flush=True)

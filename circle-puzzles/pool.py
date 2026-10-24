"""Worker pools that leave the machine usable: half the cores and below-normal priority by default.

NPROC in the environment overrides the worker count.
"""
import ctypes
import os
from multiprocessing import Pool

NPROC = int(os.environ.get('NPROC', max(1, (os.cpu_count() or 2) // 2)))
BELOW_NORMAL_PRIORITY_CLASS = 0x4000


def kernel32():
    """kernel32 with HANDLE-sized return/argument types (the int default truncates the pseudo-handle)."""
    k = ctypes.windll.kernel32
    k.GetCurrentProcess.restype = ctypes.c_void_p
    k.SetPriorityClass.argtypes = [ctypes.c_void_p, ctypes.c_ulong]
    k.GetPriorityClass.argtypes = [ctypes.c_void_p]
    k.SetProcessInformation.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_ulong]
    return k


def low_priority(extra=None):
    if os.name == 'nt':
        k = kernel32()
        k.SetPriorityClass(k.GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS)
    if extra:
        extra()


def make_pool(n=None, initializer=None):
    return Pool(n or NPROC, initializer=low_priority, initargs=(initializer,))

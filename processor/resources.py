"""
CPU and memory the container may use (cgroup v2, then v1), and the peak memory of this
process and its children (salmon).
"""
from __future__ import annotations

import math
import resource


def _read(path: str) -> str | None:
    try:
        with open(path) as fh:
            return fh.read().strip()
    except OSError:
        return None


def cgroup_cpus() -> int | None:
    """CPUs from the cgroup quota, or None without a quota."""
    v2 = _read("/sys/fs/cgroup/cpu.max")
    if v2:
        quota, period = (v2.split() + ["100000"])[:2]
        if quota != "max":
            return max(1, math.ceil(int(quota) / int(period)))
        return None
    quota, period = _read("/sys/fs/cgroup/cpu/cpu.cfs_quota_us"), _read("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    if quota and period and int(quota) > 0 and int(period) > 0:
        return max(1, math.ceil(int(quota) / int(period)))
    return None


def threads(requested: int) -> int:
    """THREADS, or (0) the cgroup CPU quota; 2 when there is no quota."""
    return requested if requested > 0 else (cgroup_cpus() or 2)


def memory_limit() -> int | None:
    """The container's memory limit in bytes, or None without one."""
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        v = _read(path)
        if v and v != "max" and int(v) < 1 << 60:
            return int(v)
    return None


def cgroup_peak_mb() -> float | None:
    """The cgroup's peak memory (includes page cache), when the kernel reports it."""
    for path in ("/sys/fs/cgroup/memory.peak", "/sys/fs/cgroup/memory/memory.max_usage_in_bytes"):
        v = _read(path)
        if v and v.isdigit():
            return round(int(v) / 2**20, 1)
    return None


def peak_rss_mb() -> dict[str, float]:
    """Peak RSS of this process and of its largest child (Linux reports KiB, macOS bytes)."""
    import sys

    unit = 2**20 if sys.platform == "darwin" else 2**10
    return {"self": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / unit, 1),
            "children": round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / unit, 1)}

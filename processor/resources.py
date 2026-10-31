"""
The memory the container may use (cgroup v2, then v1), and the peak memory of this process and
its children (salmon). Thread counts come from safe_inputs.cpu_threads().
"""
from __future__ import annotations

import resource


def _read(path: str) -> str | None:
    try:
        with open(path) as fh:
            return fh.read().strip()
    except OSError:
        return None


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

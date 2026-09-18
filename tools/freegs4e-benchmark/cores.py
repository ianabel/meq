"""PHYSICAL cores, not logical ones -- the thread count every timed arm uses.

`os.cpu_count()` returns LOGICAL CPUs, which on this machine is 16: an AMD
Ryzen 7 3800X with 8 physical cores and 2 threads per core.  Asking for 16 on
FP-heavy work hands half the "threads" to SMT siblings that share an FPU with a
thread already saturating it, and the usual result is between nothing and a
loss.

THIS IS SHARED RATHER THAN COPIED into the four harnesses that need it, for the
reason `normalisationDifferenceStep` is shared in the solver: two arms timed on
two different thread counts do not compare, and a constant written out four
times drifts.

Override with MEQ_BENCH_THREADS when deliberately sweeping the axis --
tests/performance/npc-scan.sh exists to do exactly that, and its own rule is
that OMP and MKL move TOGETHER because the element-local work is nested inside
an OpenMP region where MKL suppresses its own threading.
"""
import os


def physical_cores(default=8):
    """Unique ( socket, core ) pairs from /proc/cpuinfo, or `default`."""
    override = os.environ.get("MEQ_BENCH_THREADS", "")
    if override.strip().isdigit():
        return int(override)
    try:
        seen, socket, core = set(), None, None
        with open("/proc/cpuinfo") as handle:
            for line in handle:
                if line.startswith("physical id"):
                    socket = line.split(":")[1].strip()
                elif line.startswith("core id"):
                    core = line.split(":")[1].strip()
                elif not line.strip() and socket is not None and core is not None:
                    seen.add((socket, core))
                    socket = core = None
        if socket is not None and core is not None:
            seen.add((socket, core))
        return len(seen) or default
    except OSError:
        return default


def thread_environment(base=None):
    """`base` with both thread axes set to the physical core count.

    BOTH AXES, MATCHED.  MKL_NUM_THREADS=1 is the CTEST setting -- pinned there
    because the bit-exactness assertions between assembly modes hold at that
    value and not above it -- and is NOT the production one.  Under
    AssemblyMode::Threaded the element-local dense work is nested inside an
    active OpenMP region where MKL suppresses its own threading, so
    MKL_NUM_THREADS costs the assembly nothing while the TRACE SOLVE runs
    outside every region and takes every thread it is given.  M-78 sizes it at
    a further 7 to 9 per cent; pinning it to 1 for a timed run hands that away
    and flatters whatever the run is raced against.
    """
    n = str(physical_cores())
    out = dict(os.environ if base is None else base)
    out["OMP_NUM_THREADS"] = n
    out["MKL_NUM_THREADS"] = n
    return out

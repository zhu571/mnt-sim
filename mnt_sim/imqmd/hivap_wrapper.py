"""HIVAP Fortran wrapper for ImQMD fragment de-excitation.

Writes input.dat with (A1=Afrag-1, Z1=Zfrag-1, A2=1, Z2=1) to form
CN = (Zfrag, Afrag), runs hivap.exe, and parses SIGXPN.DAT for
evaporation branching ratios.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from functools import lru_cache

import numpy as np

_HIVAP_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "hivap_fortran"))
_HIVAP_EXE = os.path.join(_HIVAP_DIR, "hivap.exe")

# Required data files (relative to HIVAP_DIR)
_REQUIRED_FILES = [
    "hivapein.dat",
    "Mexcess95.dat",
    "mlz.dat",
    "shell.dat",
]


def _prepare_workspace(temp_dir: str) -> None:
    """Symlink or copy required data files into a temp workspace."""
    import shutil

    for name in _REQUIRED_FILES:
        src = os.path.join(_HIVAP_DIR, name)
        if os.path.exists(src):
            os.symlink(src, os.path.join(temp_dir, name))
    # Copy the binary
    dst_exe = os.path.join(temp_dir, "hivap.exe")
    if not os.path.exists(dst_exe):
        shutil.copy2(_HIVAP_EXE, dst_exe)


def _parse_sigxpn(path: str, e_star: float) -> list[tuple[int, int, float]]:
    """Parse SIGXPN.DAT and return [(final_Z, final_A, probability), ...]."""
    target_ecm = round(float(e_star), 4)
    channels: dict[tuple[int, int], float] = {}
    total_er = 0.0
    total_fiss = 0.0

    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                np_evap = int(parts[0])
                nn_evap = int(parts[1])
                ecm = float(parts[2])
                sigma_er = float(parts[4])
                sigma_fiss = float(parts[5])
            except (ValueError, IndexError):
                continue

            if abs(ecm - target_ecm) > 0.01:
                continue

            total_er += sigma_er
            total_fiss += sigma_fiss
            if sigma_er > 0:
                key = (np_evap, nn_evap)
                channels[key] = channels.get(key, 0.0) + sigma_er

    total = total_er + total_fiss
    if total <= 0:
        return []

    results: list[tuple[int, int, float]] = []
    for (np_evap, nn_evap), sigma_er in channels.items():
        results.append((np_evap, nn_evap, sigma_er / total))
    # Add fission as None marker for probability
    if total_fiss > 0:
        results.append((-1, -1, total_fiss / total))
    return results


@lru_cache(maxsize=1024)
def run_hivap(
    z_frag: int,
    a_frag: int,
    e_star: float,
    hivap_dir: str | None = None,
) -> list[tuple[int, int, float]]:
    """Run HIVAP for a single hot fragment.

    Returns list of (final_Z, final_A, probability).
    (-1, -1, prob) represents fission.

    Parameters
    ----------
    z_frag : int
        Proton number of the fragment.
    a_frag : int
        Mass number of the fragment.
    e_star : float
        Excitation energy in MeV.
    hivap_dir : str, optional
        Path to hivap_fortran/ directory.
    """
    if e_star <= 0 or a_frag <= 4:
        return [(z_frag, a_frag, 1.0)]

    work_dir = hivap_dir or _HIVAP_DIR

    with tempfile.TemporaryDirectory() as tmpdir:
        _prepare_workspace(tmpdir)

        # Write input.dat: dummy target trick
        a1 = max(a_frag - 1, 1)
        z1 = max(z_frag - 1, 0)
        input_path = os.path.join(tmpdir, "input.dat")
        with open(input_path, "w") as f:
            f.write(f"{a1} {z1} 1 1\n")
            f.write(f"{e_star:.1f}\n")
            f.write("0\n")

        # Run HIVAP
        exe = os.path.join(tmpdir, "hivap.exe")
        try:
            subprocess.run(
                [exe],
                cwd=tmpdir,
                capture_output=True,
                timeout=30,
            )
        except (subprocess.TimeoutExpired, OSError):
            return [(z_frag, a_frag, 1.0)]

        # Parse output
        sigxpn_path = os.path.join(tmpdir, "SIGXPN.DAT")
        if not os.path.exists(sigxpn_path):
            return [(z_frag, a_frag, 1.0)]

        channels = _parse_sigxpn(sigxpn_path, e_star)
        if not channels:
            return [(z_frag, a_frag, 1.0)]

        # Convert (np, nn) to (Z_final, A_final)
        results = []
        for np_evap, nn_evap, prob in channels:
            if np_evap == -1:  # fission
                results.append((-1, -1, prob))
            else:
                final_z = z_frag - np_evap
                final_a = a_frag - np_evap - nn_evap
                if final_z >= 0 and final_a >= final_z:
                    results.append((final_z, final_a, prob))
        return results


def sample_residue(
    z_frag: int,
    a_frag: int,
    e_star: float,
    rng: np.random.Generator | None = None,
) -> tuple[int, int]:
    """Sample one evaporation residue from HIVAP branching ratios."""
    rng = rng or np.random.default_rng()
    channels = run_hivap(z_frag, a_frag, e_star)

    if not channels:
        return z_frag, a_frag

    # Use all channels including fission (-1, -1)
    zs = [z for z, a, p in channels]
    a_s = [a for z, a, p in channels]
    ps = [p for z, a, p in channels]
    total_p = sum(ps)
    if total_p <= 0:
        return z_frag, a_frag
    ps = [p / total_p for p in ps]

    idx = rng.choice(len(channels), p=ps)
    zf, af = zs[idx], a_s[idx]
    if zf == -1:  # fission
        return z_frag, a_frag  # keep original fragment as placeholder
    return zf, af


__all__ = ["run_hivap", "sample_residue"]

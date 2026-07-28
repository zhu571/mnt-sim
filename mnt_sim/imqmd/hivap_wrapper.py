"""HIVAP Fortran wrapper for ImQMD fragment de-excitation.

Writes input.dat with (A1=Afrag-1, Z1=Zfrag-1, A2=1, Z2=1) to form
CN = (Zfrag, Afrag), runs hivap.exe, and parses SIGXPN.DAT for
evaporation branching ratios.

The research report's reference two-stage scheme is imQMD+GEMINI with
aden_type=-23, imf_option=2, Z_imf_min=5 and a 500 fm/c switch time.
GEMINI is not bundled in this repository, so the HIVAP path (also an
established imQMD coupling, e.g. the N=126 MNT studies) is the completed
one here; the 500 fm/c switch time is set in reaction.run_imqmd_event.
If GEMINI becomes available, mirror this module's interface
(run_*/sample_residue) so reaction.py can select it by name.
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
    (-1, -1, prob) represents fission.  Returns an EMPTY list when HIVAP
    itself failed (binary error, timeout, missing/garbled output) so that
    callers can fall back to the local evaporation chain; the trivial
    no-decay cases (E*<=0 or A<=4) return [(z_frag, a_frag, 1.0)].

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

    # Round the excitation to the 0.1 MeV written into input.dat so the
    # lru_cache actually deduplicates calls (floats otherwise never match).
    e_star = round(float(e_star), 1)
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
            return []

        # Parse output
        sigxpn_path = os.path.join(tmpdir, "SIGXPN.DAT")
        if not os.path.exists(sigxpn_path):
            return []

        channels = _parse_sigxpn(sigxpn_path, e_star)
        if not channels:
            return []

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
    """Sample one evaporation residue from HIVAP branching ratios.

    Fallback policy: if HIVAP fails to produce branching ratios (binary
    missing, timeout, unparseable output) the fragment is de-excited with
    the local Weisskopf evaporation chain instead of being returned hot.
    A fission outcome (-1, -1) keeps the parent (Z, A) as a placeholder —
    EventFragmentRecord cannot represent two fission products yet.
    """
    rng = rng or np.random.default_rng()
    # Round before the call so the lru_cache on run_hivap deduplicates.
    channels = run_hivap(int(z_frag), int(a_frag), round(float(e_star), 1))

    if not channels:
        from .decay import evaporate_full

        return evaporate_full(z_frag, a_frag, e_star, rng=rng)

    # Use all channels including fission (-1, -1)
    zs = [z for z, a, p in channels]
    a_s = [a for z, a, p in channels]
    ps = [p for z, a, p in channels]
    total_p = sum(ps)
    if total_p <= 0:
        from .decay import evaporate_full

        return evaporate_full(z_frag, a_frag, e_star, rng=rng)
    ps = [p / total_p for p in ps]

    idx = rng.choice(len(channels), p=ps)
    zf, af = zs[idx], a_s[idx]
    if zf == -1:  # fission
        return z_frag, a_frag  # keep original fragment as placeholder
    return zf, af


__all__ = ["run_hivap", "sample_residue"]

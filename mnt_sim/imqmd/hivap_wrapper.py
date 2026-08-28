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
import warnings
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

# AME mass excess (MeV) of the proton, used to turn the mass-excess table
# into binding-energy differences (S_p = B(Z,A) - B(Z-1,A-1)).
_PROTON_MASS_EXCESS = 7.289

# Myers-Swiatecki (1967) liquid-drop parameters, mirroring msben.f so the
# fallback S_p matches HIVAP's own Q-value bookkeeping when the daughter
# nuclide is missing from the experimental mass table.
_MS_A1 = 15.4941
_MS_A2 = 17.9439
_MS_A3 = 0.7053
_MS_CAY1 = 1.15303
_MS_CAY2 = 0.0
_MS_CAY3 = 200.0
_MS_GAMMA = 1.7826
_MS_MN = 8.07144
_MS_MP = 7.28899


@lru_cache(maxsize=4)
def _mass_excess_table(hivap_dir: str) -> dict[tuple[int, int], float]:
    """Load the experimental mass excess table {(A, Z): MeV} from Mexcess95.dat.

    Column 4 of Mexcess95.dat is the Audi-Wapstra mass excess in MeV
    (e.g. U-238 -> 47.304); column 3 is a smooth-model value and is not
    used here.
    """
    table: dict[tuple[int, int], float] = {}
    path = os.path.join(hivap_dir, "Mexcess95.dat")
    if not os.path.exists(path):
        return table
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                a = int(parts[0])
                z = int(parts[1])
                me = float(parts[3])
            except ValueError:
                continue
            table[(a, z)] = me
    return table


def _ld_mass_excess(z: int, a: int) -> float:
    """Myers-Swiatecki (1967) liquid-drop mass excess, matching MSBEN/SMASS."""
    zz = float(z)
    un = float(a - z)
    aa = float(a)
    a3rt = aa ** (1.0 / 3.0)
    a2rt = a3rt * a3rt
    sym = ((un - zz) / aa) ** 2
    acor = 1.0 - _MS_GAMMA * sym
    return (
        _MS_MN * un
        + _MS_MP * zz
        - _MS_A1 * acor * aa
        + _MS_A2 * acor * a2rt
        + _MS_A3 * zz * zz / a3rt
        - _MS_CAY1 * zz * zz / aa
        - _MS_CAY2 * a2rt * np.exp(-_MS_CAY3 * sym)
    )


def separation_energy(z: int, a: int, hivap_dir: str | None = None) -> float | None:
    """Proton separation energy S_p(Z,A) = B(Z,A) - B(Z-1,A-1) in MeV.

    Uses the experimental mass excesses from Mexcess95.dat (the same
    Audi-Wapstra data HIVAP reads for its fusion Q-value), with a
    Myers-Swiatecki liquid-drop fallback for nuclei absent from the table.
    Returns None when the fragment is not a bound nucleus (Z < 1).
    """
    z, a = int(z), int(a)
    if z < 1 or a < 1 or a < z:
        return None
    table = _mass_excess_table(hivap_dir or _HIVAP_DIR)
    parent = table.get((a, z))
    daughter = table.get((a - 1, z - 1))
    if parent is not None and daughter is not None:
        return daughter + _PROTON_MASS_EXCESS - parent
    return _ld_mass_excess(z - 1, a - 1) + _PROTON_MASS_EXCESS - _ld_mass_excess(z, a)


def _check_hivap_ready(hivap_dir: str | None = None) -> list[str]:
    """Return the names of missing HIVAP runtime files (empty when all present).

    Verifies the binary plus the required data files exist before handing a
    run off to the Fortran code, so callers can fail fast (and fall back to
    the local evaporation chain) instead of discovering the breakage after a
    subprocess launch.
    """
    hivap_dir = os.path.normpath(hivap_dir) if hivap_dir else _HIVAP_DIR
    missing: list[str] = []
    if not os.path.isfile(os.path.join(hivap_dir, "hivap.exe")):
        missing.append("hivap.exe")
    for name in _REQUIRED_FILES:
        if not os.path.isfile(os.path.join(hivap_dir, name)):
            missing.append(name)
    return missing


def _prepare_workspace(temp_dir: str) -> None:
    """Symlink or copy required data files into a temp workspace."""
    import shutil

    for name in _REQUIRED_FILES:
        src = os.path.join(_HIVAP_DIR, name)
        if not os.path.exists(src):
            continue
        try:
            os.symlink(src, os.path.join(temp_dir, name))
        except OSError:
            # Symlinks unavailable (e.g. restricted FS): copy instead.
            shutil.copy2(src, os.path.join(temp_dir, name))
    # Copy the binary
    dst_exe = os.path.join(temp_dir, "hivap.exe")
    if not os.path.exists(dst_exe):
        try:
            shutil.copy2(_HIVAP_EXE, dst_exe)
        except OSError:
            warnings.warn(
                f"hivap.exe not copied to workspace {temp_dir!r}; HIVAP will fall back",
                RuntimeWarning,
                stacklevel=2,
            )


def _parse_sigxpn(path: str, target_ecm: float) -> list[tuple[int, int, float]]:
    """Parse SIGXPN.DAT and return [(final_Z, final_A, probability), ...].

    ``target_ecm`` is the Ecm value written into input.dat; SIGXPN.DAT
    repeats this column verbatim, so only rows matching it are kept.
    """
    target_ecm = round(float(target_ecm), 4)
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
    no-decay cases (E*<=0, A<=4, or E* below the proton separation energy
    S_p) return [(z_frag, a_frag, 1.0)].

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

    # HIVAP treats the input "Ecm" as the compound-nucleus excitation
    # energy: EXCIT = Ecm + Q, and the dummy-target trick
    # (Z-1,A-1) + p -> (Z,A) has Q = S_p(Z,A) (the proton separation
    # energy of the compound nucleus).  Without correction the fragment is
    # de-excited at E* + S_p ~ 6-8 MeV too high, over-predicting fission
    # (audit: U-238 E*=50 -> fiss 87.7%; with S_p subtracted E*_eff~42-44
    # -> ~77%).  Subtract S_p so the written Ecm reproduces the intended
    # fragment excitation.  Below the threshold the fragment cannot decay.
    sp = separation_energy(z_frag, a_frag)
    if sp is None or e_star - sp <= 0:
        return [(z_frag, a_frag, 1.0)]
    ecm = round(e_star - sp, 1)

    work_dir = hivap_dir or _HIVAP_DIR

    with tempfile.TemporaryDirectory() as tmpdir:
        _prepare_workspace(tmpdir)

        # Write input.dat: dummy target trick
        a1 = max(a_frag - 1, 1)
        z1 = max(z_frag - 1, 0)
        input_path = os.path.join(tmpdir, "input.dat")
        with open(input_path, "w") as f:
            f.write(f"{a1} {z1} 1 1\n")
            f.write(f"{ecm:.1f}\n")
            f.write("0\n")

        # Run HIVAP
        exe = os.path.join(tmpdir, "hivap.exe")
        try:
            completed = subprocess.run(
                [exe],
                cwd=tmpdir,
                capture_output=True,
                timeout=30,
            )
        except (subprocess.TimeoutExpired, OSError):
            return []

        # Parse output.  The SIGXPN.DAT Ecm column echoes the value written
        # into input.dat, so the parse must use the corrected ecm, not e_star.
        sigxpn_path = os.path.join(tmpdir, "SIGXPN.DAT")
        if completed.returncode != 0:
            # Tolerate a non-zero exit code as long as the evaporation table
            # was actually produced (HIVAP can exit non-zero after a full,
            # usable run).  Only fall back when no output exists at all.
            if not os.path.exists(sigxpn_path):
                return []
        if not os.path.exists(sigxpn_path):
            return []

        channels = _parse_sigxpn(sigxpn_path, ecm)
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
) -> tuple[int, int] | None:
    """Sample one evaporation residue from HIVAP branching ratios.

    Returns (final_Z, final_A) of the surviving residue, or None when the
    fragment fissioned: a fission event produces two fission products, not
    an evaporation residue, and EventFragmentRecord cannot represent two
    fragments.  Callers (``_record_fragment``) skip None so fissioned
    fragments do not pollute dsigma/dZ with the un-decayed parent nucleus.

    Fallback policy: if HIVAP fails to produce branching ratios (binary
    missing, timeout, unparseable output) the fragment is de-excited with
    the local Weisskopf evaporation chain instead of being returned hot.
    """
    rng = rng or np.random.default_rng()
    # Round before the call so the lru_cache on run_hivap deduplicates.
    e_star_round = round(float(e_star), 1)
    if z_frag > 4 and e_star_round > 0:
        missing = _check_hivap_ready()
        if missing:
            warnings.warn(
                "HIVAP runtime files missing ("
                + ", ".join(missing)
                + "); falling back to the local Weisskopf evaporation chain",
                RuntimeWarning,
                stacklevel=2,
            )
    channels = run_hivap(int(z_frag), int(a_frag), e_star_round)

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
    if zf == -1:  # fission: no evaporation residue, signal with None
        return None
    return zf, af


def residue_branches(
    z_frag: int,
    a_frag: int,
    e_star: float,
    rng: np.random.Generator | None = None,
) -> list[tuple[int, int, float]]:
    """Return every HIVAP evaporation-residue branch with its probability."""
    channels = run_hivap(int(z_frag), int(a_frag), round(float(e_star), 1))
    if channels:
        return [(z, a, p) for z, a, p in channels if z >= 0 and p > 0.0]

    from .decay import evaporate_full

    zf, af = evaporate_full(z_frag, a_frag, e_star, rng=rng)
    return [(zf, af, 1.0)]


__all__ = ["residue_branches", "run_hivap", "sample_residue", "separation_energy"]

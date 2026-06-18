"""
Stopping power and range calculations for ions in gas.

Implements:
  - Bethe-Bloch formula (high energy: E > 1 MeV/u)
  - Ziegler/SRIM parameterization (low energy: E < 1 MeV/u)
  - Lindhard-Scharff (very low energy: E < 25 keV/u)
  - Energy straggling (Bohr, Tschalär)

References:
  - Ziegler, Biersack, Ziegler, SRIM-2013
  - Bethe, Ann. Phys. 5 (1930) 325
  - Bohr, K. Dan. Vidensk. Selsk. Mat. Fys. Medd. 18 (1948) 8
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional


# Gas properties: (Z, A, density [g/cm^3] at STP)
GAS_PROPERTIES = {
    'He':  (2,  4.0026, 1.786e-4),
    'H2':  (1,  1.0079, 8.988e-5),
    'N2':  (7,  14.0067, 1.251e-3),
    'Ar':  (18, 39.948,  1.784e-3),
    'Ne':  (10, 20.1797, 9.002e-4),
    'CH4': (6,  16.043,  7.168e-4),  # approx
}


def get_gas_props(gas_name: str, pressure_mbar: float = None) -> dict:
    """Get gas properties at given pressure.

    Parameters
    ----------
    gas_name : str
        Gas species ('He', 'H2', 'N2', 'Ar', etc.)
    pressure_mbar : float
        Pressure in mbar. If None, returns STP values.

    Returns
    -------
    dict with Z, A, density [g/cm^3], number_density [cm^-3]
    """
    Z, A, rho_stp = GAS_PROPERTIES.get(gas_name, (2, 4.0026, 1.786e-4))

    if pressure_mbar is not None:
        # P V = n R T => density scales with pressure
        # STP = 1013.25 mbar
        rho = rho_stp * pressure_mbar / 1013.25
    else:
        rho = rho_stp

    # Number density [atoms/cm^3]
    n_density = rho / A * 6.022e23

    return {'Z': Z, 'A': A, 'density': rho, 'number_density': n_density}


class StoppingPower:
    """Ion stopping power in gas media.

    Parameters
    ----------
    ion_Z : int
        Atomic number of the ion
    ion_A : int
        Mass number of the ion
    gas : str
        Gas species
    pressure_mbar : float
        Gas pressure in mbar
    """

    def __init__(self, ion_Z: int, ion_A: int,
                 gas: str = 'He', pressure_mbar: float = 50):
        self.ion_Z = ion_Z
        self.ion_A = ion_A
        self.gas = gas
        self.pressure = pressure_mbar

        props = get_gas_props(gas, pressure_mbar)
        self.gas_Z = props['Z']
        self.gas_A = props['A']
        self.gas_density = props['density']
        self.n_density = props['number_density']

        # Mean excitation energy [eV]
        self.I_exc = self._mean_excitation_energy(gas)

        # Atomic constants
        self.m_e = 0.511  # MeV/c^2 (electron mass)
        self.N_A = 6.022e23

    def _mean_excitation_energy(self, gas: str) -> float:
        """Mean excitation energy I [eV] for a gas."""
        table = {'He': 41.8, 'H2': 19.2, 'N2': 82.0, 'Ar': 188.0,
                 'Ne': 137.0, 'CH4': 41.7}
        return table.get(gas, 50.0)

    def _beta_gamma(self, E_MeV: float) -> Tuple[float, float]:
        """Compute β = v/c and γ from energy [MeV]."""
        M = self.ion_A * 931.494  # MeV
        T = E_MeV
        E_total = M + T
        gamma = E_total / M
        beta = np.sqrt(1 - 1 / gamma**2) if gamma >= 1 else 0
        return beta, gamma

    def bethe_bloch(self, E_MeV: float) -> float:
        """Bethe-Bloch stopping power [MeV/(mg/cm^2)] for E > 1 MeV/u."""
        if E_MeV < 0.001:
            return 0.0

        beta, gamma = self._beta_gamma(E_MeV / self.ion_A)  # per nucleon

        if beta < 0.01:
            return self._lindhard_scharff(E_MeV)

        # Bethe-Bloch formula
        K = 0.307075  # MeV cm^2/mol (4π N_A r_e^2 m_e c^2)
        z = self.ion_Z
        Z_gas = self.gas_Z
        A_gas = self.gas_A

        # Maximum energy transfer in a single collision
        W_max = 2 * self.m_e * beta**2 * gamma**2 / (
            1 + 2 * gamma * self.m_e / (self.ion_A * 931.494)
            + (self.m_e / (self.ion_A * 931.494))**2
        )

        # Density effect (simplified - negligible for gases at low P)
        delta = 0.0

        # Shell correction (simplified)
        C = 0.0

        dE_dx = (K * z**2 * Z_gas / (A_gas * beta**2)
                 * (0.5 * np.log(2 * self.m_e * beta**2 * gamma**2 * W_max
                                 / self.I_exc**2)
                    - beta**2 - delta/2 - C/Z_gas))

        # Convert from MeV/(g/cm^2) to MeV/(mg/cm^2)
        dE_dx *= 1e-3

        return max(dE_dx, 0.0)

    def _lindhard_scharff(self, E_MeV: float) -> float:
        """Low-energy stopping (Lindhard-Scharff) [MeV/(mg/cm^2)]."""
        # Electronic stopping at low energies
        epsilon = E_MeV / self.ion_A * 1e3  # keV/u
        if epsilon <= 0:
            return 0.0

        # Lindhard-Scharff reduced energy
        eps_r = (epsilon * self.ion_A * self.gas_A
                 / (self.ion_Z * self.gas_Z
                    * (self.ion_Z**(2/3) + self.gas_Z**(2/3))**(3/2)
                    * (self.ion_A + self.gas_A))
                 * 0.8853 * 0.529e-8 * 1e4
                 * 931.494 / (1.44 * 1e-13))

        # Reduced stopping power
        s_e = 0.4 * eps_r**(0.5) if eps_r < 10 else 4.0 * eps_r**(0.5) / (1 + eps_r)

        # Convert to MeV/(mg/cm^2)
        N = self.n_density
        dE_dx = (s_e * 1.44 * self.ion_Z * self.gas_Z
                 * np.sqrt(self.ion_A * self.gas_A / (self.ion_A + self.gas_A))
                 / (self.ion_Z**(2/3) + self.gas_Z**(2/3))**(3/2)
                 * 1e-3 * 931.494)

        return max(dE_dx, 0.0)

    def ziegler_stopping(self, E_MeV: float) -> float:
        """Ziegler/SRIM parameterization covering full energy range.

        Uses empirical fits for He gas stopping.
        [MeV/(mg/cm^2)]
        """
        if E_MeV <= 0:
            return 0.0

        E_per_u = E_MeV / self.ion_A  # MeV/u

        if E_per_u > 1.0:
            return self.bethe_bloch(E_MeV)
        elif E_per_u > 0.025:
            # Intermediate: empirical polynomial in sqrt(E)
            s = np.sqrt(E_per_u)
            # Ziegler-like parameterization for He
            a = np.array([0.0, 0.5, -0.02, 0.001])  # fit coefficients
            dE = (a[0] + a[1]*s + a[2]*s**2 + a[3]*s**3)
            # Scale by ion Z
            dE *= (self.ion_Z / 2)**0.7
            return max(dE * 1e-2, 0.0)
        else:
            return self._lindhard_scharff(E_MeV)

    def __call__(self, E_MeV: float) -> float:
        """Convenience: call returns stopping power [MeV/(mg/cm^2)]."""
        return self.ziegler_stopping(E_MeV)

    def range_energy(self, E0_MeV: float, n_steps: int = 500) -> Tuple[float, np.ndarray]:
        """Compute range [cm] and energy-vs-distance curve.

        Parameters
        ----------
        E0_MeV : float
            Initial energy [MeV]
        n_steps : int
            Number of integration steps

        Returns
        -------
        range_cm : float
            Total range in cm
        trajectory : ndarray (2, n_steps)
            [distance_cm, energy_MeV]
        """
        dE = E0_MeV / n_steps
        dist = 0.0
        E = E0_MeV
        trajectory = np.zeros((2, n_steps))

        for i in range(n_steps):
            dEdx = self.ziegler_stopping(E)
            if dEdx <= 0:
                dEdx = 1e-10
            dx = dE / dEdx  # mg/cm^2
            dist += dx
            E -= dE
            trajectory[0, i] = dist / (self.gas_density * 1000)  # cm -> convert density g/cm^3 to mg/cm^3? 
            # Actually: dE/dx is in MeV/(mg/cm^2). 
            # dx in mg/cm^2. To get cm: dx / (density_cm3 * 1000) = dx_mgcm2 / (rho_gcm3 * 1000 mg/g * 1 cm)
            # density is g/cm^3, so: cm = (mg/cm^2) / (rho * 1000)
            trajectory[0, i] = dx / (self.gas_density * 1000)
            trajectory[1, i] = max(E, 0)
            if i > 0:
                trajectory[0, i] += trajectory[0, i-1]

            if E <= 0:
                trajectory[0, i:] = trajectory[0, i]
                trajectory[1, i:] = 0
                break

        return trajectory[0, -1], trajectory

    def energy_straggling(self, E_MeV: float, dx_mgcm2: float) -> float:
        """Energy straggling (Bohr + Tschalär) [MeV^2].

        Returns the variance Ω² in energy loss after traversing dx.

        Parameters
        ----------
        E_MeV : float
            Initial energy [MeV]
        dx_mgcm2 : float
            Path length [mg/cm^2]

        Returns
        -------
        omega2 : float
            Energy loss variance [MeV^2]
        """
        if dx_mgcm2 <= 0 or E_MeV <= 0:
            return 0.0

        beta, gamma = self._beta_gamma(E_MeV / self.ion_A)

        # Bohr straggling
        omega2_Bohr = (4 * np.pi * (1.44e-13*1e6)**2  # (e^2)^2 in MeV^2·cm
                       * self.gas_Z * self.n_density
                       * self.ion_Z**2
                       * dx_mgcm2 / self.gas_density)  # convert back

        # Tschalär correction for lower energies
        if beta < 0.1:
            L = np.log(2 * self.m_e * beta**2 / self.I_exc * 1e6)
            omega2 = omega2_Bohr * max(0.5, L / 10)
        else:
            omega2 = omega2_Bohr

        return max(omega2, 0.0)

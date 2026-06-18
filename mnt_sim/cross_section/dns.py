"""
DNS (Di-Nuclear System) model for multi-nucleon transfer reactions.

The DNS model treats the projectile-target system as two touching
nuclei with nucleon exchange governed by a Master Equation.

Physics reference:
  - Adamian et al., Phys. Rev. C 53 (1996) 2273
  - Feng et al., At. Data Nucl. Data Tables 157 (2024) 101645
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple

from mnt_sim.data import element_to_z, typical_mass_number

__all__ = ["DNSModel", "DNSResult"]


@dataclass
class DNSResult:
    """Container for DNS cross-section calculation results."""
    Z: np.ndarray          # Proton number grid
    N: np.ndarray          # Neutron number grid
    P: np.ndarray          # Probability matrix P(Z,N)
    E_cm: float            # Center-of-mass energy [MeV]
    sigma_total: float     # Total reaction cross section [mb]
    sigma_matrix: np.ndarray  # Partial cross sections [mb]

    def to_dataframe(self):
        """Export to pandas DataFrame for analysis."""
        import pandas as pd
        Z_grid, N_grid = np.meshgrid(self.Z, self.N, indexing='ij')
        df = pd.DataFrame({
            'Z': Z_grid.ravel(),
            'N': N_grid.ravel(),
            'A': (Z_grid + N_grid).ravel(),
            'sigma_mb': self.sigma_matrix.ravel()
        })
        return df


class DNSModel:
    """Di-Nuclear System model for MNT cross sections.

    Parameters
    ----------
    projectile : str
        Projectile symbol (e.g. 'Xe', 'Ca')
    target : str
        Target symbol (e.g. 'Pb', 'U')
    E_lab : float
        Laboratory beam energy [MeV/u] (per nucleon)
    """

    # Nucleon mass [MeV/c^2]
    M_N = 931.494  # MeV/u

    def __init__(self, projectile: str, target: str,
                 E_lab: float, E_loss: float = 0.0):
        self.projectile = projectile
        self.target = target
        self.E_lab = E_lab
        self.E_loss = E_loss

        # Look up nuclear properties
        self.Zp = self._element_to_Z(projectile)
        self.Zt = self._element_to_Z(target)
        self.Ap = self._element_to_A(projectile)
        self.At = self._element_to_A(target)

        # Kinematics
        self.E_cm = self._compute_E_cm()
        self.V_C = self._compute_coulomb_barrier()
        self.E_cm_eff = self.E_cm - self.V_C  # Effective energy above barrier

    def _element_to_Z(self, symbol: str) -> int:
        """Convert element symbol to atomic number."""
        return element_to_z(symbol)

    def _element_to_A(self, symbol: str) -> int:
        """Get typical mass number for an element (most stable isotope)."""
        return typical_mass_number(symbol)

    def _compute_E_cm(self) -> float:
        """Compute center-of-mass energy from lab energy [MeV]."""
        E_total = self.Ap * self.E_lab
        mu = (self.Ap * self.At) / (self.Ap + self.At)
        return E_total * mu / self.Ap  # = E_lab * At / (Ap + At) * Ap

    def _compute_coulomb_barrier(self) -> float:
        """Coulomb barrier [MeV] using Bass potential."""
        r0 = 1.18  # fm
        R_int = r0 * (self.Ap**(1/3) + self.At**(1/3))
        return 1.44 * self.Zp * self.Zt / R_int

    def _nucleon_transfer_rate(self, Z1: int, N1: int,
                                Z2: int, N2: int,
                                direction: str = 'forward') -> float:
        """Nucleon transfer rate Λ [MeV/hbar] in the DNS model.

        Uses a simplified window formula based on the
        nuclear potential at the contact point.
        """
        # Mass numbers
        A1 = Z1 + N1
        A2 = Z2 + N2

        # Reduced mass at contact [MeV/c^2]
        mu = self.M_N * (A1 * A2) / (A1 + A2)

        # Q-value for the transfer
        # Approximate using binding energy systematics
        B1 = self._binding_energy(Z1, N1)
        B2 = self._binding_energy(Z2, N2)
        Q = (B1 + B2) - (self._binding_energy(Z1-1, N1)
                         + self._binding_energy(Z2+1, N2)) if 'proton' in direction else \
            (B1 + B2) - (self._binding_energy(Z1, N1-1)
                         + self._binding_energy(Z2, N2+1))

        # Effective excitation energy (use minimum to handle sub-barrier)
        E_eff = max(self.E_cm_eff, 0.5)
        T = np.sqrt(8.0 / (A1 + A2) * E_eff)

        # Fermi energy [MeV]
        E_F = 38.0

        # Transfer rate (window formula, simplified)
        hbar = 6.582119569e-22  # MeV·s
        omega = 5.0  # MeV, typical phonon energy

        # Rate depends on available phase space
        if Q < 0:
            # Endothermic - suppressed
            rate = omega * np.exp(Q / T)
        else:
            # Exothermic
            rate = omega * np.exp(-Q / T) if Q > 2 * T else omega

        # Mass asymmetry factor
        asym_factor = 4 * A1 * A2 / (A1 + A2)**2
        rate *= asym_factor

        return max(rate, 1e-10) * 1e21  # Scale to physical units (s^-1)

    def _binding_energy(self, Z: int, N: int) -> float:
        """Approximate binding energy [MeV] using semi-empirical mass formula."""
        A = Z + N
        if A <= 0:
            return 0.0

        # Semi-empirical mass formula coefficients [MeV]
        a_v = 15.75
        a_s = 17.80
        a_c = 0.711
        a_a = 23.70
        a_p = 11.18

        B = (a_v * A
             - a_s * A**(2/3)
             - a_c * Z * (Z - 1) / A**(1/3)
             - a_a * (A - 2*Z)**2 / A
             + a_p * (-1)**Z * (-1)**(A-Z) / A**(1/2) * (A % 2 == 0))

        return B

    def _diffusion_coefficient(self) -> float:
        """Nucleon diffusion coefficient D."""
        E_eff = max(self.E_cm_eff, 0.5)
        T = np.sqrt(8.0 * E_eff / (self.Ap + self.At))
        return 0.5 * T  # MeV, simplified

    def calculate(self, delta_Z_range: Tuple[int, int] = (-4, 4),
                  delta_N_range: Tuple[int, int] = (-6, 6),
                  n_steps: int = 200) -> DNSResult:
        """Calculate MNT cross sections using the DNS model.

        Parameters
        ----------
        delta_Z_range : (int, int)
            Range of proton transfer (min, max)
        delta_N_range : (int, int)
            Range of neutron transfer (min, max)
        n_steps : int
            Number of time steps for master equation evolution

        Returns
        -------
        DNSResult with probability and cross-section matrices
        """
        Z_min = int(delta_Z_range[0])
        Z_max = int(delta_Z_range[1])
        N_min = int(delta_N_range[0])
        N_max = int(delta_N_range[1])

        nZ = Z_max - Z_min + 1
        nN = N_max - N_min + 1
        if nZ <= 0 or nN <= 0:
            raise ValueError("delta ranges must be increasing")
        if not (Z_min <= 0 <= Z_max and N_min <= 0 <= N_max):
            raise ValueError("delta ranges must include the initial channel (0, 0)")
        if n_steps <= 0:
            raise ValueError("n_steps must be positive")

        Z_vals = np.arange(Z_min, Z_max + 1)
        N_vals = np.arange(N_min, N_max + 1)

        # Initial distribution: delta function at (0,0)
        P = np.zeros((nZ, nN))
        P[0 - Z_min, 0 - N_min] = 1.0

        # Time integration
        dt = 0.5 / n_steps  # in scaled time units

        for step in range(n_steps):
            P_new = P.copy()
            for i in range(nZ):
                for j in range(nN):
                    Z_tr = Z_vals[i]  # net protons transferred to proj
                    N_tr = N_vals[j]  # net neutrons transferred to proj

                    # Current DNS configuration
                    Z1 = self.Zp + Z_tr
                    N1 = self.Ap - self.Zp + N_tr
                    Z2 = self.Zt - Z_tr
                    N2 = self.At - self.Zt - N_tr

                    if Z1 <= 0 or Z2 <= 0 or N1 <= 0 or N2 <= 0:
                        continue

                    # Transfer rates for +1p, -1p, +1n, -1n
                    rate_p_plus = self._nucleon_transfer_rate(Z1, N1, Z2, N2, 'proton_forward')
                    rate_p_minus = self._nucleon_transfer_rate(Z2, N2, Z1, N1, 'proton_backward')
                    rate_n_plus = self._nucleon_transfer_rate(Z1, N1, Z2, N2, 'neutron_forward')
                    rate_n_minus = self._nucleon_transfer_rate(Z2, N2, Z1, N1, 'neutron_backward')

                    # Gain terms from neighboring channels
                    if i > 0 and Z1 > 0:
                        P_new[i, j] += dt * rate_p_plus * P[i-1, j]
                    if i < nZ-1 and Z2 > 0:
                        P_new[i, j] += dt * rate_p_minus * P[i+1, j]
                    if j > 0 and N1 > 0:
                        P_new[i, j] += dt * rate_n_plus * P[i, j-1]
                    if j < nN-1 and N2 > 0:
                        P_new[i, j] += dt * rate_n_minus * P[i, j+1]

                    # Loss term
                    loss = dt * (rate_p_plus + rate_p_minus
                                 + rate_n_plus + rate_n_minus) * P[i, j]
                    P_new[i, j] -= loss

            # Normalize
            total = np.sum(P_new)
            if total > 0:
                P_new /= total
            P = P_new

        # Total reaction cross section [mb]
        # Using the Wong formula for heavy-ion fusion
        hbar_omega = 4.0  # MeV, curvature of barrier
        sigma_fusion = (1.44 * self.Zp * self.Zt
                        / (self.E_cm * hbar_omega)
                        * np.log(1 + np.exp(2 * np.pi
                                            * (self.E_cm - self.V_C)
                                            / hbar_omega))
                        * 10)  # approx conversion to mb

        # Partition among transfer channels
        sigma_matrix = sigma_fusion * P

        return DNSResult(
            Z=Z_vals, N=N_vals, P=P,
            E_cm=self.E_cm,
            sigma_total=sigma_fusion,
            sigma_matrix=sigma_matrix
        )

    def angular_distribution(self, theta_lab: np.ndarray,
                             result: DNSResult) -> np.ndarray:
        """Compute angular distribution for each channel.

        Uses a simplified diffraction model. The angular spread
        depends on the transferred mass and Q-value.

        Parameters
        ----------
        theta_lab : ndarray
            Laboratory angles [degrees]
        result : DNSResult
            Calculated cross-section result

        Returns
        -------
        d2sigma_dOmega_dE : ndarray
            Double-differential cross section [mb/sr/MeV]
        """
        # Convert to radians
        theta = np.radians(theta_lab)

        # Grazing angle
        a = 1.18 * (self.Ap**(1/3) + self.At**(1/3))
        D = a + 1.44 * self.Zp * self.Zt / (2 * self.E_cm)
        theta_graz = 2 * np.arcsin(1 / (2 * self.E_cm * D
                                        / (1.44 * self.Zp * self.Zt)
                                        + 1))

        # Angular width depends on transferred angular momentum
        sigma_theta = 3.0  # degrees, typical width

        # Sigma-weighted angular distribution
        dist = np.exp(-(theta_lab - np.degrees(theta_graz))**2
                      / (2 * sigma_theta**2))
        dist_sum = np.sum(dist)
        if dist_sum > 0:
            dist /= dist_sum

        # Energy distribution: Gaussian with width from Q-value smearing
        nZ, nN = result.Z.shape[0], result.N.shape[0]
        E_step = 1.0  # MeV
        nE = 50
        E_vals = np.linspace(0.5 * result.E_cm, 1.5 * result.E_cm, nE)

        d2sigma = np.zeros((nZ, nN, len(theta), nE))
        for i in range(nZ):
            for j in range(nN):
                if result.sigma_matrix[i, j] < 1e-10:
                    continue
                # Energy spread from Q-value distribution
                sigma_E = 2.0 + 0.1 * abs(result.Z[i]) + 0.1 * abs(result.N[j])
                for k, t in enumerate(theta):
                    ang_factor = np.exp(-(np.degrees(t) - np.degrees(theta_graz)
                                          - result.Z[i] * 0.5)**2
                                        / (2 * sigma_theta**2))
                    d2sigma[i, j, k, :] = (result.sigma_matrix[i, j]
                                           * ang_factor
                                           * np.exp(-(E_vals - result.E_cm)**2
                                                    / (2 * sigma_E**2))
                                           / (np.sqrt(2*np.pi) * sigma_E))
                    norm = np.sum(d2sigma[i, j, k, :])
                    if norm > 0:
                        d2sigma[i, j, k, :] /= norm

        return d2sigma

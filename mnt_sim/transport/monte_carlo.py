"""
Monte Carlo transport simulation of recoil ions through gas cell.

Simulates the trajectory of reaction products from MNT reactions
as they slow down in the gas cell, tracking:
  - Energy loss (via stopping power)
  - Energy straggling (Bohr/Tschalär)
  - Multiple Coulomb scattering
  - Spatial deposition distribution
  - Escape probability from the gas cell

Method: Continuous Slowing Down Approximation (CSDA) + Monte Carlo steps.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
from .stopping import StoppingPower
from .gas_cell import GasCell

__all__ = ["ParticleState", "TransportMC", "TransportResult", "run_doublediff_example"]


@dataclass
class ParticleState:
    """State of a single particle during simulation."""
    x: float          # mm
    y: float          # mm
    z: float          # mm
    energy: float     # MeV
    theta: float      # deg (polar angle from z-axis)
    phi: float        # deg (azimuthal angle)
    alive: bool = True
    deposited: bool = False
    escaped: bool = False
    stopped: bool = False
    track: List[Tuple[float, float, float, float]] = field(default_factory=list)

    def __post_init__(self):
        if not self.track:
            self.track = [(self.x, self.y, self.z, self.energy)]


@dataclass
class TransportResult:
    """Result of a Monte Carlo transport simulation."""
    n_particles: int
    n_deposited: int
    n_escaped: int
    n_stopped: int
    deposition_positions: np.ndarray  # (N, 3) array [mm]
    deposition_energies: np.ndarray   # (N,) array [MeV]
    escape_energies: np.ndarray       # (N,) array [MeV]
    final_energies: np.ndarray        # (N,) array [MeV]
    initial_energies: np.ndarray      # (N,) array [MeV]

    @property
    def deposition_efficiency(self) -> float:
        """Fraction of particles that stopped inside the cell."""
        return self.n_deposited / max(1, self.n_particles)

    @property
    def escape_fraction(self) -> float:
        """Fraction that escaped from the cell."""
        return self.n_escaped / max(1, self.n_particles)

    def summary(self) -> str:
        """Print a summary of the simulation."""
        mean_dep_e = (
            np.mean(self.deposition_energies)
            if len(self.deposition_energies) > 0 else 0.0
        )
        return (
            f"Transport Simulation Summary\n"
            f"{'='*40}\n"
            f"Total particles:      {self.n_particles}\n"
            f"Deposited in cell:    {self.n_deposited} "
            f"({self.deposition_efficiency*100:.1f}%)\n"
            f"Escaped cell:         {self.n_escaped} "
            f"({self.escape_fraction*100:.1f}%)\n"
            f"Mean deposition E:    {mean_dep_e:.2f} MeV"
        )


class TransportMC:
    """Monte Carlo simulation of ion transport in gas cell.

    Parameters
    ----------
    gas_cell : GasCell
        Gas cell geometry and conditions
    ion_Z : int
        Atomic number of the ion
    ion_A : int
        Atomic number of the ion
    energy : float
        Initial ion energy [MeV]
    n_particles : int
        Number of particles to simulate
    seed : int, optional
        Random seed for reproducibility
    """

    # Radiation length for multiple scattering [g/cm^2]
    X0_VALUES = {'He': 94.3, 'H2': 63.0, 'N2': 37.9, 'Ar': 19.5}

    def __init__(self, gas_cell: GasCell,
                 ion_Z: int, ion_A: int,
                 energy: float = None,
                 n_particles: int = 1000,
                 sigma_theta: float = 3.0,  # deg, angular spread from reaction
                 seed: int = None):
        if ion_Z <= 0:
            raise ValueError("ion_Z must be positive")
        if ion_A <= 0:
            raise ValueError("ion_A must be positive")
        if n_particles <= 0:
            raise ValueError("n_particles must be positive")
        if sigma_theta < 0:
            raise ValueError("sigma_theta must be non-negative")
        self.cell = gas_cell
        self.ion_Z = ion_Z
        self.ion_A = ion_A
        self.energy = energy
        self.n_particles = n_particles
        self.sigma_theta = sigma_theta

        # Stopping power calculator
        self.stopping = StoppingPower(
            ion_Z=ion_Z, ion_A=ion_A,
            gas=gas_cell.gas,
            pressure_mbar=gas_cell.pressure_mbar
        )

        if seed is not None:
            np.random.seed(seed)

    def _multiple_scattering(self, E: float, dx_cm: float) -> Tuple[float, float]:
        """Multiple Coulomb scattering angle [deg] for step dx.

        Highland formula for rms scattering angle.
        """
        if E <= 0 or dx_cm <= 0:
            return 0.0, 0.0

        # Ion momentum
        M = self.ion_A * 931.494  # MeV
        p = np.sqrt(E * (E + 2 * M))  # MeV/c
        beta = np.sqrt(1 - (M / (M + E))**2)

        # Radiation length of gas [cm]
        rho = self.cell.gas_density()
        X0 = self.X0_VALUES.get(self.cell.gas, 94.3)
        X0_cm = X0 / rho if rho > 0 else 1e10

        # Highland formula
        theta_rms = (13.6 / (beta * p * 1e-3)  # p in MeV/c -> GeV/c
                     * self.ion_Z * np.sqrt(dx_cm / X0_cm)
                     * (1 + 0.038 * np.log(dx_cm / X0_cm)))

        theta_rms = np.degrees(theta_rms)  # convert to degrees
        theta_rms = max(theta_rms, 0.001)  # minimum spread

        theta_scat = np.random.normal(0, theta_rms)
        phi_scat = np.random.uniform(0, 360)

        return theta_scat, phi_scat

    def _step_particle(self, state: ParticleState,
                       step_size_mm: float = 0.5) -> ParticleState:
        """Advance one particle by one step."""
        if not state.alive:
            return state

        E = state.energy
        if E <= 0:
            state.alive = False
            state.stopped = True
            state.deposited = self.cell.is_within_cell(
                state.x, state.y, state.z
            )
            return state

        # Convert step to mg/cm^2
        rho = self.cell.gas_density()
        dx_mgcm2 = step_size_mm * 0.1 * rho * 1000  # mm -> cm -> mg/cm^2

        # Energy loss in this step
        dEdx = self.stopping(E)  # MeV/(mg/cm^2)
        dE = dEdx * dx_mgcm2

        # Energy straggling
        omega2 = self.stopping.energy_straggling(E, dx_mgcm2)
        dE_strag = np.random.normal(0, np.sqrt(omega2)) if omega2 > 0 else 0

        # Update energy
        new_E = E - dE + dE_strag
        new_E = max(new_E, 0.0)

        # Multiple scattering
        dx_cm = step_size_mm * 0.1
        theta_scat, phi_scat = self._multiple_scattering(E, dx_cm)

        # Update direction
        theta_rad = np.radians(state.theta)
        phi_rad = np.radians(state.phi)

        # Small-angle approximation for polar scattering
        theta_rad += np.radians(theta_scat)
        phi_rad += np.radians(phi_scat)

        # Move particle
        dx = step_size_mm * np.sin(theta_rad) * np.cos(phi_rad)
        dy = step_size_mm * np.sin(theta_rad) * np.sin(phi_rad)
        dz = step_size_mm * np.cos(theta_rad)

        state.x += dx
        state.y += dy
        state.z += dz

        # Check boundaries
        if not self.cell.is_within_cell(state.x, state.y, state.z):
            state.alive = False
            state.escaped = True
            # Record escape energy as final energy
            state.energy = new_E
        else:
            state.energy = new_E

        state.track.append((state.x, state.y, state.z, state.energy))

        # Check if stopped
        if new_E <= 0:
            state.alive = False
            state.stopped = True
            state.deposited = self.cell.is_within_cell(
                state.x, state.y, state.z
            )

        return state

    def run(self, energy: float = None,
            n_particles: int = None,
            step_size_mm: float = 0.5,
            progress: bool = False,
            entrance_pos: Tuple[float, float, float] = (0, 0, 0)) -> TransportResult:
        """Run the Monte Carlo simulation.

        Parameters
        ----------
        energy : float, optional
            Initial energy [MeV]. Overrides constructor value.
        n_particles : int, optional
            Number of particles. Overrides constructor value.
        step_size_mm : float
            Step size [mm]
        progress : bool
            Show progress bar (slow)
        entrance_pos : (float, float, float)
            Entrance position (x, y, z) [mm]

        Returns
        -------
        TransportResult
        """
        E0 = energy if energy is not None else self.energy
        N = n_particles if n_particles is not None else self.n_particles
        if E0 is None:
            raise ValueError("initial energy must be provided")
        if E0 < 0:
            raise ValueError("initial energy must be non-negative")
        if N <= 0:
            raise ValueError("n_particles must be positive")
        if step_size_mm <= 0:
            raise ValueError("step_size_mm must be positive")

        dep_positions = []
        dep_energies = []
        esc_energies = []
        final_energies = []

        for _ in range(N):
            # Initial angular distribution (from reaction kinematics)
            theta0 = abs(np.random.normal(0, self.sigma_theta))
            phi0 = np.random.uniform(0, 360)

            state = ParticleState(
                x=entrance_pos[0], y=entrance_pos[1], z=entrance_pos[2],
                energy=E0, theta=theta0, phi=phi0
            )

            while state.alive:
                self._step_particle(state, step_size_mm)

            final_energies.append(state.energy)
            if state.deposited:
                dep_positions.append([state.x, state.y, state.z])
                dep_energies.append(state.energy)
            if state.escaped:
                esc_energies.append(state.energy)

        # Convert to arrays
        dep_arr = np.array(dep_positions) if dep_positions else np.zeros((0, 3))
        dep_e_arr = np.array(dep_energies) if dep_energies else np.zeros(0)
        esc_e_arr = np.array(esc_energies) if esc_energies else np.zeros(0)
        final_arr = np.array(final_energies)

        return TransportResult(
            n_particles=N,
            n_deposited=len(dep_positions),
            n_escaped=len(esc_energies),
            n_stopped=N - len(esc_energies),
            deposition_positions=dep_arr,
            deposition_energies=dep_e_arr,
            escape_energies=esc_e_arr,
            final_energies=final_arr,
            initial_energies=np.full(N, E0)
        )

    def scan_pressure(self, pressures: List[float],
                      energy: float = None,
                      n_particles: int = 500) -> Dict[float, TransportResult]:
        """Scan gas cell performance as a function of pressure.

        Parameters
        ----------
        pressures : list of float
            Pressures [mbar] to scan
        energy : float, optional
            Ion energy [MeV]
        n_particles : int
            Particles per pressure point

        Returns
        -------
        dict: pressure -> TransportResult
        """
        results = {}
        for P in pressures:
            cell_mod = GasCell(
                gas=self.cell.gas,
                pressure_mbar=P,
                length_mm=self.cell.length_mm,
                diameter_mm=self.cell.diameter_mm,
                temperature_K=self.cell.temperature_K,
                geometry=self.cell.geometry
            )
            mc = TransportMC(
                cell_mod, self.ion_Z, self.ion_A,
                energy=energy or self.energy,
                n_particles=n_particles,
                sigma_theta=self.sigma_theta
            )
            results[P] = mc.run(progress=False)
        return results

    def scan_length(self, lengths: List[float],
                    energy: float = None,
                    n_particles: int = 500) -> Dict[float, TransportResult]:
        """Scan gas cell performance as a function of length."""
        results = {}
        for L in lengths:
            cell_mod = GasCell(
                gas=self.cell.gas,
                pressure_mbar=self.cell.pressure_mbar,
                length_mm=L,
                diameter_mm=self.cell.diameter_mm
            )
            mc = TransportMC(
                cell_mod, self.ion_Z, self.ion_A,
                energy=energy or self.energy,
                n_particles=n_particles,
                sigma_theta=self.sigma_theta
            )
            results[L] = mc.run(progress=False)
        return results


def run_doublediff_example():
    """Example: full workflow from cross section to transport."""
    print("=" * 60)
    print("MNT Cross Section → Gas Cell Transport Pipeline")
    print("=" * 60)

    # 1. Cross section
    from mnt_sim.cross_section.dns import DNSModel
    dns = DNSModel("Xe", "Pb", E_lab=8.0)
    result = dns.calculate(delta_Z_range=(-4, 4), delta_N_range=(-6, 6))
    print(f"\nTotal MNT cross section: {result.sigma_total:.1f} mb")
    df = result.to_dataframe()
    top_channels = df.nlargest(5, 'sigma_mb')
    print("Top transfer channels:")
    for _, row in top_channels.iterrows():
        print(f"  Z={int(row['Z']):+d} N={int(row['N']):+d} "
              f"A={int(row['A']):+d}  σ={row['sigma_mb']:.2f} mb")

    # 2. Double-differential -> energy-angle distributions
    theta = np.arange(0, 30, 2)
    d2sigma = dns.angular_distribution(theta, result)

    # 3. Gas cell transport for a typical channel
    from mnt_sim.transport.gas_cell import GasCell
    cell = GasCell(gas='He', pressure_mbar=50, length_mm=200, diameter_mm=40)
    print(f"\nGas cell: {cell}")

    # For the most probable channel
    top_row = top_channels.iloc[0]
    dZ, dN = int(top_row['Z']), int(top_row['N'])
    A_ion = dns.Ap + dZ + dN
    Z_ion = dns.Zp + dZ
    print(f"Simulating transport: Z={Z_ion}, A={A_ion}, E~{dns.E_cm:.1f} MeV")

    mc = TransportMC(cell, Z_ion, A_ion, energy=dns.E_cm, n_particles=2000)
    transport_result = mc.run()
    print(transport_result.summary())

    return result, transport_result

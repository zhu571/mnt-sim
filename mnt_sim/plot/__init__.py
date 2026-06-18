"""
Visualization utilities for MNT simulation results.

Provides publication-quality plotting for:
  - Cross-section matrices (Z-N heatmaps)
  - Double-differential cross sections (E-θ maps)
  - Gas cell deposition distributions
  - Energy loss curves
  - Parameter scans (pressure/length vs efficiency)
"""

import numpy as np
from typing import Optional, Dict, List, Tuple
from pathlib import Path


def plot_cross_section_matrix(result, ax=None, title: str = None,
                              save_path: str = None):
    """Plot the cross-section matrix in the Z-N plane.

    Parameters
    ----------
    result : DNSResult or GrazingResult
    ax : matplotlib.axes.Axes, optional
    title : str, optional
    save_path : str, optional
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))

    Z, N = np.meshgrid(result.Z, result.N, indexing='ij')
    A = Z + N

    # Mask zero entries for cleaner plot
    sigma = result.sigma_matrix.copy()
    sigma[sigma < 1e-10] = np.nan

    im = ax.pcolormesh(Z, N, sigma, norm=LogNorm(vmin=1e-3, vmax=sigma.max()),
                       cmap='hot', shading='auto')
    plt.colorbar(im, ax=ax, label=r'$\sigma$ (mb)')

    ax.set_xlabel('Proton transfer ΔZ')
    ax.set_ylabel('Neutron transfer ΔN')
    ax.set_title(title or 'MNT Cross-Section Matrix')
    ax.grid(True, alpha=0.3)

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()


def plot_angle_energy_map(d2sigma, theta_deg, E_MeV, channel_idx,
                          ax=None, save_path: str = None):
    """Plot double-differential cross section d²σ/dΩdE.

    Parameters
    ----------
    d2sigma : ndarray (nZ, nN, nTheta, nE)
        Double-differential cross section
    theta_deg : ndarray
        Angles [degrees]
    E_MeV : ndarray
        Energies [MeV]
    channel_idx : (int, int)
        (iZ, iN) channel index in d2sigma
    ax : matplotlib.axes.Axes, optional
    save_path : str, optional
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))

    iZ, iN = channel_idx
    data = d2sigma[iZ, iN, :, :]

    Theta, E = np.meshgrid(theta_deg, E_MeV, indexing='ij')
    im = ax.pcolormesh(Theta, E, data, shading='auto', cmap='viridis')
    plt.colorbar(im, ax=ax, label=r'$d^2\sigma/d\Omega dE$ (mb/sr/MeV)')

    ax.set_xlabel('Laboratory Angle (deg)')
    ax.set_ylabel('Energy (MeV)')
    ax.set_title(f'Double-Differential Cross Section\n'
                 f'ΔZ={d2sigma.shape[0]//2 -iZ:+d}, '
                 f'ΔN={d2sigma.shape[1]//2 - iN:+d}')

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()


def plot_deposition_distribution(result: 'TransportResult',
                                 ax=None, title: str = None,
                                 save_path: str = None,
                                 show_histograms: bool = True):
    """Plot the spatial deposition distribution in the gas cell.

    Parameters
    ----------
    result : TransportResult
    ax : matplotlib.axes.Axes, optional (3D or list of 2D)
    title : str, optional
    save_path : str, optional
    show_histograms : bool
        Show projections in x, y, z
    """
    import matplotlib.pyplot as plt

    if show_histograms:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        ax_xy = axes[0, 0]
        ax_xz = axes[0, 1]
        ax_hist_z = axes[1, 0]
        ax_hist_e = axes[1, 1]

        # XY projection
        if len(result.deposition_positions) > 0:
            xy = result.deposition_positions
            ax_xy.scatter(xy[:, 0], xy[:, 1], s=1, alpha=0.5, c='blue')
        ax_xy.set_xlabel('x (mm)')
        ax_xy.set_ylabel('y (mm)')
        ax_xy.set_title('XY Deposition')
        ax_xy.set_aspect('equal')
        ax_xy.grid(True, alpha=0.3)

        # XZ projection (side view)
        if len(result.deposition_positions) > 0:
            ax_xz.scatter(xy[:, 2], xy[:, 0], s=1, alpha=0.5, c='red')
        ax_xz.set_xlabel('z (mm)')
        ax_xz.set_ylabel('x (mm)')
        ax_xz.set_title('XZ Deposition (side view)')
        ax_xz.grid(True, alpha=0.3)

        # Z histogram
        if len(result.deposition_positions) > 0:
            ax_hist_z.hist(xy[:, 2], bins=50, alpha=0.7, color='green')
        ax_hist_z.set_xlabel('z (mm)')
        ax_hist_z.set_ylabel('Counts')
        ax_hist_z.set_title('Z Deposition Distribution')
        ax_hist_z.grid(True, alpha=0.3)

        # Energy distribution
        if len(result.deposition_energies) > 0:
            ax_hist_e.hist(result.deposition_energies, bins=50,
                          alpha=0.7, color='purple')
        ax_hist_e.set_xlabel('Final Energy (MeV)')
        ax_hist_e.set_ylabel('Counts')
        ax_hist_e.set_title('Energy Distribution at Stop')
        ax_hist_e.grid(True, alpha=0.3)

        fig.suptitle(title or 'Gas Cell Deposition Distribution')
        plt.tight_layout()
    else:
        # Simple 3D scatter (if available)
        try:
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            if len(result.deposition_positions) > 0:
                xy = result.deposition_positions
                ax.scatter(xy[:, 0], xy[:, 1], xy[:, 2],
                          c=result.deposition_energies, cmap='hot',
                          s=2, alpha=0.6)
            ax.set_xlabel('x (mm)')
            ax.set_ylabel('y (mm)')
            ax.set_zlabel('z (mm)')
            ax.set_title(title or '3D Deposition')
        except ImportError:
            print("3D plotting not available")

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()


def plot_range_curve(stopping, E0: float, save_path: str = None):
    """Plot energy vs distance curve for an ion in the gas cell."""
    import matplotlib.pyplot as plt

    range_cm, trajectory = stopping.range_energy(E0)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Energy vs distance
    ax1.plot(trajectory[0], trajectory[1], 'b-', linewidth=1.5)
    ax1.axvline(range_cm, color='r', linestyle='--', alpha=0.5,
                label=f'Range = {range_cm:.2f} cm')
    ax1.set_xlabel('Distance (cm)')
    ax1.set_ylabel('Energy (MeV)')
    ax1.set_title(f'Energy Loss in {stopping.gas.upper()}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Stopping power vs energy
    E_vals = np.logspace(-2, np.log10(E0), 200)
    dEdx_vals = np.array([stopping(E) for E in E_vals])

    ax2.loglog(E_vals, dEdx_vals, 'g-', linewidth=1.5)
    ax2.set_xlabel('Energy (MeV)')
    ax2.set_ylabel('dE/dx (MeV/(mg/cm²))')
    ax2.set_title('Stopping Power')
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()


def plot_pressure_scan(scan_results: Dict[float, 'TransportResult'],
                       save_path: str = None):
    """Plot gas cell efficiency as a function of pressure."""
    import matplotlib.pyplot as plt

    pressures = sorted(scan_results.keys())
    efficiencies = [scan_results[P].deposition_efficiency * 100 for P in pressures]
    escapes = [scan_results[P].escape_fraction * 100 for P in pressures]

    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(pressures, efficiencies, 'bo-', linewidth=2, label='Deposited')
    ax1.plot(pressures, escapes, 'rs--', linewidth=2, label='Escaped')
    ax1.set_xlabel('Pressure (mbar)')
    ax1.set_ylabel('Fraction (%)')
    ax1.set_title('Gas Cell Performance vs Pressure')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()

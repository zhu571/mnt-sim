"""
Gas cell geometry definition.

Defines the gas cell geometry used in experimental setups
at HIAF/IMP for collecting MNT reaction products.

Reference gas cell parameters:
  - Cryogenic gas cell for SHN/MNT at SHANS/IMP
  - Volume: 50-300 mm length, 20-50 mm diameter
  - Pressure: 30-100 mbar He
  - Windows: 1-5 μm Ti/Havar/degrader foils
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional
from enum import Enum


class CellGeometry(Enum):
    """Gas cell geometry types."""
    CYLINDRICAL = "cylindrical"   # Cylindrical cell
    CONICAL = "conical"           # Cone-shaped entrance
    TAPERED = "tapered"           # Tapered to nozzle


class WindowMaterial(Enum):
    """Window/degrader foil materials."""
    TITANIUM = "Ti"
    HAVAR = "Havar"
    MYLAR = "Mylar"
    NICKEL = "Ni"


@dataclass
class Foil:
    """Window/degrader foil specification."""
    material: WindowMaterial
    thickness_um: float  # μm
    diameter_mm: float   # mm

    @property
    def thickness_mgcm2(self) -> float:
        """Convert thickness to mg/cm^2."""
        densities = {
            WindowMaterial.TITANIUM: 4.51,
            WindowMaterial.HAVAR: 8.30,
            WindowMaterial.MYLAR: 1.39,
            WindowMaterial.NICKEL: 8.91,
        }
        rho = densities.get(self.material, 4.51)
        return self.thickness_um * 1e-4 * rho * 1000  # μm -> cm -> mg/cm^2


@dataclass
class GasCell:
    """Gas cell geometry and operating conditions.

    Parameters
    ----------
    gas : str
        Stopping gas ('He', 'H2', 'N2', 'Ar')
    pressure_mbar : float
        Gas pressure [mbar]
    length_mm : float
        Cell length along beam axis [mm]
    diameter_mm : float
        Cell diameter (inner) [mm]
    temperature_K : float
        Gas temperature [K]
    entrance_window : Foil, optional
        Entrance window parameters
    exit_window : Foil, optional
        Exit/collection window parameters
    geometry : CellGeometry
        Cell shape
    taper_angle_deg : float
        Taper angle if conical [degrees]
    """
    gas: str = 'He'
    pressure_mbar: float = 50.0
    length_mm: float = 200.0
    diameter_mm: float = 40.0
    temperature_K: float = 293.0
    entrance_window: Optional[Foil] = None
    exit_window: Optional[Foil] = None
    geometry: CellGeometry = CellGeometry.CYLINDRICAL
    taper_angle_deg: float = 0.0

    @property
    def effective_length_mm(self) -> float:
        """Effective stopping length [mm]."""
        return self.length_mm

    @property
    def volume_cm3(self) -> float:
        """Cell volume [cm^3]."""
        r_cm = self.diameter_mm / 20  # mm -> cm radius
        l_cm = self.length_mm / 10    # mm -> cm
        return np.pi * r_cm**2 * l_cm

    def gas_density(self) -> float:
        """Gas density [g/cm^3] at given pressure and temperature."""
        from .stopping import get_gas_props
        # Adjust density for temperature
        props = get_gas_props(self.gas, self.pressure_mbar)
        rho = props['density']
        # Ideal gas: ρ ∝ P/T
        rho_adjusted = rho * (self.pressure_mbar / 1013.25) * (273.15 / self.temperature_K)
        return rho_adjusted

    def pressure_in_atm(self) -> float:
        """Convert pressure to atm."""
        return self.pressure_mbar / 1013.25

    def is_within_cell(self, x_mm: float, y_mm: float, z_mm: float) -> bool:
        """Check if position (x, y, z) [mm] is within the cell.

        z = 0 at entrance window, z = L at exit.
        """
        in_z = 0 <= z_mm <= self.length_mm
        if not in_z:
            return False

        r = np.sqrt(x_mm**2 + y_mm**2)
        r_max = self.diameter_mm / 2

        if self.geometry == CellGeometry.CONICAL:
            # Radius changes linearly along z
            r_max = (r_max
                     + (z_mm / self.length_mm)
                     * np.tan(np.radians(self.taper_angle_deg))
                     * self.length_mm)
        elif self.geometry == CellGeometry.TAPERED:
            if z_mm > self.length_mm * 0.8:
                r_max = r_max * (1 + (z_mm - 0.8*self.length_mm)
                                 / (0.2*self.length_mm) * 0.5)

        return r <= r_max

    def __repr__(self) -> str:
        return (f"GasCell(gas={self.gas}, P={self.pressure_mbar} mbar, "
                f"L={self.length_mm} mm, D={self.diameter_mm} mm, "
                f"T={self.temperature_K} K)")

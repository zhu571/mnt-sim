"""
Gas-cell transport and deposition modules.

Components:
- stopping: SRIM-like stopping power (Bethe-Bloch, Ziegler param.)
- monte_carlo: MC simulation of ion transport through gas
- gas_cell: Gas-cell geometry and boundary conditions
"""

from .gas_cell import CellGeometry, Foil, GasCell, WindowMaterial
from .monte_carlo import ParticleState, TransportMC, TransportResult
from .stopping import GAS_PROPERTIES, StoppingPower, get_gas_props

__all__ = [
    "CellGeometry",
    "Foil",
    "GAS_PROPERTIES",
    "GasCell",
    "ParticleState",
    "StoppingPower",
    "TransportMC",
    "TransportResult",
    "WindowMaterial",
    "get_gas_props",
]

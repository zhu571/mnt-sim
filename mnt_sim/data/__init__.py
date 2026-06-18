"""
Nuclear data: masses, binding energies, and isotope properties.

Sources:
  - AME2020 atomic mass evaluation
  - NNDC Evaluated Nuclear Structure Data File (ENSDF)
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

__all__ = [
    "ELEMENT_SYMBOLS",
    "ELEMENTS",
    "MASS_EXCESS",
    "TYPICAL_MASS_NUMBERS",
    "Nuclide",
    "element_to_z",
    "get_nuclide",
    "mass_excess",
    "target_projectile_pairs",
    "typical_mass_number",
]


@dataclass
class Nuclide:
    """Nuclear species properties."""
    Z: int      # Proton number
    N: int      # Neutron number
    A: int      # Mass number
    symbol: str  # Element symbol
    mass: float  # Atomic mass [MeV/c^2]
    abundance: float  # Natural abundance (0-1)


# Element symbols keyed by atomic number.
ELEMENT_SYMBOLS = {
    0: 'n', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N',
    8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si',
    15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc',
    22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni',
    29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br',
    36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
    43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In',
    50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba',
    57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu',
    64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir',
    78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po',
    85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa',
    92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf',
    99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf',
    105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
    111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv',
    117: 'Ts', 118: 'Og',
}

# Element symbols keyed by normalized element symbol.
ELEMENTS = {symbol: Z for Z, symbol in ELEMENT_SYMBOLS.items()}

# Typical isotope mass numbers used by the simplified reaction models.
TYPICAL_MASS_NUMBERS = {
    'Ar': 40,
    'Ca': 40,
    'Ge': 74,
    'Kr': 86,
    'Ni': 58,
    'Os': 200,
    'Pt': 198,
    'Sn': 124,
    'Xe': 136,
    'Ba': 138,
    'Pb': 208,
    'Ra': 226,
    'Th': 232,
    'U': 238,
    'Cm': 248,
}


# Common stable isotopes with masses (mass excess in MeV)
# Format: (Z, N) -> mass_excess [MeV]
MASS_EXCESS = {
    # Light
    (1, 0): 7.2890,     # p (H-1)
    (1, 0): 7.2890,     # H-1
    (1, 1): 13.13572,   # H-2
    (1, 2): 14.94980,   # H-3
    (2, 2): 2.4249,     # He-4
    # Ca isotopes
    (20, 20): -34.846,  # Ca-40
    (20, 22): -38.547,  # Ca-42
    (20, 24): -41.468,  # Ca-44
    (20, 28): -44.215,  # Ca-48
    # Ni isotopes
    (28, 30): -56.373,  # Ni-58
    (28, 32): -61.155,  # Ni-60
    (28, 34): -64.471,  # Ni-62
    (28, 38): -67.097,  # Ni-64
    # Xe isotopes
    (54, 74): -71.230,  # Xe-128
    (54, 76): -72.111,  # Xe-130
    (54, 78): -73.420,  # Xe-132
    (54, 80): -74.210,  # Xe-134
    (54, 82): -74.200,  # Xe-136
    # Pb isotopes
    (82, 122): -71.100, # Pb-204
    (82, 124): -72.350, # Pb-206
    (82, 126): -73.175, # Pb-208
    # U isotopes
    (92, 142): -57.343, # U-234
    (92, 143): -57.800, # U-235
    (92, 144): -58.103, # U-236
    (92, 146): -58.164, # U-238
}


def _normalize_symbol(symbol: str) -> str:
    """Normalize an element symbol for table lookup."""
    if not isinstance(symbol, str) or not symbol.strip():
        raise ValueError("element symbol must be a non-empty string")
    stripped = symbol.strip()
    return 'n' if stripped == 'n' else stripped.capitalize()


def element_to_z(symbol: str, default: Optional[int] = None) -> int:
    """Return atomic number for an element symbol.

    Parameters
    ----------
    symbol : str
        Element symbol, case-insensitive.
    default : int, optional
        Value returned when the symbol is unknown. If omitted, a
        ValueError is raised for unknown symbols.
    """
    normalized = _normalize_symbol(symbol)
    if normalized in ELEMENTS:
        return ELEMENTS[normalized]
    if default is not None:
        return default
    raise ValueError(f"unknown element symbol: {symbol!r}")


def typical_mass_number(symbol: str, default: Optional[int] = None) -> int:
    """Return the typical mass number used by simplified models."""
    normalized = _normalize_symbol(symbol)
    if normalized in TYPICAL_MASS_NUMBERS:
        return TYPICAL_MASS_NUMBERS[normalized]
    if default is not None:
        return default
    raise ValueError(f"no typical mass number configured for {symbol!r}")


def get_nuclide(Z: int, A: int) -> Nuclide:
    """Create a Nuclide from Z and A."""
    if Z < 0:
        raise ValueError("Z must be non-negative")
    if A < Z:
        raise ValueError("A must be greater than or equal to Z")
    N = A - Z
    symbol = ELEMENT_SYMBOLS.get(Z, '?')
    mass_excess = MASS_EXCESS.get((Z, N), 0.0)
    mass = A * 931.494 + mass_excess
    return Nuclide(Z=Z, N=N, A=A, symbol=symbol, mass=mass,
                   abundance=0.0)


def mass_excess(Z: int, N: int) -> float:
    """Get mass excess [MeV] for a nuclide, interpolating when missing."""
    key = (Z, N)
    if key in MASS_EXCESS:
        return MASS_EXCESS[key]

    # Simple extrapolation using semi-empirical mass formula
    A = Z + N
    a_v, a_s, a_c, a_a, a_p = 15.75, 17.80, 0.711, 23.70, 11.18
    B = (a_v * A - a_s * A**(2/3) - a_c * Z*(Z-1)/A**(1/3)
         - a_a * (A - 2*Z)**2 / A)
    if A % 2 == 0:
        B += a_p * (-1)**Z * (-1)**(A-Z) / A**(1/2)
    return Z * 7.2890 + N * 8.0713 - B


def target_projectile_pairs() -> List[Tuple[str, str, float]]:
    """Commonly used target-projectile combinations for MNT."""
    return [
        ("Xe", "Pb", 8.0),     # Typical IMP/HIAF MNT
        ("Xe", "U", 8.0),
        ("Ca", "U", 7.0),
        ("Ca", "Pb", 7.0),
        ("Ni", "Pb", 7.5),
        ("Ni", "U", 7.5),
        ("Kr", "Pb", 8.0),
        ("Ar", "Pb", 6.5),
    ]

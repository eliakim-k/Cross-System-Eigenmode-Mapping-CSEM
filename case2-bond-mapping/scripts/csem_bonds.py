"""
Shared helpers for the bond-mapping case: parse DDEC6 bond-order files, map the
gas-phase molecule onto each adsorbed system, and compare bond orders.

All three step scripts import from here so the parsing/mapping logic lives in
one place.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

# Covalent radii (Å) for neighbor detection; bonded if r < 1.3*(r_i + r_j).
COVALENT_RADII = {"H": 0.31, "O": 0.66, "Mo": 1.54, "Cu": 1.32, "Ni": 1.24}
COV_TOL = 1.3

# Activation thresholds on the DDEC6 bond order (BO), gas -> adsorbed.
THRESHOLDS = {
    "new_bo_gas_max": 0.05,        # gas BO below this + real ads BO -> New
    "dissociated_bo_ads_max": 0.10,  # former bond with ads BO below this -> Dissociated
    "activated_delta_bo": -0.05,   # delta_bo <= this -> Activated
    "strengthened_delta_bo": 0.05,  # delta_bo >= this -> Strengthened
}


@dataclass
class DDEC6:
    """A parsed DDEC6 bond-order file (0-based atom indexing)."""
    elements: List[str]
    coords: np.ndarray                 # (n_atoms, 3) Å
    bonds: List[Tuple[int, int, float, Tuple[int, int, int]]]  # (i, j, BO, translation)

    def indices_of(self, element: str) -> List[int]:
        return [i for i, e in enumerate(self.elements) if e == element]


def read_bond_orders(path: Path) -> DDEC6:
    """Parse ``DDEC6_even_tempered_bond_orders.xyz``."""
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    n = int(lines[0].split()[0])
    elements, coords = [], []
    for i in range(n):
        p = lines[2 + i].split()
        elements.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])
    bonds, current = [], None
    for ln in lines:
        if "Printing BOs for ATOM #" in ln:
            m = re.search(r"ATOM #\s+(\d+)", ln)
            current = int(m.group(1)) - 1 if m else current   # to 0-based
        elif "Bonded to the" in ln and current is not None:
            t = re.search(r"\(\s*([\-\d]+),\s*([\-\d]+),\s*([\-\d]+)\)", ln)
            a2 = re.search(r"atom number\s+(\d+)", ln)
            bo = re.search(r"bond order\s*=\s*([\d\.\-]+)", ln)
            if t and a2 and bo:
                bonds.append((current, int(a2.group(1)) - 1, float(bo.group(1)),
                              tuple(int(x) for x in t.groups())))
    return DDEC6(elements, np.asarray(coords, dtype=float), bonds)


def _neighbors(d: DDEC6, i: int) -> List[int]:
    r1 = COVALENT_RADII.get(d.elements[i], 1.5)
    out = []
    for j, c in enumerate(d.coords):
        if j == i:
            continue
        r2 = COVALENT_RADII.get(d.elements[j], 1.5)
        if float(np.linalg.norm(d.coords[i] - c)) < COV_TOL * (r1 + r2):
            out.append(j)
    return out


def map_gas_to_adsorbed(gas: DDEC6, ads: DDEC6, h_penalty: float = 2.0) -> Dict[int, int]:
    """Hungarian atom map {gas_idx -> ads_idx} (0-based).

    Heavy atoms are matched on distance after anchoring on a unique heavy
    element; hydrogens are matched on distance plus a penalty when their bonded
    heavy neighbor does not correspond to the gas hydrogen's heavy neighbor.
    """
    common = set(gas.elements) & set(ads.elements)
    heavies = [e for e in common if e != "H"]
    unique = [e for e in heavies if gas.elements.count(e) == 1]
    anchor = unique[0] if unique else None

    trans = np.zeros(3)
    if anchor:
        gi = gas.indices_of(anchor)[0]
        ai = max(ads.indices_of(anchor), key=lambda k: ads.coords[k][2])
        trans = ads.coords[ai] - gas.coords[gi]
    gas_shift = gas.coords + trans

    mapping: Dict[int, int] = {}
    for e in heavies:
        gl, al = gas.indices_of(e), ads.indices_of(e)
        if not gl or not al:
            continue
        cost = np.array([[np.linalg.norm(gas_shift[g] - ads.coords[a]) for a in al] for g in gl])
        r, c = linear_sum_assignment(cost)
        for ri, ci in zip(r, c):
            mapping[gl[ri]] = al[ci]

    gh, ah = gas.indices_of("H"), ads.indices_of("H")
    if gh and ah:
        cost = np.zeros((len(gh), len(ah)))
        for i, g in enumerate(gh):
            ng = [mapping[k] for k in _neighbors(gas, g) if k in mapping]
            for j, a in enumerate(ah):
                dist = float(np.linalg.norm(gas_shift[g] - ads.coords[a]))
                na = _neighbors(ads, a)
                pen = 0.0 if (any(m in na for m in ng) or not ng) else h_penalty
                cost[i, j] = dist + pen
        r, c = linear_sum_assignment(cost)
        for ri, ci in zip(r, c):
            mapping[gh[ri]] = ah[ci]
    return mapping


def _index_intra_bonds(d: DDEC6, min_bo: float) -> Dict[Tuple[int, int], float]:
    """{(i, j) sorted 0-based: BO} for intra-cell bonds above min_bo."""
    out: Dict[Tuple[int, int], float] = {}
    for i, j, bo, tr in d.bonds:
        if tr != (0, 0, 0) or bo < min_bo:
            continue
        k = (i, j) if i <= j else (j, i)
        if k not in out or bo > out[k]:
            out[k] = bo
    return out


def classify(bo_gas: float, bo_ads: float) -> str:
    """Activation label for one bond from its gas and adsorbed bond orders."""
    th = THRESHOLDS
    gas_missing = np.isnan(bo_gas) or bo_gas <= th["new_bo_gas_max"]
    ads_missing = np.isnan(bo_ads)
    if gas_missing and not ads_missing and bo_ads >= th["new_bo_gas_max"]:
        return "New"
    if (not np.isnan(bo_gas)) and bo_gas > th["new_bo_gas_max"] and (
        ads_missing or bo_ads <= th["dissociated_bo_ads_max"]
    ):
        return "Dissociated"
    if np.isnan(bo_gas) or np.isnan(bo_ads):
        return "Unchanged"
    delta = bo_ads - bo_gas
    if delta <= th["activated_delta_bo"]:
        return "Activated"
    if delta >= th["strengthened_delta_bo"]:
        return "Strengthened"
    return "Unchanged"


def distance(d: DDEC6, i: int, j: int) -> float:
    return float(np.linalg.norm(d.coords[i] - d.coords[j]))

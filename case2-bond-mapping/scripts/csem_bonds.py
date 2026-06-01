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
    "new_bo_gas_max": 0.05,
    "dissociated_bo_ads_max": 0.10,
    "activated_delta_bo": -0.05,
    "strengthened_delta_bo": 0.05,
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
            current = int(m.group(1)) - 1 if m else current
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


def _kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Proper rotation R with (P - mean P) @ R.T best matching (Q - mean Q)."""
    H = (P - P.mean(0)).T @ (Q - Q.mean(0))
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def map_gas_to_adsorbed(gas: DDEC6, ads: DDEC6, h_penalty: float = 2.0,
                        max_iter: int = 30) -> Dict[int, int]:
    """Rigid (rotation + translation) atom map {gas_idx -> ads_idx}, 0-based.

    The gas molecule is superposed onto the adsorbed molecule by an iterative
    Kabsch (rotation) + per-element Hungarian loop, minimizing RMSD:

    1. anchor on the unique heavy element (Mo) for an initial translation;
    2. assign atoms per element by minimum distance (heavy atoms by distance;
       hydrogens by distance plus a chemical-ancestry penalty that discourages,
       e.g., an O-H hydrogen from mapping onto a metal-H hydride);
    3. refit the rotation from the current correspondence and repeat.

    Only adsorbed atoms whose element also occurs in the gas phase (i.e. the
    adsorbate, not the metal surface) are candidates. The best-RMSD mapping
    encountered is returned.
    """
    common = set(gas.elements) & set(ads.elements)
    heavies = [e for e in common if e != "H"]
    unique = [e for e in heavies if gas.elements.count(e) == 1]
    anchor = unique[0] if unique else None
    ads_mol = [i for i, e in enumerate(ads.elements) if e in common]

    # Initial transform: translate so the anchor heavy atoms coincide.
    R = np.eye(3)
    if anchor:
        g0 = gas.indices_of(anchor)[0]
        anchors = [a for a in ads.indices_of(anchor) if a in ads_mol] or ads.indices_of(anchor)
        a0 = max(anchors, key=lambda k: ads.coords[k][2])  # topmost = adsorbate
        gas_com, ads_com = gas.coords[g0], ads.coords[a0]
    else:
        gas_com = gas.coords.mean(0)
        ads_com = ads.coords[ads_mol].mean(0)

    def assign(R, gas_com, ads_com) -> Dict[int, int]:
        T = (gas.coords - gas_com) @ R.T + ads_com   # gas atoms in the ads frame
        mapping: Dict[int, int] = {}
        for e in heavies:
            gl = gas.indices_of(e)
            al = [a for a in ads.indices_of(e) if a in ads_mol]
            cost = np.array([[np.linalg.norm(T[g] - ads.coords[a]) for a in al] for g in gl])
            r, c = linear_sum_assignment(cost)
            for ri, ci in zip(r, c):
                mapping[gl[ri]] = al[ci]
        gh = gas.indices_of("H")
        ah = [a for a in ads.indices_of("H") if a in ads_mol]
        if gh and ah:
            cost = np.zeros((len(gh), len(ah)))
            for i, g in enumerate(gh):
                ng = [mapping[k] for k in _neighbors(gas, g) if k in mapping]
                for j, a in enumerate(ah):
                    dist = float(np.linalg.norm(T[g] - ads.coords[a]))
                    na = _neighbors(ads, a)
                    pen = 0.0 if (any(m in na for m in ng) or not ng) else h_penalty
                    cost[i, j] = dist + pen
            r, c = linear_sum_assignment(cost)
            for ri, ci in zip(r, c):
                mapping[gh[ri]] = ah[ci]
        return mapping

    mapping = assign(R, gas_com, ads_com)
    best = (np.inf, mapping)
    for _ in range(max_iter):
        gl = sorted(mapping)
        P = gas.coords[gl]
        Q = ads.coords[[mapping[g] for g in gl]]
        gas_com, ads_com = P.mean(0), Q.mean(0)
        R = _kabsch_rotation(P, Q)
        rmsd = float(np.sqrt(np.mean(((P - gas_com) @ R.T + ads_com - Q) ** 2)))
        new = assign(R, gas_com, ads_com)
        if rmsd < best[0]:
            best = (rmsd, new)
        if new == mapping:
            break
        mapping = new
    return best[1]


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

#!/usr/bin/env python3
"""
Step 2 - Compare bond orders across gas, Cu(111), and Ni(111).

Using the atom maps from step 1, every intramolecular bond of the gas-phase
molecule is located in each adsorbed system and its DDEC6 bond order recorded.
The change relative to the gas phase, delta_BO = BO_ads - BO_gas, and an
activation label (Activated / Strengthened / Dissociated / Unchanged) follow
directly. New adsorbate-surface bonds (one mapped atom + one metal atom) are
reported separately per surface.

Inputs : data/.../bond_orders.xyz, results/atom_map_<S>.json
Outputs: results/bond_match_<S>.csv          (per-surface, incl. surface bonds)
         results/intramolecular_comparison.csv (gas bonds x gas/Cu/Ni)
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from csem_bonds import (read_bond_orders, _index_intra_bonds, classify,
                        distance, COVALENT_RADII)

ROOT = Path(__file__).resolve().parent.parent
SURFACES = ["Cu", "Ni"]
MIN_BO = 0.05


def per_surface(gas, ads, mapping, S):
    """Full per-surface bond match (intramolecular mapped bonds + surface bonds)."""
    g2a = mapping
    mapped_ads = set(g2a.values())
    gi = _index_intra_bonds(gas, MIN_BO)
    ai = _index_intra_bonds(ads, MIN_BO)
    rows, seen = [], set()
    for (g1, g2), bo_g in gi.items():
        if g1 not in g2a or g2 not in g2a:
            continue
        a1, a2 = g2a[g1], g2a[g2]
        key = (a1, a2) if a1 <= a2 else (a2, a1)
        seen.add(key)
        bo_a = ai.get(key, np.nan)
        rows.append(dict(elem1=gas.elements[g1], elem2=gas.elements[g2],
                         gas1=g1, gas2=g2, ads1=a1, ads2=a2,
                         bo_gas=bo_g, bo_ads=bo_a,
                         delta_bo=(bo_a - bo_g) if not np.isnan(bo_a) else np.nan,
                         length_gas=distance(gas, g1, g2),
                         length_ads=distance(ads, a1, a2) if not np.isnan(bo_a) else np.nan,
                         kind="intramolecular", label=classify(bo_g, bo_a)))
    for (a1, a2), bo_a in ai.items():
        in1, in2 = a1 in mapped_ads, a2 in mapped_ads
        key = (a1, a2)
        if in1 and in2 and key not in seen:               # new intramolecular bond
            rows.append(dict(elem1=ads.elements[a1], elem2=ads.elements[a2],
                             gas1=-1, gas2=-1, ads1=a1, ads2=a2,
                             bo_gas=np.nan, bo_ads=bo_a, delta_bo=np.nan,
                             length_gas=np.nan, length_ads=distance(ads, a1, a2),
                             kind="intramolecular", label=classify(np.nan, bo_a)))
        elif in1 ^ in2:                                    # adsorbate-surface bond
            mol, surf = (a1, a2) if in1 else (a2, a1)
            rows.append(dict(elem1=ads.elements[mol], elem2=ads.elements[surf],
                             gas1=-1, gas2=-1, ads1=mol, ads2=surf,
                             bo_gas=np.nan, bo_ads=bo_a, delta_bo=np.nan,
                             length_gas=np.nan, length_ads=distance(ads, mol, surf),
                             kind="adsorbate-surface", label="New"))
    return pd.DataFrame(rows)


def main():
    gas = read_bond_orders(ROOT / "data" / "gas" / "bond_orders.xyz")
    gi = _index_intra_bonds(gas, MIN_BO)
    comp = {(g1, g2): dict(bond=f"{gas.elements[g1]}{g1}-{gas.elements[g2]}{g2}",
                           elem_pair=f"{gas.elements[g1]}-{gas.elements[g2]}",
                           bo_gas=bo) for (g1, g2), bo in gi.items()}

    for S in SURFACES:
        ads = read_bond_orders(ROOT / "data" / S / "bond_orders.xyz")
        mapping = {int(k): int(v) for k, v in
                   json.loads((ROOT / "results" / f"atom_map_{S}.json").read_text())
                   ["mapping_gas_to_ads"].items()}
        df = per_surface(gas, ads, mapping, S)
        df.to_csv(ROOT / "results" / f"bond_match_{S}.csv", index=False)
        summary = df["label"].value_counts().to_dict()
        print(f"[{S}] {summary}")
        ai = _index_intra_bonds(ads, MIN_BO)
        for (g1, g2), rec in comp.items():
            a1, a2 = mapping[g1], mapping[g2]
            key = (a1, a2) if a1 <= a2 else (a2, a1)
            bo_a = ai.get(key, np.nan)
            rec[f"bo_{S}"] = bo_a
            rec[f"dBO_{S}"] = (bo_a - rec["bo_gas"]) if not np.isnan(bo_a) else np.nan
            rec[f"label_{S}"] = classify(rec["bo_gas"], bo_a)

    out = pd.DataFrame(list(comp.values())).sort_values("bo_gas", ascending=False)
    out.to_csv(ROOT / "results" / "intramolecular_comparison.csv", index=False)
    print("Wrote results/bond_match_Cu.csv, bond_match_Ni.csv, intramolecular_comparison.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Step 1 - Map the gas-phase molecule onto each adsorbed system.

To compare a bond across systems we must know which adsorbed atom is which
gas-phase atom. This step solves that atom correspondence for gas -> Cu and gas -> Ni with an
iterative Kabsch rotation + per-element Hungarian assignment that minimizes RMSD
(anchor on Mo; heavy atoms by distance; hydrogens by distance plus a
chemical-ancestry penalty), reading the geometries straight from the DDEC6
bond-order files.

Input : data/gas/bond_orders.xyz, data/<S>/bond_orders.xyz   (S = Cu, Ni)
Output: results/atom_map_<S>.json   ({gas_index: ads_index}, 0-based)
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from csem_bonds import read_bond_orders, map_gas_to_adsorbed

ROOT = Path(__file__).resolve().parent.parent
SURFACES = ["Cu", "Ni"]


def main():
    (ROOT / "results").mkdir(exist_ok=True)
    gas = read_bond_orders(ROOT / "data" / "gas" / "bond_orders.xyz")
    for S in SURFACES:
        ads = read_bond_orders(ROOT / "data" / S / "bond_orders.xyz")
        mapping = map_gas_to_adsorbed(gas, ads)
        if len(mapping) != len(gas.elements):
            raise SystemExit(f"[{S}] incomplete map: {len(mapping)}/{len(gas.elements)}")
        pairs = {f"{g} ({gas.elements[g]})": f"{a} ({ads.elements[a]})"
                 for g, a in sorted(mapping.items())}
        out = ROOT / "results" / f"atom_map_{S}.json"
        out.write_text(json.dumps({
            "surface": S,
            "mapping_gas_to_ads": {int(g): int(a) for g, a in mapping.items()},
            "readable": pairs,
        }, indent=2))
        print(f"[{S}] mapped {len(mapping)} molecule atoms -> {out.name}")


if __name__ == "__main__":
    main()

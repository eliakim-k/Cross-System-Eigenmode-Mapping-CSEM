# Case 2 — Bond mapping (gas phase vs Cu(111), Ni(111))

Track how each intramolecular bond of an adsorbed molecule changes relative to
the gas phase, using DDEC6 bond orders.

## The problem

Adsorption strengthens some bonds and weakens (activates) others. To see this we
compare each bond's order on the surface with its gas-phase value — but the same
bond carries different atom indices in each calculation, and new molecule–metal
bonds appear on the surface. The bonds must be put into correspondence first.

## The approach

First map the gas-phase molecule onto the adsorbed system (which adsorbed atom is
which gas atom), then compare bond orders through that map.

**Atom map (gas → adsorbed).** The gas molecule is superposed onto the adsorbed
molecule by an iterative **Kabsch rotation + per-element Hungarian** assignment
that minimizes RMSD:

1. anchor on the unique heavy element (Mo) for an initial translation;
2. assign atoms per element by nearest distance — heavy atoms directly,
   hydrogens with a chemical-ancestry penalty (a gas O–H hydrogen is discouraged
   from mapping onto an adsorbed Mo–H hydride unless the geometry demands it);
3. refit the rotation from the current correspondence and repeat until it
   converges, keeping the best-RMSD mapping.

**Bond comparison.** Each gas-phase bond is located on the surface through the
map and its DDEC6 bond order recorded, giving

```
ΔBO = BO_ads − BO_gas
```

and an activation label from `ΔBO` (with floors for appearing/disappearing
bonds):

| label | condition |
|-------|-----------|
| Activated | `ΔBO ≤ −0.05` (bond weakened) |
| Strengthened | `ΔBO ≥ +0.05` |
| Dissociated | present in gas, essentially absent on the surface |
| Unchanged | `|ΔBO| < 0.05` |
| New | absent in gas, present on the surface (e.g. molecule–metal) |

New adsorbate–surface bonds (one mapped molecule atom + one metal atom) are
reported separately per surface.

## Application: H3MoOH on Cu(111) and Ni(111)

Bond orders are from a DDEC6 (Chargemol) analysis of periodic plane-wave DFT
(**VASP**) electron densities, for the gas-phase molecule and the two adsorbed
systems. The molecule is H3MoOH (Mo, O, 4 H).

<p align="center">
  <img src="results/bond_order_comparison.png" width="620" alt="Intramolecular bond orders: gas vs Cu vs Ni">
</p>

Each intramolecular bond's DDEC6 bond order is shown for the gas phase and the
two surfaces. `results/intramolecular_comparison.csv` tabulates `bo_gas`,
`bo_Cu`, `bo_Ni`, `ΔBO`, and the activation label per surface;
`results/bond_match_<S>.csv` gives the full per-surface match including the new
adsorbate–surface bonds.

## Data

`data/{gas,Cu,Ni}/bond_orders.xyz` — DDEC6 `even_tempered` bond-order files
(`Chargemol`). Each file is self-contained: an atom block (element, x, y, z, sum
of bond orders) followed by per-atom bonded-pair listings with bond orders, so it
provides both the geometry (for the atom map) and the bond orders (for the
comparison).

## Usage

```bash
pip install -r requirements.txt
python scripts/step1_map_atoms.py        # -> results/atom_map_Cu.json, atom_map_Ni.json
python scripts/step2_match_bonds.py      # -> results/bond_match_<S>.csv, intramolecular_comparison.csv
python scripts/step3_plot_bond_changes.py  # -> results/bond_order_comparison.png
```

Surfaces, the bond-order floor, and the activation thresholds are set in the
configuration blocks at the top of `scripts/csem_bonds.py` and the step scripts.

## License

MIT (see the repository root `LICENSE`).

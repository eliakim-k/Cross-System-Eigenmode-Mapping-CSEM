# Data

Vibrational data for **H3MoOH** adsorbed on Cu(111) and Ni(111), computed in this
work with periodic plane-wave DFT (**VASP**). The vibrational analysis is a
finite-displacement run (`IBRION = 5`); frequencies and displacement eigenvectors
are read from the resulting `OUTCAR`.

The raw `OUTCAR` files are large and are *not* shipped. Instead, for each surface
the eigenvectors are projected onto the six adsorbate atoms (Mo, O, 4 H) and
L2-normalized to unit length, then written to CSV. The mapping pipeline
(`scripts/step1_match_modes.py`) re-normalizes defensively, so the result is
unchanged by the scaling of the input vectors.

Per surface, `<surface>` is `Cu` or `Ni`.

## `POSCAR_<surface>` — slab geometry (VASP POSCAR)

Standard VASP POSCAR for the adsorbed slab. The adsorbate atoms are listed
**first** (element order `H O Mo`), so the six leading rows are the adsorbate.
Used only to align Ni onto Cu (per-element Hungarian permutation + Kabsch
rotation); the slab substrate atoms are not read by the pipeline.

## `eigenvectors_<surface>.csv` — adsorbate-projected modes (analysis input)

One row per vibrational mode (10 modes per surface).

| Column                      | Meaning                                                              |
|-----------------------------|----------------------------------------------------------------------|
| `mode_index`                | mode number (0-based), ordered by descending frequency               |
| `freq_cm-1`                 | mode frequency (cm⁻¹)                                                 |
| `dx_a0,dy_a0,dz_a0` … `_a5` | unit-normalized displacement of each adsorbate atom (Cartesian, Å)   |

The six atom blocks `a0…a5` follow the POSCAR adsorbate order (`H H H H O Mo`).
Each row is L2-normalized over all 18 displacement components, so
`|v| = 1` (verified). The overlap `S_ij = |v_i(Cu)·v_j(Ni)|` is therefore a true
cosine in [0, 1].

## `freqs_<surface>.csv` — frequency list (convenience)

| Column       | Meaning                          |
|--------------|----------------------------------|
| `mode_index` | mode number (matches the CSV)    |
| `freq_cm-1`  | mode frequency (cm⁻¹)            |

A standalone copy of the `mode_index, freq_cm-1` columns for quick reference; the
mapping reads its frequencies from `eigenvectors_<surface>.csv`.

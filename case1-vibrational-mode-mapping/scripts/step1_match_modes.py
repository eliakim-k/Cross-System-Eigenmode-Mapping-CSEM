#!/usr/bin/env python3
"""
Step 1 - Match vibrational eigenmodes between the Cu(111)- and Ni(111)-adsorbed
states of the same molecule.

A vibrational mode is a unit displacement vector over the adsorbate atoms. Two
modes "are the same motion" to the extent their displacement vectors point the
same way, i.e. their absolute dot product (cosine similarity) is close to 1.

To compare modes the two adsorbates must live in the same atom order and the
same Cartesian frame, so the moved system (Ni) is first rigidly aligned onto the
reference (Cu): a per-species Hungarian assignment fixes the atom permutation and
Kabsch fixes the rotation. That rotation/permutation is then applied to the moved
system's mode vectors, and the overlap matrix S[i, j] = |v_i(ref) . v_j(mov)| is
formed. A greedy one-to-one pass pairs each reference mode with a distinct moved
mode in descending overlap.

The mode vectors are the adsorbate-projected displacement eigenvectors read from
the VASP OUTCAR of the finite-displacement run; they are L2-normalized again here
so the overlap is a cosine regardless of how the inputs are scaled.

Inputs  : data/POSCAR_<ref>, data/POSCAR_<mov>            (geometry; adsorbate first)
          data/eigenvectors_<ref>.csv, _<mov>.csv         (per-mode displacements)
Outputs : results/matches.csv          (one row per matched mode pair)
          results/overlap_matrix.csv    (full |ref . mov| matrix)
          results/mode_alignment.json   (rotation, permutation, RMSD)
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent.parent
REF, MOV = "Cu", "Ni"                      # reference / moved system labels
ADSORBATE_ELEMENTS = {"H", "O", "Mo"}      # everything else is substrate


# --------------------------------------------------------------------------- #
# Minimal POSCAR reader (Cartesian Å, adsorbate = atoms whose element is in
# ADSORBATE_ELEMENTS, which by convention are listed first)
# --------------------------------------------------------------------------- #
def read_adsorbate(poscar: Path):
    lines = [l.strip() for l in poscar.read_text().splitlines() if l.strip()]
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    elements = [s for s, c in zip(species, counts) for _ in range(c)]
    idx = 7
    if lines[idx][0].lower() == "s":         # selective dynamics
        idx += 1
    direct = lines[idx][0].lower() in ("d", "l")
    idx += 1
    n_ads = sum(c for s, c in zip(species, counts) if s in ADSORBATE_ELEMENTS)
    coords = np.array([[float(x) for x in lines[idx + i].split()[:3]] for i in range(n_ads)])
    coords = coords @ lat if direct else coords * scale
    return elements[:n_ads], coords


def kabsch_rotation(P, Q):
    """Proper rotation taking centered P onto centered Q."""
    H = (P - P.mean(0)).T @ (Q - Q.mean(0))
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def align_permutation(P_ref, P_mov, symbols, max_iter=20):
    """Iterative per-species Hungarian + Kabsch. Returns best (R, perm)."""
    n = len(symbols)
    perm = list(range(n))
    best = (np.inf, np.eye(3), perm[:])
    cur = P_mov.copy()
    for _ in range(max_iter):
        R = kabsch_rotation(cur, P_ref)
        diff = (P_ref - P_ref.mean(0)) - (cur - cur.mean(0)) @ R.T
        rmsd = float(np.sqrt((diff ** 2).mean()))
        if rmsd < best[0]:
            best = (rmsd, R, perm[:])
        if rmsd < 1e-3:
            break
        P_rot = (P_mov - P_mov.mean(0)) @ R.T + P_ref.mean(0)
        new_perm = [0] * n
        for s in sorted(set(symbols)):
            ix = [i for i, x in enumerate(symbols) if x == s]
            cost = np.sum((P_ref[ix][:, None, :] - P_rot[ix][None, :, :]) ** 2, axis=-1)
            r, c = linear_sum_assignment(cost)
            for a, b in zip(r, c):
                new_perm[ix[a]] = ix[b]
        if new_perm == perm:
            break
        perm, cur = new_perm, P_mov[new_perm]
    return best[1], best[2], best[0]


def load_modes(csv: Path):
    df = pd.read_csv(csv)
    disp = [c for c in df.columns if c.startswith(("dx", "dy", "dz"))]
    vecs = df[disp].to_numpy(float)
    return df["freq_cm-1"].to_numpy(float), vecs


def greedy_one_to_one(S):
    nA, nB = S.shape
    order = np.argsort(-S.ravel(), kind="stable")
    used_i, used_j, pairs = np.zeros(nA, bool), np.zeros(nB, bool), []
    for idx in order:
        i, j = int(idx // nB), int(idx % nB)
        if used_i[i] or used_j[j]:
            continue
        used_i[i] = used_j[j] = True
        pairs.append((i, j, float(S[i, j])))
        if len(pairs) >= min(nA, nB):
            break
    pairs.sort()
    return pairs


def main():
    d = ROOT / "data"
    sym_ref, P_ref = read_adsorbate(d / f"POSCAR_{REF}")
    sym_mov, P_mov = read_adsorbate(d / f"POSCAR_{MOV}")
    if sorted(sym_ref) != sorted(sym_mov):
        raise SystemExit(f"Adsorbate composition differs: {sym_ref} vs {sym_mov}")
    n_ads = len(sym_ref)

    R, perm, rmsd = align_permutation(P_ref, P_mov, sym_ref)

    f_ref, V_ref = load_modes(d / f"eigenvectors_{REF}.csv")
    f_mov, V_mov = load_modes(d / f"eigenvectors_{MOV}.csv")
    nmodes = V_mov.shape[0]

    # Apply alignment (atom permutation + rotation) to the moved mode vectors.
    Vm = V_mov.reshape(nmodes, n_ads, 3)[:, perm, :] @ R.T
    Vm = Vm.reshape(nmodes, n_ads * 3)
    V_ref = V_ref / np.linalg.norm(V_ref, axis=1, keepdims=True)
    Vm = Vm / np.linalg.norm(Vm, axis=1, keepdims=True)

    S = np.abs(V_ref @ Vm.T)
    pairs = greedy_one_to_one(S)

    (ROOT / "results").mkdir(exist_ok=True)
    rows = [{
        f"{REF}_mode": i, f"freq_{REF}_cm-1": float(f_ref[i]),
        f"{MOV}_mode": j, f"freq_{MOV}_cm-1": float(f_mov[j]),
        "delta_cm-1": float(f_mov[j] - f_ref[i]), "overlap": s,
    } for i, j, s in pairs]
    pd.DataFrame(rows).to_csv(ROOT / "results" / "matches.csv", index=False)

    pd.DataFrame(
        S, index=[f"{REF}_{i}" for i in range(S.shape[0])],
        columns=[f"{MOV}_{j}" for j in range(S.shape[1])],
    ).to_csv(ROOT / "results" / "overlap_matrix.csv")

    (ROOT / "results" / "mode_alignment.json").write_text(json.dumps({
        "reference": REF, "moved": MOV, "n_adsorbate_atoms": n_ads,
        "rotation": R.tolist(), "permutation": list(map(int, perm)),
        "rmsd_A": rmsd, "symbols": sym_ref,
    }, indent=2))

    print(f"Aligned {MOV} onto {REF}: RMSD = {rmsd:.3f} Å, perm = {perm}")
    print(f"Matched {len(pairs)} modes; mean overlap = {np.mean([s for *_ , s in pairs]):.3f}")
    print(f"Wrote results/matches.csv, results/overlap_matrix.csv, results/mode_alignment.json")


if __name__ == "__main__":
    main()

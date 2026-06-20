#!/usr/bin/env python
"""
00_prepare_receptor.py — Prepare the CDK2 receptor (PDB 2R3I) for docking and MD.

Pipeline:
  1. Download the crystal structure (if absent).
  2. Extract the co-crystal inhibitor (SCF) -> reference for the docking box & redock control.
  3. Clean the protein with PDBFixer:
       - revert the oxidized cysteine CSD177 -> CYS (crystallization artifact);
       - rebuild the disordered beta3/alphaC loop (~res 46-52, 'LDTETEG'), missing in the
         crystal, so the MD topology has no chain break (loop is distal to the ATP site);
       - drop heteroatoms/waters; add missing heavy atoms and hydrogens at pH 7.0.
  4. Define the docking box from the inhibitor's heavy-atom envelope.
  5. Write a Vina-ready rigid receptor PDBQT (OpenBabel).

Outputs (under <workdir>):
  receptor/cdk2_prepared.pdb     cleaned, protonated protein
  receptor/ref_ligand_xtal.pdb   crystal pose of SCF (heavy atoms; box + redock reference)
  receptor/receptor.pdbqt        rigid receptor for AutoDock Vina
  config/box.json                docking box center & size
"""
import argparse
import json
import subprocess
import urllib.request
from pathlib import Path

import numpy as np


def log(msg):
    print(f"[prep] {msg}", flush=True)


def download_pdb(pdb_id, dest):
    if dest.exists():
        log(f"{dest.name} already present, skipping download")
        return
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    log(f"downloading {url}")
    urllib.request.urlretrieve(url, dest)


def extract_ligand(pdb_path, resname, chain, out_pdb):
    """Pull the HETATM block for `resname`/`chain` and return its heavy-atom coordinates."""
    rows = []
    with open(pdb_path) as fh:
        for ln in fh:
            if ln.startswith("HETATM") and ln[17:20].strip() == resname and ln[21] == chain:
                if ln[76:78].strip() != "H" and not ln[12:16].strip().startswith("H"):
                    rows.append(ln)
    if not rows:
        raise RuntimeError(f"ligand {resname}/{chain} not found in {pdb_path}")
    with open(out_pdb, "w") as fh:
        fh.writelines(rows)
        fh.write("END\n")
    xyz = np.array([[float(r[30:38]), float(r[38:46]), float(r[46:54])] for r in rows])
    log(f"extracted {len(rows)} heavy atoms of {resname} -> {out_pdb.name}")
    return xyz


def define_box(xyz, buffer_a, min_a):
    center = xyz.mean(axis=0)
    extent = xyz.max(axis=0) - xyz.min(axis=0)
    size = np.maximum(extent + buffer_a, min_a)
    return center, size, extent


def clean_protein(raw_pdb, chain, out_pdb):
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(raw_pdb))

    # Identify gaps, then keep only INTERNAL missing loops (do not extend termini).
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        chain_len = len(list(chains[key[0]].residues()))
        if key[1] == 0 or key[1] == chain_len:
            del fixer.missingResidues[key]
    log(f"internal missing residues to model: {fixer.missingResidues}")

    # Revert nonstandard residues (CSD177 -> CYS) BEFORE removing heteroatoms,
    # otherwise the residue would be deleted and leave a gap.
    fixer.findNonstandardResidues()
    log(f"nonstandard residues: {[(r.name, repl) for r, repl in fixer.nonstandardResidues]}")
    fixer.replaceNonstandardResidues()

    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    n_res = sum(1 for _ in fixer.topology.residues())
    n_atom = sum(1 for _ in fixer.topology.atoms())
    log(f"wrote {out_pdb.name}: {n_res} residues, {n_atom} atoms")


def receptor_to_pdbqt(protein_pdb, out_pdbqt):
    # Rigid receptor for Vina: OpenBabel assigns AutoDock atom types + Gasteiger charges
    # and merges nonpolar hydrogens (-xr).
    subprocess.run(["obabel", str(protein_pdb), "-O", str(out_pdbqt), "-xr"], check=True)
    n = sum(1 for ln in open(out_pdbqt) if ln.startswith(("ATOM", "HETATM")))
    log(f"wrote {out_pdbqt.name} via OpenBabel ({n} receptor atoms)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", default="2R3I")
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--ligand", default="SCF")
    ap.add_argument("--box-buffer", type=float, default=10.0,
                    help="Angstrom added to each box dimension beyond the ligand envelope")
    ap.add_argument("--box-min", type=float, default=20.0,
                    help="minimum box edge length (Angstrom)")
    args = ap.parse_args()

    wd = Path(args.workdir).resolve()
    (wd / "receptor").mkdir(parents=True, exist_ok=True)
    (wd / "config").mkdir(parents=True, exist_ok=True)

    raw = wd / "receptor" / f"{args.pdb}.pdb"
    download_pdb(args.pdb, raw)

    ref_pdb = wd / "receptor" / "ref_ligand_xtal.pdb"
    xyz = extract_ligand(raw, args.ligand, args.chain, ref_pdb)
    center, size, extent = define_box(xyz, args.box_buffer, args.box_min)
    box = {
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x": round(float(size[0]), 3),
        "size_y": round(float(size[1]), 3),
        "size_z": round(float(size[2]), 3),
        "source": f"{args.pdb} ligand {args.ligand}/{args.chain}",
        "ligand_extent": [round(float(e), 2) for e in extent],
    }
    (wd / "config" / "box.json").write_text(json.dumps(box, indent=2))
    log(f"box center=({box['center_x']}, {box['center_y']}, {box['center_z']}) "
        f"size=({box['size_x']}, {box['size_y']}, {box['size_z']})")

    protein_pdb = wd / "receptor" / "cdk2_prepared.pdb"
    clean_protein(raw, args.chain, protein_pdb)
    receptor_to_pdbqt(protein_pdb, wd / "receptor" / "receptor.pdbqt")
    log("receptor preparation complete")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
02_dock.py — Dock the CDK2 screening library with AutoDock Vina.

- Rigid-receptor docking of every ligand in ligands/pdbqt/ into the ATP-site box (config/box.json).
- Ligands are distributed across CPU workers; each worker computes the receptor affinity maps
  once and reuses them (Vina Python API).
- Positive control: the native co-crystal inhibitor SCF is re-docked and the best pose is
  compared to its crystal coordinates (redocking RMSD) — a standard sanity check on the setup.

Outputs (under <workdir>):
  docking/docking_scores.csv     chembl_id,label,pchembl,mw,clogp,vina_score,status (sorted)
  docking/poses/<id>.pdbqt       best pose per ligand (gitignored)
  docking/redock_control.json    SCF redock score + RMSD-to-crystal
"""
import argparse
import csv
import json
import subprocess
import tempfile
from multiprocessing import Pool
from pathlib import Path


def log(msg):
    print(f"[dock] {msg}", flush=True)


def load_manifest(path):
    meta = {}
    if not path.exists():
        return meta
    with open(path) as fh:
        for row in csv.DictReader(fh):
            meta[row["chembl_id"]] = row
    return meta


def dock_worker(payload):
    chunk, receptor, center, size, ex, seed, poses_dir = payload
    from vina import Vina
    v = Vina(sf_name="vina", cpu=1, seed=seed, verbosity=0)
    v.set_receptor(receptor)
    v.compute_vina_maps(center=center, box_size=size)
    out = []
    for cid, path in chunk:
        try:
            v.set_ligand_from_file(path)
            v.dock(exhaustiveness=ex, n_poses=20)
            score = float(v.energies(n_poses=1)[0][0])
            v.write_poses(str(Path(poses_dir) / f"{cid}.pdbqt"), n_poses=1, overwrite=True)
            out.append((cid, score, "ok"))
        except Exception as exc:
            out.append((cid, None, f"fail:{exc.__class__.__name__}"))
    return out


def chunkify(items, n):
    chunks = [[] for _ in range(n)]
    for i, it in enumerate(items):
        chunks[i % n].append(it)
    return [c for c in chunks if c]


# ---------- redocking positive control ----------

def fetch_ligand_smiles(comp_id):
    import requests
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}", timeout=60)
    r.raise_for_status()
    desc = r.json().get("rcsb_chem_comp_descriptor", {})
    return desc.get("smiles_stereo") or desc.get("smiles") or desc.get("SMILES")


def prepare_ligand_pdbqt(smiles, out_pdbqt, seed=0xC0FFEE):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = seed
    AllChem.EmbedMolecule(mol, p)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    s = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mol)[0])
    txt = s[0] if isinstance(s, tuple) else s
    Path(out_pdbqt).write_text(txt)


def best_rms_to_crystal(smiles, crystal_pdb, docked_pdbqt):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    tmpl = Chem.MolFromSmiles(smiles)
    cryst = Chem.MolFromPDBFile(str(crystal_pdb), sanitize=False, removeHs=True)
    cryst = AllChem.AssignBondOrdersFromTemplate(tmpl, cryst)
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tf:
        sdf = tf.name
    subprocess.run(["obabel", str(docked_pdbqt), "-O", sdf], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    docked = Chem.MolFromMolFile(sdf, removeHs=True)
    docked = AllChem.AssignBondOrdersFromTemplate(tmpl, docked)
    return float(AllChem.GetBestRMS(Chem.RemoveHs(docked), Chem.RemoveHs(cryst)))


def redock_control(receptor, center, size, ex, seed, wd, control_id):
    from vina import Vina
    res = {"control_id": control_id}
    try:
        smiles = fetch_ligand_smiles(control_id)
        res["smiles"] = smiles
        ctrl_pdbqt = wd / "docking" / f"control_{control_id}.pdbqt"
        prepare_ligand_pdbqt(smiles, ctrl_pdbqt)
        v = Vina(sf_name="vina", cpu=4, seed=seed, verbosity=0)
        v.set_receptor(receptor)
        v.compute_vina_maps(center=center, box_size=size)
        v.set_ligand_from_file(str(ctrl_pdbqt))
        v.dock(exhaustiveness=ex, n_poses=20)
        res["vina_score"] = round(float(v.energies(n_poses=1)[0][0]), 2)
        pose = wd / "docking" / f"control_{control_id}_pose.pdbqt"
        v.write_poses(str(pose), n_poses=1, overwrite=True)
        res["redock_rmsd"] = round(
            best_rms_to_crystal(smiles, wd / "receptor" / "ref_ligand_xtal.pdb", pose), 2)
        res["status"] = "ok"
    except Exception as exc:
        res["status"] = f"fail:{exc.__class__.__name__}:{exc}"
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--cpus", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--control-id", default="SCF")
    args = ap.parse_args()

    wd = Path(args.workdir).resolve()
    box = json.loads((wd / "config" / "box.json").read_text())
    center = [box["center_x"], box["center_y"], box["center_z"]]
    size = [box["size_x"], box["size_y"], box["size_z"]]
    receptor = str(wd / "receptor" / "receptor.pdbqt")
    poses_dir = wd / "docking" / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)

    ligands = sorted((lig.stem, str(lig)) for lig in (wd / "ligands" / "pdbqt").glob("*.pdbqt"))
    log(f"docking {len(ligands)} ligands, exhaustiveness={args.exhaustiveness}, cpus={args.cpus}")

    chunks = chunkify(ligands, args.cpus)
    payloads = [(c, receptor, center, size, args.seed, args.seed, str(poses_dir)) for c in chunks]
    results = []
    with Pool(len(payloads)) as pool:
        for part in pool.map(dock_worker, payloads):
            results.extend(part)
    scores = {cid: (sc, st) for cid, sc, st in results}
    n_ok = sum(1 for sc, st in scores.values() if st == "ok")
    log(f"docking complete: {n_ok}/{len(ligands)} ok")

    # redocking positive control
    log(f"redocking native control {args.control_id} ...")
    ctrl = redock_control(receptor, center, size, args.exhaustiveness, args.seed, wd, args.control_id)
    (wd / "docking" / "redock_control.json").write_text(json.dumps(ctrl, indent=2))
    log(f"control: {ctrl}")

    # merge with manifest, sort by score
    meta = load_manifest(wd / "ligands" / "ligands.csv")
    rows = []
    for cid, (sc, st) in scores.items():
        m = meta.get(cid, {})
        rows.append([cid, m.get("label", "?"), m.get("pchembl", ""), m.get("mw", ""),
                     m.get("clogp", ""), "" if sc is None else round(sc, 2), st])
    rows.sort(key=lambda r: (r[5] == "", r[5] if r[5] != "" else 0.0))
    out = wd / "docking" / "docking_scores.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chembl_id", "label", "pchembl", "mw", "clogp", "vina_score", "status"])
        w.writerows(rows)
    log(f"wrote {out}")


if __name__ == "__main__":
    main()

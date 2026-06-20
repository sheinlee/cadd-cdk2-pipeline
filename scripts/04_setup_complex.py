#!/usr/bin/env python
"""
04_setup_complex.py — Parametrize one docked hit and build an Amber-ready complex for MD/MM-GBSA.

For a given ligand (chembl_id), starting from its best Vina pose:
  1. Rebuild a chemically-correct 3D ligand: PDBQT -> SDF (OpenBabel) -> assign bond orders from
     the manifest SMILES (RDKit) -> add hydrogens on the docked coordinates.
  2. GAFF2 parameters + AM1-BCC charges via antechamber/parmchk2.
  3. tleap: ff14SB protein + GAFF2 ligand -> dry complex/receptor/ligand prmtops (for MM-GBSA)
     and a TIP3P-solvated, neutralized complex (for MD).

Outputs (under <workdir>/md/<cid>/):
  ligand.mol2, ligand.frcmod
  complex_dry.prmtop, receptor_dry.prmtop, ligand_dry.prmtop   (MM-GBSA endpoints)
  complex_solv.prmtop, complex_solv.inpcrd, complex_solv.pdb    (MD input)
"""
import argparse
import csv
import subprocess
from pathlib import Path


def log(msg):
    print(f"[setup] {msg}", flush=True)


def run(cmd, cwd=None):
    log("$ " + " ".join(str(c) for c in cmd))
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed (rc={r.returncode})\nSTDOUT:\n{r.stdout[-2000:]}\n"
                           f"STDERR:\n{r.stderr[-2000:]}")
    return r


def smiles_for(manifest, cid):
    with open(manifest) as fh:
        for row in csv.DictReader(fh):
            if row["chembl_id"] == cid:
                return row["smiles"]
    raise KeyError(f"{cid} not in {manifest}")


def rebuild_ligand(pose_pdbqt, smiles, out_mol2, work):
    """Docked PDBQT pose -> correct-topology, H-added ligand mol2. Returns net formal charge."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    raw_sdf = work / "pose_raw.sdf"
    run(["obabel", str(pose_pdbqt), "-O", str(raw_sdf)])
    tmpl = Chem.MolFromSmiles(smiles)
    pose = Chem.MolFromMolFile(str(raw_sdf), removeHs=True, sanitize=False)
    pose = AllChem.AssignBondOrdersFromTemplate(tmpl, pose)
    pose = Chem.AddHs(pose, addCoords=True)
    net_charge = Chem.GetFormalCharge(pose)
    h_sdf = work / "pose_H.sdf"
    Chem.MolToMolFile(pose, str(h_sdf))
    run(["obabel", str(h_sdf), "-O", str(out_mol2)])
    return net_charge


def write_leap(work, ligand_mol2, ligand_frcmod, receptor_pdb, solvate_buffer):
    leap = f"""source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams {ligand_frcmod.name}
LIG = loadmol2 {ligand_mol2.name}
saveamberparm LIG ligand_dry.prmtop ligand_dry.inpcrd
REC = loadpdb {receptor_pdb.name}
saveamberparm REC receptor_dry.prmtop receptor_dry.inpcrd
COM = combine {{ REC LIG }}
saveamberparm COM complex_dry.prmtop complex_dry.inpcrd
solvateOct COM TIP3PBOX {solvate_buffer}
addIonsRand COM Na+ 0
addIonsRand COM Cl- 0
saveamberparm COM complex_solv.prmtop complex_solv.inpcrd
savepdb COM complex_solv.pdb
quit
"""
    (work / "leap.in").write_text(leap)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--cid", required=True)
    ap.add_argument("--pose", default=None, help="docked pose pdbqt (default docking/poses/<cid>.pdbqt)")
    ap.add_argument("--solvate-buffer", type=float, default=10.0)
    args = ap.parse_args()

    wd = Path(args.workdir).resolve()
    work = wd / "md" / args.cid
    work.mkdir(parents=True, exist_ok=True)
    pose = Path(args.pose) if args.pose else wd / "docking" / "poses" / f"{args.cid}.pdbqt"
    smiles = smiles_for(wd / "ligands" / "ligands.csv", args.cid)
    log(f"hit {args.cid}  pose={pose.name}  smiles={smiles}")

    ligand_mol2 = work / "ligand.mol2"
    net_charge = rebuild_ligand(pose, smiles, work / "ligand_in.mol2", work)
    log(f"net formal charge = {net_charge}")

    # antechamber (GAFF2 + AM1-BCC) then parmchk2
    run(["antechamber", "-i", "ligand_in.mol2", "-fi", "mol2", "-o", "ligand.mol2", "-fo", "mol2",
         "-c", "bcc", "-nc", net_charge, "-at", "gaff2", "-rn", "LIG", "-pf", "y"], cwd=work)
    run(["parmchk2", "-i", "ligand.mol2", "-f", "mol2", "-o", "ligand.frcmod", "-s", "gaff2"], cwd=work)

    # clean protein for tleap (re-add H consistently)
    run(["pdb4amber", "-i", str(wd / "receptor" / "cdk2_prepared.pdb"),
         "-o", str(work / "receptor_amber.pdb"), "-y"], cwd=work)

    write_leap(work, ligand_mol2, work / "ligand.frcmod", work / "receptor_amber.pdb",
               args.solvate_buffer)
    r = run(["tleap", "-f", "leap.in"], cwd=work)
    tail = r.stdout[-1500:]
    log("tleap done; checking outputs")
    for f in ["complex_solv.prmtop", "complex_solv.inpcrd", "complex_dry.prmtop",
              "receptor_dry.prmtop", "ligand_dry.prmtop"]:
        ok = (work / f).exists()
        log(f"  {'OK ' if ok else 'MISSING'} {f}")
        if not ok:
            raise RuntimeError(f"tleap did not produce {f}\nleap tail:\n{tail}")
    log(f"complex ready in {work}")


if __name__ == "__main__":
    main()

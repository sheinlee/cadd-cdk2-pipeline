#!/usr/bin/env python
"""
01_fetch_ligands.py — Build the CDK2 docking library of actives + property-matched decoys.

Two interchangeable sources (--source):
  dude   (default) : the standard DUD-E CDK2 benchmark — curated actives + property-matched
                     decoys, served as static files (robust, widely recognized).
  chembl           : live ChEMBL — actives are molecules with a measured CDK2 (CHEMBL301)
                     potency pChEMBL >= --active-pchembl (IC50/Ki/Kd), decoys are random
                     ChEMBL small molecules (not CDK2 binders). Demonstrates live-data handling
                     but depends on the EBI API being healthy.

In both cases decoys are nearest-(MW, cLogP) matched to the selected actives locally, so
enrichment cannot be won on trivial size/lipophilicity differences. Every ligand is then
standardized (largest fragment) -> 3D embedded (ETKDGv3) -> MMFF-optimized -> written as an
AutoDock PDBQT via Meeko.

Outputs (under <workdir>):
  ligands/ligands.csv        manifest: chembl_id,label,pchembl,mw,clogp,prep_status,smiles
  ligands/pdbqt/<id>.pdbqt   docking inputs (gitignored; regenerable)
"""
import argparse
import csv
import random
import time
from pathlib import Path

import numpy as np
import requests

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
ORGANIC = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "B"}


def log(msg):
    print(f"[ligands] {msg}", flush=True)


def chembl_get(session, endpoint, params, retries=4):
    url = f"{CHEMBL}/{endpoint}.json"
    params = dict(params)
    for attempt in range(retries):
        try:
            r = session.get(url, params=params, timeout=120)
            r.raise_for_status()
            return r.json()
        except Exception as exc:
            if attempt == retries - 1:
                raise
            log(f"  request retry {attempt + 1} ({exc.__class__.__name__})")
            time.sleep(3 * (attempt + 1))


def fetch_actives(session, target, pchembl_min, types, max_records):
    collected = {}
    offset, limit, total = 0, 1000, None
    while True:
        d = chembl_get(session, "activity", {
            "target_chembl_id": target, "pchembl_value__gte": pchembl_min,
            "standard_type__in": ",".join(types), "limit": limit, "offset": offset})
        if total is None:
            total = d["page_meta"]["total_count"]
            log(f"actives: {total} activity records for {target}")
        for a in d["activities"]:
            smi, cid, pv = a.get("canonical_smiles"), a.get("molecule_chembl_id"), a.get("pchembl_value")
            if not smi or not cid or pv is None:
                continue
            pv = float(pv)
            if cid not in collected or pv > collected[cid][1]:
                collected[cid] = (smi, pv)
        offset += limit
        if offset >= min(total, max_records):
            break
    log(f"actives: {len(collected)} unique molecules")
    return collected


def fetch_dude(session, lig_dir, n_decoy_candidates, seed,
               base="https://dude.docking.org/targets/cdk2"):
    """Download (cache) the DUD-E CDK2 actives/decoys SMILES sets -> raw {id: smiles}.
    DUD-E decoys are property-matched to the actives by construction; we subsample a candidate
    pool of decoys (cheaper than standardizing all ~28k) for the local nearest-property match."""
    lig_dir.mkdir(parents=True, exist_ok=True)
    files = {"active": lig_dir / "dude_actives.ism", "decoy": lig_dir / "dude_decoys.ism"}
    urls = {"active": f"{base}/actives_final.ism", "decoy": f"{base}/decoys_final.ism"}
    for k in files:
        if not files[k].exists() or files[k].stat().st_size == 0:
            log(f"downloading DUD-E {k}: {urls[k]}")
            r = session.get(urls[k], timeout=120)
            r.raise_for_status()
            files[k].write_text(r.text)
    def parse(path, id_col):  # line: "SMILES <id> [<chembl_id>]"
        out = {}
        for line in path.read_text().splitlines():
            p = line.split()
            if len(p) > id_col:
                out[p[id_col]] = p[0]
        return out
    actives = parse(files["active"], 2)   # SMILES dude_id CHEMBL_id
    decoys = parse(files["decoy"], 1)     # SMILES decoy_id
    log(f"DUD-E CDK2: {len(actives)} actives, {len(decoys)} decoys available")
    dec_ids = list(decoys)
    random.Random(seed).shuffle(dec_ids)
    return actives, {i: decoys[i] for i in dec_ids[:n_decoy_candidates]}


def fetch_decoy_pool(session, n_candidates, seed, page=200):
    """Randomly sample small molecules from ChEMBL (no server-side property-range filter, which
    triggers 500s); property selection is done locally. One bad page is skipped, not fatal."""
    rng = random.Random(seed)
    base = {"molecule_type": "Small molecule"}
    total = chembl_get(session, "molecule", {**base, "limit": 1})["page_meta"]["total_count"]
    n_pages = max(1, n_candidates // page)
    log(f"small-molecule universe: {total}; sampling ~{n_candidates} candidates over {n_pages} pages")
    pool = {}
    max_off = max(0, total - page)
    for _ in range(n_pages):
        off = rng.randint(0, max_off)
        try:
            d = chembl_get(session, "molecule", {**base, "limit": page, "offset": off})
        except Exception as exc:
            log(f"  skip decoy page @offset {off} ({exc.__class__.__name__})")
            continue
        for m in d["molecules"]:
            ms = m.get("molecule_structures") or {}
            smi, cid = ms.get("canonical_smiles"), m.get("molecule_chembl_id")
            if smi and cid:
                pool[cid] = smi
    log(f"decoy candidate pool: {len(pool)} unique molecules")
    return pool


def standardize(smiles):
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        return None
    return mol, Descriptors.MolWt(mol), Crippen.MolLogP(mol)


def druglike(mol, mw, logp, mw_lo=250, mw_hi=600, logp_lo=-1.0, logp_hi=6.0):
    if not (mw_lo <= mw <= mw_hi) or not (logp_lo <= logp <= logp_hi):
        return False
    return all(a.GetSymbol() in ORGANIC for a in mol.GetAtoms())


def build_records(raw, label, exclude_keys=None):
    """raw: {cid: smiles or (smiles, pchembl)} -> list of standardized, drug-like record dicts."""
    from rdkit import Chem
    out, seen_keys = [], set()
    for cid, val in raw.items():
        smi, pchembl = (val if isinstance(val, tuple) else (val, None))
        std = standardize(smi)
        if std is None:
            continue
        mol, mw, logp = std
        if not druglike(mol, mw, logp):
            continue
        key = Chem.MolToInchiKey(mol)
        if not key or key in seen_keys or (exclude_keys and key in exclude_keys):
            continue
        seen_keys.add(key)
        out.append({"cid": cid, "smiles": Chem.MolToSmiles(mol), "mol": mol,
                    "mw": mw, "logp": logp, "pchembl": pchembl, "label": label, "inchikey": key})
    return out


def match_decoys(actives, decoys, ratio, seed):
    rng = random.Random(seed)
    mw_sd = np.std([d["mw"] for d in decoys]) or 1.0
    lp_sd = np.std([d["logp"] for d in decoys]) or 1.0
    targets = [a for a in actives for _ in range(ratio)]
    rng.shuffle(targets)
    used, chosen = set(), []
    for a in targets:
        best_i, best_d = None, 1e18
        for i, d in enumerate(decoys):
            if i in used:
                continue
            dist = abs(a["mw"] - d["mw"]) / mw_sd + abs(a["logp"] - d["logp"]) / lp_sd
            if dist < best_d:
                best_d, best_i = dist, i
        if best_i is None:
            break
        used.add(best_i)
        chosen.append(decoys[best_i])
    return chosen


def prepare_pdbqt(mol, out_path, seed):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    m = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(m, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(m, params) != 0:
            return "embed_fail"
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        pass
    try:
        setups = MoleculePreparation().prepare(m)
        s = PDBQTWriterLegacy.write_string(setups[0])
        txt = s[0] if isinstance(s, tuple) else s
    except Exception as exc:
        return f"meeko_fail:{exc.__class__.__name__}"
    if not txt or "ATOM" not in txt:
        return "pdbqt_empty"
    out_path.write_text(txt)
    return "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--source", choices=["dude", "chembl"], default="dude",
                    help="dude = standard DUD-E CDK2 benchmark (robust); chembl = live ChEMBL")
    ap.add_argument("--target", default="CHEMBL301")
    ap.add_argument("--active-pchembl", type=float, default=6.5)
    ap.add_argument("--n-active", type=int, default=60)
    ap.add_argument("--ratio", type=int, default=4, help="decoys per active")
    ap.add_argument("--decoy-pages", type=int, default=10)
    ap.add_argument("--max-active-records", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    wd = Path(args.workdir).resolve()
    lig_dir = wd / "ligands"
    pdbqt_dir = lig_dir / "pdbqt"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"Accept": "application/json"})

    # --- actives ---
    rng = random.Random(args.seed)
    if args.source == "dude":
        raw_act, raw_dec = fetch_dude(session, lig_dir, args.n_active * args.ratio * 5, args.seed)
        actives = build_records(raw_act, "active")
        rng.shuffle(actives)
        actives = actives[:args.n_active]
    else:
        raw_act = fetch_actives(session, args.target, args.active_pchembl,
                                ["IC50", "Ki", "Kd"], args.max_active_records)
        actives = build_records(raw_act, "active")
        rng.shuffle(actives)
        actives = sorted(actives[:args.n_active], key=lambda r: -(r["pchembl"] or 0))
    active_keys = {a["inchikey"] for a in actives}
    log(f"actives selected: {len(actives)} (source={args.source})")

    # --- decoys: property-matched on (MW, cLogP), excluding actives ---
    if args.source == "dude":
        decoys_all = build_records(raw_dec, "decoy", exclude_keys=active_keys)
    else:
        mws = np.array([a["mw"] for a in actives])
        lps = np.array([a["logp"] for a in actives])
        log(f"actives MW {mws.min():.0f}-{mws.max():.0f}, cLogP {lps.min():.1f}-{lps.max():.1f}")
        pool_raw = fetch_decoy_pool(session, n_candidates=max(6000, args.n_active * args.ratio * 25),
                                    seed=args.seed)
        decoys_all = build_records(pool_raw, "decoy", exclude_keys=active_keys)
    log(f"decoy candidates drug-like (excl. actives): {len(decoys_all)}")
    decoys = match_decoys(actives, decoys_all, args.ratio, args.seed)
    log(f"decoys matched: {len(decoys)} (target {len(actives) * args.ratio})")

    # --- 3D + PDBQT prep + manifest ---
    ligands = actives + decoys
    manifest = lig_dir / "ligands.csv"
    n_ok = 0
    with open(manifest, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chembl_id", "label", "pchembl", "mw", "clogp", "prep_status", "smiles"])
        for i, r in enumerate(ligands):
            status = prepare_pdbqt(r["mol"], pdbqt_dir / f"{r['cid']}.pdbqt", seed=args.seed + i)
            n_ok += status == "ok"
            w.writerow([r["cid"], r["label"], f"{r['pchembl']:.2f}" if r["pchembl"] else "",
                        f"{r['mw']:.1f}", f"{r['logp']:.2f}", status, r["smiles"]])
            if (i + 1) % 50 == 0:
                log(f"  prepared {i + 1}/{len(ligands)} ...")
    n_act = sum(1 for r in ligands if r["label"] == "active")
    log(f"DONE: {len(ligands)} ligands ({n_act} active / {len(ligands) - n_act} decoy), "
        f"{n_ok} pdbqt OK -> {manifest}")


if __name__ == "__main__":
    main()

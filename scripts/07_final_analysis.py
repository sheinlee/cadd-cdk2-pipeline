#!/usr/bin/env python
"""
07_final_analysis.py — Combine docking + MD + MM-GBSA into the final ranking and figures.

For the hits advanced to MD it merges:
  * Vina docking score
  * MD pose stability (ligand heavy-atom RMSD-to-start)
  * MM-GBSA binding free energy (dG, with per-residue decomposition for the best hit)

and produces the headline comparison: does the physics-based MM-GBSA re-rank the docking hits,
and which poses are stable?

Outputs (under <workdir>/results/):
  final_ranking.csv          merged table (sorted by MM-GBSA dG)
  docking_vs_mmgbsa.png      Vina score vs MM-GBSA dG (re-ranking)
  rmsd_stability.png         ligand RMSD vs time for each MD hit
  mmgbsa_decomposition.png   top per-residue contributions for the best stable hit
  final_summary.json
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def log(msg):
    print(f"[final] {msg}", flush=True)


def load_docking(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["status"] == "ok" and r["vina_score"]:
                out[r["chembl_id"]] = {"label": r["label"], "vina": float(r["vina_score"])}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=".")
    args = ap.parse_args()
    wd = Path(args.workdir).resolve()
    res = wd / "results"
    res.mkdir(parents=True, exist_ok=True)

    dock = load_docking(wd / "docking" / "docking_scores.csv")
    rows = []
    for mj in sorted((wd / "md").glob("*/mmgbsa.json")):
        cid = mj.parent.name
        mg = json.loads(mj.read_text())
        sj = mj.parent / "md_summary.json"
        md = json.loads(sj.read_text()) if sj.exists() else {}
        rows.append({
            "cid": cid,
            "label": dock.get(cid, {}).get("label", "?"),
            "vina": dock.get(cid, {}).get("vina"),
            "dG_mmgbsa": mg.get("dG_bind_kcal_mol"),
            "dG_std": mg.get("std"),
            "lig_rmsd_mean": md.get("ligand_rmsd_mean"),
            "lig_rmsd_final": md.get("ligand_rmsd_final"),
            "stable": md.get("stable"),
        })
    if not rows:
        log("no MM-GBSA results found yet — run the MD pipeline first")
        return

    rows.sort(key=lambda r: (r["dG_mmgbsa"] is None, r["dG_mmgbsa"] if r["dG_mmgbsa"] is not None else 0))
    out_csv = res / "final_ranking.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {out_csv} ({len(rows)} hits)")

    # --- docking vs MM-GBSA ---
    have = [r for r in rows if r["vina"] is not None and r["dG_mmgbsa"] is not None]
    if have:
        plt.figure(figsize=(5, 4.2))
        for r in have:
            stable = r["stable"]
            plt.scatter(r["vina"], r["dG_mmgbsa"],
                        c=("#2ca02c" if stable else "#d62728"),
                        marker=("o" if r["label"] == "active" else "s"), s=60,
                        edgecolors="k", linewidths=0.5)
            plt.annotate(r["cid"][-4:], (r["vina"], r["dG_mmgbsa"]), fontsize=6,
                         xytext=(3, 3), textcoords="offset points")
        if len(have) > 2:
            x = np.array([r["vina"] for r in have]); y = np.array([r["dG_mmgbsa"] for r in have])
            rho = np.corrcoef(x, y)[0, 1]
            plt.title(f"Docking vs MM-GBSA (Pearson r={rho:.2f})")
        plt.xlabel("Vina docking score (kcal/mol)")
        plt.ylabel("MM-GBSA dG_bind (kcal/mol)")
        plt.scatter([], [], c="#2ca02c", edgecolors="k", label="stable pose")
        plt.scatter([], [], c="#d62728", edgecolors="k", label="drifted pose")
        plt.legend(fontsize=8); plt.tight_layout()
        plt.savefig(res / "docking_vs_mmgbsa.png", dpi=150); plt.close()

    # --- RMSD stability over time ---
    plt.figure(figsize=(5.5, 4))
    for r in rows:
        f = wd / "md" / r["cid"] / "rmsd.csv"
        if f.exists():
            d = np.loadtxt(f, delimiter=",", skiprows=1)
            if d.ndim == 2 and len(d):
                plt.plot(d[:, 0], d[:, 1], lw=1.2,
                         label=f"{r['cid'][-5:]} ({'stable' if r['stable'] else 'drift'})")
    plt.axhline(3.0, color="grey", ls="--", lw=1)
    plt.xlabel("time (ns)"); plt.ylabel("ligand RMSD to start (A)")
    plt.title("Pose stability during MD")
    plt.legend(fontsize=7, ncol=2); plt.tight_layout()
    plt.savefig(res / "rmsd_stability.png", dpi=150); plt.close()

    # --- per-residue decomposition for the best stable hit ---
    best = next((r for r in rows if r["stable"] and r["dG_mmgbsa"] is not None), None)
    if best:
        mg = json.loads((wd / "md" / best["cid"] / "mmgbsa.json").read_text())
        tr = mg.get("top_residues", [])
        if tr:
            names = [t["residue"] for t in tr][::-1]
            vals = [t["kcal_mol"] for t in tr][::-1]
            plt.figure(figsize=(5, 4))
            plt.barh(names, vals, color="#1f77b4")
            plt.xlabel("per-residue contribution (kcal/mol)")
            plt.title(f"MM-GBSA hotspots: {best['cid']}")
            plt.tight_layout()
            plt.savefig(res / "mmgbsa_decomposition.png", dpi=150); plt.close()

    summary = {
        "n_hits_md": len(rows),
        "n_stable": sum(1 for r in rows if r["stable"]),
        "best_hit": best["cid"] if best else None,
        "best_dG_mmgbsa": best["dG_mmgbsa"] if best else None,
        "ranking": rows,
    }
    (res / "final_summary.json").write_text(json.dumps(summary, indent=2))
    log(f"final summary: best={summary['best_hit']} dG={summary['best_dG_mmgbsa']} "
        f"stable={summary['n_stable']}/{summary['n_hits_md']}")


if __name__ == "__main__":
    main()

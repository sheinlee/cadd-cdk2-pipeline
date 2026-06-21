#!/usr/bin/env python
"""
08_fep_analysis.py — Analyze relative binding FEP results vs experiment (CDK2 / JACS benchmark).

Runs `openfe gather --report ddg` over the completed transformations to get per-edge predicted
DDG (= DG_j - DG_i) with uncertainty, computes the experimental DDG from the benchmark IC50s
(ligands.yml), and produces the predicted-vs-experimental correlation (MUE / RMSE / Pearson).

Outputs (under <fepdir>/results/):
  fep_ddg.tsv          raw openfe gather output (per edge)
  fep_ddg.csv          merged predicted vs experimental DDG
  fep_correlation.png  predicted vs experimental DDG
  fep_summary.json     MUE, RMSE, Pearson, per-edge table
"""
import argparse
import csv
import json
import math
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

RT = 0.5961  # kcal/mol at 298 K


def log(m):
    print(f"[fep] {m}", flush=True)


def experimental_dG(ligands_yml):
    d = yaml.safe_load(open(ligands_yml))
    factor = {"m": 1.0, "mm": 1e-3, "um": 1e-6, "nm": 1e-9, "pm": 1e-12}
    out = {}
    for name, rec in d.items():
        m = rec["measurement"]
        out[name] = RT * math.log(float(m["value"]) * factor[m["unit"].lower()])
    return out


def parse_gather_tsv(tsv):
    """Return list of (lig_i, lig_j, ddg_pred, unc)."""
    rows = []
    with open(tsv) as fh:
        reader = list(csv.reader(fh, delimiter="\t"))
    header = [h.strip().lower() for h in reader[0]]

    def col(*subs):
        for i, h in enumerate(header):
            if any(s in h for s in subs):
                return i
        return None

    ci, cj = col("ligand_i", "ligand_a"), col("ligand_j", "ligand_b")
    cd = col("ddg")
    cu = col("uncertainty", "error", "unc")
    for r in reader[1:]:
        if not r or len(r) <= max(x for x in (ci, cj, cd) if x is not None):
            continue
        try:
            rows.append((r[ci].strip(), r[cj].strip(), float(r[cd]),
                         float(r[cu]) if cu is not None and r[cu] not in ("", "nan") else float("nan")))
        except (ValueError, IndexError):
            continue
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fepdir", default=".")
    ap.add_argument("--ligands-yml",
                    default="plb/data/cdk2/00_data/ligands.yml")
    args = ap.parse_args()
    fepdir = Path(args.fepdir).resolve()
    res = fepdir / "results"
    res.mkdir(parents=True, exist_ok=True)

    tsv = res / "fep_ddg.tsv"
    log("running openfe gather --report ddg ...")
    subprocess.run(["openfe", "gather", str(res), "--report", "ddg", "--tsv",
                    "-o", str(tsv), "--allow-partial"], check=True)

    pred = parse_gather_tsv(tsv)
    exp_dg = experimental_dG(fepdir / args.ligands_yml)
    log(f"{len(pred)} edges from gather; {len(exp_dg)} experimental dG values")

    merged = []
    for li, lj, ddg, unc in pred:
        if li in exp_dg and lj in exp_dg:
            ddg_exp = exp_dg[lj] - exp_dg[li]
            merged.append({"edge": f"{li}->{lj}", "ddg_pred": round(ddg, 2),
                           "unc": round(unc, 2) if not math.isnan(unc) else None,
                           "ddg_exp": round(ddg_exp, 2)})
    if not merged:
        log("no edges with both prediction and experimental data yet")
        return

    with open(res / "fep_ddg.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["edge", "ddg_pred", "unc", "ddg_exp"])
        w.writeheader()
        w.writerows(merged)

    p = np.array([m["ddg_pred"] for m in merged])
    e = np.array([m["ddg_exp"] for m in merged])
    u = np.array([m["unc"] if m["unc"] is not None else 0.0 for m in merged])
    mue = float(np.mean(np.abs(p - e)))
    rmse = float(np.sqrt(np.mean((p - e) ** 2)))
    pearson = float(np.corrcoef(p, e)[0, 1]) if len(p) > 1 else float("nan")

    lim = max(2.5, np.max(np.abs(np.concatenate([p, e]))) + 0.5)
    plt.figure(figsize=(4.8, 4.8))
    plt.fill_between([-lim, lim], [-lim - 1, lim - 1], [-lim + 1, lim + 1],
                     color="grey", alpha=0.15, label="±1 kcal/mol")
    plt.plot([-lim, lim], [-lim, lim], "k--", lw=1)
    plt.errorbar(e, p, yerr=u, fmt="o", ms=8, capsize=3, color="#1f77b4")
    for m in merged:
        plt.annotate(m["edge"].replace("lig_", ""), (m["ddg_exp"], m["ddg_pred"]),
                     fontsize=6, xytext=(4, 4), textcoords="offset points")
    plt.xlim(-lim, lim); plt.ylim(-lim, lim)
    plt.xlabel("experimental ΔΔG (kcal/mol)")
    plt.ylabel("FEP predicted ΔΔG (kcal/mol)")
    plt.title(f"CDK2 relative FEP vs experiment\nMUE={mue:.2f}  RMSE={rmse:.2f}  r={pearson:.2f}  (n={len(p)})")
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(res / "fep_correlation.png", dpi=150); plt.close()

    summary = {"n_edges": len(merged), "MUE": round(mue, 2), "RMSE": round(rmse, 2),
               "pearson_r": round(pearson, 2) if not math.isnan(pearson) else None,
               "edges": merged}
    (res / "fep_summary.json").write_text(json.dumps(summary, indent=2))
    log(f"summary: {summary}")


if __name__ == "__main__":
    main()

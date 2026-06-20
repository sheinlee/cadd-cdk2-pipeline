#!/usr/bin/env python
"""
03_analyze_docking.py — Validate the docking screen against known CDK2 actives.

Treats labelled actives (ChEMBL potent binders) as positives and property-matched decoys as
negatives, and asks how well the Vina score ranks actives above decoys:

  * ROC-AUC
  * Enrichment Factor at 1 / 5 / 10 %
  * score distributions (actives vs decoys)

Outputs (under <workdir>):
  results/docking_metrics.json      AUC, EF, counts, redocking RMSD
  results/roc_curve.png
  results/score_distribution.png
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
    print(f"[analyze] {msg}", flush=True)


def load_scores(path):
    actives, decoys, rows = [], [], []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["status"] != "ok" or r["vina_score"] == "":
                continue
            score = float(r["vina_score"])
            rows.append((r["chembl_id"], r["label"], score))
            (actives if r["label"] == "active" else decoys).append(score)
    return rows, np.array(actives), np.array(decoys)


def enrichment_factor(rows, frac):
    # better (more negative) score ranks first
    ordered = sorted(rows, key=lambda x: x[2])
    n = len(ordered)
    n_top = max(1, int(round(frac * n)))
    n_act = sum(1 for _, lab, _ in ordered if lab == "active")
    if n_act == 0:
        return float("nan")
    top_act = sum(1 for _, lab, _ in ordered[:n_top] if lab == "active")
    return (top_act / n_top) / (n_act / n)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=".")
    args = ap.parse_args()
    wd = Path(args.workdir).resolve()
    res = wd / "results"
    res.mkdir(parents=True, exist_ok=True)

    rows, act, dec = load_scores(wd / "docking" / "docking_scores.csv")
    log(f"actives={len(act)} decoys={len(dec)}")

    from sklearn.metrics import roc_auc_score, roc_curve
    y_true = np.array([1] * len(act) + [0] * len(dec))
    y_score = np.concatenate([-act, -dec])      # higher = better -> negate Vina score
    auc = float(roc_auc_score(y_true, y_score))
    fpr, tpr, _ = roc_curve(y_true, y_score)

    efs = {f"EF{int(f*100)}": round(enrichment_factor(rows, f), 2) for f in (0.01, 0.05, 0.10)}

    control = {}
    ctrl_path = wd / "docking" / "redock_control.json"
    if ctrl_path.exists():
        control = json.loads(ctrl_path.read_text())

    metrics = {
        "n_active": len(act), "n_decoy": len(dec),
        "roc_auc": round(auc, 3), **efs,
        "active_score_mean": round(float(act.mean()), 2),
        "decoy_score_mean": round(float(dec.mean()), 2),
        "active_score_best": round(float(act.min()), 2),
        "redock_control": {k: control.get(k) for k in ("control_id", "vina_score", "redock_rmsd", "status")},
    }
    (res / "docking_metrics.json").write_text(json.dumps(metrics, indent=2))
    log(f"metrics: {metrics}")

    # ROC curve
    plt.figure(figsize=(4.5, 4.5))
    plt.plot(fpr, tpr, lw=2, label=f"Vina (AUC={auc:.2f})")
    plt.plot([0, 1], [0, 1], "k--", lw=1, label="random")
    plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
    plt.title("CDK2 docking enrichment (ROC)")
    plt.legend(loc="lower right"); plt.tight_layout()
    plt.savefig(res / "roc_curve.png", dpi=150); plt.close()

    # score distribution
    plt.figure(figsize=(5, 4))
    bins = np.linspace(min(act.min(), dec.min()), max(act.max(), dec.max()), 30)
    plt.hist(dec, bins=bins, alpha=0.6, label=f"decoys (n={len(dec)})", color="#888")
    plt.hist(act, bins=bins, alpha=0.6, label=f"actives (n={len(act)})", color="#d62728")
    if control.get("vina_score") is not None:
        plt.axvline(control["vina_score"], color="navy", ls="--",
                    label=f"native {control.get('control_id','ctrl')} ({control['vina_score']})")
    plt.xlabel("Vina score (kcal/mol)"); plt.ylabel("count")
    plt.title("Docking score: actives vs decoys")
    plt.legend(); plt.tight_layout()
    plt.savefig(res / "score_distribution.png", dpi=150); plt.close()
    log(f"wrote figures + metrics to {res}")


if __name__ == "__main__":
    main()

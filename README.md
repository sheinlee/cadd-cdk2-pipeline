# CDK2 Structure-Based Virtual Screening — Docking → MD → MM-GBSA

> **Self-initiated practice project.** A complete, scriptable, reproducible computational
> pipeline for structure-based virtual screening against **cyclin-dependent kinase 2 (CDK2)**.
> The goal is to demonstrate an end-to-end CADD workflow with *honest method validation* —
> **not** to claim discovery of novel inhibitors.

---

## TL;DR

- **Target:** CDK2, ATP-competitive site. Crystal structure **PDB `2R3I`** (1.28 Å, monomeric).
- **Pipeline:** rigid-receptor docking (**AutoDock Vina**) → physical filtering by short MD
  (**OpenMM**) → endpoint free-energy rescoring (**MM-GBSA**, AmberTools) → ranking & analysis.
- **Validation:** enrichment of known CDK2 actives (ChEMBL) against property-matched decoys;
  redocking of the co-crystal ligand (SCF) as a positive control.
- **Compute:** a single NVIDIA RTX 4090 (MD) + CPU (docking, MM-GBSA).

The scientific question is deliberately modest and *checkable*: **how much does each added layer
of physics (MD pose filtering, then MM-GBSA) change the ranking, and does it recover known
actives better than docking alone?**

---

## Why CDK2

CDK2 is a textbook structure-based design target and a good choice for a clean methods demo:

- **Well-defined, deep ATP pocket** at the kinase hinge — docking behaves well here, unlike
  shallow/solvent-exposed sites.
- **Non-covalent, ATP-competitive inhibitors dominate** — this matches the assumptions of
  docking and MM-GBSA (no covalent-warhead confounder).
- **Abundant, consistent public bioactivity data** (ChEMBL) for enrichment validation.
- **Monomeric ~300-residue system** — small enough for short explicit-solvent MD on one GPU.

Structure `2R3I` was chosen for its ultra-high resolution (1.28 Å) and a drug-like
pyrazolo[1,5-a]pyrimidine inhibitor (`SCF`) suitable as a redocking control. (An oxidized
cysteine `CSD`, a crystallization artifact, is reverted to `CYS` during preparation.)

---

## Pipeline

```
PDB 2R3I ─┐
          ├─▶ 0. Receptor prep (clean chain A, CSD→CYS, protonate, define box)
ChEMBL  ──┘        │
                   ▼
          1. Ligand library (actives + decoys → 3D → pdbqt)
                   │
                   ▼
          2. Docking (AutoDock Vina)  ──▶  enrichment / ROC vs known actives
                   │   top-N hits + redock control
                   ▼
          3. Short MD (OpenMM, explicit solvent)  ──▶  RMSD stability filter
                   │   stable poses only
                   ▼
          4. MM-GBSA rescoring (AmberTools MMPBSA.py) + per-residue decomposition
                   │
                   ▼
          5. Analysis: docking vs MM-GBSA ranking, interaction fingerprints, final shortlist
```

---

## Repository layout

```
cadd-cdk2-pipeline/
├── README.md
├── LICENSE
├── env/                 # conda environment specification
├── config/              # target & run configuration (PDB id, box, MD params)
├── scripts/             # one script per pipeline stage (00_… … 07_…)
├── slurm/               # SLURM submission scripts for the GPU/MD steps
├── results/             # small result tables and figures (committed)
└── notebooks/           # summary notebook (figures + conclusions)
```

Large artifacts (MD trajectories, checkpoints, full ligand sets) are **not** committed; see
`.gitignore`. Every result figure/table is regenerable from the scripts.

---

## Reproduce

> Full commands are filled in as each stage lands. High level:

```bash
# 1. environment (Linux, conda/mamba)
conda env create -f env/cadd_env.yml   # or: mamba env create -f env/cadd_env.yml
conda activate cadd

# 2. run the pipeline (scripts/00_… onward)
```

---

## Results

*Populated as the pipeline runs.* Will include:

- Docking-score distribution and enrichment (ROC-AUC / EF) of known actives vs decoys.
- Per-pose RMSD stability over the short MD (which docked poses survive).
- Docking vs MM-GBSA ranking comparison and the final top-N shortlist with ΔG_bind.
- Per-residue MM-GBSA decomposition and MD interaction fingerprints at the hinge.

---

## Limitations & next steps

- Endpoint MM-GBSA is an approximate free-energy method; it improves on docking but is **not**
  a rigorous binding free energy. The natural next step is **relative binding FEP** on a
  congeneric subset (CDK2 is a classic FEP benchmark) — scoped out here for compute-time reasons.
- Rigid-receptor docking; protein flexibility is only sampled in the post-docking MD.
- No experimental validation — predictions are relative rankings, not measured affinities.

---

## Tools & environment

RDKit · AutoDock Vina · Meeko · OpenBabel · PDBFixer · AmberTools (antechamber, tleap,
MMPBSA.py, cpptraj) · OpenMM · OpenFF/openmmforcefields (GAFF) · MDAnalysis · scikit-learn ·
matplotlib. Hardware: NVIDIA RTX 4090.

## License

MIT — see [LICENSE](LICENSE).

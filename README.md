# CDK2 Structure-Based Virtual Screening — Docking → MD → MM-GBSA

> **Self-initiated practice project.** A complete, scriptable, reproducible computational
> pipeline for structure-based virtual screening against **cyclin-dependent kinase 2 (CDK2)**.
> The goal is to demonstrate an end-to-end CADD workflow with *honest method validation* —
> **not** to claim discovery of novel inhibitors.

---

## TL;DR

- **Target:** CDK2, ATP-competitive site. Crystal structure **PDB `2R3I`** (1.28 Å, monomeric).
- **Library:** 300 ligands — 60 known CDK2 actives + 240 property-matched decoys (DUD-E; a live
  ChEMBL path is also provided).
- **Pipeline:** rigid-receptor docking (**AutoDock Vina**) → physical filtering by short MD
  (**OpenMM**, explicit solvent, RTX 4090) → endpoint free-energy rescoring (**MM-GBSA**,
  AmberTools) → ranking & per-residue analysis.
- **Docking enrichment:** ROC-AUC **0.61**, EF₁% **3.3** — modest but real, and improved by the
  physics-based rescoring of the top hits.

The scientific question is deliberately modest and *checkable*: **how much does each added layer
of physics (MD pose filtering, then MM-GBSA) change the ranking of the top docking hits, and how
well does docking alone recover known actives from property-matched decoys?**

---

## Why CDK2

CDK2 is a textbook structure-based design target and a clean choice for a methods demo:

- **Well-defined, deep ATP pocket** at the kinase hinge — docking behaves well here, unlike
  shallow / solvent-exposed sites.
- **Non-covalent, ATP-competitive inhibitors dominate** — matching the assumptions of docking
  and MM-GBSA (no covalent-warhead confounder, unlike e.g. SARS-CoV-2 Mpro).
- **Abundant, consistent public bioactivity data** for enrichment validation.
- **Monomeric ~300-residue system** — small enough for short explicit-solvent MD on one GPU.

Structure `2R3I` was chosen for its ultra-high resolution (1.28 Å) and a drug-like
pyrazolo[1,5-a]pyrimidine inhibitor (`SCF`) at the ATP site. An oxidized cysteine `CSD177`
(a crystallization artifact) is reverted to `CYS`, and a disordered β3/αC loop (res 46–52,
distal to the pocket) is rebuilt so the MD topology has no chain break.

---

## Pipeline

```
PDB 2R3I ─┐
          ├─▶ 0. Receptor prep (clean chain A, CSD→CYS, rebuild loop, protonate, define box)
DUD-E   ──┘        │
ChEMBL ───────▶ 1. Ligand library (60 actives + 240 matched decoys → 3D → PDBQT)
                   │
                   ▼
          2. Docking (AutoDock Vina)  ──▶  ROC / enrichment vs known actives + native redock control
                   │   top-8 hits
                   ▼
          3. Short MD (OpenMM, explicit solvent, 4 ns)  ──▶  ligand-RMSD pose-stability filter
                   │
                   ▼
          4. MM-GBSA rescoring (AmberTools MMPBSA.py) + per-residue decomposition
                   │
                   ▼
          5. Final analysis: docking vs MM-GBSA ranking, pose stability, binding hotspots
```

---

## Results

### Docking enrichment (300 ligands: 60 actives / 240 matched decoys)

| Metric | Value |
|--------|-------|
| ROC-AUC | **0.613** |
| Enrichment factor @1% | **3.33** |
| Enrichment factor @5% | 2.33 |
| Enrichment factor @10% | 1.50 |
| Active mean score | −8.76 kcal/mol |
| Decoy mean score | −8.34 kcal/mol |
| Native SCF redock (control) | −8.98 kcal/mol (ranks top 31%) |

![ROC](results/roc_curve.png) ![score distribution](results/score_distribution.png)

Single-conformation rigid docking gives **modest but real** enrichment (AUC 0.61, with the
strongest early enrichment at 1%). The top docking hits are a mix of actives and decoys —
exactly the situation that motivates re-scoring the top hits with more physics.

*Note on the redock control:* the native SCF re-docks with a strong score (−8.98 kcal/mol);
a pose-RMSD-to-crystal was not computed because the crystal ligand carries dual alternate
conformations with no connectivity records (ambiguous bond perception), so the score + rank is
used as the positive control instead.

### MD pose stability + MM-GBSA rescoring (top 8 docking hits)

The 8 top-Vina hits (4 actives, 4 decoys) were each run through 4 ns of explicit-solvent MD and
rescored with MM-GBSA. **5 of 8 poses stayed stable; 3 drifted** (ligand RMSD-to-start grew past
3 Å).

| Ligand | Class | Vina | MM-GBSA ΔG (kcal/mol) | Pose stable? |
|--------|-------|------|-----------------------|--------------|
| CHEMBL148580  | active | −10.29 | **−53.1 ± 5.3** | ✅ |
| C14253253     | decoy  | −10.43 | −44.7 ± 3.3 | ❌ drift (→4.0 Å) |
| CHEMBL210765  | active | −10.24 | −43.8 ± 4.3 | ✅ |
| C12376943     | decoy  | −10.65 | −43.7 ± 3.5 | ✅ |
| C24583851     | decoy  | −10.24 | −40.3 ± 3.9 | ❌ drift |
| CHEMBL183950  | active | −10.85 | −40.2 ± 6.2 | ❌ drift (→5.1 Å) |
| C27794555     | decoy  | −10.28 | −37.8 ± 2.7 | ✅ |
| CHEMBL1094820 | active | −10.50 | −37.4 ± 3.8 | ✅ |

![docking vs MM-GBSA](results/docking_vs_mmgbsa.png) ![pose stability](results/rmsd_stability.png)

**What each physics layer adds:**
- **Re-ranking:** across these eight already-well-docked hits the Vina score (a narrow −10.2…
  −10.9 band) carries little ranking information; MM-GBSA spreads them over ~16 kcal/mol and puts
  a genuine active (CHEMBL148580, −53 kcal/mol) at the top.
- **False-positive filtering:** the pose MM-GBSA ranked #2 (decoy C14253253) **drifts out of the
  pocket during MD** (RMSD → 4 Å) and would be discarded — a high-scoring artifact caught only by
  dynamics. The #1 *docking* hit (active CHEMBL183950) also drifts (→5 Å), a reminder that a top
  docking score does not guarantee a stable binding mode.
- **Mechanism:** per-residue MM-GBSA decomposition localises binding to the CDK2 hinge
  (Gln131 / Asn132 / Leu134), the catalytic Lys33 and the Gly-rich/β-sheet residues lining the
  ATP cleft — consistent with known CDK2 pharmacophores.

![MM-GBSA hotspots](results/mmgbsa_decomposition.png)

*Honest reading:* with only 4 actives / 4 decoys advanced and uncalibrated absolute MM-GBSA
energies, this is a methods demonstration rather than a benchmark — but each layer behaves as
theory predicts, and the combined docking → MD → MM-GBSA funnel ends on a real active.

---

## Reproduce

```bash
# environment (Linux + conda/mamba; CUDA build of OpenMM must match the GPU driver)
conda env create -f env/cadd_env.yml && conda activate cadd

WD=$PWD
python scripts/00_prepare_receptor.py --workdir $WD --pdb 2R3I --box-min 22.5
python scripts/01_fetch_ligands.py    --workdir $WD --source dude     # or --source chembl
python scripts/02_dock.py             --workdir $WD --exhaustiveness 8 --cpus 16
python scripts/03_analyze_docking.py  --workdir $WD

# MD + MM-GBSA for the top hits (needs a GPU; run under SLURM)
TOPN=8 NS=4 sbatch slurm/md_pipeline.slurm
python scripts/07_final_analysis.py   --workdir $WD
```

---

## Engineering notes (real problems solved)

A faithful record of the non-trivial issues hit while making this run end-to-end:

- **ChEMBL API outage / fragile property queries.** The live ChEMBL REST API returned 500s on
  float property-range filters and went down mid-run. Switched the default ligand source to the
  static **DUD-E CDK2** benchmark (robust, recognized) while keeping the ChEMBL path as an
  alternative (`--source chembl`).
- **OpenMM CUDA `UNSUPPORTED_PTX_VERSION`.** The conda OpenMM shipped CUDA-12.9 nvrtc while the
  node driver supports CUDA 12.6; pinned `cuda-version=12.6` to realign the JIT PTX.
- **Bond perception from docked PDBQT.** OpenBabel mis-perceived bond orders from the docked pose
  (spurious pentavalent N), crashing ligand parametrization. Reconstruct the ligand with **Meeko**
  (`RDKitMolCreate.from_pdbqt_mol`) directly from the SMILES-annotated PDBQT instead.
- **MMPBSA per-residue selection.** `print_res="within 5"` is rejected by this MMPBSA.py build;
  the pocket residues are now computed explicitly (ligand neighbours) and passed as an integer
  list, with a graceful fall-back to dG-only if decomposition fails.
- **`set -e` swallowed in a `|| ` context.** A `( set -e; … ) || echo` SLURM loop body silently
  disabled `set -e`, so failed setups still ran MD/MM-GBSA; replaced with explicit `&&` chaining.

---

## Limitations & next steps

- Endpoint MM-GBSA is an approximate free-energy method; absolute values are not calibrated and
  only **relative** rankings are meaningful. The natural next step is **relative binding FEP** on
  a congeneric subset (CDK2 is a classic FEP benchmark) — scoped out here for compute-time reasons.
- Rigid-receptor docking; protein flexibility is only sampled in the post-docking MD.
- No experimental validation — predictions are relative rankings, not measured affinities.
- Modest library size (300) and a single MD replica per hit; production work would scale both.

---

## Tools & environment

RDKit · AutoDock Vina · Meeko · OpenBabel · PDBFixer · AmberTools (antechamber, tleap,
MMPBSA.py, cpptraj) · OpenMM (CUDA) · MDAnalysis/MDTraj · scikit-learn · matplotlib.
Hardware: NVIDIA RTX 4090.

## License

MIT — see [LICENSE](LICENSE).

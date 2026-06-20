#!/bin/bash
# MD + MM-GBSA pipeline for the top docking hits on kepler8 (1x RTX 4090, via SLURM).
# Usage:  TOPN=8 NS=4 sbatch slurm/md_pipeline.slurm
#   TOPN  number of top-Vina-score hits to advance (default 8)
#   NS    production MD length in ns per hit (default 4)
#   WD    working dir (default $HOME/cadd_cdk2)
#SBATCH -J cdk2_md
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=28000M
#SBATCH -t 1-00:00:00
#SBATCH -o cdk2_md.%j.out

set -u
WD=${WD:-$HOME/cadd_cdk2}
TOPN=${TOPN:-8}
NS=${NS:-4}

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate cadd          # sets AMBERHOME for antechamber/tleap/MMPBSA.py
cd "$WD"
echo "host=$(hostname) gpu=$CUDA_VISIBLE_DEVICES TOPN=$TOPN NS=$NS WD=$WD"

# top-N hits by Vina score (col6), status==ok (col7); most negative = best
HITS=$(tr -d '\r' < docking/docking_scores.csv \
       | awk -F, 'NR>1 && $7=="ok" && $6!="" {print $6","$1}' \
       | sort -t, -k1,1g | head -n "$TOPN" | cut -d, -f2)
echo "top hits: $HITS"

for cid in $HITS; do
  echo "######## $cid ########"
  (
    set -e
    python scripts/04_setup_complex.py --workdir "$WD" --cid "$cid"
    python scripts/05_run_md.py        --workdir "$WD" --cid "$cid" --ns "$NS"
    python scripts/06_mmgbsa.py        --workdir "$WD" --cid "$cid"
  ) || echo "!!!! $cid FAILED — continuing with next hit"
done
echo "MD + MM-GBSA pipeline complete"

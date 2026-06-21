#!/bin/bash
# Relative binding FEP (openfe quickrun) for selected transformations on kepler8 (1x RTX 4090).
# Each "transformation" is one leg (complex or solvent) of one edge; an edge = 2 transformations.
# Usage:  TRANSFORMS="rbfe_ligA_ligB_complex rbfe_ligA_ligB_solvent" sbatch slurm/fep_run.slurm
#SBATCH -J cdk2_fep
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=28000M
#SBATCH -t 2-00:00:00
#SBATCH -o cdk2_fep.%j.out

set -u
WD=${WD:-$HOME/cadd_cdk2/fep}
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate fep
cd "$WD"
mkdir -p results work
echo "host=$(hostname) gpu=$CUDA_VISIBLE_DEVICES TRANSFORMS=$TRANSFORMS"

for t in $TRANSFORMS; do
  echo "######## $t ########"
  date
  openfe quickrun "network_setup/transformations/${t}.json" \
        -o "results/${t}.json" -d "work/${t}" \
    && echo "OK $t" \
    || echo "!!!! $t FAILED"
done
echo "FEP run complete"

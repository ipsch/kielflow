#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=6-00:00:00

make bin/kielflow
./bin/kielflow
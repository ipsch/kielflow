#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=6-00:00:00
#SBATCH --mail-type=END,FAIL                          # notifications for job done & fail
#SBATCH --mail-user=schnell@theo-physik.uni-kiel.de   # send-to address

make bin/kielflow
./bin/kielflow
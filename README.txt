## Required module
module add R/3.6.1-gcb03
module add gcc/gcc/7.1.0-fasrc01

## Previously: module add gcc/7.3.0-gcb01

## fido setup (for labraduck)
# to build .tar.gz
R CMD build fido --no-build-vignettes
# to install locally
install.packages("fido", repos = NULL, type = "source")

## Run job with multiple cores
sbatch --cpus-per-task=4 job.slurm

# ASV-level jobs have been attempted but none have finished yet.
# These take at least 32GB RAM and will probably take hours.

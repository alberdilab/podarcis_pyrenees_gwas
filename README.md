# Podarcis Pyrenees GWAS
GWAS analysis of Podarcis lizards in an environmental transect in the Pyrenees

## Create plink bed files from VCF

The GWAS analysis tool we will employ require genotype information to be in bed format. This can be achieved using plink2. Plink also enables us to include sex information and to assign unique codes to all variants, which is necessary for downstream analyses.

```sh
#Create the batch file
cat <<EOF > 1_input.sh
#!/bin/bash
#SBATCH --job-name=1_input
#SBATCH --nodes=1
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=128gb
#SBATCH --time=1:00:00

# Load required modules in the server
module load openblas/0.3.24 plink/2.0.0

# Run format conversion
plink2 --vcf all.filtered.vcf --max-alleles 2 --set-all-var-ids @_# --allow-extra-chr --make-bed --out all.filtered2
EOF

#Launch the batch file
sbatch 1_input.sh
```

## Run GWAS with GridLMM

[GridLMM](https://github.com/deruncie/GridLMM) is a package for fitting linear mixed models (LMMs) with multiple random effects, which is useful for GWAS analyses because it contains a fitting process optimised for repeated evaluation of the random effect model with different sets of fixed effects. As we have (slightly) different dominance values for the different time points, we run 7 different GWAS analyses associating dominance indices in each time point with the mouse genotypes. This approach provides an overview of the signal variation through time. We also include a GWAS analysis with the mean dominance value of each mouse. The R code required for this analysis can be found in 'code/gwas.R'. Output data is stored as a Rdata file in the 'results' directory.

```sh
#Create the batch file
cat <<EOF > 2_gwas.sh
#!/bin/bash
#SBATCH --job-name=2_gwas
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=500gb
#SBATCH --time=24:00:00

# Load R (with all dependencies already installed)
module load gcc R/4.2.1

# Run R scripts for mean dominance and time-specific dominance metrics
Rscript code/gwas.R -i "all.filtered2" -m "metadata.tsv" -o "gwas_results.Rdata"
EOF

#Launch the batch file
sbatch 2_gwas.sh
```

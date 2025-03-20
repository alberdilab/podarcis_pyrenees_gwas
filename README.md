# Podarcis Pyrenees GWAS
GWAS analysis of Podarcis lizards in an environmental transect in the Pyrenees

## Run GWAS using plink

The GWAS analysis tool we will employ require genotype information to be in bed format. This can be achieved using plink2. Plink also enables us to include sex information and to assign unique codes to all variants, which is necessary for downstream analyses.

```sh
#Create the batch file
cat <<EOF > gwas.sh
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
plink2 --vcf all.filtered.vcf --pheno metadata.tsv --glm allow-no-covars --max-alleles 2 --set-all-var-ids @_# --allow-extra-chr --make-bed --out all.filtered2
EOF

#Launch the batch file
sbatch gwas.sh
```

## Filter and summarise GWAS

```py
import pandas as pd

# Load PLINK --glm results
file_path = "all.filtered2.elevation.glm.linear"  # Update with actual file path
df = pd.read_csv(file_path, sep="\t")

# Convert numerical columns (handling 'NA' as missing values)
numeric_cols = ["OBS_CT", "BETA", "SE", "T_STAT", "P"]
df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")

# Step 1: Remove Variants with Low Sample Counts
min_sample_count = 5  # Change as needed
df_filtered = df[df["OBS_CT"] >= min_sample_count]

# Step 2: Remove Variants with 'NA' BETA (usually constant or monomorphic alleles)
df_filtered = df_filtered.dropna(subset=["BETA"])

# Step 3: Identify Significant Variants (p < 0.05)
significance_threshold = 0.05  # Adjust for multiple testing correction if needed
df_significant = df_filtered[df_filtered["P"] < significance_threshold]

# Step 4: Highlight Variants with Large Effect Sizes (absolute BETA > threshold)
large_effect_threshold = 50  # Change as needed
df_large_effect = df_filtered[abs(df_filtered["BETA"]) > large_effect_threshold]

# Display the filtered data
import ace_tools as tools

tools.display_dataframe_to_user(name="Filtered PLINK --glm Results", dataframe=df_filtered)

# Save filtered results
df_filtered.to_csv("filtered_glm_results.txt", sep="\t", index=False)
```

## Manhattan plot
```py
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load PLINK --glm results
file_path = "all.filtered2.elevation.glm.linea"  # Update with actual file path
df = pd.read_csv(file_path, delim_whitespace=True)

# Convert numeric columns and handle missing values
df["P"] = pd.to_numeric(df["P"], errors="coerce")
df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

# Drop rows with missing P-values
df = df.dropna(subset=["P"])

# Convert chromosome names if necessary (e.g., removing "NC_" prefix)
df["CHROM"] = df["#CHROM"].astype(str).str.replace("NC_", "").str.split(".").str[0]
df["CHROM"] = pd.to_numeric(df["CHROM"], errors="coerce")

# Sort dataframe by chromosome and position
df = df.sort_values(by=["CHROM", "POS"])

# Assign index for plotting
df["index"] = range(len(df))

# Define colors for chromosomes
chromosomes = df["CHROM"].unique()
colors = ["blue", "red"] * (len(chromosomes) // 2 + 1)

# Create Manhattan plot
plt.figure(figsize=(12, 6))

for i, chrom in enumerate(chromosomes):
    subset = df[df["CHROM"] == chrom]
    plt.scatter(subset["index"], -np.log10(subset["P"]), 
                color=colors[i % 2], s=10)

# Add genome-wide significance threshold line (Bonferroni correction)
significance_threshold = 0.05 / len(df)
plt.axhline(-np.log10(significance_threshold), color='black', linestyle='dashed', linewidth=1)

# Labels and formatting
plt.xlabel("")
plt.ylabel("-log10(P-value)")
plt.title("Manhattan Plot of GWAS Results")

# Remove legend and X-axis ticks for large datasets
plt.xticks([])  # Hide X-axis values
plt.legend().remove()  # Remove legend

# Save the figure to a file
output_file = "manhattan_plot.png"
plt.savefig(output_file, dpi=300, bbox_inches="tight")
```

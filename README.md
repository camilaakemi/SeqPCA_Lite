# SeqPCA-Lite: A Simplified Sequence PCA and Clustering Pipeline

**SeqPCA-Lite** is a lightweight, Python-based pipeline designed for the rapid filtering and sequence-space exploration of enzyme candidates (e.g., esterases, α/β-hydrolases) mined from large databases like UniProt/TrEMBL.

Instead of relying on computationally intensive multi-sequence alignments and complex evolutionary models, this pipeline employs a straightforward **One-Hot Encoding** followed by **Principal Component Analysis (PCA)** and **K-means clustering**. It is designed to rapidly identify canonical structural representatives from large sequence datasets using straightforward structural and physicochemical rules.

## ⚠️ Scope and Limitations
It is important to note the boundaries of this tool:
* **Linear Correlation Only:** This pipeline uses a One-Hot Encoding + PCA approach designed for simplified, fast sequence correlation and clustering. It is a linear mathematical model.
* **No Evolutionary Covariance:** It **does not** account for complex evolutionary covariance, epistasis, or advanced conservation scoring metrics (e.g., `scorecons`). It simply treats each amino acid at each position as an independent categorical feature.
* **Best Use Case:** This tool is best used as a "first-pass" filter to reduce redundancy and group thousands of distant homologous sequences into manageable, representative sub-families before moving on to rigorous phylogenetic analysis, 3D modeling, or molecular dynamics.

## Requirements
It is recommended to use `mamba` or `conda` to create an isolated environment to run these scripts:

```bash
mamba create -n seqpca-lite python=3.11 biopython pandas numpy scikit-learn seaborn matplotlib plotly -y
mamba activate seqpca-lite
```

## Available Scripts
1. filter_motif.py
Filters a combined FASTA file to retain only sequences containing a specific functional motif. By default, it searches for G[A-Z]S[A-Z]G, a classic nucleophilic elbow motif found in α/β-hydrolases.

Usage:

```bash

python filter_motif.py -i input.fasta -o filtered_motif.fasta -m "G[A-Z]S[A-Z]G"
```

2. filter_pI_length.py
Performs a statistical and physicochemical refinement to ensure candidates are viable for specific environmental applications (e.g., soil bioremediation).

Length Filter: Analyzes sequence length distribution via the Shapiro-Wilk test. If normal, it selects sequences within ±1 standard deviation. If non-normal, it uses the 15th-85th percentiles to exclude truncated fragments or massive fusions.

pI Filter: Discards sequences with an Isoelectric Point above a defined threshold. This is particularly useful for avoiding enzymes that would immobilize by adsorbing to negatively charged clay particles in soil environments.

Outputs colorblind-friendly histograms and scatter plots for selection justification.

Usage:

```bash
# Note: Adjust the --max_pi parameter according to your environmental target.
# Example: Use 5.5 for highly acidic enzymes suitable for soil, or change to 7.0/8.0 for other contexts.
python filter_pI_length.py -i filtered_motif.fasta --save_fasta candidates_stats.fasta --max_pi 5.5
```

3. pca_clustering_pipeline.py
The core unified script of the SeqPCA-Lite pipeline. It executes the sequence-space calculation in a single run:

One-Hot Encoding: Converts the alignment into a multi-dimensional binary matrix, without assuming any prior relationship between amino acids.

PCA Calculation: Calculates Principal Components and feature loadings.

K-Means Clustering: Groups sequences mathematically and identifies the centroid (canonical representative) of each cluster.

Signature Identification: Calculates which specific amino acid variations drive the formation of each super-cluster.

Visualization: Outputs an interactive .html Plotly 2D scatter plot of the PC space.

Usage:

```bash
# Ensure your sequences are aligned first (e.g., using MAFFT).
# Note: Adjust the -k parameter to the number of super-clusters you expect or wish to evaluate.
python pca_clustering_pipeline.py -i alignment.aln.fasta -k 3 -o pca_results
```

Outputs
pca_results_clusters.csv: Sequence mapping to their respective clusters.

pca_results_signatures.csv: The list of amino acid defining features (loadings x cluster centers).

pca_results_centroids.fasta: The FASTA files containing only the canonical models.

pca_results_plot.html: Interactive 2D scatter plot for visual inspection.


Qualquer outra dúvida para estruturar o seu código no GitHub, é só avisar!

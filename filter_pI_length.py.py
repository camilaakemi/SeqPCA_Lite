import argparse
import sys
from pathlib import Path
from typing import Tuple, Dict

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# --- Configuration: Colorblind Palette ---
# Manually defining accessible colors to ensure consistency across plots.
# These hex codes are based on the Seaborn/Matplotlib "colorblind" palette.
CB_GREY = "#7f7f7f"    # Grey (Histogram Bars / Discarded points)
CB_BLUE = "#377eb8"    # Blue (Mean line)
CB_ORANGE = "#ff7f00"  # Orange (Selected points highlight)
CB_GREEN = "#4daf4a"   # Green (Selected Region shading)
CB_PURPLE = "#984ea3"  # Purple (Median line)
CB_RED = "#e41a1c"     # Red (Cutoff lines)
CB_DARK = ".15"        # Dark Grey (KDE line)

def parse_fasta_data(file_path: Path) -> Tuple[pd.DataFrame, Dict[str, SeqRecord]]:
    """
    Parses a FASTA file, calculates physicochemical properties (pI/Length), 
    and stores original BioPython records for later retrieval.
    
    Args:
        file_path (Path): Path to the input FASTA file.
        
    Returns:
        pd.DataFrame: DataFrame containing 'ID', 'Length', and 'pI'.
        Dict[str, SeqRecord]: Dictionary mapping IDs to original SeqRecord objects.
    """
    print(f"[*] Loading sequences from: {file_path}")
    data = []
    record_dict = {}

    try:
        for record in SeqIO.parse(file_path, "fasta"):
            # Clean sequence: remove stop codons, ensure uppercase
            seq_str = str(record.seq).upper().replace("*", "")
            
            # Skip empty sequences
            if len(seq_str) == 0: 
                continue

            try:
                analyser = ProteinAnalysis(seq_str)
                # Store numerical data for analysis
                data.append({
                    "ID": record.id,
                    "Length": len(seq_str),
                    "pI": analyser.isoelectric_point()
                })
                # Store original object for saving
                record_dict[record.id] = record
            except ValueError:
                # Skip sequences with ambiguous amino acids (X, B, Z, J)
                continue

    except FileNotFoundError:
        print(f"[x] Error: File '{file_path}' not found.")
        sys.exit(1)

    if not data:
        print("[x] Error: No valid sequences found in the file.")
        sys.exit(1)

    return pd.DataFrame(data), record_dict

def calculate_statistics(df: pd.DataFrame) -> Tuple[float, float, str]:
    """
    Performs Shapiro-Wilk normality test and calculates dynamic length cutoffs.
    
    Returns:
        lower (float): Lower length cutoff.
        upper (float): Upper length cutoff.
        method_name (str): Description of the statistical method used.
    """
    # Shapiro-Wilk Test
    # H0: Distribution is Normal. If p > 0.05, we cannot reject H0.
    shapiro_stat, shapiro_p = stats.shapiro(df['Length'])
    is_normal = shapiro_p > 0.05
    
    lower, upper = 0.0, 0.0
    method_name = ""

    if is_normal:
        method_name = "Normal Distribution (Mean +/- 1 Sigma)"
        mean = df['Length'].mean()
        std = df['Length'].std()
        lower = mean - std
        upper = mean + std
        print(f"[*] Stats: Distribution is NORMAL (p={shapiro_p:.4f}). Using Parametric approach.")
    else:
        method_name = "Non-Normal Distribution (Percentiles 15%-85%)"
        lower = df['Length'].quantile(0.15)
        upper = df['Length'].quantile(0.85)
        print(f"[*] Stats: Distribution is NON-NORMAL (p={shapiro_p:.4f}). Using Robust approach.")

    print(f"[*] Cutoffs calculated: {lower:.1f} - {upper:.1f} aa")
    return lower, upper, method_name

def plot_histogram_justification(df: pd.DataFrame, lower: float, upper: float, method_name: str, output_path: Path):
    """
    Generates a Histogram with Mean, Median, and Cutoff regions using a Colorblind palette.
    This plot serves as justification for the chosen statistical method.
    """
    sns.set_style("whitegrid")
    plt.figure(figsize=(10, 6))

    # 1. Histogram + KDE (Kernel Density Estimate)
    sns.histplot(
        df['Length'], 
        kde=True, 
        color=CB_GREY, 
        alpha=0.5,
        line_kws={'color': CB_DARK, 'linewidth': 2},
        label='Sequences'
    )

    # 2. Statistical Indicators
    mean_val = df['Length'].mean()
    median_val = df['Length'].median()

    plt.axvline(mean_val, color=CB_BLUE, linestyle='--', linewidth=2, label=f'Mean ({mean_val:.0f})')
    plt.axvline(median_val, color=CB_PURPLE, linestyle=':', linewidth=3, label=f'Median ({median_val:.0f})')

    # 3. Selection Region
    plt.axvline(lower, color=CB_RED, linestyle='-', linewidth=2, label='Cutoffs')
    plt.axvline(upper, color=CB_RED, linestyle='-', linewidth=2)
    # Shaded region indicating selected range
    plt.axvspan(lower, upper, color=CB_GREEN, alpha=0.2, label='Selected Region')

    # Formatting
    plt.title(f"Length Distribution & Selection Criteria\n{method_name}", fontsize=14)
    plt.xlabel("Sequence Length (aa)", fontsize=12)
    plt.ylabel("Count", fontsize=12)
    plt.legend(loc='upper right')
    plt.tight_layout()

    plt.savefig(output_path, dpi=300)
    print(f"[*] Histogram saved to '{output_path}'")

def plot_scatter_selection(df: pd.DataFrame, method_name: str, output_path: Path):
    """
    Generates a pI vs Length scatter plot using a Colorblind palette.
    Highlights selected candidates in Orange over a Grey background.
    """
    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 8))

    # Map boolean status to Colorblind colors
    palette = {False: CB_GREY, True: CB_ORANGE}
    
    sns.scatterplot(
        data=df,
        x="Length",
        y="pI",
        hue="Selected",
        palette=palette,
        hue_order=[False, True], # Ensures 'False' (Grey) is plotted first as background
        s=15,                    # Small marker size to reduce overplotting
        alpha=0.4,               # Transparency to show density
        edgecolor=None
    )

    # Custom Legend
    handles, _ = plt.gca().get_legend_handles_labels()
    plt.legend(handles, ["Not selected", "Selected"], title=None, loc='upper right')

    plt.title(f"Candidate Selection: pI vs. Length", fontsize=15)
    plt.xlabel("Sequence Length (aa)", fontsize=13)
    plt.ylabel("Isoelectric Point (pI)", fontsize=13)
    plt.tight_layout()

    plt.savefig(output_path, dpi=300)
    print(f"[*] Scatter plot saved to '{output_path}'")

def main():
    # Argument Parsing
    parser = argparse.ArgumentParser(
        description="Comprehensive Statistical Analysis & Filtering Pipeline (Colorblind Friendly)."
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to input FASTA file.")
    
    # Output file arguments
    parser.add_argument("--out_hist", type=Path, default=Path("plot_distribution_cb.png"), help="Filename for Histogram plot.")
    parser.add_argument("--out_scatter", type=Path, default=Path("plot_selection_cb.png"), help="Filename for Scatter plot.")
    parser.add_argument("--save_fasta", type=Path, default=None, help="Filename to save the filtered FASTA sequences.")
    
    # Filter parameter
    parser.add_argument("--max_pi", type=float, default=5.5, help="Maximum pI threshold (default: 5.5).")
    
    args = parser.parse_args()

    # 1. Load Data
    df, record_dict = parse_fasta_data(args.input)

    # 2. Calculate Stats & Cutoffs
    low, high, method = calculate_statistics(df)

    # 3. Apply Filters (Vectorized)
    mask_length = df['Length'].between(low, high)
    mask_pi = df['pI'] < args.max_pi
    df['Selected'] = mask_length & mask_pi

    # 4. Generate Plots
    plot_histogram_justification(df, low, high, method, args.out_hist)
    plot_scatter_selection(df, method, args.out_scatter)

    # 5. Save Selected Sequences
    if args.save_fasta:
        # Filter IDs based on the 'Selected' column
        selected_ids = df[df['Selected']]['ID'].tolist()
        # Retrieve original records
        selected_records = [record_dict[id_] for id_ in selected_ids]
        
        if selected_records:
            SeqIO.write(selected_records, args.save_fasta, "fasta")
            print(f"[*] SUCCESS: Saved {len(selected_records)} sequences to '{args.save_fasta}'")
        else:
            print("[!] WARNING: No sequences passed the filters. Nothing was saved.")

if __name__ == "__main__":
    main()
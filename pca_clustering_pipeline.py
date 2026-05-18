#!/usr/bin/env python3
"""
pca_clustering_pipeline.py

Description:
A unified pipeline that:
1. One-hot encodes a multiple sequence alignment (MSA).
2. Performs Principal Component Analysis (PCA).
3. Clusters the sequences in the PC space using K-Means.
4. Identifies canonical representatives (centroids).
5. Calculates defining amino acid signatures for each cluster.
6. Generates an interactive Plotly visualization.
"""

import argparse
import sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from Bio import AlignIO, SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import OneHotEncoder
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

# Standard amino acids for One-Hot Encoding
AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Unified PCA and Clustering pipeline for MSA.")
    parser.add_argument("-i", "--input", required=True, help="Input MSA in FASTA format.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output files.")
    parser.add_argument("-n", "--n_components", type=int, default=10, help="Number of PCA components.")
    parser.add_argument("-k", "--n_clusters", type=int, required=True, help="Number of k-means clusters.")
    return parser.parse_args()

def process_alignment(filename):
    print(f"[*] Loading alignment: {filename}")
    try:
        alignment = AlignIO.read(filename, "fasta")
        records_dict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    except Exception as e:
        print(f"[x] Error reading alignment: {e}", file=sys.stderr)
        sys.exit(1)
        
    seq_ids = [record.id for record in alignment]
    seq_matrix = np.array([list(record.seq) for record in alignment])
    print(f"    -> {seq_matrix.shape[0]} sequences, {seq_matrix.shape[1]} positions.")
    return seq_matrix, seq_ids, records_dict

def run_pca_workflow(seq_matrix, n_components):
    n_sequences, n_positions = seq_matrix.shape
    encoded_matrix = np.zeros((n_sequences, n_positions * len(AMINO_ACIDS)), dtype=np.float32)
    encoder = OneHotEncoder(categories=[AMINO_ACIDS], handle_unknown='ignore', sparse_output=False, dtype=np.float32)
    
    print("[*] Performing One-Hot Encoding...")
    for i in range(n_positions):
        column = seq_matrix[:, i].reshape(-1, 1)
        start_idx = i * len(AMINO_ACIDS)
        end_idx = (i + 1) * len(AMINO_ACIDS)
        encoded_matrix[:, start_idx:end_idx] = encoder.fit_transform(column)
        
    print(f"[*] Running PCA (n_components={n_components})...")
    pca = PCA(n_components=min(n_components, n_sequences, encoded_matrix.shape[1]))
    scores = pca.fit_transform(encoded_matrix)
    
    # Format PC columns and Loadings
    pc_cols = [f"PC{i+1}" for i in range(pca.n_components_)]
    df_scores = pd.DataFrame(scores, columns=pc_cols)
    
    positions, amino_acids = [], []
    for i in range(n_positions):
        for aa in AMINO_ACIDS:
            positions.append(i + 1)
            amino_acids.append(aa)
            
    df_loadings = pd.DataFrame(pca.components_.T, columns=pc_cols)
    df_loadings.insert(0, "Position", positions)
    df_loadings.insert(1, "AminoAcid", amino_acids)
    
    return df_scores, df_loadings, pca.explained_variance_ratio_, pc_cols

def run_clustering_and_signatures(df_scores, df_loadings, pc_cols, n_clusters):
    print(f"[*] Running K-Means (k={n_clusters})...")
    X = df_scores[pc_cols].values
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    df_scores['Cluster'] = kmeans.fit_predict(X)
    
    closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, X)
    rep_indices = closest.tolist()
    
    # Calculate Signatures
    cluster_centers = df_scores.groupby('Cluster')[pc_cols].mean()
    loadings_matrix = df_loadings[pc_cols].values
    centers_matrix = cluster_centers.T.values
    signature_scores = np.dot(loadings_matrix, centers_matrix)
    
    df_signatures = df_loadings.copy()
    for i, cluster_id in enumerate(cluster_centers.index):
        df_signatures[f"Signature_Cluster_{cluster_id}"] = signature_scores[:, i]
        
    return df_scores, df_signatures, rep_indices, cluster_centers

def export_results(args, df_scores, df_signatures, rep_indices, seq_ids, records_dict, variance, pc_cols, cluster_centers):
    # Save Canonical FASTA
    rep_ids = [seq_ids[idx] for idx in rep_indices]
    SeqIO.write([records_dict[rid] for rid in rep_ids], f"{args.output_prefix}_centroids.fasta", "fasta")
    
    print("\n--- Canonical Cluster Representatives (Centroids) ---")
    for i, seq_id in enumerate(rep_ids):
        print(f"Cluster {i} Centroid: {seq_id}")
        
    # Print Top Signatures
    print("\n--- Top 5 Defining Features per Cluster ---")
    for cluster_id in cluster_centers.index:
        col_name = f"Signature_Cluster_{cluster_id}"
        top_features = df_signatures.nlargest(5, col_name)
        print(f"\nCluster {cluster_id}:")
        for _, row in top_features.iterrows():
            print(f"  Pos: {int(row['Position']):<4} AA: {row['AminoAcid']:<2} Score: {row[col_name]:.4f}")
            
    # Plotly Visualization (PC1 vs PC2)
    fig = go.Figure()
    for cluster_id in sorted(df_scores['Cluster'].unique()):
        df_cluster = df_scores[df_scores['Cluster'] == cluster_id]
        fig.add_trace(go.Scatter(
            x=df_cluster['PC1'], y=df_cluster['PC2'], mode='markers',
            name=f'Cluster {cluster_id}', text=[seq_ids[i] for i in df_cluster.index], hoverinfo='text'
        ))
        
    fig.update_layout(
        title="PCA & Clustering Space",
        xaxis_title=f"PC1 ({variance[0]*100:.2f}%)",
        yaxis_title=f"PC2 ({variance[1]*100:.2f}%)",
        template='simple_white'
    )
    fig.write_html(f"{args.output_prefix}_plot.html")
    
    # Save CSVs
    df_scores['SequenceID'] = seq_ids
    df_scores.to_csv(f"{args.output_prefix}_clusters.csv", index=False)
    df_signatures.to_csv(f"{args.output_prefix}_signatures.csv", index=False)
    print(f"\n[*] All results saved with prefix: {args.output_prefix}")

def main():
    args = parse_arguments()
    seq_matrix, seq_ids, records_dict = process_alignment(args.input)
    df_scores, df_loadings, variance, pc_cols = run_pca_workflow(seq_matrix, args.n_components)
    df_scores, df_signatures, rep_indices, cluster_centers = run_clustering_and_signatures(df_scores, df_loadings, pc_cols, args.n_clusters)
    export_results(args, df_scores, df_signatures, rep_indices, seq_ids, records_dict, variance, pc_cols, cluster_centers)

if __name__ == "__main__":
    main()

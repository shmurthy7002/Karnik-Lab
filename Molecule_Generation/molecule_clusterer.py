import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
import numpy as np
from sklearn.cluster import AffinityPropagation
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import csv
import sys
import umap
import os

print(f"RDKit version: {rdkit.__version__}")
 
def get_smiles_from_csv(csv_file):
    smiles_list = []
    with open(os.path.join("csv_output", csv_file)) as file_obj:
        reader_obj = csv.reader(file_obj) 
        
        for row in reader_obj:
            smiles_list.append(row[0])
    return smiles_list

def smiles_to_fingerprints(smiles_list):
    fingerprints = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fingerprints.append(fingerprint)
    return fingerprints

def evaluate_clusters(similarity_matrix, labels):
    distance_matrix = 1 - similarity_matrix
    score = silhouette_score(distance_matrix, labels, metric='precomputed')
    print(f"Silhouette Score: {score}")

def compute_similarity_matrix(fingerprints):
    num_mols = len(fingerprints)
    similarity_matrix = np.zeros((num_mols, num_mols))
    for i in range(num_mols):
        for j in range(num_mols):
            if i <= j:
                similarity = TanimotoSimilarity(fingerprints[i], fingerprints[j])
                similarity_matrix[i, j] = similarity
                similarity_matrix[j, i] = similarity
    return similarity_matrix

def cluster_molecules(similarity_matrix):
    af = AffinityPropagation(affinity='precomputed', random_state=42)
    af.fit(similarity_matrix)
    return af.labels_, af.cluster_centers_indices_

def visualize_clusters(similarity_matrix, labels, perplexity=None, learning_rate=200):
        
        # Visualizing using T-SNE
        n_samples = similarity_matrix.shape[0]
        if perplexity is None or perplexity >= n_samples:
            perplexity = min(30, n_samples - 1)
        tsne = TSNE(n_components=2, metric='precomputed', perplexity=perplexity, learning_rate=learning_rate, init='random', random_state=42)
        tsne_result = tsne.fit_transform(1 - similarity_matrix)
        plt.figure(figsize=(12, 8))
        scatter = plt.scatter(tsne_result[:, 0], tsne_result[:, 1], c=labels, cmap='viridis', alpha=0.7)
        plt.colorbar(scatter)
        plt.title("Molecular Clusters")
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")
        plt.show()

        # Visualizing using UMAP
        umap_model = umap.UMAP(metric='precomputed', random_state=42)
        umap_result = umap_model.fit_transform(1 - similarity_matrix)
        plt.figure(figsize=(12, 8))
        scatterplt = plt.scatter(umap_result[:, 0], umap_result[:, 1], c=labels, cmap='viridis', alpha=0.7)
        plt.colorbar(scatterplt)
        plt.title("Molecular Clusters (UMAP)")
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
        plt.show()

def main(smiles_list):
    fingerprints = smiles_to_fingerprints(smiles_list)
    
    similarity_matrix = compute_similarity_matrix(fingerprints)
    
    labels, cluster_centers = cluster_molecules(similarity_matrix)
    
    print("Number of clusters:", len(cluster_centers))
    for cluster_id in np.unique(labels):
        print(f"Cluster {cluster_id}:")
        for i, label in enumerate(labels):
            if label == cluster_id:
                print(f"  - Molecule {i}: {smiles_list[i]}")

    evaluate_clusters(similarity_matrix, labels)
    
    visualize_clusters(similarity_matrix, labels)

if __name__ == "__main__":
    main(get_smiles_from_csv(sys.argv[1]))

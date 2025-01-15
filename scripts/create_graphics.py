import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from ete3 import Tree, TreeStyle, NodeStyle
import os

os.makedirs("figures", exist_ok=True)
#########################################################
#########################################################
def plot_hmm_evalue_distribution(hmm_results_file):


    plt.hist(hmm_results["E-value"], bins=50, color="skyblue", edgecolor="black")
    plt.xscale("log")
    plt.title("Distribution of HMM Hit Scores (E-values)")
    plt.xlabel("E-value")
    plt.ylabel("Frequency")
    plt.savefig("figures/hmm_evalue_distribution.png")
    plt.close()

#########################################################
#########################################################

def plot_cluster_distributions(tree_clusters_file):
    with open(tree_clusters_file, "r") as f:
        clusters = [line.split()[1] for line in f]

    cluster_counts = pd.Series(clusters).value_counts()

    cluster_counts.plot(kind="bar", figsize=(10, 6), color="orange", edgecolor="black")
    plt.title("Number of Sequences per Cluster")
    plt.xlabel("Cluster")
    plt.ylabel("Number of Sequences")
    plt.tight_layout()
    plt.savefig("figures/cluster_bar_chart.png")

  
    plt.close()

    cluster_counts.plot(kind="hist", bins=30, color="purple", edgecolor="black")
    plt.title("Distribution of Cluster Sizes")
    plt.xlabel("Number of Sequences")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig("figures/cluster_size_distribution.png")

  
    plt.close()




#########################################################
#########################################################
def plot_representative_lengths(clustered_sequences_file):
    lengths = [len(record.seq) for record in SeqIO.parse(clustered_sequences_file, "fasta")]

    plt.hist(lengths, bins=30, color="green", edgecolor="black")
    plt.title("Distribution of Cluster Representative Lengths")
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig("figures/representative_lengths_distribution.png")
    plt.close()





#########################################################
#########################################################
def plot_phylogenetic_tree(phylo_tree_file, tree_clusters_file):
    tree = Tree(phylo_tree_file)
    cluster_assignments = {}

    with open(tree_clusters_file, "r") as f:
        for line in f:
            seq_id, cluster_id = line.strip().split()
            cluster_assignments[seq_id] = cluster_id

    for leaf in tree:
        cluster = cluster_assignments.get(leaf.name, "NA")
        style = NodeStyle()
         #singleton
        if cluster == "-1":
            style["fgcolor"] = "red"
            style["size"] = 8
        else:
            style["fgcolor"] = "blue"
            style["size"] = 10
        leaf.set_style(style)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.scale = 20
    tree.render("figures/phylogenetic_tree_with_clusters.png", tree_style=ts)

def main():
    hmm_results_file = "hmmsearch_results.txt"
    tree_clusters_file = "tree_clusters.txt"
    clustered_sequences_file = "clustered_sequences.fasta"
    phylo_tree_file = "phylogenetic_tree.nwk"

    plot_hmm_evalue_distribution(hmm_results_file)
    plot_cluster_distributions(tree_clusters_file)
    plot_representative_lengths(clustered_sequences_file)
    plot_phylogenetic_tree(phylo_tree_file, tree_clusters_file)


if __name__ == "__main__":
    main()

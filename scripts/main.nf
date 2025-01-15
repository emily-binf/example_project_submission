#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process HMMSearch {
    input:
    path protein_db
    path hmm_model
    output:
    path "hmmsearch_results.txt"

    script:
    """
    hmmsearch --tblout hmmsearch_results.txt $hmm_model $protein_db
    """
}

process FilterHits {
    input:
    path "hmmsearch_results.txt"
    path protein_db
    output:
    path "filtered_hits.fasta"

    script:
    """

    grep -v "^#" hmmsearch_results.txt | awk '\$5 < 1e-3 {print \$1}' > hits_ids.txt

    if [ -s hits_ids.txt ]; then
        seqtk subseq ${protein_db} hits_ids.txt > filtered_hits.fasta
    else
        echo "No sequences passed e-value threshold" > filtered_hits.fasta
    fi





    """
}



process UsearchCluster {
    input:
    path "filtered_hits.fasta"
    output:
    path "clustered_sequences.fasta"

    script:
    """
    if [ -s filtered_hits.fasta ]; then
        usearch -cluster_fast filtered_hits.fasta -id 0.9 -centroids clustered_sequences.fasta
    else
        echo "0 sequence from clustering" > clustered_sequences.fasta
    fi
    """
}

process BuildTree {
    input:
    path "clustered_sequences.fasta"
    output:
    path "phylogenetic_tree.nwk"

    script:
    """
    if [ -s clustered_sequences.fasta ]; then
        mafft --auto clustered_sequences.fasta > aligned_sequences.fasta
        FastTree -nt aligned_sequences.fasta > phylogenetic_tree.nwk
    else
        echo "0 sequence from tree building" > phylogenetic_tree.nwk
    fi
    """
}



process CalculateThreshold {
    input:
    path "phylogenetic_tree.nwk" 

    output:
    path "calculated_threshold.txt" 

    script:
    """
    python - <<EOF
from ete3 import Tree
import numpy as np

# Load the Newick tree
tree = Tree("phylogenetic_tree.nwk")

# Calculate all pairwise distances
distances = [tree.get_distance(leaf1, leaf2)
             for leaf1 in tree for leaf2 in tree if leaf1 != leaf2]
print(distances)
#median
threshold = np.median(distances)
print(threshold)

max_threshold = 0.5
min_threshold = 0.2

if threshold > max_threshold:
        threshold = max_threshold
if threshold < min_threshold:
        threshold = min_threshold

with open("calculated_threshold.txt", "w") as f:
    f.write(str(threshold) + "\\n")
EOF
    """
}




process TreeCluster {
    input:
    path "calculated_threshold.txt"
    path "phylogenetic_tree.nwk"
    output:
    path "tree_clusters.txt"

    script:
    """
    calculated_threshold=\$(cat calculated_threshold.txt)
    #TreeCluster.py -i phylogenetic_tree.nwk -m max_clade -t 0.5 > tree_clusters.txt
    TreeCluster.py -i phylogenetic_tree.nwk -m avg_clade -t \$calculated_threshold > tree_clusters.txt


    """
}

process SelectRepresentatives {
    input:
    path "tree_clusters.txt"
    path "aligned_sequences.fasta"
    output:
    path "tree_clusters_reps.txt"

    script:
    """
    python - <<EOF
import random
from collections import defaultdict

clusters = defaultdict(list)
with open('tree_clusters.txt', 'r') as file:
    for line in file:
        seq_id, cluster_id = line.strip().split()
        clusters[cluster_id].append(seq_id)

with open('tree_clusters_reps.txt', 'w') as output:
    for cluster_id, seq_ids in clusters.items():
        if cluster_id == "-1":
        
            for seq_id in seq_ids:
                output.write(f"{seq_id} {cluster_id}\\n")
        else:
            # Randomly select one representative for other clusters
            representative = random.choice(seq_ids)
            output.write(f"{representative} {cluster_id}\\n")
EOF
    """

}



process GenerateSummary {
    input:
    path "hmmsearch_results.txt"         
    path "filtered_hits.fasta"   
    path "clustered_sequences.fasta"  // Uclust 
    path "tree_clusters.txt"

    output:
    path "summary.txt"

    script:
    """
  
    total_hmm_hits=\$(grep -c -v '^#' hmmsearch_results.txt)


    filtered_hits=\$(grep -c '^>' filtered_hits.fasta)


    representative_sequences=\$(grep -c '^>' clustered_sequences.fasta)
    
    non_singleton_clusters=\$(awk '{print \$2}' tree_clusters.txt | grep -v -c '^-1')
    singletons=\$(awk '{print \$2}' tree_clusters.txt | grep -c '^-1')


    echo "Summary of Analysis:" > summary.txt
    echo "--------------------" >> summary.txt
    echo "Total number of HMM hits: \$total_hmm_hits" >> summary.txt
    echo "Number of HMM hits after filtering: \$filtered_hits" >> summary.txt
    echo "Number of representative sequences after uclust: \$representative_sequences" >> summary.txt
    echo "Number of non-singleton clusters: \$non_singleton_clusters" >> summary.txt
    echo "Number of total singletons (-1): \$singletons" >> summary.txt
    """
}





workflow {
    Channel
        .fromPath(params.protein_db)
        .set { protein_db }

    Channel
        .fromPath(params.hmm_model)
        .set { hmm_model }

    HMMSearch(protein_db, hmm_model)
    FilterHits(HMMSearch.out, protein_db)
    UsearchCluster(FilterHits.out)
    BuildTree(UsearchCluster.out)
    CalculateThreshold(BuildTree.out)
    TreeCluster(CalculateThreshold.out,BuildTree.out)
    SelectRepresentatives(TreeCluster.out, UsearchCluster.out)
    GenerateSummary(HMMSearch.out, FilterHits.out, UsearchCluster.out, TreeCluster.out)
}

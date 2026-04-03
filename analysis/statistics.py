# analysis/statistics.py
# This module provides simple statistical functions for Gene Ontology analysis.
# These are plain functions (not a class) that work with an AnnotationManager.
#
# Functions:
#   build_gene_term_matrix()   - create a binary matrix of genes vs GO terms
#   namespace_distribution()   - count annotations per namespace (BP, MF, CC)
#   evidence_distribution()    - count annotations per evidence code

import numpy as np
import pandas as pd


# ============================================================
# build_gene_term_matrix()
# ============================================================
def build_gene_term_matrix(annotation_manager, top_n=100):
    """
    Build a binary matrix where rows are genes and columns are GO terms.

    A cell is 1 if that gene is annotated with that GO term, 0 otherwise.

    We only include the top N most-annotated genes to keep the matrix
    a reasonable size.

    Parameters
    ----------
    annotation_manager : AnnotationManager
        The annotation manager that holds all the data.
    top_n : int
        How many of the most-annotated genes to include (default 100).

    Returns
    -------
    matrix : numpy.ndarray
        A 2D array of shape (num_genes, num_terms) with 0s and 1s.
    gene_list : list
        The gene symbols, in the same order as the matrix rows.
    term_list : list
        The GO term IDs, in the same order as the matrix columns.
    """
    print("Building gene-term matrix...")

    # Step 1: Get the top N most-annotated genes
    top_genes = annotation_manager.get_top_genes(n=top_n)

    # Make a list of just the gene symbols (without the counts)
    gene_list = []
    for gene_symbol, count in top_genes:
        gene_list.append(gene_symbol)

    print("  Selected " + str(len(gene_list)) + " top genes")

    # Step 2: For each gene, get its GO terms
    gene_terms_dict = {}
    all_terms_set = set()

    for gene in gene_list:
        # Get the GO terms for this gene
        terms = annotation_manager.get_terms_for_gene(gene)

        # Convert to a set for fast lookup later
        terms_set = set(terms)

        # Store it
        gene_terms_dict[gene] = terms_set

        # Add these terms to our collection of all terms
        for term in terms_set:
            all_terms_set.add(term)

    # Convert the set of all terms to a sorted list
    
    term_list = sorted(list(all_terms_set))

    print("  Found " + str(len(term_list)) + " unique GO terms")

    # Step 3: Create the matrix
    # Rows = genes, Columns = GO terms
    num_genes = len(gene_list)
    num_terms = len(term_list)

    # Start with a matrix of all zeros
    matrix = np.zeros((num_genes, num_terms), dtype=int)

    # Fill in 1s where a gene has a GO term
    for i in range(num_genes):
        gene = gene_list[i]
        gene_terms = gene_terms_dict[gene]

        for j in range(num_terms):
            term = term_list[j]

            # Check if this gene has this GO term
            if term in gene_terms:
                matrix[i][j] = 1

        # Print progress every 25 genes
        if (i + 1) % 25 == 0:
            print("  Processed " + str(i + 1) + " / " + str(num_genes) + " genes...")

    print("  Matrix shape: " + str(num_genes) + " genes x " + str(num_terms) + " terms")
    print("Building gene-term matrix... done!")

    return matrix, gene_list, term_list


# ============================================================
# namespace_distribution()
# ============================================================
def namespace_distribution(annotation_manager):
    """
    Count how many annotations fall into each GO namespace.

    The three namespaces are:
      BP = Biological Process   (aspect code "P")
      MF = Molecular Function   (aspect code "F")
      CC = Cellular Component   (aspect code "C")

    Parameters
    ----------
    annotation_manager : AnnotationManager
        The annotation manager that holds all the data.

    Returns
    -------
    dist : dict
        A dictionary like {"BP": 50000, "MF": 20000, "CC": 15000}
    """
    # Get the summary from the annotation manager
    summary = annotation_manager.summary()

    aspect_to_namespace = {
        "P": "BP",
        "F": "MF",
        "C": "CC"
    }

    # Build the result dictionary
    dist = {}

    # Get the by_namespace counts from the summary
    by_namespace = summary["by_namespace"]

    # Loop through the mapping and convert
    for aspect_code, namespace_name in aspect_to_namespace.items():
        # Check if this aspect code exists in the data
        if aspect_code in by_namespace:
            dist[namespace_name] = by_namespace[aspect_code]
        else:
            # If the aspect code is not in the data, set count to 0
            dist[namespace_name] = 0

    return dist


# ============================================================
# evidence_distribution()
# ============================================================
def evidence_distribution(annotation_manager):
    """
    Count how many annotations use each evidence code.

    Parameters
    ----------
    annotation_manager : AnnotationManager
        The annotation manager that holds all the data.

    Returns
    -------
    dist : dict
        A dictionary like {"IEA": 400000, "IDA": 50000, "TAS": 30000, ...}
    """
    # Get the summary from the annotation manager
    summary = annotation_manager.summary()

    # The by_evidence field already has what we need
    dist = summary["by_evidence"]

    return dist


# ============================================================
# matrix_summary()
# ============================================================
def matrix_summary(annotation_manager, top_n=20):
    """
    Build the gene-term binary matrix and return summary statistics.

    This function uses NumPy to compute properties of the matrix,
    such as density (how many cells are non-zero), average annotations
    per gene, and average genes per term.

    Parameters
    ----------
    annotation_manager : AnnotationManager
        The annotation manager that holds all the data.
    top_n : int
        How many of the most-annotated genes to include (default 20).

    Returns
    -------
    info : dict
        A dictionary with keys:
        - num_genes: number of rows in the matrix
        - num_terms: number of columns in the matrix
        - total_ones: total number of 1s (gene-term links)
        - density: fraction of cells that are 1 (between 0 and 1)
        - avg_annotations_per_gene: average number of GO terms per gene
        - avg_genes_per_term: average number of genes per GO term
        - gene_list: list of gene symbols (same order as matrix rows)
    """
    # Build the matrix using the function defined above
    matrix, gene_list, term_list = build_gene_term_matrix(annotation_manager, top_n)

    # Count the total number of 1s in the matrix using NumPy
    total_ones = int(np.sum(matrix))

    # Get the matrix dimensions
    num_genes = matrix.shape[0]
    num_terms = matrix.shape[1]

    # Density is the fraction of cells that contain a 1
    total_cells = num_genes * num_terms
    if total_cells > 0:
        density = total_ones / total_cells
    else:
        density = 0.0

    # Average annotations per gene = mean of row sums
    row_sums = np.sum(matrix, axis=1)
    avg_per_gene = float(np.mean(row_sums))

    # Average genes per term = mean of column sums
    col_sums = np.sum(matrix, axis=0)
    avg_per_term = float(np.mean(col_sums))

    info = {
        "num_genes": num_genes,
        "num_terms": num_terms,
        "total_ones": total_ones,
        "density": round(density, 4),
        "avg_annotations_per_gene": round(avg_per_gene, 1),
        "avg_genes_per_term": round(avg_per_term, 1),
        "gene_list": gene_list,
    }

    return info


# ============================================================
# gene_cosine_similarity()
# ============================================================
def gene_cosine_similarity(annotation_manager, top_n=10):
    """
    Compute pairwise cosine similarity between the top N genes
    using the binary gene-term matrix and NumPy vector operations.

    Cosine similarity measures how similar two gene annotation
    profiles are by computing:
        cos(A, B) = (A . B) / (||A|| * ||B||)

    where A and B are binary vectors of GO term annotations.

    Parameters
    ----------
    annotation_manager : AnnotationManager
        The annotation manager that holds all the data.
    top_n : int
        How many of the most-annotated genes to compare (default 10).

    Returns
    -------
    pairs : list of tuple
        Each tuple is (gene_a, gene_b, score) sorted by score
        descending, excluding self-pairs.
    """
    # Build the matrix
    matrix, gene_list, term_list = build_gene_term_matrix(annotation_manager, top_n)

    num_genes = len(gene_list)
    pairs = []

    # Compare every pair of genes
    for i in range(num_genes):
        for j in range(i + 1, num_genes):
            vec_a = matrix[i].astype(float)
            vec_b = matrix[j].astype(float)

            # Compute the dot product using NumPy
            dot_product = np.dot(vec_a, vec_b)

            # Compute the magnitude (L2 norm) of each vector
            norm_a = np.sqrt(np.sum(vec_a ** 2))
            norm_b = np.sqrt(np.sum(vec_b ** 2))

            # Compute cosine similarity
            if norm_a > 0 and norm_b > 0:
                score = dot_product / (norm_a * norm_b)
            else:
                score = 0.0

            pairs.append((gene_list[i], gene_list[j], round(float(score), 4)))

    # Sort by similarity score, highest first
    pairs.sort(key=lambda x: x[2], reverse=True)

    return pairs

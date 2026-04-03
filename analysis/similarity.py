# analysis/similarity.py
# This module defines classes that measure how similar two genes are.
# We compare genes by looking at which GO terms they share.
#
# Classes:
#   SimilarityMeasure   - base class (returns 0.0 by default)
#   JaccardSimilarity   - uses Jaccard index  (intersection / union)
#   OverlapSimilarity   - uses overlap coeff  (intersection / min size)


# ============================================================
# SimilarityMeasure  (base class)
# ============================================================
class SimilarityMeasure:
    """
    Base class for all similarity measures.
    Subclasses override the calculate() method with their own formula.
    """

    def __init__(self, annotation_manager):
        # Store the annotation manager so we can look up gene annotations
        self._annotation_manager = annotation_manager

    def calculate(self, gene_a, gene_b):
        """
        Calculate the similarity between two genes.
        The base implementation just returns 0.0.
        Subclasses will override this with a real formula.
        """
        return 0.0

    def get_name(self):
        """Return the name of this similarity measure."""
        return "Base Similarity"

    def get_description(self):
        """Return a short description of what this measure does."""
        return "Base similarity measure"


# ============================================================
# JaccardSimilarity
# ============================================================
class JaccardSimilarity(SimilarityMeasure):
    """
    Jaccard similarity between two genes.

    The Jaccard index is calculated as:
        |intersection| / |union|

    In plain English: count how many GO terms the two genes share,
    then divide by the total number of distinct GO terms across both genes.

    The result is a number between 0.0 (nothing in common) and
    1.0 (exactly the same GO terms).
    """

    def calculate(self, gene_a, gene_b):
        """
        Calculate the Jaccard similarity between gene_a and gene_b.

        Steps:
        1. Get the GO terms for gene_a
        2. Get the GO terms for gene_b
        3. Convert both lists to sets
        4. Compute intersection and union
        5. Return |intersection| / |union|
        """
        # Step 1: get GO terms for gene A from the annotation manager
        terms_a = self._annotation_manager.get_terms_for_gene(gene_a)

        # Step 2: get GO terms for gene B
        terms_b = self._annotation_manager.get_terms_for_gene(gene_b)

        # Step 3: convert both to sets so we can use set operations
        set_a = set(terms_a)
        set_b = set(terms_b)

        # If both sets are empty, the genes have no annotations at all
        if len(set_a) == 0 and len(set_b) == 0:
            return 0.0

        # Step 4: find the intersection (GO terms they share)
        intersection = set_a.intersection(set_b)

        # Find the union (all GO terms from both genes combined)
        union = set_a.union(set_b)

        # Step 5: calculate the Jaccard index
        jaccard_score = len(intersection) / len(union)

        return jaccard_score

    def get_name(self):
        """Return the name of this similarity measure."""
        return "Jaccard Similarity"

    def get_description(self):
        """Return a short description of what this measure does."""
        return "Measures similarity as the size of the intersection divided by the size of the union of GO term sets."


# ============================================================
# OverlapSimilarity
# ============================================================
class OverlapSimilarity(SimilarityMeasure):
    """
    Overlap coefficient between two genes.

    The overlap coefficient is calculated as:
        |intersection| / min(|set_a|, |set_b|)

    This is useful when one gene has many more annotations than the other.
    If the smaller set is completely contained in the larger set,
    the overlap coefficient will be 1.0.
    """

    def calculate(self, gene_a, gene_b):
        """
        Calculate the overlap coefficient between gene_a and gene_b.

        Steps:
        1. Get the GO terms for each gene
        2. Convert to sets
        3. Compute intersection
        4. Divide by the size of the smaller set
        """
        # Step 1: get GO terms for both genes
        terms_a = self._annotation_manager.get_terms_for_gene(gene_a)
        terms_b = self._annotation_manager.get_terms_for_gene(gene_b)

        # Step 2: convert to sets
        set_a = set(terms_a)
        set_b = set(terms_b)

        # If either set is empty, there is no overlap
        if len(set_a) == 0 or len(set_b) == 0:
            return 0.0

        # Step 3: find the intersection
        intersection = set_a.intersection(set_b)

        # Step 4: find the size of the smaller set
        smaller_size = min(len(set_a), len(set_b))

        # Calculate the overlap coefficient
        overlap_score = len(intersection) / smaller_size

        return overlap_score

    def get_name(self):
        """Return the name of this similarity measure."""
        return "Overlap Coefficient"

    def get_description(self):
        """Return a short description of what this measure does."""
        return "Measures similarity as the size of the intersection divided by the size of the smaller GO term set."

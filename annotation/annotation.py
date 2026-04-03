# annotation/annotation.py
# This module defines classes for Gene Ontology annotations.
# An annotation links a gene to a GO term with evidence about why.
#
# Classes:
#   Annotation             - base class for any annotation
#   ExperimentalAnnotation - for annotations backed by lab experiments
#   ElectronicAnnotation   - for annotations made by computers (IEA)
#   AnnotationManager      - holds all annotations and lets you query them
#
# Helper function:
#   create_annotation()    - factory that picks the right class automatically

import pandas as pd


# ============================================================
# Evidence codes that count as experimental
# These are the codes assigned when a real experiment was done
# ============================================================
EXPERIMENTAL_CODES = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"]


# ============================================================
# Annotation  (base class)
# ============================================================
class Annotation:
    """
    Base class that represents one gene-to-GO-term annotation.

    Each annotation says: "this gene is associated with this GO term,
    and here is the evidence code that tells us why we believe it."
    """

    def __init__(self, gene_symbol, go_id, evidence_code, aspect, date):
        # Store everything as private attributes (leading underscore)
        self._gene_symbol = gene_symbol      # e.g. "TP53"
        self._go_id = go_id                  # e.g. "GO:0006915"
        self._evidence_code = evidence_code  # e.g. "IDA"
        self._aspect = aspect                # "P", "F", or "C"
        self._date = date                    # date string from GAF file

    # --- properties so outsiders can read but not change the values ---

    @property
    def gene_symbol(self):
        """Return the gene symbol, like TP53."""
        return self._gene_symbol

    @property
    def go_id(self):
        """Return the GO term id, like GO:0006915."""
        return self._go_id

    @property
    def evidence_code(self):
        """Return the evidence code, like IDA or IEA."""
        return self._evidence_code

    @property
    def aspect(self):
        """Return the aspect: P (process), F (function), or C (component)."""
        return self._aspect

    @property
    def date(self):
        """Return the annotation date string."""
        return self._date

    # --- methods ---

    def describe(self):
        """
        Return a short human-readable description of this annotation.
        Format: "gene_symbol -- go_id (evidence_code)"
        """
        return self._gene_symbol + " -- " + self._go_id + " (" + self._evidence_code + ")"

    def reliability_score(self):
        """
        Return a number between 0 and 1 that says how reliable
        this annotation is.  The base class just returns 0.5.
        Subclasses override this to give better scores.
        """
        return 0.5

    def __str__(self):
        """When you print() an annotation, show the describe() text."""
        return self.describe()

    def __repr__(self):
        """Show a helpful representation in the debugger / REPL."""
        return "Annotation(" + self.describe() + ")"


# ============================================================
# ExperimentalAnnotation  (for lab-verified evidence)
# ============================================================
class ExperimentalAnnotation(Annotation):
    """
    An annotation that comes from a real laboratory experiment.
    Evidence codes: EXP, IDA, IPI, IMP, IGI, IEP.
    These are the most trustworthy annotations.
    """

    def describe(self):
        """Add [Experimental] tag in front of the base description."""
        # First get the base description from the parent class
        base_text = super().describe()
        # Then add the prefix
        return "[Experimental] " + base_text

    def reliability_score(self):
        """
        Experimental annotations are the most reliable.
        Return 1.0 (maximum reliability).
        """
        return 1.0


# ============================================================
# ElectronicAnnotation  (for computer-predicted evidence)
# ============================================================
class ElectronicAnnotation(Annotation):
    """
    An annotation that was created automatically by a computer program.
    Evidence code: IEA (Inferred from Electronic Annotation).
    These are less reliable than experimental annotations.
    """

    def describe(self):
        """Add [Electronic] tag in front of the base description."""
        # First get the base description from the parent class
        base_text = super().describe()
        # Then add the prefix
        return "[Electronic] " + base_text

    def reliability_score(self):
        """
        Electronic annotations are less reliable.
        Return 0.5 (moderate reliability).
        """
        return 0.5


# ============================================================
# create_annotation()  -  helper / factory function
# ============================================================
def create_annotation(gene_symbol, go_id, evidence_code, aspect, date):
    """
    Create the right kind of Annotation object based on the evidence code.

    - If the evidence code is experimental (EXP, IDA, etc.)
      we return an ExperimentalAnnotation.
    - If the evidence code is IEA (electronic)
      we return an ElectronicAnnotation.
    - For anything else we return a plain Annotation.
    """
    # Check if it is an experimental evidence code
    if evidence_code in EXPERIMENTAL_CODES:
        return ExperimentalAnnotation(gene_symbol, go_id, evidence_code, aspect, date)

    # Check if it is the electronic evidence code
    elif evidence_code == "IEA":
        return ElectronicAnnotation(gene_symbol, go_id, evidence_code, aspect, date)

    # For everything else, just use the base class
    else:
        return Annotation(gene_symbol, go_id, evidence_code, aspect, date)


# ============================================================
# AnnotationManager
# ============================================================
class AnnotationManager:
    """
    Holds all the annotations and provides ways to query them.

    This class keeps:
      - A list of Annotation objects (built from a sample of rows)
      - Lookup dicts so we can quickly find annotations by gene or GO term
      - The full original DataFrame for statistical queries

    We only build Annotation objects for the first 50,000 rows because
    the GAF file can have 600K+ rows and creating that many objects
    uses too much memory.  But we keep the full DataFrame so that
    statistical queries (counts, top genes, etc.) are accurate.
    """

    # How many rows to turn into Annotation objects
    SAMPLE_SIZE = 50000

    def __init__(self):
        # List that holds all Annotation objects we created
        self._annotations = []

        # Dictionary: gene_symbol -> list of Annotation objects
        self._by_gene = {}

        # Dictionary: go_id -> list of Annotation objects
        self._by_term = {}

        # The full original DataFrame (all rows, not just the sample)
        self._df = None

    def build_from_data(self, annotations_df):
        """
        Build the annotation collection from a pandas DataFrame.

        Parameters
        ----------
        annotations_df : pd.DataFrame
            The DataFrame that came from parsing the GAF file.
            Expected columns: db_object_symbol, go_id, evidence_code,
                              aspect, date  (plus others we don't need here).
        """
        # Keep the full dataframe for statistical queries later
        self._df = annotations_df

        print("Building annotations...")

        # Figure out how many rows to turn into objects
        total_rows = len(annotations_df)
        rows_to_process = min(self.SAMPLE_SIZE, total_rows)

        print("  Total rows in DataFrame: " + str(total_rows))
        print("  Building Annotation objects for first " + str(rows_to_process) + " rows...")

        # We only loop through the first SAMPLE_SIZE rows
        sample_df = annotations_df.head(rows_to_process)

        # Keep a counter so we can print progress
        counter = 0

        # Loop through each row in the sample
        for index, row in sample_df.iterrows():
            # Get the values we need from this row
            gene_symbol = row["db_object_symbol"]
            go_id = row["go_id"]
            evidence_code = row["evidence_code"]
            aspect = row["aspect"]
            date = row["date"]

            # Use the factory function to create the right kind of annotation
            annotation = create_annotation(gene_symbol, go_id, evidence_code, aspect, date)

            # Add it to our main list
            self._annotations.append(annotation)

            # --- Group by gene symbol ---
            # If we have not seen this gene before, make a new empty list
            if gene_symbol not in self._by_gene:
                self._by_gene[gene_symbol] = []
            # Add the annotation to that gene's list
            self._by_gene[gene_symbol].append(annotation)

            # --- Group by GO term ---
            # Same idea: make a new list if needed, then append
            if go_id not in self._by_term:
                self._by_term[go_id] = []
            self._by_term[go_id].append(annotation)

            # Print progress every 10000 rows
            counter = counter + 1
            if counter % 10000 == 0:
                print("  Processed " + str(counter) + " / " + str(rows_to_process) + " rows...")

        # Done!
        print("Building annotations... " + str(len(self._annotations)) + " annotation objects loaded")
        print("  Unique genes (in objects): " + str(len(self._by_gene)))
        print("  Unique GO terms (in objects): " + str(len(self._by_term)))

    # --------------------------------------------------------
    # Query methods that use the Annotation objects
    # --------------------------------------------------------

    def get_by_gene(self, symbol):
        """
        Get all Annotation objects for a given gene symbol.
        Returns an empty list if the gene is not found.
        """
        if symbol in self._by_gene:
            return self._by_gene[symbol]
        else:
            return []

    def get_by_term(self, go_id):
        """
        Get all Annotation objects for a given GO term.
        Returns an empty list if the term is not found.
        """
        if go_id in self._by_term:
            return self._by_term[go_id]
        else:
            return []

    def get_by_evidence(self, code):
        """
        Get all Annotation objects that have a specific evidence code.
        This loops through all objects and filters them.
        """
        
        results = []

    
        for annotation in self._annotations:
            
            if annotation.evidence_code == code:
                results.append(annotation)

        return results

    # --------------------------------------------------------
    # Query methods that use the FULL DataFrame
    # These give accurate counts over all 600K+ rows.
    # --------------------------------------------------------

    def get_terms_for_gene(self, symbol):
        """
        Get all GO term IDs that a gene is annotated with.
        Uses the FULL DataFrame so we don't miss anything.
        Returns a list of GO ID strings.
        """
        
        if self._df is None:
            return []

        # Filter the dataframe to rows for this gene
        gene_rows = self._df[self._df["db_object_symbol"] == symbol]

        # Get the unique GO IDs from those rows
        term_list = gene_rows["go_id"].unique().tolist()

        return term_list

    def get_genes_for_term(self, go_id):
        """
        Get all gene symbols that are annotated with a given GO term.
        Uses the FULL DataFrame so we don't miss anything.
        Returns a list of gene symbol strings.
        """
        # Make sure we have the dataframe
        if self._df is None:
            return []

        # Filter the dataframe to rows for this GO term
        term_rows = self._df[self._df["go_id"] == go_id]

        # Get the unique gene symbols from those rows
        gene_list = term_rows["db_object_symbol"].unique().tolist()

        return gene_list

    # --------------------------------------------------------
    # Summary / statistics methods  (use FULL DataFrame)
    # --------------------------------------------------------

    def summary(self):
        """
        Return a dictionary with summary statistics about all annotations.
        Uses the FULL DataFrame for accurate counts.

        Returns something like:
        {
            "total": 612345,
            "by_namespace": {"P": 300000, "F": 150000, "C": 162345},
            "by_evidence": {"IEA": 400000, "IDA": 50000, ...}
        }
        """
        # Start with an empty result dictionary
        result = {}

        if self._df is None:
            result["total"] = 0
            result["by_namespace"] = {}
            result["by_evidence"] = {}
            return result

        # Total number of annotations
        result["total"] = len(self._df)

        # Count by namespace (the "aspect" column: P, F, or C)
        namespace_counts = self._df["aspect"].value_counts()
        # Convert to a regular dictionary
        by_namespace = {}
        for aspect_code, count in namespace_counts.items():
            by_namespace[aspect_code] = int(count)
        result["by_namespace"] = by_namespace

        # Count by evidence code
        evidence_counts = self._df["evidence_code"].value_counts()
        # Convert to a regular dictionary
        by_evidence = {}
        for evidence_code, count in evidence_counts.items():
            by_evidence[evidence_code] = int(count)
        result["by_evidence"] = by_evidence

        return result

    def get_all_genes(self):
        """
        Return a sorted list of all unique gene symbols.
        Uses the FULL DataFrame.
        """
        if self._df is None:
            return []

        # Get unique gene symbols
        all_genes = self._df["db_object_symbol"].unique().tolist()

        # Sort them alphabetically
        all_genes.sort()

        return all_genes

    def get_top_genes(self, n=10):
        """
        Return the top N most-annotated genes.
        Each item is a tuple: (gene_symbol, annotation_count).
        Uses the FULL DataFrame.
        """
        if self._df is None:
            return []

        # Count how many annotations each gene has
        gene_counts = self._df["db_object_symbol"].value_counts()

        # Get the top N
        top = gene_counts.head(n)

        # Convert to a list of tuples
        result = []
        for gene_symbol, count in top.items():
            result.append((gene_symbol, int(count)))

        return result

    def get_top_terms(self, n=10):
        """
        Return the top N most-used GO terms.
        Each item is a tuple: (go_id, annotation_count).
        Uses the FULL DataFrame.
        """
        if self._df is None:
            return []

        # Count how many annotations each GO term has
        term_counts = self._df["go_id"].value_counts()

        # Get the top N
        top = term_counts.head(n)

        # Convert to a list of tuples
        result = []
        for go_id, count in top.items():
            result.append((go_id, int(count)))

        return result

    # --------------------------------------------------------
    # Detailed annotation queries  (use FULL DataFrame)
    # These return lists of dictionaries with annotation info.
    # --------------------------------------------------------

    def get_annotations_for_term(self, go_id, limit=50):
        """
        Get annotation details for a specific GO term.
        Uses the FULL DataFrame so we don't miss any annotations.

        Parameters
        ----------
        go_id : str
            The GO term id (e.g. "GO:0008150").
        limit : int
            Maximum number of annotations to return (default 50).

        Returns
        -------
        list of dict
            Each dictionary has keys: "gene_symbol", "evidence_code".
        """
        if self._df is None:
            return []

        # Filter the dataframe to rows for this GO term
        term_rows = self._df[self._df["go_id"] == go_id]

        # Build a list of dictionaries, limited to the first N rows
        results = []
        for _, row in term_rows.head(limit).iterrows():
            results.append({
                "gene_symbol": row["db_object_symbol"],
                "evidence_code": row["evidence_code"],
            })

        return results

    def get_annotations_for_gene(self, symbol):
        """
        Get annotation details for a specific gene symbol.
        Uses the FULL DataFrame so we don't miss any annotations.
        Tries the exact symbol first, then uppercase if not found.

        Parameters
        ----------
        symbol : str
            The gene symbol (e.g. "TP53").

        Returns
        -------
        list of dict
            Each dictionary has keys: "go_id", "evidence_code", "aspect".
        """
        if self._df is None:
            return []

        # Try to find rows for this gene symbol
        gene_rows = self._df[self._df["db_object_symbol"] == symbol]

        # If nothing found, try uppercase version
        if len(gene_rows) == 0:
            gene_upper = symbol.upper()
            if gene_upper != symbol:
                gene_rows = self._df[self._df["db_object_symbol"] == gene_upper]

        # Build a list of dictionaries with the annotation details
        results = []
        for _, row in gene_rows.iterrows():
            results.append({
                "go_id": row["go_id"],
                "evidence_code": row["evidence_code"],
                "aspect": row.get("aspect", ""),
            })

        return results

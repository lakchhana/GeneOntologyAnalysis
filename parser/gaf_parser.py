# gaf_parser.py
# This module reads Gene Annotation Format (GAF) files.
# GAF files contain gene-to-GO-term annotations.
# They are tab-separated and may be gzipped (.gaf.gz).

import pandas as pd


class GAFParser:
    """
    A parser for GAF (Gene Annotation Format) files.

    GAF files are tab-separated with 17 columns and no header row.
    Lines starting with '!' are comments and should be skipped.
    The file might be gzipped — pandas handles that automatically.
    """

    def parse(self, filepath):
        """
        Parse a GAF file and return the annotations as a DataFrame.

        Parameters:
            filepath - path to the .gaf or .gaf.gz file

        Returns:
            A pandas DataFrame with all the annotation rows
        """

        print("Parsing GAF file...")

        # These are the 17 standard column names for a GAF 2.1 file.
        # We set them manually because the file has no header row.
        column_names = [
            "db",                    # database contributing the annotation (e.g. UniProtKB)
            "db_object_id",          # unique id in the database (e.g. P12345)
            "db_object_symbol",      # gene or protein symbol (e.g. TP53)
            "qualifier",             # relationship qualifier (e.g. NOT, contributes_to)
            "go_id",                 # the GO term id (e.g. GO:0008150)
            "db_reference",          # literature reference (e.g. PMID:12345)
            "evidence_code",         # evidence code (e.g. IDA, IEA)
            "with_from",             # supporting evidence ids
            "aspect",                # which ontology: P, F, or C
            "db_object_name",        # full name of the gene/protein
            "db_object_synonym",     # alternative names, pipe-separated
            "db_object_type",        # type of entity (e.g. protein, gene)
            "taxon",                 # species taxon id (e.g. taxon:9606)
            "date",                  # annotation date (YYYYMMDD)
            "assigned_by",           # who made the annotation
            "annotation_extension",  # extra relationship info (optional column)
            "gene_product_form_id",  # specific isoform id (optional column)
        ]

        # Use pandas read_csv to load the file.
        annotations = pd.read_csv(
            filepath,
            comment="!",
            sep="\t",
            header=None,
            names=column_names,
            dtype=str,           # read everything as strings to be safe
            low_memory=False,    # avoid mixed type warnings
        )

        # Print some info about what we loaded
        print("Loaded " + str(len(annotations)) + " annotations")
        print("Unique GO terms: " + str(annotations["go_id"].nunique()))
        print("Unique genes: " + str(annotations["db_object_symbol"].nunique()))

        # Return the dataframe
        return annotations

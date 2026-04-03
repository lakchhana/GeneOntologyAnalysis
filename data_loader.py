"""
data_loader.py
==============
This file has one main function: load_all_data().
It reads the OBO and GAF data files, builds the ontology and annotations,
and creates the similarity calculators. Then it returns everything
in a dictionary so the Flask app can use it.

This keeps the loading logic separate from the web app, which makes
the code easier to understand and test.
"""

import os

# Import the parsers that read the raw data files
from parser.obo_parser import OBOParser
from parser.gaf_parser import GAFParser

# Import the ontology class that stores GO terms and their relationships
from ontology.ontology import Ontology

# Import the annotation manager that stores gene-to-GO-term links
from annotation.annotation import AnnotationManager

# Import the similarity calculators
from analysis.similarity import JaccardSimilarity, OverlapSimilarity


def load_all_data():
    """
    Load and build all data objects needed by the Flask app.

    This function does the following steps:
    1. Check that the data files exist on disk
    2. Parse the OBO file to get GO terms and edges
    3. Parse the GAF file to get gene annotations
    4. Build the Ontology from terms and edges
    5. Build the AnnotationManager from annotation rows
    6. Create JaccardSimilarity and OverlapSimilarity calculators
    7. Return everything in a dictionary

    Returns:
        A dictionary with keys:
            "ontology"            -> Ontology object
            "annotation_manager"  -> AnnotationManager object
            "jaccard"             -> JaccardSimilarity object
            "overlap"             -> OverlapSimilarity object

    Raises:
        FileNotFoundError if the data files are missing.
    """

    # ---------------------------------------------------------------
    # Figure out where the data files should be
    # ---------------------------------------------------------------
    # The data folder is inside the same directory as this script
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")

    # These are the two files we need
    obo_file = os.path.join(data_dir, "go-basic.obo")
    gaf_file = os.path.join(data_dir, "goa_human.gaf.gz")

    # ---------------------------------------------------------------
    # Step 1: Check that both data files exist
    # ---------------------------------------------------------------
    print("=" * 60)
    print("Loading Gene Ontology data...")
    print("=" * 60)

    # Check for the OBO file
    if not os.path.exists(obo_file):
        raise FileNotFoundError(
            "Could not find the OBO file at: " + obo_file + "\n"
            "Please download go-basic.obo and place it in the data/ folder.\n"
            "You can get it from: http://purl.obolibrary.org/obo/go/go-basic.obo"
        )

    # Check for the GAF file
    if not os.path.exists(gaf_file):
        raise FileNotFoundError(
            "Could not find the GAF file at: " + gaf_file + "\n"
            "Please download goa_human.gaf.gz and place it in the data/ folder.\n"
            "You can get it from: https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
        )

    print("Data files found!")
    print("  OBO file: " + obo_file)
    print("  GAF file: " + gaf_file)
    print()

    # ---------------------------------------------------------------
    # Step 2: Parse the OBO file to get GO terms and edges
    # ---------------------------------------------------------------
    try:
        print("--- Parsing OBO file ---")
        obo_parser = OBOParser()
        terms_df, edges = obo_parser.parse(obo_file)
        print("OBO parsing complete!")
        print()
    except Exception as e:
        print("ERROR: Something went wrong while parsing the OBO file.")
        print("Details: " + str(e))
        raise

    # ---------------------------------------------------------------
    # Step 3: Parse the GAF file to get gene annotations
    # ---------------------------------------------------------------
    try:
        print("--- Parsing GAF file ---")
        gaf_parser = GAFParser()
        annotations_df = gaf_parser.parse(gaf_file)
        print("GAF parsing complete!")
        print()
    except Exception as e:
        print("ERROR: Something went wrong while parsing the GAF file.")
        print("Details: " + str(e))
        raise

    # ---------------------------------------------------------------
    # Step 4: Build the Ontology from terms and edges
    # ---------------------------------------------------------------
    try:
        print("--- Building Ontology ---")
        ontology = Ontology()
        ontology.build_from_data(terms_df, edges)
        print("Ontology built! Total terms: " + str(ontology.total_terms()))
        print()
    except Exception as e:
        print("ERROR: Something went wrong while building the ontology.")
        print("Details: " + str(e))
        raise

    # ---------------------------------------------------------------
    # Step 5: Build the AnnotationManager from annotation rows
    # ---------------------------------------------------------------
    try:
        print("--- Building Annotation Manager ---")
        am = AnnotationManager()
        am.build_from_data(annotations_df)
        summary = am.summary()
        print("Annotation manager built!")
        print("  Total annotations: " + str(summary["total"]))
        print()
    except Exception as e:
        print("ERROR: Something went wrong while building annotations.")
        print("Details: " + str(e))
        raise

    # ---------------------------------------------------------------
    # Step 6: Create the similarity calculators
    # ---------------------------------------------------------------
    print("--- Creating similarity calculators ---")
    jaccard = JaccardSimilarity(am)
    overlap = OverlapSimilarity(am)
    print("Similarity calculators ready!")
    print()

    # ---------------------------------------------------------------
    # Step 7: Return everything in a dictionary
    # ---------------------------------------------------------------
    print("=" * 60)
    print("All data loaded successfully!")
    print("=" * 60)
    print()

    result = {
        "ontology": ontology,
        "annotation_manager": am,
        "jaccard": jaccard,
        "overlap": overlap,
    }

    return result

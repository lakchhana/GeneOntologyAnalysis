# tests/test_parser.py
# Simple tests for the OBO and GAF parsers

import sys
import os

# Add the project root to the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parser.obo_parser import OBOParser
from parser.gaf_parser import GAFParser


def test_obo_parser():
    """Test that the OBO parser reads the file correctly."""
    parser = OBOParser()
    data_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data", "go-basic.obo"
    )

    if not os.path.exists(data_path):
        print("  SKIPPED - data file not found")
        return

    terms_df, edges = parser.parse(data_path)

    # Check that we got some terms
    assert len(terms_df) > 40000, f"Expected >40000 terms, got {len(terms_df)}"

    # Check that the DataFrame has the right columns
    expected_columns = ["go_id", "name", "namespace", "definition", "is_obsolete"]
    for col in expected_columns:
        assert col in terms_df.columns, f"Missing column: {col}"

    # Check that edges were found
    assert len(edges) > 50000, f"Expected >50000 edges, got {len(edges)}"

    # Check that an edge is a tuple of two GO IDs
    first_edge = edges[0]
    assert len(first_edge) == 2
    assert first_edge[0].startswith("GO:")
    assert first_edge[1].startswith("GO:")

    print("  test_obo_parser PASSED")


def test_gaf_parser():
    """Test that the GAF parser reads the file correctly."""
    parser = GAFParser()
    data_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data", "goa_human.gaf.gz"
    )

    if not os.path.exists(data_path):
        print("  SKIPPED - data file not found")
        return

    annotations_df = parser.parse(data_path)

    # Check that we got lots of annotations
    assert len(annotations_df) > 500000, f"Expected >500000 annotations, got {len(annotations_df)}"

    # Check that important columns exist
    assert "db_object_symbol" in annotations_df.columns
    assert "go_id" in annotations_df.columns
    assert "evidence_code" in annotations_df.columns

    print("  test_gaf_parser PASSED")


# Run all tests
if __name__ == "__main__":
    print("Running parser tests...")
    test_obo_parser()
    test_gaf_parser()
    print("\nAll parser tests passed!")

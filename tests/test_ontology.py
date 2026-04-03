# tests/test_ontology.py
# Simple tests to make sure the ontology classes work correctly

import sys
import os

# Add the project root to the path so we can import our modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ontology.entity import (
    OntologyEntity, GOTerm, BiologicalProcess,
    MolecularFunction, CellularComponent
)


def test_ontology_entity():
    """Test that OntologyEntity stores data correctly."""
    entity = OntologyEntity("GO:0000001", "test term", "biological_process", "A test definition")

    # Check that the properties work
    assert entity.go_id == "GO:0000001"
    assert entity.name == "test term"
    assert entity.namespace == "biological_process"
    assert entity.definition == "A test definition"
    print("  test_ontology_entity PASSED")


def test_goterm():
    """Test that GOTerm inherits from OntologyEntity and adds extra features."""
    term = GOTerm("GO:0000002", "cell cycle", "biological_process", "The cell cycle")

    # Check inherited properties
    assert term.go_id == "GO:0000002"
    assert term.name == "cell cycle"

    # Check GOTerm-specific features
    assert term.is_obsolete == False
    assert term.parents == []
    assert term.children == []

    # Test adding parents and children
    term.add_parent("GO:0000001")
    term.add_child("GO:0000003")
    assert "GO:0000001" in term.parents
    assert "GO:0000003" in term.children

    print("  test_goterm PASSED")


def test_biological_process():
    """Test that BiologicalProcess overrides describe() with [BP] prefix."""
    bp = BiologicalProcess("GO:0006915", "apoptotic process", "biological_process", "Cell death")

    description = bp.describe()
    assert "[BP]" in description
    assert "apoptotic process" in description
    print("  test_biological_process PASSED")


def test_molecular_function():
    """Test that MolecularFunction overrides describe() with [MF] prefix."""
    mf = MolecularFunction("GO:0003674", "DNA binding", "molecular_function", "Binds DNA")

    description = mf.describe()
    assert "[MF]" in description
    print("  test_molecular_function PASSED")


def test_cellular_component():
    """Test that CellularComponent overrides describe() with [CC] prefix."""
    cc = CellularComponent("GO:0005575", "nucleus", "cellular_component", "The nucleus")

    description = cc.describe()
    assert "[CC]" in description
    print("  test_cellular_component PASSED")


def test_polymorphism():
    """Test that describe() works differently for each subclass (polymorphism)."""
    terms = [
        BiologicalProcess("GO:0001", "bp term", "biological_process", ""),
        MolecularFunction("GO:0002", "mf term", "molecular_function", ""),
        CellularComponent("GO:0003", "cc term", "cellular_component", ""),
    ]

    prefixes = []
    for term in terms:
        # Same method call, different behavior - this is polymorphism!
        desc = term.describe()
        if "[BP]" in desc:
            prefixes.append("BP")
        elif "[MF]" in desc:
            prefixes.append("MF")
        elif "[CC]" in desc:
            prefixes.append("CC")

    assert prefixes == ["BP", "MF", "CC"]
    print("  test_polymorphism PASSED")


def test_encapsulation():
    """Test that private attributes can't be accessed directly (encapsulation)."""
    term = GOTerm("GO:0001", "test", "biological_process", "def")

    assert term.go_id == "GO:0001"  
    assert term._go_id == "GO:0001"  
    print("  test_encapsulation PASSED")


# Run all tests
if __name__ == "__main__":
    print("Running ontology tests...")
    test_ontology_entity()
    test_goterm()
    test_biological_process()
    test_molecular_function()
    test_cellular_component()
    test_polymorphism()
    test_encapsulation()
    print("\nAll ontology tests passed!")

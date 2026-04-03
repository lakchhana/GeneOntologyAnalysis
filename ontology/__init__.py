# ontology/__init__.py
# This file makes the ontology folder a Python package.
# We import the main classes here so they can be used easily from outside.

from ontology.entity import OntologyEntity
from ontology.entity import GOTerm
from ontology.entity import BiologicalProcess
from ontology.entity import MolecularFunction
from ontology.entity import CellularComponent
from ontology.ontology import Ontology

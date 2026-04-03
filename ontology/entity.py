# ontology/entity.py
# This file defines the class hierarchy for Gene Ontology entities.
# We have a base class (OntologyEntity) and several subclasses that
# represent different types of GO terms.


# ============================================================
# Base class: OntologyEntity
# This is the parent class for all ontology objects.
# It stores the basic information that every GO entity has.
# ============================================================

class OntologyEntity:
    """
    Base class for all Gene Ontology entities.
    Stores the GO ID, name, namespace, and definition.
    """

    def __init__(self, go_id, name, namespace, definition):
         
        self._go_id = go_id
        self._name = name
        self._namespace = namespace
        self._definition = definition

    @property
    def go_id(self):
        # Returns the GO ID string, like "GO:0008150"
        return self._go_id

    @property
    def name(self):
        # Returns the human-readable name of this entity
        return self._name

    @property
    def namespace(self):
        # Returns which namespace this belongs to
        return self._namespace

    @property
    def definition(self):
        # Returns the text definition of this entity
        return self._definition

    # --- Regular methods ---

    def describe(self):
        # This function returns a short description string.
        return self._go_id + " - " + self._name

    def get_info(self):
        # This function returns a dictionary with all the basic attributes.
        info = {
            "go_id": self._go_id,
            "name": self._name,
            "namespace": self._namespace,
            "definition": self._definition,
        }
        return info

    def __str__(self):
         
        return self.describe()

    def __repr__(self):
        
     return "OntologyEntity(" + self._go_id + ")"


# ============================================================
# GOTerm class — inherits from OntologyEntity
# This adds parent/child relationships and obsolete flag.
# Most GO terms in the ontology will be this type or a subclass.
# ============================================================

class GOTerm(OntologyEntity):
    """
    Represents a single Gene Ontology term.
    Inherits basic info from OntologyEntity and adds
    parent/child relationships and an obsolete flag.
    """

    def __init__(self, go_id, name, namespace, definition, is_obsolete=False):
        
        super().__init__(go_id, name, namespace, definition)

        self._parents = []   # list of parent GO IDs
        self._children = []  # list of child GO IDs
        self._is_obsolete = is_obsolete

    # --- Property methods for the new attributes ---

    @property
    def parents(self):
        # Returns the list of parent GO IDs
        return self._parents

    @property
    def children(self):
        # Returns the list of child GO IDs
        return self._children

    @property
    def is_obsolete(self):
        # Returns True if this term is obsolete, False otherwise
        return self._is_obsolete

    # --- Methods to build relationships ---

    def add_parent(self, parent_id):
        # This function adds a parent GO ID to this term's parent list.
        if parent_id not in self._parents:
            self._parents.append(parent_id)

    def add_child(self, child_id):
        # This function adds a child GO ID to this term's children list.
        if child_id not in self._children:
            self._children.append(child_id)

    # --- Override methods from parent class ---

    def describe(self):
        
        description = self._go_id + " - " + self._name
        description = description + " (" + self._namespace + ")"
        description = description + " [parents: " + str(len(self._parents))
        description = description + ", children: " + str(len(self._children)) + "]"
        if self._is_obsolete:
            description = description + " [OBSOLETE]"
        return description

    def get_info(self):
        info = super().get_info()

        # Now add the new fields to the dictionary.
        info["parents"] = self._parents
        info["children"] = self._children
        info["is_obsolete"] = self._is_obsolete

        return info

    def __repr__(self):
        # Override repr to show this is a GOTerm, not just OntologyEntity.
        return "GOTerm(" + self._go_id + ")"


# ============================================================
# BiologicalProcess — inherits from GOTerm
# Represents terms in the "biological_process" namespace.
# ============================================================

class BiologicalProcess(GOTerm):
    """
    A GO term that belongs to the biological_process namespace.
    Examples: cell division, DNA repair, signal transduction.
    """

    def __init__(self, go_id, name, definition, is_obsolete=False):
        
        super().__init__(go_id, name, "biological_process", definition, is_obsolete)

    def describe(self):
        base_description = super().describe()
        return "[BP] " + base_description

    def __repr__(self):
        # Show the class name in repr output
        return "BiologicalProcess(" + self._go_id + ")"


# ============================================================
# MolecularFunction — inherits from GOTerm
# Represents terms in the "molecular_function" namespace.
# ============================================================

class MolecularFunction(GOTerm):
    """
    A GO term that belongs to the molecular_function namespace.
    Examples: kinase activity, DNA binding, transporter activity.
    """

    def __init__(self, go_id, name, definition, is_obsolete=False):
        # The namespace is always "molecular_function" for this class.
        super().__init__(go_id, name, "molecular_function", definition, is_obsolete)

    def describe(self):
        base_description = super().describe()
        return "[MF] " + base_description

    def __repr__(self):
        # Show the class name in repr output
        return "MolecularFunction(" + self._go_id + ")"


# ============================================================
# CellularComponent — inherits from GOTerm
# Represents terms in the "cellular_component" namespace.
# ============================================================

class CellularComponent(GOTerm):
    """
    A GO term that belongs to the cellular_component namespace.
    Examples: nucleus, mitochondrion, plasma membrane.
    """

    def __init__(self, go_id, name, definition, is_obsolete=False):
        # The namespace is always "cellular_component" for this class.
        super().__init__(go_id, name, "cellular_component", definition, is_obsolete)

    def describe(self):
        
        base_description = super().describe()
        return "[CC] " + base_description

    def __repr__(self):
        # Show the class name in repr output
        return "CellularComponent(" + self._go_id + ")"

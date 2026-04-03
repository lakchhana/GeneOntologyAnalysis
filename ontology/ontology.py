# ontology/ontology.py
# This file defines the Ontology class which is the main container
# for all GO terms and their relationships. It uses a NetworkX
# directed graph to store the parent-child edges between terms.

import networkx as nx

# Import our custom classes from entity.py
from ontology.entity import BiologicalProcess
from ontology.entity import MolecularFunction
from ontology.entity import CellularComponent


class Ontology:
    """
    The Ontology class holds all GO terms and their relationships.
    It uses a NetworkX DiGraph (directed graph) to represent the
    hierarchy. Edges go from child to parent (child -> parent).
    """

    def __init__(self):
        self._graph = nx.DiGraph()
        self._terms = {}

    def build_from_data(self, terms_df, edges):
        """
        Build the ontology from a DataFrame of terms and a list of edges.

        Parameters:
        - terms_df: a pandas DataFrame with columns like
          'go_id', 'name', 'namespace', 'definition'
        - edges: a list of tuples like (child_id, parent_id)
          representing parent-child relationships
        """

        # Step 1: Create GOTerm objects from the DataFrame rows. 
        print("Building ontology...")

        # Keep a counter so we can show progress
        count = 0

        for index, row in terms_df.iterrows():
            # Get the values from this row
            go_id = row["go_id"]
            name = row["name"]
            namespace = row["namespace"]
            definition = row["definition"]

            # Check if there is an is_obsolete column
            if "is_obsolete" in row:
                is_obsolete = row["is_obsolete"]
            else:
                is_obsolete = False

            # Create the right type of object based on the namespace.
            if namespace == "biological_process":
                term = BiologicalProcess(go_id, name, definition, is_obsolete)
            elif namespace == "molecular_function":
                term = MolecularFunction(go_id, name, definition, is_obsolete)
            elif namespace == "cellular_component":
                term = CellularComponent(go_id, name, definition, is_obsolete)
            else:
                 print("Warning: unknown namespace '" + str(namespace) + "' for " + str(go_id))
                continue

            # Store the term object in our dictionary
            self._terms[go_id] = term

            # Also add this GO ID as a node in the graph
            self._graph.add_node(go_id)

            # Update the counter
            count = count + 1

        print("  " + str(count) + " terms loaded.")

        # Step 2: Add all edges to the NetworkX graph.
        # Each edge is a tuple (child_id, parent_id).
        # The direction is child -> parent.
        edge_count = 0

        for edge in edges:
            child_id = edge[0]
            parent_id = edge[1]

            # Only add the edge if both terms exist in our dictionary
            if child_id in self._terms and parent_id in self._terms:
                self._graph.add_edge(child_id, parent_id)
                edge_count = edge_count + 1

        print("  " + str(edge_count) + " edges loaded.")

        # Step 3: Set the parents and children lists on each GOTerm object.
         for child_id, parent_id in self._graph.edges():
            # Get the actual GOTerm objects
            child_term = self._terms.get(child_id)
            parent_term = self._terms.get(parent_id)

            # Add the relationship to both term objects
            if child_term is not None and parent_term is not None:
                child_term.add_parent(parent_id)
                parent_term.add_child(child_id)

        print("Building ontology... done!")

    def get_term(self, go_id):
        """
        This function returns the GOTerm object for a given GO ID.
        If the GO ID is not found, it returns None.
        """
        # Use .get() so we don't get a KeyError if the ID doesn't exist
        result = self._terms.get(go_id, None)
        return result

    def get_parents(self, go_id):
        """
        This function returns a list of parent GOTerm objects
        for the given GO ID. Parents are one step up in the hierarchy.
        """
        # First check if the term exists
        term = self._terms.get(go_id)
        if term is None:
            return []

        # Build a list of parent GOTerm objects
        result = []
        for parent_id in term.parents:
            parent_term = self._terms.get(parent_id)
            if parent_term is not None:
                result.append(parent_term)

        return result

    def get_children(self, go_id):
        """
        This function returns a list of child GOTerm objects
        for the given GO ID. Children are one step down in the hierarchy.
        """
        # First check if the term exists
        term = self._terms.get(go_id)
        if term is None:
            return []

        # Build a list of child GOTerm objects
        result = []
        for child_id in term.children:
            child_term = self._terms.get(child_id)
            if child_term is not None:
                result.append(child_term)

        return result

    def get_ancestors(self, go_id):
        """
        This function returns a set of all ancestor GOTerm objects.
        Ancestors are all terms above this one in the hierarchy
        (parents, grandparents, etc.). We use NetworkX to walk
        up the graph from the given node.
        """
        # Check if the node exists in the graph
        if go_id not in self._graph:
            return set()

        # Use networkx to find all nodes reachable from go_id
        ancestor_ids = nx.descendants(self._graph, go_id)

        # Convert the GO IDs to actual GOTerm objects
        result = set()
        for ancestor_id in ancestor_ids:
            ancestor_term = self._terms.get(ancestor_id)
            if ancestor_term is not None:
                result.add(ancestor_term)

        return result

    def get_descendants(self, go_id):
        """
        This function returns a set of all descendant GOTerm objects.
        Descendants are all terms below this one in the hierarchy
        (children, grandchildren, etc.). We use NetworkX to walk
        down the graph from the given node.
        """
        # Check if the node exists in the graph
        if go_id not in self._graph:
            return set()

        # Use networkx ancestors() because in our graph edges go
        descendant_ids = nx.ancestors(self._graph, go_id)

        # Convert the GO IDs to actual GOTerm objects
        result = set()
        for descendant_id in descendant_ids:
            descendant_term = self._terms.get(descendant_id)
            if descendant_term is not None:
                result.add(descendant_term)

        return result

    def get_subgraph(self, go_id, depth=2):
        """
        This function returns a list of GOTerm objects that are
        within N steps of the given GO ID. It looks both up (parents)
        and down (children) in the hierarchy.
        """
        # Check if the node exists in the graph
        if go_id not in self._graph:
            return []

        # visited keeps track of which nodes we have already seen
        visited = set()
        visited.add(go_id)

        # current_level holds the nodes we are currently exploring
        current_level = [go_id]

        # Loop for each depth level
        for step in range(depth):
            # next_level will hold the nodes for the next step
            next_level = []

            for node_id in current_level:

                # Get predecessors (children in the ontology)
                for predecessor in self._graph.predecessors(node_id):
                    if predecessor not in visited:
                        visited.add(predecessor)
                        next_level.append(predecessor)

                # Get successors (parents in the ontology)
                for successor in self._graph.successors(node_id):
                    if successor not in visited:
                        visited.add(successor)
                        next_level.append(successor)

            # Move to the next level
            current_level = next_level

        # Convert all visited GO IDs to GOTerm objects
        result = []
        for node_id in visited:
            term = self._terms.get(node_id)
            if term is not None:
                result.append(term)

        return result

    def search(self, keyword):
        """
        This function searches for GO terms whose name contains
        the given keyword. The search is case insensitive.
        Returns a list of matching GOTerm objects.
        """
        # Convert the keyword to lowercase for case-insensitive matching
        keyword_lower = keyword.lower()

        # Loop through all terms and check if the name contains the keyword
        result = []
        for go_id in self._terms:
            term = self._terms[go_id]
            # Convert the term name to lowercase and check for the keyword
            if keyword_lower in term.name.lower():
                result.append(term)

        return result

    def get_all_terms(self):
        """
        This function returns a list of all GOTerm objects
        stored in the ontology.
        """
        # Get all the values from the dictionary and put them in a list
        result = []
        for go_id in self._terms:
            result.append(self._terms[go_id])
        return result

    def get_terms_by_namespace(self, namespace):
        """
        This function returns a list of GOTerm objects that belong
        to the specified namespace. For example, you can pass
        "biological_process" to get only biological process terms.
        """
        # Loop through all terms and filter by namespace
        result = []
        for go_id in self._terms:
            term = self._terms[go_id]
            if term.namespace == namespace:
                result.append(term)
        return result

    def total_terms(self):
        """
        This function returns the total number of terms in the ontology.
        It just returns the length of the terms dictionary.
        """
        return len(self._terms)

# obo_parser.py
# This module reads and parses Gene Ontology OBO files.
# OBO is a plain-text format used to describe ontology terms and their relationships.
# We parse it line by line ourselves (no external OBO library needed).

import pandas as pd


class OBOParser:
    """
    A parser for OBO-format Gene Ontology files.

    The OBO file contains blocks like [Term] and [Typedef].
    We only care about [Term] blocks. Each term has fields like
    id, name, namespace, def, is_a, and is_obsolete.
    """

    def parse(self, filepath):
        """
        Parse an OBO file and return the GO terms and edges.

        Parameters:
            filepath - path to the .obo file on disk

        Returns:
            A tuple of (terms_dataframe, edges_list)
            - terms_dataframe is a pandas DataFrame with columns:
              [go_id, name, namespace, definition, is_obsolete]
            - edges_list is a list of (child_id, parent_id) tuples
              built from the is_a relationships
        """

        print("Parsing OBO file...")

        # This list will hold one dictionary per GO term we find
        all_terms = []

        # This list will hold (child, parent) tuples for the is_a edges
        all_edges = []

        # Open the file and read all lines into a list
        with open(filepath, "r") as f:
            lines = f.readlines()

        # We need to track whether we are currently inside a [Term] block
        inside_term_block = False

        # We also need to skip [Typedef] blocks
        inside_typedef_block = False

        # Variables to hold the current term's data as we read it
        current_id = None
        current_name = None
        current_namespace = None
        current_definition = None
        current_is_obsolete = False
        current_parents = []  # a term can have multiple is_a parents

        # Loop through each line in the file
        for line in lines:

            # Remove the newline character at the end
            line = line.strip()

            if line == "[Term]":
                
                if inside_term_block and current_id is not None:
                    
                    term_dict = {
                        "go_id": current_id,
                        "name": current_name,
                        "namespace": current_namespace,
                        "definition": current_definition,
                        "is_obsolete": current_is_obsolete,
                    }
                    all_terms.append(term_dict)

                    # Add all is_a edges for this term
                    for parent_id in current_parents:
                        all_edges.append((current_id, parent_id))

                # Reset everything for the new term
                inside_term_block = True
                inside_typedef_block = False
                current_id = None
                current_name = None
                current_namespace = None
                current_definition = None
                current_is_obsolete = False
                current_parents = []
                continue

            # Check if we hit a [Typedef] block — we want to skip these
            if line == "[Typedef]":
                # Save the previous term if we were reading one
                if inside_term_block and current_id is not None:
                    term_dict = {
                        "go_id": current_id,
                        "name": current_name,
                        "namespace": current_namespace,
                        "definition": current_definition,
                        "is_obsolete": current_is_obsolete,
                    }
                    all_terms.append(term_dict)

                    for parent_id in current_parents:
                        all_edges.append((current_id, parent_id))

                # Mark that we are now in a typedef block (skip it)
                inside_term_block = False
                inside_typedef_block = True
                current_id = None
                current_name = None
                current_namespace = None
                current_definition = None
                current_is_obsolete = False
                current_parents = []
                continue

            # If we are inside a Typedef block, skip the line
            if inside_typedef_block:
                continue

            # If we are not inside a Term block, skip the line
            if not inside_term_block:
                continue

            # Now we know we are inside a [Term] block.
            # Parse the different fields we care about.

            # Parse the GO id (e.g., "id: GO:0008150")
            if line.startswith("id:"):
                # Get everything after "id: "
                current_id = line[len("id:"):].strip()

            # Parse the name (e.g., "name: biological_process")
            elif line.startswith("name:"):
                current_name = line[len("name:"):].strip()

            # Parse the namespace (e.g., "namespace: biological_process")
            elif line.startswith("namespace:"):
                current_namespace = line[len("namespace:"):].strip()

            elif line.startswith("def:"):
                # Find the first quote and the second quote
                first_quote = line.find('"')
                second_quote = line.find('"', first_quote + 1)

                if first_quote != -1 and second_quote != -1:
                    # Extract the text between the quotes
                    current_definition = line[first_quote + 1 : second_quote]
                else:
                    # If we can't find quotes, just grab whatever is after "def:"
                    current_definition = line[len("def:"):].strip()

            elif line.startswith("is_a:"):
                # Get the part after "is_a: "
                rest = line[len("is_a:"):].strip()
                parent_id = rest.split("!")[0].strip()

                # Add to the list of parents for this term
                current_parents.append(parent_id)

            # Check if the term is marked as obsolete
            elif line.startswith("is_obsolete:"):
                value = line[len("is_obsolete:"):].strip()
                if value == "true":
                    current_is_obsolete = True

        if inside_term_block and current_id is not None:
            term_dict = {
                "go_id": current_id,
                "name": current_name,
                "namespace": current_namespace,
                "definition": current_definition,
                "is_obsolete": current_is_obsolete,
            }
            all_terms.append(term_dict)

            for parent_id in current_parents:
                all_edges.append((current_id, parent_id))

        terms_dataframe = pd.DataFrame(
            all_terms,
            columns=["go_id", "name", "namespace", "definition", "is_obsolete"],
        )

        # Print out some info so the user knows what happened
        print("Found " + str(len(all_terms)) + " terms")
        print("Found " + str(len(all_edges)) + " is_a edges")

        # Return the dataframe and edges as a tuple
        return (terms_dataframe, all_edges)

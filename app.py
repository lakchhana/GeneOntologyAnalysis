"""
Gene Ontology Analysis System - Flask Web Application
=====================================================
This is the main web application file. It sets up all the routes
(pages) for our Gene Ontology analysis tool.

It uses the data_loader module to load all the real data at startup,
then each route uses the real ontology, annotation, and similarity
objects to answer user requests.
"""

import os
import json
from flask import (
    Flask, render_template, request, redirect,
    url_for, flash, jsonify
)

# Import our data loader and statistics functions
from data_loader import load_all_data
from analysis.statistics import (
    namespace_distribution, evidence_distribution,
    matrix_summary, gene_cosine_similarity,
)

# ---------------------------------------------------------------------------
# App configuration
# ---------------------------------------------------------------------------

app = Flask(__name__)

# Secret key is needed for flash messages to work
app.secret_key = os.environ.get(
    "GO_ANALYSIS_SECRET_KEY",
    "go-analysis-secret-key-change-in-production",
)

# Flask 3 rejects requests for unknown Host headers unless they are trusted.
default_trusted_hosts = "127.0.0.1,localhost,go-project.eshagh.dev"
app.config["TRUSTED_HOSTS"] = [
    host.strip()
    for host in os.environ.get(
        "GO_ANALYSIS_TRUSTED_HOSTS",
        default_trusted_hosts,
    ).split(",")
    if host.strip()
]

# Configure file uploads
UPLOAD_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024  # 50 MB max upload

# Make sure the upload folder exists
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

# Allowed file extensions for upload
ALLOWED_EXTENSIONS = {"obo", "gaf", "gz", "txt"}


def allowed_file(filename):
    """Check if the uploaded file has an allowed extension."""
    # Handle double extensions like .gaf.gz
    lower_name = filename.lower()
    if lower_name.endswith(".gaf.gz"):
        return True
    if "." in lower_name:
        ext = lower_name.rsplit(".", 1)[1]
        return ext in ALLOWED_EXTENSIONS
    return False


# ---------------------------------------------------------------------------
# Load all data when the app starts
# ---------------------------------------------------------------------------
# This runs ONCE when the server starts up. It parses the OBO and GAF files,
# builds the ontology and annotation manager, and creates the similarity
# calculators. Everything is stored in the "data" dictionary.

print("Starting up the Gene Ontology Analysis app...")
print("This may take a minute while we load the data files...")
print()

# Call load_all_data() to get all our objects
data = load_all_data()

# Pull out the individual objects so they are easy to use in routes
ontology = data["ontology"]
annotation_manager = data["annotation_manager"]
jaccard_calculator = data["jaccard"]
overlap_calculator = data["overlap"]

print("App is ready to go!")
print()


# ---------------------------------------------------------------------------
# Helper: map aspect codes to friendly namespace names
# ---------------------------------------------------------------------------
# The GAF file uses single-letter codes: "P", "F", "C"
# We want to show nicer labels like "BP", "MF", "CC" in the web interface.

ASPECT_MAP = {
    "P": "BP",
    "F": "MF",
    "C": "CC",
}


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.route("/")
def index():
    """Homepage -- shows a welcome message and quick stats."""

    # Get the total number of GO terms from the ontology
    total_terms = ontology.total_terms()

    # Get the total number of annotations from the annotation manager
    summary = annotation_manager.summary()
    total_annotations = summary["total"]

    # Get the total number of unique genes
    total_genes = len(annotation_manager.get_all_genes())

    # Render the homepage template with the real numbers
    return render_template(
        "index.html",
        total_terms=total_terms,
        total_annotations=total_annotations,
        total_genes=total_genes,
    )


@app.route("/search/term")
def search_term():
    """Search for GO terms by name or ID."""

    # Get the search query from the URL parameters (?q=something)
    query = request.args.get("q", "").strip()
    results = []

    if query:
        # Check if the query looks like a GO ID (starts with "GO:")
        # If so, try to look it up directly first
        if query.upper().startswith("GO:"):
            # Try to get the exact term by its GO ID
            term_obj = ontology.get_term(query.upper())

            if term_obj is not None:
                # We found an exact match! Put it in the results list
                results.append({
                    "go_id": term_obj.go_id,
                    "name": term_obj.name,
                    "namespace": term_obj.namespace,
                })
            else:
                # The exact GO ID was not found, so try a keyword search
                # (maybe the user typed a partial GO ID)
                matching_terms = ontology.search(query)

                # Convert each GOTerm object to a simple dictionary
                for t in matching_terms:
                    results.append({
                        "go_id": t.go_id,
                        "name": t.name,
                        "namespace": t.namespace,
                    })
        else:
            # It does not look like a GO ID, so search by keyword in name
            matching_terms = ontology.search(query)

            # Convert each GOTerm object to a simple dictionary
            for t in matching_terms:
                results.append({
                    "go_id": t.go_id,
                    "name": t.name,
                    "namespace": t.namespace,
                })

        # Limit the results to 50 so the page does not get too big
        results = results[:50]

    return render_template("search_term.html", query=query, results=results)


@app.route("/term/<go_id>")
def term_detail(go_id):
    """Show details for a specific GO term."""

    # Try to look up the term by its GO ID
    term_obj = ontology.get_term(go_id)

    # If the term was not found, show an error and go back to search
    if term_obj is None:
        flash("Term '" + go_id + "' was not found in the ontology.", "warning")
        return redirect(url_for("search_term"))

    # --- Get parent terms ---
    # ontology.get_parents() returns a list of GOTerm objects
    parent_objects = ontology.get_parents(go_id)

    # Convert parent objects to simple dictionaries for the template
    parents = []
    for p in parent_objects:
        parents.append({
            "go_id": p.go_id,
            "name": p.name,
        })

    # --- Get child terms ---
    # ontology.get_children() returns a list of GOTerm objects
    child_objects = ontology.get_children(go_id)

    # Convert child objects to simple dictionaries for the template
    children = []
    for c in child_objects:
        children.append({
            "go_id": c.go_id,
            "name": c.name,
        })

    # --- Get annotations for this term ---
    # We use a public method on the AnnotationManager to find all genes
    # annotated with this term. We limit to 50 so the page stays fast.
    annotations = annotation_manager.get_annotations_for_term(go_id, limit=50)

    # --- Get neighbourhood (nearby terms within 2 steps) ---
    # ontology.get_subgraph() returns a list of GOTerm objects
    neighbour_objects = ontology.get_subgraph(go_id, depth=2)

    # Convert to simple dictionaries for the template
    neighbourhood = []
    for n in neighbour_objects:
        # Skip the term itself -- we only want neighbors
        if n.go_id != go_id:
            neighbourhood.append({
                "go_id": n.go_id,
                "name": n.name,
            })

    # --- Build graph data for the DAG visualisation ---
    # We create lists of nodes and edges that vis.js can understand.
    # The current term is highlighted in a different colour.
    graph_nodes = []
    graph_edges = []

    # Add the current term as the central node (highlighted)
    graph_nodes.append({
        "id": term_obj.go_id,
        "label": term_obj.name[:30],
        "color": "#0d6efd",
        "font": {"color": "#ffffff"},
    })

    # Add parent nodes and edges (current -> parent)
    for p in parent_objects:
        graph_nodes.append({
            "id": p.go_id,
            "label": p.name[:30],
            "color": "#198754",
        })
        graph_edges.append({
            "from": term_obj.go_id,
            "to": p.go_id,
            "arrows": "to",
            "label": "is_a",
        })

    # Add child nodes and edges (child -> current)
    for c in child_objects:
        graph_nodes.append({
            "id": c.go_id,
            "label": c.name[:30],
            "color": "#fd7e14",
        })
        graph_edges.append({
            "from": c.go_id,
            "to": term_obj.go_id,
            "arrows": "to",
            "label": "is_a",
        })

    # Also add grandparents (parents of parents) for a richer graph
    for p in parent_objects:
        grandparents = ontology.get_parents(p.go_id)
        for gp in grandparents:
            # Only add if not already in nodes
            if not any(n["id"] == gp.go_id for n in graph_nodes):
                graph_nodes.append({
                    "id": gp.go_id,
                    "label": gp.name[:30],
                    "color": "#6c757d",
                })
            graph_edges.append({
                "from": p.go_id,
                "to": gp.go_id,
                "arrows": "to",
                "label": "is_a",
            })

    # Convert to JSON strings for the template
    graph_nodes_json = json.dumps(graph_nodes)
    graph_edges_json = json.dumps(graph_edges)

    # --- Build the term dictionary for the template ---
    term = {
        "go_id": term_obj.go_id,
        "name": term_obj.name,
        "namespace": term_obj.namespace,
        "definition": term_obj.definition if term_obj.definition else "No definition available.",
        "parents": parents,
        "children": children,
        "annotations": annotations,
        "neighbourhood": neighbourhood,
    }

    return render_template(
        "term_detail.html",
        term=term,
        graph_nodes=graph_nodes_json,
        graph_edges=graph_edges_json,
    )


@app.route("/search/gene")
def search_gene():
    """Search for genes by symbol."""

    # Get the search query from the URL parameters
    query = request.args.get("q", "").strip()
    results = []

    if query:
        # Get ALL gene symbols from the annotation manager
        all_genes = annotation_manager.get_all_genes()

        # Filter: keep genes where the query appears in the gene name
        # We convert both to uppercase so the search is case-insensitive
        query_upper = query.upper()
        for gene in all_genes:
            if query_upper in gene.upper():
                results.append(gene)

        # Limit to 50 results so the page stays manageable
        results = results[:50]

    return render_template("search_gene.html", query=query, results=results)


@app.route("/gene/<symbol>")
def gene_detail(symbol):
    """Show all GO annotations for a specific gene."""

    # Get optional namespace filter from the URL parameter
    namespace_filter = request.args.get("namespace", "All")

    # Use a public method on the AnnotationManager to get ALL annotations
    # for this gene. We pass the uppercase symbol because gene symbols
    # are conventionally uppercase. The method also tries uppercase
    # internally if the exact match is not found.
    gene_upper = symbol.upper()
    raw_annotations = annotation_manager.get_annotations_for_gene(gene_upper)

    # Build a list of dictionaries with extra details for the template.
    # We look up the GO term name from the ontology and map the aspect
    # code to a nicer label.
    annotations = []
    for ann in raw_annotations:
        go_id = ann["go_id"]
        evidence_code = ann["evidence_code"]
        aspect = ann["aspect"]

        # Try to get the GO term name from the ontology
        term_obj = ontology.get_term(go_id)

        # If the term exists in the ontology, use its name
        if term_obj is not None:
            term_name = term_obj.name
        else:
            term_name = go_id

        # Map the single-letter aspect code to a nicer label
        aspect_label = ASPECT_MAP.get(aspect, aspect)

        annotations.append({
            "go_id": go_id,
            "term_name": term_name,
            "evidence_code": evidence_code,
            "aspect": aspect_label,
        })

    # Apply the namespace filter if the user picked one
    # (e.g., only show "BP" annotations)
    if namespace_filter != "All":
        annotations = [a for a in annotations if a["aspect"] == namespace_filter]

    # Use the uppercase symbol for display (genes are usually shown uppercase)
    display_symbol = symbol.upper()

    return render_template(
        "gene_detail.html",
        symbol=display_symbol,
        annotations=annotations,
        namespace_filter=namespace_filter,
    )


@app.route("/similarity", methods=["GET", "POST"])
def similarity():
    """Similarity calculator -- compare two genes using different methods."""

    result = None
    gene_a = ""
    gene_b = ""
    method = "jaccard"

    if request.method == "POST":
        # Get form data from the submitted form
        gene_a = request.form.get("gene_a", "").strip()
        gene_b = request.form.get("gene_b", "").strip()
        method = request.form.get("method", "jaccard")

        # Basic validation: make sure both fields are filled in
        if not gene_a or not gene_b:
            flash("Please enter both gene symbols.", "warning")
        else:
            # Convert gene names to uppercase (standard for gene symbols)
            gene_a_upper = gene_a.upper()
            gene_b_upper = gene_b.upper()

            # Check if both genes exist in the annotation manager
            # We check by trying to get their annotations
            terms_a = annotation_manager.get_terms_for_gene(gene_a_upper)
            terms_b = annotation_manager.get_terms_for_gene(gene_b_upper)

            if len(terms_a) == 0:
                flash("Gene '" + gene_a_upper + "' was not found in the annotations.", "warning")
            elif len(terms_b) == 0:
                flash("Gene '" + gene_b_upper + "' was not found in the annotations.", "warning")
            else:
                # Both genes exist! Calculate the similarity score
                if method == "jaccard":
                    score = jaccard_calculator.calculate(gene_a_upper, gene_b_upper)
                else:
                    score = overlap_calculator.calculate(gene_a_upper, gene_b_upper)

                # Build the result dictionary to show on the page
                result = {
                    "gene_a": gene_a_upper,
                    "gene_b": gene_b_upper,
                    "method": method,
                    "score": round(score, 4),
                }

    return render_template(
        "similarity.html",
        result=result,
        gene_a=gene_a,
        gene_b=gene_b,
        method=method,
    )


@app.route("/statistics")
def statistics():
    """Statistics dashboard -- shows charts and tables."""

    # --- Namespace distribution ---
    # This gives us something like: {"BP": 50000, "MF": 20000, "CC": 15000}
    ns_dist = namespace_distribution(annotation_manager)

    # --- Evidence code distribution ---
    # This gives us something like: {"IEA": 400000, "IDA": 50000, ...}
    ev_dist = evidence_distribution(annotation_manager)

    # --- Top genes by annotation count ---
    # get_top_genes returns a list of tuples: [(gene, count), ...]
    top_genes_raw = annotation_manager.get_top_genes(10)

    # Convert to a list of dictionaries for the template
    top_genes = []
    for gene_symbol, count in top_genes_raw:
        top_genes.append({
            "symbol": gene_symbol,
            "count": count,
        })

    # --- Top terms by annotation count ---
    # get_top_terms returns a list of tuples: [(go_id, count), ...]
    top_terms_raw = annotation_manager.get_top_terms(10)

    # Convert to a list of dictionaries, and look up the term name
    top_terms = []
    for go_id, count in top_terms_raw:
        # Try to get the term name from the ontology
        term_obj = ontology.get_term(go_id)
        if term_obj is not None:
            term_name = term_obj.name
        else:
            term_name = go_id  # fallback: just show the GO ID

        top_terms.append({
            "go_id": go_id,
            "name": term_name,
            "count": count,
        })

    # --- Gene-term matrix summary (uses NumPy) ---
    # Build the binary matrix for the top 20 genes and get statistics
    mat_info = matrix_summary(annotation_manager, top_n=20)

    # --- Cosine similarity between top genes (uses NumPy) ---
    # Compute pairwise cosine similarity for the top 10 genes
    cosine_pairs = gene_cosine_similarity(annotation_manager, top_n=10)

    # Take only the top 10 most similar pairs for the template
    top_cosine_pairs = []
    for gene_a, gene_b, score in cosine_pairs[:10]:
        top_cosine_pairs.append({
            "gene_a": gene_a,
            "gene_b": gene_b,
            "score": score,
        })

    # --- Build the stats dictionary that the template expects ---
    stats = {
        "namespace_counts": ns_dist,
        "evidence_counts": ev_dist,
        "top_genes": top_genes,
        "top_terms": top_terms,
        "matrix": mat_info,
        "cosine_pairs": top_cosine_pairs,
    }

    # Convert data to JSON strings so JavaScript (Chart.js) can use them
    namespace_labels = json.dumps(list(ns_dist.keys()))
    namespace_values = json.dumps(list(ns_dist.values()))
    evidence_labels = json.dumps(list(ev_dist.keys()))
    evidence_values = json.dumps(list(ev_dist.values()))

    return render_template(
        "statistics.html",
        stats=stats,
        namespace_labels=namespace_labels,
        namespace_values=namespace_values,
        evidence_labels=evidence_labels,
        evidence_values=evidence_values,
    )


@app.route("/upload", methods=["GET", "POST"])
def upload():
    """Upload OBO or GAF files and optionally reload the data."""

    # We need to declare these as global so we can reassign them
    # when the user clicks the "Reload Data" button.
    global ontology, annotation_manager, jaccard_calculator, overlap_calculator

    if request.method == "POST":

        # Check if the user clicked the "Reload Data" button
        # (the reload button has a hidden field named "action" with value "reload")
        action = request.form.get("action", "")

        if action == "reload":
            # Reload all data from the files in the data/ folder
            try:
                new_data = load_all_data()
                ontology = new_data["ontology"]
                annotation_manager = new_data["annotation_manager"]
                jaccard_calculator = new_data["jaccard"]
                overlap_calculator = new_data["overlap"]
                flash("Data reloaded successfully! The system is now using the new files.", "success")
            except FileNotFoundError as e:
                flash("Could not reload: " + str(e), "danger")
            except Exception as e:
                flash("Error while reloading data: " + str(e), "danger")
            return redirect(url_for("upload"))

        # Otherwise, the user is uploading a file
        # Check if a file was submitted
        if "file" not in request.files:
            flash("No file selected.", "danger")
            return redirect(url_for("upload"))

        file = request.files["file"]

        # Check if the user actually selected a file
        if file.filename == "":
            flash("No file selected.", "danger")
            return redirect(url_for("upload"))

        # Check if the file type is allowed
        if file and allowed_file(file.filename):
            # Save the file to the data/ folder
            filename = file.filename
            save_path = os.path.join(app.config["UPLOAD_FOLDER"], filename)
            file.save(save_path)
            flash("File '" + filename + "' uploaded successfully! "
                  "Click 'Reload Data' below to start using the new file.", "success")
            return redirect(url_for("upload"))
        else:
            flash("Invalid file type. Please upload .obo, .gaf, or .txt files.", "danger")
            return redirect(url_for("upload"))

    # GET request -- just show the upload form
    return render_template("upload.html")


# ---------------------------------------------------------------------------
# Run the app
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Debug mode is ON during development so we see errors in the browser
    app.run(debug=True, port=5000)

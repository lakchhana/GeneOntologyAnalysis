# Gene Ontology Analysis System

A Python application for exploring and analysing the Gene Ontology (GO) and human gene annotations. Built as a project for the Advanced Programming course (2025/2026).

## What This Project Does

- **Parses** the Gene Ontology file (OBO format) and human annotation file (GAF format)
- **Represents** GO terms as objects using an inheritance-based class hierarchy
- **Analyses** gene similarity using Jaccard and Overlap coefficients
- **Displays** everything through a web interface built with Flask

## Requirements

- Python 3.10 or higher
- pip (Python package manager)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/lakchhana/GeneOntologyAnalysis.git
cd GeneOntologyAnalysis
```

2. Install the required packages:
```bash
pip install -r requirements.txt
```

3. Download the data files into the `data/` folder:

- **GO Ontology file**: Download `go-basic.obo` from
  https://current.geneontology.org/ontology/go-basic.obo

- **Human Annotation file**: Download `goa_human.gaf.gz` from
  https://current.geneontology.org/annotations/goa_human.gaf.gz

You can use these commands:
```bash
mkdir -p data
curl -L -o data/go-basic.obo https://current.geneontology.org/ontology/go-basic.obo
curl -L -o data/goa_human.gaf.gz https://current.geneontology.org/annotations/goa_human.gaf.gz
```

## How to Run

Start the Flask web application:
```bash
python3 app.py
```

The first time you run it, it will take about 1 minute to load all the data. After that, open your browser and go to:

```
http://localhost:5000
```

## Features

- **Search GO Terms**: Search by GO ID (e.g., `GO:0006915`) or keyword (e.g., `apoptosis`)
- **View Term Details**: See term definition, parent/child terms, and associated genes
- **Search Genes**: Search by gene symbol (e.g., `TP53`, `BRCA1`)
- **View Gene Annotations**: See all GO annotations for a specific gene
- **Similarity Calculator**: Compare two genes using Jaccard or Overlap similarity
- **Statistics Dashboard**: Bar charts and tables showing annotation distributions
- **File Upload**: Upload custom OBO or GAF files

## Project Structure

```
GeneOntologyAnalysis/
├── app.py                  # Flask web application (main entry point)
├── data_loader.py          # Loads and initializes all data
├── requirements.txt        # Python dependencies
├── data/                   # Data files (OBO and GAF)
├── parser/                 # Module for reading OBO and GAF files
│   ├── obo_parser.py       # Hand-written OBO file parser
│   └── gaf_parser.py       # GAF file parser using pandas
├── ontology/               # Module for GO term classes and graph
│   ├── entity.py           # OntologyEntity, GOTerm, BP, MF, CC classes
│   └── ontology.py         # Ontology class with NetworkX graph
├── annotation/             # Module for gene annotation classes
│   └── annotation.py       # Annotation classes and AnnotationManager
├── analysis/               # Module for similarity and statistics
│   ├── similarity.py       # Jaccard and Overlap similarity classes
│   └── statistics.py       # NumPy matrix and distribution functions
├── templates/              # HTML templates for Flask
├── static/                 # CSS files
└── report/                 # Project report (PDF)
```

## Technologies Used

| Library | Purpose |
|---------|---------|
| Flask | Web framework for the user interface |
| Pandas | Reading and organizing tabular data |
| NumPy | Numerical operations and matrices |
| NetworkX | Graph structure for the GO hierarchy |
| Bootstrap 5 | CSS framework for styling the web pages |
| Chart.js | JavaScript library for charts on the statistics page |

## Data

- **Gene Ontology (go-basic.obo)**: Contains ~48,000 GO terms organized as a directed acyclic graph (DAG)
- **Human Annotations (goa_human.gaf.gz)**: Contains ~937,000 annotations linking human genes to GO terms

## Author
Lakchhana Gurung
Sumeyye Uzun 
Pegah Mohebbi
Advanced Programming 2025/2026
University of Bologna

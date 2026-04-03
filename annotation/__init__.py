# annotation/__init__.py
# This file makes the annotation folder a Python package.
# We import the main classes so they can be used easily from outside.

from annotation.annotation import Annotation
from annotation.annotation import ExperimentalAnnotation
from annotation.annotation import ElectronicAnnotation
from annotation.annotation import AnnotationManager
from annotation.annotation import create_annotation

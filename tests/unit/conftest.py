import importlib.util
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN = os.path.join(REPO, "bin")
FIX = os.path.join(REPO, "tests", "fixtures")
SCHEMA = os.path.join(REPO, "assets", "schema", "cohort.schema.json")
sys.path.insert(0, BIN)


def load(name):
    """Import a bin/ script as a module."""
    spec = importlib.util.spec_from_file_location(name, os.path.join(BIN, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

from pathlib import Path

# Repository root
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Project directories
FILTERS_DIR = PROJECT_ROOT / "filters"
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# Stellar atmosphere model grid
MODEL_GRID_DIR = Path(
    "/home/luis/pysynphot_models/trds/grid/ck04models"
)

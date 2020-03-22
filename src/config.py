"""Defines various global variables, such as the project root file path and a dictionary of amino acids."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODEL_DIR = PROJECT_ROOT / "models/models"
TENSORBOARD_DIR = PROJECT_ROOT / "models/tensorboard_logs"

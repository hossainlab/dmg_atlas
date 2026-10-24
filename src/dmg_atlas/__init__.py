"""dmg_atlas — reproducible code for the human H3K27M DMG single-cell atlas."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
CONFIGS = ROOT / "configs"
LOGS = ROOT / "logs"

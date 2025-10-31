from pathlib import Path

import seaborn as sns

ASSIGNMENT_DIR = Path(__file__).resolve().parent
REPORT_DIR = ASSIGNMENT_DIR / "report"
OUTPUT_DIR = REPORT_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

EXPORT_DPI = 600
PAPER_WIDTH_IN = 8.0

def setup_plotting():
    sns.set_theme(style="whitegrid", context="paper")

"""q4_ftd_verify.py -- run the declining-phase 2-D dynamo's own eight verification checks, unchanged.

The solver (sc24_final/ftd_core.py) and its checks (sc24_final/ftd_verify.py) are imported read-only from the
declining-phase project; only the log is redirected here so nothing in that project is overwritten."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
FTD_DIR = os.environ.get("FTD_DIR", r"C:\Users\Shivansh Mishra\OneDrive\Desktop\declining phase paper\sc24_final")
sys.path.insert(0, FTD_DIR)
import ftd_verify as V                                   # noqa: E402

V.OUT = os.path.join(HERE, "out", "ftd_verify")
sys.exit(0 if V.main() else 1)

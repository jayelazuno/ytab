#!/usr/bin/env python3
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[2] / "src"))
from ytab.screen_coverage.runner import main

if __name__ == "__main__":
    main()

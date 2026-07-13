#!/usr/bin/env python3
"""Benchmark entry point for the PSI evidence gate; nonzero is intentional today."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def main() -> None:
    result = subprocess.run([sys.executable, str(HERE / "audit.py")], check=False)
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()

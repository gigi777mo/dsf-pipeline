"""Console entry point: dsf-run"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> None:
    """Run scripts/run_dsf.py from the installed or source tree."""
    # Prefer packaged sibling scripts if present; else repo scripts/
    here = Path(__file__).resolve()
    candidates = [
        here.parents[2] / "scripts" / "run_dsf.py",  # src/dsf_pipeline -> repo root
        here.parents[1] / "scripts" / "run_dsf.py",
        Path.cwd() / "scripts" / "run_dsf.py",
    ]
    script = next((p for p in candidates if p.exists()), None)
    if script is None:
        print(
            "dsf-run: could not find scripts/run_dsf.py.\n"
            "Clone the repo and run from the project root, or:\n"
            "  python scripts/run_dsf.py --melt your.csv --out results/",
            file=sys.stderr,
        )
        sys.exit(1)
    # Ensure scripts/ is on path for io_melt / tm_calc imports
    sys.path.insert(0, str(script.parent))
    if argv is not None:
        sys.argv = [str(script), *argv]
    else:
        sys.argv[0] = str(script)
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()

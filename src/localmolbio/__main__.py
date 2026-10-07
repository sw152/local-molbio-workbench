"""python -m localmolbio uses the same explicit entry point as the molbio command."""
from .cli import main

raise SystemExit(main())

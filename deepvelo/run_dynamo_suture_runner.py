#!/usr/bin/env python3
from pathlib import Path
import run_dynamo_analysis as rda

rda.OUTPUT_DIR = rda.ROOT / "deepvelo"
rda.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

input_path = rda.ROOT / "scvelo" / "E17_suture_velocity.h5ad"
if not input_path.exists():
    print(f"Input file not found: {input_path}")
else:
    rda.process_file(input_path)

# SutureRNAVelocity

RNA velocity analysis of mouse cranial suture scRNA-seq data (E17) from [GSE163693](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163693).

## Data

- **Source**: GEO GSE163693 / ENA
- **Samples**: E17 mouse cranial suture, 3 batches (SRR13288596–SRR13288598)
- **Annotation**: Gencode vM25 (GRCm38/mm10)

## Setup

```bash
# Install dependencies with uv
uv sync

# Install velocyto separately (requires gcc for OpenMP support on macOS)
uv pip install cython pysam
CC=/opt/homebrew/bin/gcc-15 uv pip install --no-build-isolation velocyto
```

## Workflow

See [RNA_velocity_workflow.md](RNA_velocity_workflow.md) for the full step-by-step guide.

1. **Download BAM files** from ENA (`aria2c_download.txt`)
2. **Generate loom files** with velocyto (`run_velocyto.sh`)
3. **Run RNA velocity analysis** with scVelo (`run_scvelo.py`)

```bash
# Step 1: Download BAMs
aria2c --input-file=aria2c_download.txt --max-concurrent-downloads=3 --split=8 --max-connection-per-server=8 --continue=true

# Step 2: Run velocyto (requires BAM files + gencode.vM25.nochr.gtf)
bash run_velocyto.sh

# Step 3: Run scVelo analysis
uv run python run_scvelo.py
```

## Output

- `E17_suture_velocity.h5ad` — Annotated data with velocity estimates
- `figures/` — UMAP, velocity stream, latent time plots
- `velocity_driver_genes.csv` — Top velocity driver genes per cluster

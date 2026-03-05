#!/bin/bash
set -euo pipefail

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

ANNOTATION="gencode.vM25.nochr.gtf"
BARCODES="E17_barcodes.tsv"
OUTDIR="loom_output"

mkdir -p "$OUTDIR"

echo "=== [$(date)] Skipping batch_1 (already done) ==="

echo "=== [$(date)] Starting velocyto batch_2 ==="
velocyto run \
  -b "$BARCODES" \
  -o "$OUTDIR" \
  --samtools-threads 4 \
  downloads/SRR13288597/E17_batch_2_possorted_genome_bam.bam \
  "$ANNOTATION"
echo "=== [$(date)] batch_2 done ==="

echo "=== [$(date)] Starting velocyto batch_3 ==="
velocyto run \
  -b "$BARCODES" \
  -o "$OUTDIR" \
  --samtools-threads 4 \
  downloads/SRR13288598/E17_batch_3_possorted_genome_bam.bam \
  "$ANNOTATION"
echo "=== [$(date)] batch_3 done ==="

echo "=== [$(date)] All velocyto runs complete ==="

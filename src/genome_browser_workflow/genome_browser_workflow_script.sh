#!/bin/bash
#SBATCH --job-name=genome_browser_brdu
#SBATCH --output=genome_browser_brdu_%j.out
#SBATCH --error=genome_browser_brdu_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=normal

module load python/3.13.7
module load samtools
module load modkit

# Directory structure
# workflows/
# ├── data/
# │   ├── bam/
# │   ├── bed/
# │   ├── bedgraph/
# │   ├── sorted_bam/
# │   ├── index_sorted_bam_bai/
# │   └── liftover_chains/
# ├── src/
# │   ├── rainplot_workflow/
# │   ├── genome_browser_workflow/
# │   └── utils/
# ├── tools/
# └── results/
#     ├── rainplot_results/
#     └── genome_browser_results/

WORKFLOW_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

BAM_DIR="$WORKFLOW_ROOT/data/bam"
BEDGRAPH_DIR="$WORKFLOW_ROOT/data/bedgraph"
SORTED_BAM_DIR="$WORKFLOW_ROOT/data/sorted_bam"
INDEX_DIR="$WORKFLOW_ROOT/data/index_sorted_bam_bai"
SRC_DIR="$WORKFLOW_ROOT/src/genome_browser_workflow"
RESULTS_DIR="$WORKFLOW_ROOT/results/genome_browser_results"

mkdir -p "$BAM_DIR" "$BEDGRAPH_DIR" "$SORTED_BAM_DIR" "$INDEX_DIR" "$RESULTS_DIR"

# Usage: bash src/genome_browser_workflow/genome_browser_workflow_script.sh BAM [output_prefix]
BAM_FILE=$1
OUTPUT_PREFIX=$2

if [ -z "$BAM_FILE" ]; then
    echo "[ERROR] Missing required BAM file argument."
    echo "Usage: bash src/genome_browser_workflow/genome_browser_workflow_script.sh BAM [output_prefix]"
    exit 1
fi

BAM="$BAM_DIR/$BAM_FILE"

if [ ! -f "$BAM" ]; then
    echo "[ERROR] BAM file not found: $BAM"
    exit 1
fi

# Virtual environment setup
VENV_DIR="$SRC_DIR/.genome_browser_env"

if [ ! -d "$VENV_DIR" ]; then
    echo "[INFO] Virtual environment not found. Creating $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
    if [ $? -ne 0 ]; then
        echo "[ERROR] Failed to create virtual environment. Exiting."
        exit 1
    fi
    echo "[INFO] Virtual environment created."
fi

echo "[INFO] Activating virtual environment..."
source "$VENV_DIR/bin/activate"
if [ $? -ne 0 ]; then
    echo "[ERROR] Failed to activate virtual environment. Exiting."
    exit 1
fi

echo "[INFO] Installing required libraries from requirements.txt..."
pip install -r "$SRC_DIR/requirements.txt" --quiet
if [ $? -ne 0 ]; then
    echo "[ERROR] pip install failed. Exiting."
    deactivate
    exit 1
fi
echo "[INFO] Libraries installed successfully."

# Raw data extraction
CMD="python $SRC_DIR/raw_data_extraction_on_bam.py $BAM"

if [ -n "$OUTPUT_PREFIX" ]; then
    CMD="$CMD -o $OUTPUT_PREFIX"
fi

echo "[INFO] Running genome browser raw data extraction..."
eval "$CMD"

# Check expected output files
if [ -n "$OUTPUT_PREFIX" ]; then
    POSITIVE_OUTPUT="$BEDGRAPH_DIR/${OUTPUT_PREFIX}_positive.bedgraph"
    NEGATIVE_OUTPUT="$BEDGRAPH_DIR/${OUTPUT_PREFIX}_negative.bedgraph"
else
    BAM_BASENAME="$(basename "$BAM_FILE" .bam)"
    POSITIVE_OUTPUT="$BEDGRAPH_DIR/${BAM_BASENAME}_positive.bedgraph"
    NEGATIVE_OUTPUT="$BEDGRAPH_DIR/${BAM_BASENAME}_negative.bedgraph"
fi

if [ ! -s "$POSITIVE_OUTPUT" ]; then
    echo "[ERROR] Positive bedgraph file is empty or was not created: $POSITIVE_OUTPUT"
    deactivate
    exit 1
fi

if [ ! -s "$NEGATIVE_OUTPUT" ]; then
    echo "[ERROR] Negative bedgraph file is empty or was not created: $NEGATIVE_OUTPUT"
    deactivate
    exit 1
fi

echo "[INFO] Positive strand bedgraph created: $POSITIVE_OUTPUT"
echo "[INFO] Negative strand bedgraph created: $NEGATIVE_OUTPUT"
echo "[INFO] Genome browser raw extraction complete."

deactivate
echo "[INFO] Virtual environment deactivated."

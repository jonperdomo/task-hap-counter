#!/bin/bash
#SBATCH --job-name=SNVTestData
#SBATCH --mem=50G
#SBATCH --cpus-per-task=6
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%j-%x.out

set -e  # Exit on error

# Define paths
if [ -n "${SLURM_SUBMIT_DIR}" ]; then
    SCRIPT_DIR=$SLURM_SUBMIT_DIR
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
TEST_DATA_DIR="$SCRIPT_DIR/test_data"
BAM_FILE="$TEST_DATA_DIR/giab_2023.05.hg002.haplotagged.chr16_28000000_29000000.processed.30x.bam"
VCF_FILE="$TEST_DATA_DIR/giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.vcf.gz"
CLI_PY="$SCRIPT_DIR/src/hap_counter/cli.py"
OUTFILE="$TEST_DATA_DIR/giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.snv_counts.tsv"

# Run cli.py with BAM and VCF from test_data
echo "[debug] Running cli.py with BAM and VCF files"
python3 "$CLI_PY" --bam "$BAM_FILE" --vcf "$VCF_FILE" --output "$OUTFILE"
echo "[debug] Output saved to $OUTFILE"

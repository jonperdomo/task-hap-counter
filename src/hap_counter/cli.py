import argparse
import sys
from hap_counter import count_snvs, save_counts_to_tsv


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(description="Haplotype counter CLI")
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--vcf', required=True, help='Input VCF.GZ file')
    parser.add_argument('--output', default='snv_counts.tsv', help='Output TSV file')
    args = parser.parse_args()

    print("[debug] Calling count_snvs with BAM and VCF files")
    counts = count_snvs(args.bam, args.vcf)

    output_file = args.output
    save_counts_to_tsv(counts, output_file)


if __name__ == "__main__":
    main()

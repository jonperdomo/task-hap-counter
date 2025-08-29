import sys

def save_counts_to_tsv(counts, output_file):
    """
    Save SNV counts to a TSV file.
    Format:
    chrom    pos    h1_REF  h2_REF  h1_ALT  h2_ALT  REF ALT
    where h1 and h2 are the haplotypes (1-based) and ALT/REF are the alleles.
    h<H>_<A> = Number of supporting reads for the given haplotype and allele.
    """
    with open(output_file, 'w') as f:
        f.write("chrom\tpos\th1_REF\th1_ALT\th2_REF\th2_ALT\tREF\tALT\n")
        for chrom in counts:
            for key, value in counts[chrom].items():
                pos1, ref, alt = key
                ref1_count = value["ref1"]
                alt1_count = value["alt1"]
                ref2_count = value["ref2"]
                alt2_count = value["alt2"]

                f.write(f"{chrom}\t{pos1}\t{ref1_count}\t{alt1_count}\t{ref2_count}\t{alt2_count}\t{ref}\t{alt}\n")

    print(f"[debug] SNV counts saved to {output_file}", file=sys.stderr)

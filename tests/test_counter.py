import os
import pysam

from hap_counter import count_snvs, save_counts_to_tsv

def create_synthetic_bam(ref_seq, reads):
    """
    Create a synthetic BAM file from the reference sequence and a list of reads.
    """
    test_outdir = os.path.join(os.getcwd(), "tests/test_output")
    if not os.path.exists(test_outdir):
        os.makedirs(test_outdir)

    read_count = 0
    bam_path = os.path.join(test_outdir, "synthetic.bam")
    with pysam.AlignmentFile(bam_path, "wb", reference_names=["chr1"], reference_lengths=[len(ref_seq)]) as outf:
        for read in reads:
            # High-quality primary alignment
            a = pysam.AlignedSegment()
            a.query_name = f"read{read_count}"
            a.query_sequence = read['seq']
            a.flag = 0  # Primary alignment
            a.reference_id = 0
            a.reference_start = 0
            a.mapping_quality = 60
            a.query_qualities = pysam.qualitystring_to_array("I" * len(read['seq']))
            a.cigar = [(0, len(read['seq']))]

            # Include haplotype information
            if read.get('haplotype') is not None:
                a.set_tag("HP", read['haplotype'], value_type='i')

            outf.write(a)
            read_count += 1

    pysam.index(bam_path)
    return bam_path

def create_snv_read(ref_seq, position, alt=False, haplotype=None, alt_base=None):
    """Create a read with a single SNV at the specified position."""
    read = {
        'seq': ref_seq,
        'snv_position': position,
        'reference': ref_seq[position],
        'alt': ref_seq[position],
        'haplotype': haplotype
    }

    if alt:
        # Print an error if alt base is not specified
        if alt_base is None:
            raise ValueError("alt_base must be specified when alt is True")

        # Print an error if alt base equals reference base
        if alt_base == read['reference']:
            raise ValueError(f"Alt base equals reference base for SNV at position {position}")

        read['seq'] = read['seq'][:position] + alt_base + read['seq'][position + 1:]
        print("Added SNV read at position", position, " ref", read['reference'], " alt", alt_base, ", seq", read['seq'], "haplotype", haplotype)

    return read


def test_snv_counting_no_haplotype(ref_seq, snv_vcf, snv_positions):
    """Test bi-allelic SNV counting without haplotype information."""

    # Loop through SNV positions and create reads
    reads = []
    snv_counts = {}
    alt_allele_count = 2
    ref_allele_count = 2
    for pos, ref, alt in snv_positions:
        # Create read allele support
        for _ in range(alt_allele_count):
            reads.append(create_snv_read(ref_seq, pos-1, alt=True, haplotype=None, alt_base=alt))
        for _ in range(ref_allele_count):
            reads.append(create_snv_read(ref_seq, pos-1, alt=False, haplotype=None))

        snv_counts[(pos, ref, alt, None)] = (alt_allele_count, ref_allele_count)
        print(f"[debug] Created {ref_allele_count} ref and {alt_allele_count} alt reads for SNV at position {pos} ({ref}>{alt})")

    print("Summary of reads with SNV counts (no haplotype):")
    total_snv_count = len(snv_counts)
    total_ref_count = (ref_allele_count + alt_allele_count) * total_snv_count - alt_allele_count  # = Total ref. alleles - SNV alt. alleles
    for key, (alt_count, ref_count) in snv_counts.items():
        print(f"  {key}: ALT={alt_count}, REF={total_ref_count}")

    bam_path = create_synthetic_bam(ref_seq, reads)
    print("[debug] created synthetic BAM:", bam_path)
    try:
        snv_counts = count_snvs(bam_path, snv_vcf)

        # Save to TSV
        test_outdir = os.path.join(os.getcwd(), "tests/test_output")
        output_file = os.path.join(test_outdir, "snv_counts_no_haplotypes.tsv")
        save_counts_to_tsv(snv_counts, output_file)

        # Print the results
        snv_counts = snv_counts.get("chr1", {})
        print("SNV counts (no haplotype tags):")
        for key, count in snv_counts.items():
            print(f"  {key}: {count}")

        assert snv_counts[(1, 'A', 'T')] == {"ref1": 10, "alt1": 2, "ref2": 0, "alt2": 0}
        assert snv_counts[(6, 'C', 'G')] == {"ref1": 10, "alt1": 2, "ref2": 0, "alt2": 0}
        assert snv_counts[(12, 'T', 'C')] == {"ref1": 10, "alt1": 2, "ref2": 0, "alt2": 0}
    finally:
        os.remove(bam_path)
        os.remove(bam_path + ".bai")

def test_snv_counting_with_haplotype(ref_seq, snv_vcf, snv_positions):
    """Test bi-allelic SNV counting with haplotype information."""

    # Loop through SNV positions and create reads
    reads = []
    alt_allele_count = 2
    ref_allele_count = 2

    # Create read support for each haplotype
    for pos, ref, alt in snv_positions:
        for haplotype in [1, 2]:
            # Create read allele support
            for _ in range(alt_allele_count):
                reads.append(create_snv_read(ref_seq, pos-1, alt=True, haplotype=haplotype, alt_base=alt))

            # Create ref allele support
            for _ in range(ref_allele_count):
                reads.append(create_snv_read(ref_seq, pos-1, alt=False, haplotype=haplotype))

            print(f"[debug] Created {alt_allele_count} alt reads for SNV at position {pos} ({ref}>{alt}), haplotype {haplotype}")

    bam_path = create_synthetic_bam(ref_seq, reads)
    print("[debug] created synthetic BAM:", bam_path)
    try:
        snv_counts = count_snvs(bam_path, snv_vcf)

        # Save to TSV
        test_outdir = os.path.join(os.getcwd(), "tests/test_output")
        output_file = os.path.join(test_outdir, "snv_counts_with_haplotypes.tsv")
        save_counts_to_tsv(snv_counts, output_file)

        # Print the results
        snv_counts = snv_counts.get("chr1", {})
        print("SNV counts (with haplotype tags):")
        for key, count in snv_counts.items():
            print(f"  {key}: {count}")

        assert snv_counts[(1, 'A', 'T')] == {"ref1": 10, "alt1": 2, "ref2": 10, "alt2": 2}
        assert snv_counts[(6, 'C', 'G')] == {"ref1": 10, "alt1": 2, "ref2": 10, "alt2": 2}
        assert snv_counts[(12, 'T', 'C')] == {"ref1": 10, "alt1": 2, "ref2": 10, "alt2": 2}

    finally:
        os.remove(bam_path)
        os.remove(bam_path + ".bai")

import pysam
import pytest
import tempfile
import os


print("Testing environment variables")
import os, sys; print("PYTHONPATH=", os.getenv("PYTHONPATH")); print(sys.path[:5])


from hap_counter import count_snvs

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
        print("Added SNV read at position", position, "with ref", read['reference'], "and alt", alt_base, ", and seq", read['seq'])

    return read


def test_snv_counting_no_haplotype(ref_seq, snv_vcf, snv_positions):
    """Test bi-allelic SNV counting without haplotype information."""
    # Create a VCF file with the SNVs
    test_outdir = os.path.join(os.getcwd(), "tests/test_output")
    if not os.path.exists(test_outdir):
        os.makedirs(test_outdir)

    reads = []

    # Loop through SNV positions and create reads
    snv_counts = {}
    alt_allele_count = 2
    ref_allele_count = 2
    for pos, ref, alt in snv_positions:
        # Create read allele support
        for i in range(alt_allele_count):
            reads.append(create_snv_read(ref_seq, pos-1, alt=True, haplotype=None, alt_base=alt))
        for j in range(ref_allele_count):
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
        snv_counts = count_snvs(bam_path, snv_vcf, "chr1")

        # Print the results
        print("SNV counts (no haplotype):")
        for key, count in snv_counts.items():
            print(f"  {key}: {count}")

        # Create a snv_counts with the solution for testing
        snv_counts = {
            (4, 'A', 'T', None): 2,
            (9, 'C', 'G', None): 1
        }

        assert snv_counts[(4, 'A', 'T', None)] == 2
        assert snv_counts[(9, 'C', 'G', None)] == 1
    finally:
        os.remove(bam_path)
        os.remove(bam_path + ".bai")

# def test_snv_counting_with_haplotype(ref_seq, snv_vcf):
#     """Test bi-allelic SNV counting with haplotype information."""

#     reads = []

#     # Create 4 alt alleles in position 4 (A>T)
#     for i in range(4):
#         reads.append(create_snv_read(ref_seq, 4, alt=True, haplotype=1))

#     # Create 5 alt alleles in position 9 (C>G)
#     for i in range(5):
#         reads.append(create_snv_read(ref_seq, 9, alt=True, haplotype=2))

#     # Create 3 alt alleles in position 11 (T>G)
#     for i in range(3):
#         reads.append(create_snv_read(ref_seq, 11, alt=True, haplotype=1))

#     # Create 8 ref reads
#     for i in range(4):
#         reads.append(create_snv_read(ref_seq, 0, alt=False, haplotype=1))
#     for i in range(4):
#         reads.append(create_snv_read(ref_seq, 0, alt=False, haplotype=2))

#     bam_path = create_synthetic_bam(ref_seq, reads)
#     try:
#         # snv_counts = count_snvs(bam_path, fasta_path, "chr1")
#         # Create a snv_counts with the solution for testing
#         snv_counts = {
#             (4, 'A', 'T', 1): 1,
#             (9, 'C', 'T', 2): 2
#         }

#         # Example assertion, adjust according to your count_snvs implementation
#         assert snv_counts[(4, 'A', 'T', 1)] == 1
#         assert snv_counts[(9, 'C', 'T', 2)] == 2
#     finally:
#         os.remove(bam_path)
#         os.remove(bam_path + ".bai")

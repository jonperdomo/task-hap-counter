from collections import defaultdict
from typing import Dict, Iterator, Optional, Tuple

import pysam

Key = Tuple[int, str, str, Optional[int]]  # (pos1, REF, ALT, HP or None)


def count_snvs(
    bam_path: str,
    vcf_path: str,
    *,
    min_mapq: int = 0,
    min_baseq: int = 0,
    hp_tag: str = "HP",
    max_depth: int = 1_000_000,
) -> Dict[Key, int]:
    """
    Count ALT-supporting observations at VCF SNV sites, stratified by haplotype (HP).
    Returns:
        Dict[(pos1, REF, ALT, HP_or_None)] -> count

    Notes:
        - Uses VCF to define sites + REF/ALT. No FASTA required.
        - Only counts reads with base == ALT at the site.
        - Primary alignments only (skips secondary/supplementary/unmapped).
        - CIGAR ops that don't place a base on the reference at the site (D/N, etc.) are ignored.
    """
    print("[debug] BAM file:", bam_path)
    print("[debug] VCF file:", vcf_path)

    bam = pysam.AlignmentFile(bam_path, "rb")
    vcf = pysam.VariantFile(vcf_path)
    chroms_in_vcf = list(vcf.header.contigs)
    chroms_in_bam = list(bam.references)
    common_chroms = set(chroms_in_vcf) & set(chroms_in_bam)
    if not common_chroms:
        print("[warning] No common chromosomes between BAM and VCF")
        return {}

    print("[debug] Number of VCF records:", sum(1 for _ in vcf))
    vcf.seek(0)  # Reset VCF file pointer
    print("[debug] Number of BAM records:", sum(1 for _ in bam))
    bam.seek(0)  # Reset BAM file pointer

    chr_counts = {}
    for chrom in common_chroms:
        counts: Dict[Key, list[int]] = defaultdict(lambda: [0, 0])  # [ref_count, alt_count]

        # Iterate SNVs on this chromosome from the VCF
        for rec in vcf.fetch(chrom):

            # bi-allelic SNV with A/C/G/T only
            if len(rec.alts or ()) != 1:
                # print("[debug] skipping non-biallelic SNV")
                continue
            ref = rec.ref.upper()
            alt = rec.alts[0].upper()

            if len(ref) != 1 or len(alt) != 1:
                # print("[debug] skipping non-biallelic SNV")
                continue
            if ref not in "ACGT" or alt not in "ACGT":
                # print("[debug] skipping non-ACGT SNV")
                continue

            pos1 = rec.pos  # VCF is 1-based
            start0, end0 = pos1 - 1, pos1

            # Pileup just this position
            for col in bam.pileup(
                chrom,
                start0,
                end0,
                truncate=True,
                stepper="samtools",
                max_depth=max_depth,
            ):
                if col.reference_pos != start0:
                    # print("[debug] skipping non-matching pileup")
                    continue

                for pr in col.pileups:
                    rd = pr.alignment

                    # Only primary, mapped reads
                    if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                        # print("[debug] skipping non-primary read")
                        continue
                    if rd.mapping_quality < min_mapq:
                        # print("[debug] skipping low mapping quality read")
                        continue
                    if pr.is_del or pr.is_refskip:
                        # print("[debug] skipping deletion or reference skip")
                        continue  # no base here

                    qpos = pr.query_position  # 1-based
                    if qpos is None:
                        # print("[debug] skipping unmapped read")
                        continue

                    # Base/qual filters
                    qb = rd.query_sequence[qpos].upper()  # Sequence is 0-based
                    if min_baseq and rd.query_qualities is not None:
                        if rd.query_qualities[qpos] < min_baseq:
                            # print("[debug] skipping low base quality read")
                            continue

                    hp = None
                    if rd.has_tag(hp_tag):
                        try:
                            hp = int(rd.get_tag(hp_tag))
                            # print("[debug] Read haplotype tag:", hp, " for pos", pos1, "with ref", ref, "and read base", qb)
                        except Exception:
                            hp = None

                    key: Key = (pos1, ref, alt, hp)
                    if qb == ref:
                        counts[key][0] += 1  # Increment REF count for haplotype
                    elif qb == alt:
                        counts[key][1] += 1  # Increment ALT count for haplotype

        chr_counts[chrom] = counts

    vcf.close()
    bam.close()

    # Convert keys from (pos1, ref, alt, hp) to (pos1, ref, alt) with values as
    # counts for ref1, alt1, ref2, alt2
    chr_final_counts = {}
    for chrom, counts in chr_counts.items():
        final_counts = defaultdict(lambda: {"ref1": 0, "alt1": 0, "ref2": 0, "alt2": 0})
        for (pos1, ref, alt, hp), (ref_count, alt_count) in counts.items():
            if hp == 1:
                final_counts[(pos1, ref, alt)]["ref1"] += ref_count
                final_counts[(pos1, ref, alt)]["alt1"] += alt_count
            elif hp == 2:
                final_counts[(pos1, ref, alt)]["ref2"] += ref_count
                final_counts[(pos1, ref, alt)]["alt2"] += alt_count
            elif hp is None:
                final_counts[(pos1, ref, alt)]["ref1"] += ref_count
                final_counts[(pos1, ref, alt)]["alt1"] += alt_count
                final_counts[(pos1, ref, alt)]["ref2"] = 0
                final_counts[(pos1, ref, alt)]["alt2"] = 0
        chr_final_counts[chrom] = final_counts

    return chr_final_counts

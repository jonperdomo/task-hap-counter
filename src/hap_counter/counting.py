from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterator, Optional, Tuple

import pysam

Key = Tuple[int, str, str, Optional[int]]  # (pos1, REF, ALT, HP or None)


def count_snvs(
    bam_path: str,
    vcf_path: str,
    chrom: str,
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
    bam = pysam.AlignmentFile(bam_path, "rb")
    vcf = pysam.VariantFile(vcf_path)

    counts: Dict[Key, int] = defaultdict(int)

    def iter_alt_obs() -> Iterator[Key]:
        # Iterate SNVs on this chromosome from the VCF
        for rec in vcf.fetch(chrom):
            # bi-allelic SNV with A/C/G/T only
            if len(rec.alts or ()) != 1:
                print("[debug] skipping non-biallelic SNV")
                continue
            ref = rec.ref.upper()
            alt = rec.alts[0].upper()

            # print("[debug] ref:", ref, ", alt:", alt, ", pos:", rec.pos)

            if len(ref) != 1 or len(alt) != 1:
                print("[debug] skipping non-biallelic SNV")
                continue
            if ref not in "ACGT" or alt not in "ACGT":
                print("[debug] skipping non-ACGT SNV")
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
                    print("[debug] skipping non-matching pileup")
                    continue

                for pr in col.pileups:
                    rd = pr.alignment

                    # Only primary, mapped reads
                    if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                        print("[debug] skipping non-primary read")
                        continue
                    if rd.mapping_quality < min_mapq:
                        print("[debug] skipping low mapping quality read")
                        continue
                    if pr.is_del or pr.is_refskip:
                        print("[debug] skipping deletion or reference skip")
                        continue  # no base here

                    qpos = pr.query_position  # 1-based
                    if qpos is None:
                        print("[debug] skipping unmapped read")
                        continue

                    # Base/qual filters
                    qb = rd.query_sequence[qpos].upper()  # Sequence is 0-based
                    if min_baseq and rd.query_qualities is not None:
                        if rd.query_qualities[qpos] < min_baseq:
                            print("[debug] skipping low base quality read")
                            continue

                    # Count only ALT-supporting observations
                    if qb == alt or qb == ref:
                        hp: Optional[int] = None
                        if rd.has_tag(hp_tag):
                            try:
                                hp = int(rd.get_tag(hp_tag))
                            except Exception:
                                hp = None
                        yield (pos1, ref, qb, hp)
                        

    for key in iter_alt_obs():
        counts[key] += 1

    vcf.close()
    bam.close()
    return dict(counts)

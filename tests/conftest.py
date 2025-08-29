# tests/conftest.py
import os
import pysam
import pytest
from pathlib import Path

@pytest.fixture(scope="session")
def ref_seq():
    # chr1: 1..12  A C G T A C G T A  C  G  T
    return "ACGTACGTACGT"

@pytest.fixture(scope="session")
def snv_positions():
    return [(1, 'A', 'T'), (6, 'C', 'G'), (12, 'T', 'C')]

@pytest.fixture(scope="session")
def snv_vcf(tmp_path_factory, snv_positions):
    """
    Create a single VCF of bi-allelic SNVs once per test session (1-based):
      chr1:1  A>T
      chr1:6  C>G
      chr1:12 T>C
    Returns path to bgzipped VCF (*.vcf.gz).
    """
    test_outdir = os.path.join(os.getcwd(), "tests/test_output")
    os.makedirs(test_outdir, exist_ok=True)
    vcf_txt = os.path.join(test_outdir, "snvs.vcf")
    print("[debug] Creating synthetic VCF:", vcf_txt)
    with open(vcf_txt, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos, ref, alt in snv_positions:
            f.write(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    # bgzip + tabix index (creates snvs.vcf.gz and snvs.vcf.gz.tbi)
    pysam.tabix_index(str(vcf_txt), preset="vcf", force=True, keep_original=True)
    vcfgz = str(vcf_txt) + ".gz"
    return vcfgz

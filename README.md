# task-hap-counter
A solution for test-task-hap-counter (https://github.com/oxygen311/test-task-hap-counter/tree/main):

> Produce a program in Python programming language, which given read alignments (in BAM format with haplotype tags) and a set of variants (in phased VCF format) computes support for ALT and REF alleles across reads assigned to individual haplotypes.

## Setup
- Download test data
- Define dependencies
  - conda env create -f [environment.yml](environment.yml)
- Outline project structure

```
hap-counter
├─ README.md
├─ environment.yml
├─ src/
│  └─ hap-counter/
│     ├─ __init__.py
│     ├─ cli.py                  # Typer CLI: parses args, calls runner
│     ├─ runner.py               # orchestration: streams VCF → counts → TSV
│     ├─ io.py                   # BAM/VCF open helpers, region parsing, bgzip/index checks
│     ├─ filters.py              # read/base quality filters; primary-only logic
│     ├─ haplotag.py             # HP tag parsing (1/2), fallbacks, validation
│     ├─ counter.py              # core: per-SNV pileup → classify base vs REF/ALT by HP
│     ├─ models.py               # dataclasses (SNV, CountsRow), constants, enums
│     ├─ tsv_writer.py           # buffered TSV writer (+ optional extra columns)
│     └─ regions.py              # optional: restrict to --region BED or chr:start-end
├─ tests/
│  ├─ conftest.py
│  ├─ test_counter.py            # unit tests for SNV counting on synthetic reads
│  ├─ test_filters.py
│  ├─ test_haplotag.py
│  └─ data/                      # tiny synthetic BAM/VCF for tests
│     ├─ toy.bam
│     ├─ toy.bam.bai
│     ├─ toy.vcf.gz
│     └─ toy.vcf.gz.tbi
├─ scripts/
│  ├─ make_tiny_test_data.sh     # optional: generate synthetic BAM/VCF with samtools/bcftools
│  └─ bench.sh                   # optional: quick perf smoke test
└─ examples/
   └─ run_example.md             # “how to run” with example commands
```

- Define unit tests

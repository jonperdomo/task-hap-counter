# task-hap-counter
A solution for test-task-hap-counter (https://github.com/oxygen311/test-task-hap-counter/tree/main):

> Produce a program in Python programming language, which given read alignments (in BAM format with haplotype tags) and a set of variants (in phased VCF format) computes support for ALT and REF alleles across reads assigned to individual haplotypes.

## Setup
- Define dependencies
  - conda env create -f [environment.yml](environment.yml)
- Outline project structure

```
hap-counter
├─ README.md              # overview
├─ environment.yml        # dependencies
├─ src/
│  └─ hap-counter/
│     ├─ __init__.py
│     ├─ cli.py           # entry point
│     ├─ counting.py      # counter for SNV read support per haplotype
│     ├─ save_to_tsv.py   # saves results in TSV format
├─ tests/
│  ├─ conftest.py         # creates VCF for unit tests
│  └─ test_counter.py     # unit tests for SNV counting on synthetic reads
```

- Define unit tests
- Run test data

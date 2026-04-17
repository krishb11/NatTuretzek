# NatTuretzek

NatTuretzek is an R package for early-stage comparative RNA-seq workflows across species and developmental conditions.

## What was modernized

This package has been updated toward a more conventional **R package style**:

- clearer function signatures (fewer interactive prompts, more explicit parameters),
- stronger input validation and errors,
- roxygen documentation for onboarding,
- helper for package-shipped data directories.

## Installation

```r
# from local clone
install.packages(".", repos = NULL, type = "source")
```

## Quick onboarding

### 1) Kallisto quantification with `MyMuse`

```r
library(NatTuretzek)

muse <- MyMuse(
  cdsfile = "transcripts.fa",
  fastqfiles_dir = "fastq/",
  replicates = 3,
  paired_end = TRUE,
  bootstrap_samples = 100,
  kallisto_bin = "kallisto"
)

head(muse$est_counts)
head(muse$tpms)
```

For single-end libraries, provide fragment stats:

```r
muse <- MyMuse(
  cdsfile = "transcripts.fa",
  fastqfiles_dir = "fastq_single/",
  replicates = 3,
  paired_end = FALSE,
  fragment_length = 200,
  sd = 20
)
```

### 2) BLAST + Exonerate prep with `MyWitness`

```r
wit <- MyWitness(
  blasttype = "blastn",
  pathtogenome = "genome.fa",
  pathtoqueryfile = "queries.fa",
  species = "my_species"
)
```

### 3) Parse GFF and recover transcript FASTA with `MyDeep`

```r
deep <- MyDeep(
  file = "exonerate_output.gff",
  species = "my_species",
  genomefile = "genome.fa"
)
```


### Configuring BLAST/exonerate binaries and command options

If BLAST or exonerate are not on your PATH, pass explicit executable paths:

```r
wit <- MyWitness(
  blasttype = "blastn",
  pathtogenome = "genome.fa",
  pathtoqueryfile = "queries.fa",
  species = "my_species",
  blast_path = "/opt/ncbi/bin/blastn",
  exonerate_path = "/opt/exonerate/bin/exonerate",
  blast_params = list(
    outfmt = 6,
    max_target_seqs = 5,
    num_threads = 16,
    split = 1000,
    additional_args = "-evalue 1e-5"
  ),
  exonerate_params = list(
    model = "e2g",
    bestn = 3,
    minintron = 10,
    maxintron = 50000,
    gfffile = "custom_exonerate.gff",
    additional_args = "--percent 50"
  )
)
```

The constructor validates executable discovery and throws a clear error if tools are missing.

## Test data location for future FASTA/FASTQ files

A dedicated package data path is included for future kallisto test files:

```r
kallisto_test_data_dir()
```

Place your test FASTA/FASTQ files under `inst/extdata/kallisto/` in the source tree.

## Notes

- External tools (`kallisto`, BLAST, `exonerate`) must be installed and available on your PATH.
- This package currently focuses on workflow automation; downstream statistical interpretation is left to the analyst.

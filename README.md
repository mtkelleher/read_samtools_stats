# read_samtools_stats

Simple R functions for parsing statistics output from **samtools stats** and **bcftools stats** into tidy data frames.

## Overview

`samtools stats` produces a multi-section text report describing the alignment characteristics of a BAM/CRAM file. Each line is prefixed with a two- or three-letter tag identifying its section (e.g. `SN` for summary numbers, `MAPQ` for mapping qualities, `COV` for coverage). `bcftools stats` produces a multi-section text report describing the characteristics of a VCF/BCF file. Both the `read_samtools_stats()` and `read_bcftools_stats()` functions extract one section at a time and returns a tibble with meaningful column names, ready for downstream analysis in the tidyverse.

## Requirements

* R (≥ 4.0 recommended)
* [tidyverse](https://www.tidyverse.org/)
* samtools (version tested 1.19.2)
* bcftools (version tested 1.19)

## Installation

This is a single-file utility — there's nothing to install. Source the functions directly from GitHub:

```r
library(tidyverse)

source("https://raw.githubusercontent.com/mtkelleher/read_samtools_stats/main/read_samtools_stats.R")
source("https://raw.githubusercontent.com/mtkelleher/read_samtools_stats/main/read_bcftools_stats.R")
```

Or clone the repo and source the functions locally:

```r
library(tidyverse)

source("path/to/read_samtools_stats.R")
source("path/to/read_bcftools_stats.R")
```

## Usage

### `read_samtools_stats()`

First, generate a stats file with samtools:

```bash
samtools stats aligned.bam > aligned.stats.txt
```

Then, in R:

```r
library(tidyverse)

source("~/read_samtools_stats/read_samtools_stats.R")

# Summary numbers (one-row tibble, one column per metric)
SN <- read_samtools_stats("aligned.stats.txt", section = "SN")

# Mapping quality distribution
MAPQ <- read_samtools_stats("aligned.stats.txt", section = "MAPQ")

# For multiple samples start SN
SN <- NULL

# For sample_name in the imputed assembly alignment samtools stats files
for (sample_name in sample_names){

  file_path <- paste0(path_to_stats_files, sample_name, ".stats")

  # Run the read_samtools_stats function
  SN_temp <- read_samtools_stats(file_path, "SN")

  # Add the sample_name to the dataframe
  SN_temp$sample_name <- sample_name

  # Add the alignment description
  SN_temp$alignment_type <- "Reference"

  # Bind SN_temp rows to SN
  SN <- bind_rows(SN, SN_temp)

}
```

### `read_samtools_stats()` Arguments

| Argument          | Description                                                                       |
| ----------------- | --------------------------------------------------------------------------------- |
| `stats_file_path` | Path to a `samtools stats` output file.                                           |
| `section`         | Two- or three-letter string identifying which section to parse (see table below). |

### Supported `samtools stats` sections

| Section | Description                                                  |
| ------- | ------------------------------------------------------------ |
| `SN`    | Summary Numbers                                              |
| `FFQ`   | First Fragment Qualities                                     |
| `LFQ`   | Last Fragment Qualities                                      |
| `GCF`   | GC Content of first fragments                                |
| `GCL`   | GC Content of last fragments                                 |
| `GCC`   | ACGT content per cycle                                       |
| `GCT`   | ACGT content per cycle, read oriented                        |
| `FBC`   | ACGT content per cycle for first fragments                   |
| `FTC`   | ACGT raw counters for first fragments                        |
| `LBC`   | ACGT content per cycle for last fragments                    |
| `LTC`   | ACGT raw counters for last fragments                         |
| `IS`    | Insert sizes                                                 |
| `RL`    | Read lengths                                                 |
| `FRL`   | Read lengths - first fragments                               |
| `LRL`   | Read lengths - last fragments                                |
| `MAPQ`  | Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP) |
| `ID`    | Indel distribution                                           |
| `IC`    | Indels per cycle                                             |
| `COV`   | Coverage distribution                                        |
| `GCD`   | GC-depth                                                     |

Passing any other value to `section` will raise an error.

Sections returned without column renaming retain their original `V1`, `V2`, ... column names. Refer to the [samtools stats documentation](http://www.htslib.org/doc/samtools-stats.html) for the meaning of each column.

### Notes on the `SN` section

The `SN` (summary numbers) section is the most commonly used. It is pivoted wide so that each metric becomes its own column, with metric names cleaned up along the way:

* whitespace trimmed
* spaces and dashes replaced with underscores
* characters `( ) % :` removed
* `1st` replaced with `first` (to produce valid R names)

This makes it easy to `bind_rows()` or `map_dfr()` across samples.

### `read_bcftools_stats()`

First, generate a stats file with bcftools:

```bash
bcftools stats variants.vcf.gz > variants.stats
```

Then, in R:

```r
library(tidyverse)

source("~/read_samtools_stats/read_bcftools_stats.R")

SN <- read_bcftools_stats("~/variants.stats",  section = "SN")

# Using with aws.s3
library(aws.s3)
Sys.setenv(AWS_DEFAULT_REGION = s3_region)

# Get object
obj <- get_object(
  object = s3_filepath,
  bucket = s3_bucket,
  as = "text"
)

# Start a text connection
s3_text_connection <- textConnection(obj)

# The text connection will only work for one file read, but can be run again
DP <- read_bcftools_stats(s3_text_connection, section = "DP")
```

The function returns the parsed BCFtools statistics as a tidy data frame.

## See also

* [`samtools stats` documentation](http://www.htslib.org/doc/samtools-stats.html)
* [`bcftools stats` documentation](https://samtools.github.io/bcftools/bcftools.html#stats)
* [SAMtools/BCFtools/HTSlib](https://www.htslib.org/)


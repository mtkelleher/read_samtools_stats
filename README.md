# read_samtools_stats
A simple R function that parses info from a file generated with the samtools stats tool into a tidy data frame.

## Overview
 
`samtools stats` produces a multi-section text report describing the alignment characteristics of a BAM/CRAM file. Each line is prefixed with a two- or three-letter tag identifying its section (e.g. `SN` for summary numbers, `MAPQ` for mapping qualities, `COV` for coverage). This function extracts one section at a time and returns it as a tibble with meaningful column names, ready for downstream analysis in the tidyverse.
 
## Requirements
 
- R (≥ 4.0 recommended)
- [tidyverse](https://www.tidyverse.org/)
## Installation
 
This is a single-file utility — there's nothing to install. Source it directly from GitHub:
 
```r
library(tidyverse)
source("https://raw.githubusercontent.com/<user>/read_samtools_stats/main/read_samtools_stats.R")
```
 
Or clone the repo and source it locally:
 
```r
library(tidyverse)
source("path/to/read_samtools_stats.R")
```
 
## Usage
 
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
 
## Arguments
 
| Argument | Description |
|---|---|
| `stats_file_path` | Path to a `samtools stats` output file. |
| `section` | Two- or three-letter string identifying which section to parse (see table below). |
 
## Supported sections
 
| Section | Description |
|---------|-------------|
| `SN`    | Summary Numbers |
| `FFQ`   | First Fragment Qualities |
| `LFQ`   | Last Fragment Qualities |
| `GCF`   | GC Content of first fragments |
| `GCL`   | GC Content of last fragments |
| `GCC`   | ACGT content per cycle |
| `GCT`   | ACGT content per cycle, read oriented |
| `FBC`   | ACGT content per cycle for first fragments |
| `FTC`   | ACGT raw counters for first fragments |
| `LBC`   | ACGT content per cycle for last fragments |
| `LTC`   | ACGT raw counters for last fragments |
| `IS`    | Insert sizes |
| `RL`    | Read lengths |
| `FRL`   | Read lengths - first fragments |
| `LRL`   | Read lengths - last fragments |
| `MAPQ`  | Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP) |
| `ID`    | Indel distribution |
| `IC`    | Indels per cycle |
| `COV`   | Coverage distribution |
| `GCD`   | GC-depth |


 
Passing any other value to `section` will raise an error.
 
Sections marked "no" under column renaming are returned with their original `V1, V2, ...` columns — refer to the [samtools stats documentation](http://www.htslib.org/doc/samtools-stats.html) for the meaning of each column.
 
## Example
 
```r
library(tidyverse)
source("read_samtools_stats.R")
 
stats_file <- "sample.stats.txt"
 
# Pull summary numbers into a one-row tibble
sn <- read_samtools_stats(stats_file, "SN")
sn$reads_mapped
sn$error_rate
 
# Plot the mapping quality distribution
read_samtools_stats(stats_file, "MAPQ") %>%
  ggplot(aes(x = mapq, y = count)) +
  geom_col() +
  labs(x = "Mapping quality", y = "Read count")
 
# Combine summary stats across many samples
stats_files <- list.files("stats/", pattern = "\\.stats\\.txt$", full.names = TRUE)
 
all_sn <- stats_files %>%
  set_names(basename) %>%
  map_dfr(read_samtools_stats, section = "SN", .id = "sample")
```
 
## Notes on the `SN` section
 
The `SN` (summary numbers) section is the most commonly used. It is pivoted wide so that each metric becomes its own column, with metric names cleaned up along the way:
 
- whitespace trimmed
- spaces and dashes replaced with underscores
- characters `( ) % :` removed
- `1st` replaced with `first` (to produce valid R names)
This makes it easy to `bind_rows()` or `map_dfr()` across samples.
 
## See also
 
- [`samtools stats` documentation](http://www.htslib.org/doc/samtools-stats.html)
- [samtools](https://www.htslib.org/)
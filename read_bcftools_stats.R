# read_bcftools_stats.R
# Description: Parser for bcftools stats output files.
# Required packages: tidyverse

#### Function to read in bcftools stats file ####
read_bcftools_stats <- function(stats_file_path, section) {
  # Read the lines from the stats file
  lines <- readLines(stats_file_path)
  
  raw <- read.delim(
    text = lines[str_detect(lines, paste0("^", section))],
    header = FALSE,
    sep = "\t"
  )
  
  ## Summary Numbers ##
  if (section == "SN") {
    # Clean the raw df
    SN <- raw %>%
      mutate(
        id = V2,
        metric = V3,
        value = V4
      ) %>%
      # Select relevant columns
      select(id, metric, value) %>%
      # Trim the whitespace
      mutate(metric = str_trim(metric)) %>%
      # Replace spaces with underscore
      mutate(metric = str_replace_all(metric, " ", "_")) %>%
      # Replace dashes with underscore
      mutate(metric = str_replace_all(metric, "-", "_")) %>%
      # Mutate the values to numeric
      mutate(value = as.numeric(value)) %>%
      # Pivot metric to column names
      pivot_wider(names_from = metric, values_from = value)
    
    return(SN)
  }
  
  ## transitions/transversions ##
  else if (section == "TSTV") {
    # Clean the raw df
    TSTV <- raw %>%
      mutate(
        id = V2,
        ts = V3,
        tv = V4,
        ts_tv_ratio	= V5,
        ts_1st_ALT = V6,
        tv_1st_ALT = V7,
        ts_tv_ratio_1st_ALT = V8
      ) %>%
      # Select relevant columns
      select(-V1, -V2, -V3, -V4, -V5, -V6, -V7, -V8)
    
    return(TSTV)
  }
  
  ## Singleton stats ##
  else if (section == "SiS") {
    # Clean the raw df
    SiS <- raw %>%
      mutate(
        id = V2,
        allele_count = V3,
        number_SNPs = V4,
        number_transitions = V5,
        number_transversions = V6,
        number_indels = V7,
        repeat_consistent = V8,
        repeat_inconsistent = V9,
        not_applicable = V10
      ) %>%
      # Select relevant columns
      select(-V1, -V2, -V3, -V4, -V5, -V6, -V7, -V8, -V9, -V10) %>%
      # Mutate numeric columns
      mutate(across(-id, as.numeric))
    
    return(SiS)
  }
  
  ## Stats by non-reference allele frequency ##
  else if (section == "AF") {
    # Clean the raw df
    AF <- raw %>%
      mutate(
        id = V2,
        allele_frequency = V3,
        number_SNPs = V4,
        number_transitions = V5,
        number_transversions = V6,
        number_indels = V7,
        repeat_consistent = V8,
        repeat_inconsistent = V9,
        not_applicable = V10
      ) %>%
      # Select relevant columns
      select(-V1, -V2, -V3, -V4, -V5, -V6, -V7, -V8, -V9, -V10) %>%
      # Mutate numeric columns
      mutate(across(-id, as.numeric))
    
    return(AF)
  }
  
  ## Stats by quality ##
  else if (section == "QUAL") {
    # Clean the raw df
    QUAL <- raw %>%
      mutate(
        id = V2,
        quality = V3,
        number_SNPs = V4,
        number_transitions_1st_ALT = V5,
        number_transversions_1st_ALT = V6,
        number_indels = V7
      ) %>%
      # Select relevant columns
      select(-V1, -V2, -V3, -V4, -V5, -V6, -V7) %>%
      # Mutate numeric columns
      mutate(across(-id, as.numeric))
    
    return(QUAL)
  }
  
  ## InDel distribution ##
  else if (section == "IDD") {
    # Clean the raw df
    IDD <- raw %>%
      mutate(
        id = V2,
        length = V3,
        number_sites = V4,
        number_genotypes = V5,
        mean_VAF = V6
      ) %>%
      # Select relevant columns
      select(-V1, -V2, -V3, -V4, -V5, -V6) %>%
      # Mutate numeric columns
      mutate(across(-id, as.numeric))
    
    return(IDD)
  }
  
  ## Substitution types ##
  else if (section == "ST") {
    # Clean the raw df
    ST <- raw %>%
      mutate(
        id = V2,
        type = V3,
        count = V4
      ) %>%
      # Select relevant columns
      select(-V1, -V2, -V3, -V4) %>%
      # Mutate numeric columns
      mutate(count = as.numeric(count))
    
    return(ST)
  }
  
  ## Depth distribution ##
  else if (section == "DP") {
    # Clean the raw df
    DP <- raw %>%
      mutate(
        id = V2,
        bin = V3,
        number_genotypes = V4,
        fraction_genotypes = V5,
        number_sites = V6,
        fraction_sites = V7
      ) %>%
      # Select relevant columns
      select(-V1, -V2, -V3, -V4, -V5, -V6, -V7) %>%
      # Mutate numeric columns
      mutate(across(-id, as.numeric))
    
    return(DP)
  }
  
  else {
    stop("Invalid section provided")
  }
}




# read_samtools_stats.R
# Description: Parser for samtools stats output files.
# Required packages: tidyverse

#### Function to read in samtools stats file ####
read_samtools_stats <- function(stats_file_path, section) {
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
        metric = V2,
        value = V3,
        comment = V4
      ) %>%
      # Select relevant columns
      select(metric, value) %>%
      # Replace 1st with first
      mutate(metric = str_replace_all(metric, "1st", "first")) %>%
      # Remove all characters in (%):
      mutate(metric = str_remove_all(metric, "[(%):]")) %>%
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
  
  ## First Fragment Qualities ##
  else if (section == "FFQ") {
    FFQ <- raw
    
    return(FFQ)
  }

  ## Last Fragment Qualities ##
  else if (section == "LFQ") {
    LFQ <- raw
    
    return(LFQ)
  }
  
  ## GC Content of first fragments ##
  else if (section == "GCF") {
    GCF <- raw
    
    return(GCF)
  }
  
  ## GC Content of last fragments ##
  else if (section == "GCL") {
    GCL <- raw
    
    return(GCL)
  }
  
  ## ACGT content per cycle ##
  else if (section == "GCC") {
    GCC <- raw %>%
      mutate(
        cycle = V2,
        A_percentage_of_bases = V3,
        C_percentage_of_bases = V4,
        G_percentage_of_bases = V5,
        T_percentage_of_bases = V6,
        N_percentage_of_bases = V7,
        zero_percentage_of_bases = V8
      ) %>%
      select(
        cycle, A_percentage_of_bases, C_percentage_of_bases, G_percentage_of_bases,
        T_percentage_of_bases, N_percentage_of_bases, zero_percentage_of_bases)
    
    return(GCC)
  }
  
  ## ACGT content per cycle, read oriented ##
  else if (section == "GCT") {
    GCT <- raw %>%
      mutate(
        cycle = V2,
        A_percentage_of_bases = V3,
        C_percentage_of_bases = V4,
        G_percentage_of_bases = V5,
        T_percentage_of_bases = V6,
        N_percentage_of_bases = V7,
        zero_percentage_of_bases = V8
      ) %>%
      select(
        cycle, A_percentage_of_bases, C_percentage_of_bases, G_percentage_of_bases,
        T_percentage_of_bases, N_percentage_of_bases, zero_percentage_of_bases)
    
    return(GCT)
  }
  
  ## ACGT content per cycle for first fragments ##
  else if (section == "FBC") {
    FBC <- raw
    
    return(FBC)
  }
  
  ## ACGT raw counters for first fragments ##
  else if (section == "FTC") {
    FTC <- raw
    
    return(FTC)
  }
  
  ## ACGT content per cycle for last fragments ##
  else if (section == "LBC") {
    LBC <- raw
    
    return(LBC)
  }
  
  ## ACGT raw counters for last fragments ##
  else if (section == "LTC") {
    LTC <- raw
    
    return(LTC)
  }
  
  ## Insert sizes ##
  else if (section == "IS") {
    IS <- raw
    
    return(IS)
  }
  
  ## Read lengths ##
  else if (section == "RL") {
    RL <- raw
    
    return(RL)
  }
  
  ## Read lengths - first fragments ##
  else if (section == "FRL") {
    FRL <- raw
    
    return(FRL)
  }
  
  ## Read lengths - last fragments ##
  else if (section == "LRL") {
    LRL <- raw
    
    return(LRL)
  }
  
  ## Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP) ##
  else if (section == "MAPQ") {
    MAPQ <- raw %>%
      mutate(mapq = V2, count = V3) %>%
      select(mapq, count)
    
    return(MAPQ)
  }
  
  ## Indel distribution ##
  else if (section == "ID") {
    ID <- raw %>%
      mutate(
        length = V2,
        number_of_insertions = V3,
        number_of_deletions = V4
      ) %>%
      select(length, number_of_insertions, number_of_deletions)
    
    return(ID)
  }
  
  ## Indels per cycle ##
  else if (section == "IC") {
    IC <- raw %>%
      mutate(
        cycle = V2,
        number_of_insertions_fwd = V3,
        number_of_insertions_rev = V4,
        number_of_deletions_fwd = V5,
        number_of_deletions_rev = V6
      ) %>%
      select(
        cycle, number_of_insertions_fwd, number_of_insertions_rev,
        number_of_deletions_fwd, number_of_deletions_rev
        )
    
    return(IC)
  }
  
  ## Coverage distribution ##
  else if (section == "COV") {
    COV <- raw %>%
      mutate(
        inclusive_coverage_range_min_max = V2,
        max_portion_of_depth_range = V3,
        reference_sites = V4,
      ) %>%
      select(
        inclusive_coverage_range_min_max, max_portion_of_depth_range,
        reference_sites
      )
    
    return(COV)
  }
  
  ## GC-depth ##
  else if (section == "GCD") {
    GCD <- raw %>%
      mutate(
        GC_percentage = V2,
        unique_sequence_percentiles = V3,
        depth_percentile_10th = V4,
        depth_percentile_25th = V5,
        depth_percentile_50th = V6,
        depth_percentile_75th = V7,
        depth_percentile_90th = V8
      ) %>%
      select(
        GC_percentage, unique_sequence_percentiles, depth_percentile_10th,
        depth_percentile_25th, depth_percentile_50th, depth_percentile_75th,
        depth_percentile_90th
      )
    
    return(GCD)
  }
    
  else {
    stop("Invalid section provided")
  }
}




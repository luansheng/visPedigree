# Load necessary libraries
devtools::load_all()
library(data.table)

# Function to analyze a single pedigree file
analyze_pedigree <- function(filepath, name) {
  cat(sprintf("\n\n======================================================\n"))
  cat(sprintf("ANALYSIS REPORT FOR: %s\n", name))
  cat(sprintf("======================================================\n"))
  
  # Load data
  ped_dt <- fread(filepath)
  cat(sprintf("Loaded %d rows from %s\n", nrow(ped_dt), name))
  print(head(ped_dt))
  
  # Preprocessing column names
  if ("AnimalID" %in% names(ped_dt)) setnames(ped_dt, "AnimalID", "Ind")
  if ("SireID" %in% names(ped_dt)) setnames(ped_dt, "SireID", "Sire")
  if ("DamID" %in% names(ped_dt)) setnames(ped_dt, "DamID", "Dam")
  if ("SexID" %in% names(ped_dt)) setnames(ped_dt, "SexID", "Sex")
  
  # Specific fix for sex column in first file if needed
  if ("Sex" %in% names(ped_dt)) {
     ped_dt[, Sex := tolower(Sex)]
  }
  
  # 1. Tidy Pedigree
  cat("\n--- 1. Tidy Pedigree ---\n")
  tryCatch({
    # Use tidyped to clean and structure the data
    # We add generation numbers (addgen) and numeric IDs (addnum)
    tped <- tidyped(ped_dt, addgen = TRUE, addnum = TRUE)
    cat("Tidy pedigree created successfully.\n")
    print(tped)
    
    # 2. Basic Statistics
    cat("\n--- 2. Pedigree Statistics (pedstats) ---\n")
    # We use the new accelerated pedstats
    # For the first file with 'Year', we use it as timevar
    time_col <- if ("Year" %in% names(tped)) "Year" else NULL
    
    stats <- pedstats(tped, timevar = time_col, calc_ecg = TRUE, calc_genint = TRUE)
    print(stats)
    
    # 3. Inbreeding Statistics
    cat("\n--- 3. Inbreeding Analysis ---\n")
    # Calculate inbreeding if not already done by tidyped/pedstats
    if (!"f" %in% names(tped)) {
      tped <- inbreed(tped)
    }
    
    cat("Summary of Inbreeding Coefficient (F):\n")
    print(summary(tped$f))
    
    cat("\nInbreeding Classes:\n")
    print(pedinbreed_class(tped))
    
    # 4. Effective Population Size (Ne)
    cat("\n--- 4. Effective Population Size (Ne) ---\n")
    if (!is.null(time_col)) {
      ne_res <- pedne(tped, by = time_col)
      print(ne_res)
    } else {
      cat("No time variable found, calculating Ne by Generation...\n")
      ne_res <- pedne(tped, by = "Gen")
      print(ne_res)
    }
    
    # 5. Founder Contributions
    cat("\n--- 5. Founder Contributions (Top 10) ---\n")
    # Define candidates as the last generation
    max_gen <- max(tped$Gen)
    cand <- tped[Gen == max_gen, Ind]
    cat(sprintf("Analyzing contributions to %d individuals in Gen %d\n", length(cand), max_gen))
    
    contrib <- pedcontrib(tped, cand = cand, mode = "founder", top = 10)
    print(contrib)
    
  }, error = function(e) {
    cat(sprintf("Error during analysis: %s\n", e$message))
  })
}

# Run analysis for both files
ped1_path <- system.file("extdata", "G09_prunedPed.csv", package = "visPedigree")
if (ped1_path == "") ped1_path <- "inst/extdata/G09_prunedPed.csv" # Fallback for dev mode
analyze_pedigree(ped1_path, "G09_prunedPed.csv")

ped2_path <- system.file("extdata", "g2016G00to2025G09Ped.csv", package = "visPedigree")
if (ped2_path == "") ped2_path <- "inst/extdata/g2016G00to2025G09Ped.csv" # Fallback for dev mode
analyze_pedigree(ped2_path, "g2016G00to2025G09Ped.csv")

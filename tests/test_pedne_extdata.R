
library(data.table)
devtools::load_all(".")

# Function to run tests on a CSV
run_pedne_test <- function(csv_path) {
  cat("\n================================================================================\n")
  cat("Testing file:", basename(csv_path), "\n")
  cat("================================================================================\n")
  
  # Read CSV
  raw_ped <- fread(csv_path)
  
  # Tidy
  # tidyped expects first 3 cols to be Ind, Sire, Dam
  # Let's ensure duplicate IDs are removed or handled if needed, usually tidyped handles it
  tped <- tryCatch({
    tidyped(raw_ped, addgen = TRUE)
  }, error = function(e) {
    cat("Error in tidyped():", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(tped)) return()
  
  # Determine 'by' column
  if ("Year" %in% names(tped)) {
    by_col <- "Year"
  } else if ("BirthYear" %in% names(tped)) {
    by_col <- "BirthYear"
  } else {
    by_col <- "Gen"
    cat("No Year/BirthYear column found. Using 'Gen' as grouping variable.\n")
  }
  
  cat("Using grouping column:", by_col, "\n")
  
  # Test 1: Demographic Ne
  cat("\n--- Method: Demographic ---\n")
  tryCatch({
    ne_demo <- pedne(tped, method = "demographic", by = by_col)
    print(ne_demo)
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
  
  # Test 2: Coancestry Ne
  # This uses the new C++ function I just fixed
  cat("\n--- Method: Coancestry (Sampled) ---\n")
  tryCatch({
    # Use smaller sample size for speed in test
    ne_co <- pedne(tped, method = "coancestry", by = by_col, nsamples = 200)
    print(ne_co)
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
  
  # Test 3: Inbreeding Ne
  # inbreeding method requires 'f' and 'ECG'
  # pedne calcs them internally if missing, but let's see
  cat("\n--- Method: Inbreeding ---\n")
  tryCatch({
    ne_inb <- pedne(tped, method = "inbreeding", by = by_col)
    print(ne_inb)
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

# Path to extdata
extdata_path <- "inst/extdata"
files <- list.files(extdata_path, pattern = "\\.csv$", full.names = TRUE)

for (f in files) {
  run_pedne_test(f)
}

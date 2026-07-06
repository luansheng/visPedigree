required <- c("DESCRIPTION", "renv.lock", file.path("renv", "settings.json"))
missing <- required[!file.exists(required)]

if (length(missing)) {
  stop(
    "Run bootstrap from the visPedigree project root. Missing: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

message("Restoring the renv library from renv.lock...")
renv::restore(prompt = FALSE)

message("Checking runtime and development dependencies...")
status <- renv::status(dev = TRUE, library = renv::paths$library())

if (!isTRUE(status$synchronized)) {
  stop("renv restore completed, but the project is not synchronized.", call. = FALSE)
}

message("visPedigree development environment is ready.")

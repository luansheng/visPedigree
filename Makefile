R_SCRIPT ?= Rscript
TEST_FILTER ?=

.PHONY: bootstrap status test check docs

bootstrap:
	$(R_SCRIPT) scripts/bootstrap.R

status:
	$(R_SCRIPT) -e 'status <- renv::status(dev = TRUE, library = renv::paths$$library()); quit(status = if (isTRUE(status[["synchronized"]])) 0L else 1L)'

test:
	$(R_SCRIPT) -e 'filter <- "$(TEST_FILTER)"; devtools::test(filter = if (nzchar(filter)) filter else NULL, stop_on_failure = TRUE)'

check:
	$(R_SCRIPT) -e 'devtools::check(cran = TRUE, manual = TRUE, error_on = "never")'

docs:
	$(R_SCRIPT) -e 'pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)'

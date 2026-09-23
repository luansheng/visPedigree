# renv is a local development tool here: renv.lock is not tracked (see
# .gitignore), so on CI there is no lockfile for it to restore and nothing for
# it to manage.
#
# Activating it on CI is not merely useless, it breaks the pkgdown workflow.
# renv masks install.packages(), so the pak that setup-r-dependencies installs
# lands in renv/library as a symlink into renv's cache, while the packages pak
# then installs are real directories. The workflow caches renv/library but not
# renv's cache, so once that cache is restored the symlink dangles and
# library(pak) fails with "there is no package called 'pak'".
#
# CI is set to true by GitHub Actions, so renv activates locally and is skipped
# on the runner, where setup-r-dependencies falls back to its usual
# R_LIBS_SITE path and caches real directories.
if (!nzchar(Sys.getenv("CI"))) source("renv/activate.R")

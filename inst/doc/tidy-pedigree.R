## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installdevtools,eval=TRUE-------------------------------------------
suppressPackageStartupMessages(is_installed <- require(devtools))
if (!is_installed) {
  install.packages("devtools")
  suppressPackageStartupMessages(require(devtools))
}

## ----installvisPed, eval=TRUE--------------------------------------------
suppressPackageStartupMessages(is_installed <- require(visPedigree))
if (!is_installed) {
  install_github("luansheng/visPedigree")  
  suppressPackageStartupMessages(require(visPedigree))
}

## ----gettingdataset,eval=FALSE-------------------------------------------
#  data(package="visPedigree")

## ----simpleped-----------------------------------------------------------
head(simple_ped)
tail(simple_ped)
# The number of individuals in the pedigree dataset
nrow(simple_ped)
# Individual records with missing parents
simple_ped[Sire %in% c("0", "*", "NA", NA) |
             Dam %in% c("0", "*", "NA", NA)]

## ------------------------------------------------------------------------
x <- data.table::copy(simple_ped)
x[ID == "J2F588", Sire := "J0Z167"]
y <- tidyped(x)

## ----tidyped-------------------------------------------------------------
tidy_simple_ped <- tidyped(simple_ped)
head(tidy_simple_ped)
tail(tidy_simple_ped)
nrow(tidy_simple_ped)

## ------------------------------------------------------------------------
tidy_simple_ped_no_gen_num <-
  tidyped(simple_ped, addgen = FALSE, addnum = FALSE)
  head(tidy_simple_ped_no_gen_num)

## ----writeped,eval=FALSE-------------------------------------------------
#  saved_ped <- data.table::copy(tidy_simple_ped)
#  saved_ped[is.na(Sire), Sire := "0"]
#  saved_ped[is.na(Dam), Dam := "0"]
#  data.table::fwrite(
#    x = saved_ped,
#    file = "tidysimpleped.csv",
#    sep = ",",
#    quote = FALSE
#  )

## ------------------------------------------------------------------------
tidy_simple_ped_J5X804_ancestors <-
  tidyped(ped = tidy_simple_ped_no_gen_num, cand = "J5X804")
  tail(tidy_simple_ped_J5X804_ancestors)

## ------------------------------------------------------------------------
tidy_simple_ped_J5X804_ancestors_2 <-
  tidyped(ped = tidy_simple_ped_no_gen_num,
  cand = "J5X804",
  tracegen = 2)
  print(tidy_simple_ped_J5X804_ancestors_2)

## ------------------------------------------------------------------------
tidy_simple_ped_J0Z990_offspring <-
  tidyped(ped = tidy_simple_ped_no_gen_num, cand = "J0Z990", trace = "down")
  print(tidy_simple_ped_J0Z990_offspring)


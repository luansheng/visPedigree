## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installdevtools, include=FALSE, eval=TRUE---------------------------
suppressPackageStartupMessages(is_installed <- require(devtools))
if (!is_installed) {
  install.packages("devtools")
  suppressPackageStartupMessages(require(devtools))
} 

## ----installvisPed, include=FALSE, eval=TRUE-----------------------------
suppressPackageStartupMessages(is_installed <- require(visPedigree))
if (!is_installed) {
  install_github("luansheng/visPedigree")  
  suppressPackageStartupMessages(require(visPedigree))
} 

## ----smallped, fig.width=6.5, fig.height=6.5, fig.show="hold"------------
tidy_small_ped <-
  tidyped(ped = small_ped,
          cand = c("Y", "Z1", "Z2"))
visped(tidy_small_ped, compact = TRUE, file = "smallped.pdf")

## ----vissimpleped,  fig.width=6.5, fig.height=6.5, fig.show="hold"-------
tidy_simple_ped <- tidyped(simple_ped)
visped(tidy_simple_ped)

## ------------------------------------------------------------------------
suppressMessages(visped(tidy_simple_ped, showgraph = FALSE, file="simpleped.pdf"))

## ----deepped, eval=FALSE-------------------------------------------------
#  cand_J11_labels <- deep_ped[(substr(Ind, 1, 3) == "K11"), Ind]
#  visped(tidyped(deep_ped, cand = cand_J11_labels, tracegen = 3))

## ----reduceped1, fig.width=6.5, fig.height=6.5---------------------------
cand_J11_labels <- deep_ped[(substr(Ind,1,3) == "K11"),Ind]
visped(
  tidyped(
    deep_ped,
    cand = cand_J11_labels,
    trace = "up",
    tracegen = 3
  ),
  compact = TRUE,
  showgraph = TRUE,
  file = "deepped1.pdf"
)

## ----reduceped2, fig.width=6.5, fig.height=6.5---------------------------
visped(
  tidyped(
    deep_ped,
    cand = cand_J11_labels,
    trace = "up",
    tracegen = 3
  ),
  compact = TRUE,
  cex = 0.83,
  showgraph = FALSE,
  file = "deepped2.pdf"
)

## ----reduceped3, fig.width=6.5, fig.height=6.5---------------------------
suppressMessages(visped(
  tidyped(
    deep_ped,
    cand = cand_J11_labels, 
    tracegen = 3),
  compact = TRUE,
  outline = TRUE,
  showgraph = TRUE,
  file = "deepped3.pdf"
))

## ----pedofoneind, fig.width=6.5, fig.height=6.5--------------------------
suppressWarnings(J110550G_ped <-
                   tidyped(deep_ped, cand = "K110550H"))
suppressMessages(visped(J110550G_ped, showgraph = TRUE, file = "K110550HGped.pdf"))

## ----optiMate, fig.width=6.5, fig.height=6.5-----------------------------
cand_2007_G8_labels <-
  big_family_size_ped[(Year == 2007) & (substr(Ind, 1, 2) == "G8"), Ind]
suppressWarnings(
  cand_2007_G8_tidy_ped_ancestor_2 <-
    tidyped(
      big_family_size_ped,
      cand = cand_2007_G8_labels,
      trace = "up",
      tracegen = 2
    )
)
sire_label <-
  unique(cand_2007_G8_tidy_ped_ancestor_2[Ind %in% cand_2007_G8_labels,
                                          Sire])
dam_label <-
  unique(cand_2007_G8_tidy_ped_ancestor_2[Ind %in% cand_2007_G8_labels,
                                          Dam])
sire_dam_label <- unique(c(sire_label, dam_label))
sire_dam_label <- sire_dam_label[!is.na(sire_dam_label)]
sire_dam_ped <-
  cand_2007_G8_tidy_ped_ancestor_2[Ind %in% sire_dam_label]
sire_dam_ped <-
  sire_dam_ped[, FamilyID := paste(Sire, Dam, sep = "")]
family_size <- sire_dam_ped[, .N, by = c("FamilyID")]
fullsib_family_label <- unique(sire_dam_ped$FamilyID)
suppressMessages(
  visped(
    cand_2007_G8_tidy_ped_ancestor_2,
    compact = TRUE,
    outline = TRUE,
    showgraph = TRUE
  )
)


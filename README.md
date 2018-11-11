# visPedigree
Tidying and visualization for animal pedigree. 
This package takes a pedigree, checks duplicated, bisexual individuals, detects pedigree loop, adds missing founders, sorts parents before offspring, and traces the pedigree of the candidates. It outputs a hierarchical graph for all individuals in the pedigree. It can draw the graph of a very large pedigree (> 10,000 individuals per generation) by compacting the full-sib individuals. It is especially effective for drawing the pedigree of aquatic animal, which usually including many full-sib families per generation in the nucleus breeding population.
## To obtain visPedigree:
 * From GitHub:
   * Install devtools from CRAN
   ```R
   install.packages("devtools")
   ```
   * install the latest development version directly in R using the[ `devtools`](https://github.com/hadley/devtools) package:
   ```R
   library(devtools)
   install_github("luansheng/visPedigree")
   ```
## Vignette
Simplified Chinese link: https://luansheng.netlify.com/post/the-first-package-vispedigree-0-1/     
English link: https://luansheng.netlify.com/post/vispedigree-use-guide/


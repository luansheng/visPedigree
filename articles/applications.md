# Research using visPedigree

Since its GitHub release in 2018, `visPedigree` has been used in
livestock breeding, conservation genetics, human genealogy, and
large-scale pedigree research. This page includes only applications
documented in a paper’s methods, supplementary material, or source code;
citation alone is not treated as evidence of use.

## Featured applications

Conservation genetics

### Eastern black rhinoceros

*Natural dispersal is better than translocation for reducing risks of
inbreeding depression in eastern black rhinoceros*

PNAS, 2025

**Documented use:** Constructed observational pedigrees for each
geographic location and used those pedigrees to classify ancestry
cohorts.

[Article](https://doi.org/10.1073/pnas.2414412122) · [Methods
evidence](https://eprints.gla.ac.uk/352795/4/352795.pdf)

Human genealogy

### Historical human populations

*The long-lasting legacy of reproduction: lifetime reproductive success
shapes expected genetic contributions of humans after 10 generations*

Proceedings of the Royal Society B, 2023

**Documented use:** Traced deceased descendants and quantified whether
they continued the lineage, did not reproduce in the population, or
emigrated.

[Article](https://doi.org/10.1098/rspb.2023.0287) · [Methods
evidence](https://pure.rug.nl/ws/portalfiles/portal/762936487/rspb.2023.0287.pdf)

Large-scale pedigree QC

### Labrador Retriever

1,486,764

pedigree records

*randPedPCA: rapid approximation of principal components from large
pedigrees*

Genetics Selection Evolution, 2025

**Documented use:** Located the two individuals responsible for pedigree
loops; both loops resulted from incorrect sire assignments.

[Article](https://doi.org/10.1186/s12711-025-00994-y) · [Results
evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC12392600/)

## Explore the corresponding workflows

- [Pedigree
  visualization](https://luansheng.github.io/visPedigree/articles/draw-pedigree.md)
- [Pedigree preparation, loop detection, and
  tracing](https://luansheng.github.io/visPedigree/articles/tidy-pedigree.md)
- [Pedigree analysis and genetic
  diversity](https://luansheng.github.io/visPedigree/articles/pedigree-analysis.md)
- [Relationship matrices and matrix-free
  computation](https://luansheng.github.io/visPedigree/articles/relationship-matrix.md)

## Verified publications using visPedigree

The entries below are ordered by year. Each use statement is based on an
explicit methods or results statement in the linked source.

### 2025

**Mellya, R. V. K., et al.**  
[Natural dispersal is better than translocation for reducing risks of
inbreeding depression in eastern black rhinoceros (*Diceros bicornis
michaeli*)](https://doi.org/10.1073/pnas.2414412122).  
*Proceedings of the National Academy of Sciences*, 122(23), e2414412122.

**Verified use:** Construction of observational pedigrees for each
geographic location; the pedigrees were then used with known migration
and translocation histories to classify ancestry cohorts ([Methods,
p. 9](https://eprints.gla.ac.uk/352795/4/352795.pdf)).

------------------------------------------------------------------------

**Lee, H., Craddock, R. F., Gorjanc, G., & Becher, H.**  
[randPedPCA: rapid approximation of principal components from large
pedigrees](https://doi.org/10.1186/s12711-025-00994-y).  
*Genetics Selection Evolution*, 57, 46.

**Verified use:** Location of the two individuals responsible for
pedigree loops in 1,486,764 Labrador Retriever records. The paper
reports that both loops were caused by incorrect sire assignments
([Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC12392600/)).

------------------------------------------------------------------------

**Ferdosi, M. H., & Johnston, D.**  
[Identification of relationships and pedigree simplification using graph
theory](https://www.aaabg.org/aaabghome/AAABG26papers/A%2097-Ferdosi%20250424%20Final.pdf).  
*Proceedings of the Association for the Advancement of Animal Breeding
and Genetics*, 26, 303–306.

**Verified use:** Plotting of the 15-individual test pedigree shown in
Figure 1. The graph-theory analysis itself used `igraph`, not
`visPedigree` ([Methods,
p. 304](https://www.aaabg.org/aaabghome/AAABG26papers/A%2097-Ferdosi%20250424%20Final.pdf)).

### 2024

**Arias, K. D., et al.**  
[Population dynamics of potentially harmful haplotypes: a pedigree
analysis](https://doi.org/10.1186/s12864-024-10407-x).  
*BMC Genomics*, 25, 487.

**Verified use:** Visualization of a pedigree containing 471 Gochu
Asturcelta pigs from 51 families ([Methods, “Pedigree, cohorts and
genotyping”](https://link.springer.com/article/10.1186/s12864-024-10407-x)).

### 2023

**Young, E. A., Chesterton, E., Lummaa, V., Postma, E., & Dugdale, H.
L.**  
[The long-lasting legacy of reproduction: lifetime reproductive success
shapes expected genetic contributions of humans after 10
generations](https://doi.org/10.1098/rspb.2023.0287).  
*Proceedings of the Royal Society B*, 290(1998), 20230287.

**Verified use:** Tracing of deceased descendants to classify lineage
continuation, non-reproduction, and emigration ([Methods, section
2(d)](https://pure.rug.nl/ws/portalfiles/portal/762936487/rspb.2023.0287.pdf)).

------------------------------------------------------------------------

**Jones, T. B., Manseau, M., Merriell, B., Pittoello, G., Hervieux, D.,
& Wilson, P. J.**  
[Novel multilayer network analysis to assess variation in the spatial
co-occurrences of close kin in wild caribou
populations](https://doi.org/10.1016/j.gecco.2023.e02688).  
*Global Ecology and Conservation*, 47, e02688.

**Verified use:** Construction and visualization of pedigrees after
parentage inference with COLONY. Those pedigrees supported
classification of first-, second-, and third-order relationships for
downstream network analysis ([Methods, section
2.2](https://www.cclmportal.ca/sites/default/files/2023-10/1-s2.0-S2351989423003232-main.pdf)).

------------------------------------------------------------------------

**Vander Jagt, C. J., et al.**  
[Investigating the genetic cause of wry face in Australian Jersey
cattle](https://www.aaabg.org/aaabghome/AAABG25papers/58VanderJagt25238.pdf).  
*Proceedings of the Association for the Advancement of Animal Breeding
and Genetics*, 25, 238–241.

**Verified use:** Visualization of the pedigree used in an investigation
of common ancestors among Jersey cattle affected by wry face. The paper
attributes common-ancestor identification to manual investigation and
custom scripts, not to `visPedigree` ([Methods,
p. 239](https://www.aaabg.org/aaabghome/AAABG25papers/58VanderJagt25238.pdf)).

### 2022

**Arias, K. D., et al.**  
[Understanding Mendelian errors in SNP arrays data using a Gochu
Asturcelta pig pedigree: genomic alterations, family size and calling
errors](https://doi.org/10.1038/s41598-022-24340-0).  
*Scientific Reports*, 12, 19686.

**Verified use:** Visualization of the Gochu Asturcelta pedigree
presented in Figure 1 of the Mendelian-error study ([Methods, “Basic
statistical analyses and visualization of
data”](https://www.nature.com/articles/s41598-022-24340-0)).

### 2021

**León Rubio, J., da Silva Faria, R., Infante Gonzalez, J., Rincón
Lozano, Y., & Dominguez-Castaño, P.**  
[Genealogical analyses in open population of Silla Argentino horses
belonging to the national police of
Colombia](https://doi.org/10.18548/aspe/0009.22).  
*SPERMOVA*, 11(2), 166–172.

**Verified use:** Preparation of pedigree data for a total population of
1,861 Silla Argentino horses. Population parameters and
genetic-diversity statistics were subsequently calculated with ENDOG
([Methods, full
text](https://www.researchgate.net/publication/357713082_Genealogical_analyses_in_open_population_of_Silla_Argentino_horses_belonging_to_the_national_police_of_Colombia)).

## Software integrations

### MoBPS

The MoBPS function `get.pedigree.visual()` directly calls
[`visPedigree::tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
and
[`visPedigree::visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
to prepare and visualize pedigrees.

**Evidence:** [MoBPS source code at a fixed repository
revision](https://github.com/cran/MoBPS/blob/133a18bddfa50848ebd3a81b466bf314851dba69/R/get.pedigree.visual.R).

## Ecosystem recognition

### CRAN Task View: Agricultural Science

The [CRAN Task View: Agricultural
Science](https://cran.r-project.org/web/views/Agriculture.html) lists
`visPedigree` in its Animal Science section as a package for visualizing
complex animal pedigrees. This is ecosystem recognition, not a
publication or CRAN endorsement.

## Using visPedigree in your research?

If `visPedigree` contributes to published work, please cite the package.
Run the following command for the citation associated with the installed
version:

``` r

citation("visPedigree")
```

For version 1.9.0:

Luan, S. (2026). *visPedigree: Tidying, Analysis, and Fast Visualization
of Animal and Plant Pedigrees*. R package version 1.9.0.
<https://github.com/luansheng/visPedigree>.

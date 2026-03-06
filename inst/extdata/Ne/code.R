library(data.table)
library(magrittr)

library(ggplot2)
library(optiSel)

# method 1
ped <- fread("inst/extdata/Ne/2025PVBPGK1ped.csv")
ind_names_v <- ped[AnimalID %like% "^GK"]$AnimalID
tidy_ped <- tidyped(ped = ped)
coancestry_ne <- pedne(tidy_ped, cand=ind_names_v, by="Gen", method = "coancestry", nsamples = 1000, ncores = 4)
inbreeding_ne <- pedne(tidy_ped, cand=ind_names_v, by="Gen", method = "inbreeding")
demographic_ne <- pedne(tidy_ped, cand=ind_names_v, by="Gen", method = "demographic")

# method 2

source("inst/extdata/Ne/functions.R")
tidy_ped <- tidyped(ped = ped,cand = ind_names_v)
setDF(tidy_ped)
pedigree.pre <- prePed(tidy_ped)
current_gen <- max(unique(pedigree.pre$Gen))
keep <- pedigree.pre$Indiv[pedigree.pre$Gen %in% current_gen]    
keep.completeness <-
  summary(pedigree.pre, keep.only = keep, d = (current_gen-1)) 

cat(
  paste(
    "等价完整世代数（the number of equivalent complete generations，equiGen）是评价系谱完整度的一个重要参数。它指的是个体已知的祖先个体，在可追溯的每个世代所占的比例之和。",
    "育种群体",
    "世代家系个体的",
    "equiGen描述性统计参数如下所述。其中，",
    DescStatChinese(
      keep.completeness$equiGen,
      scale.unit = "个",
      digit.num = 2
    ),
    "。",
    "\n\n",
    sep = ""
  )
)
cat(
  paste(
    "\n\n系谱完整度指数（Index of pedigree completeness，PCI）是评价近交系数，估计育种值的一个重要参考指标。因为当个体的双亲有一个缺失时，其PCI指数就会为0，对于系谱缺失非常敏感。",
    "育种群体",
    "世代家系个体的",
    "PCI描述统计参数如下所述。其中，",
    DescStatChinese(keep.completeness$PCI, digit.num = 2),
    "。",
    "PCI越接近1，系谱越完整。",
    "\n\n",
    sep = ""
  )
)
cat(
  paste(
    "\n\n",
    "育种群体",
    "世代家系",
    "双亲均可追溯世代数（包括双亲）（Number of fully traced generations,fullGen）的描述统计参数如下所述。其中，",
    DescStatChinese(
      keep.completeness$fullGen,
      scale.unit = "个",
      digit.num = 2
    ),
    "。\n\n",
    sep = ""
  )
)
cat(
  paste(
    "\n\n",
    "育种群体",
    "世代家系",
    "最大可追溯世代数（Number of maximum generations traced,maxGen）的描述统计参数如下所述。其中，",
    DescStatChinese(
      keep.completeness$maxGen,
      scale.unit = "个",
      digit.num = 2
    ),
    "。\n\n",
    sep = ""
  )
)
keep.kin <- pedIBD(pedigree.pre, keep.only = keep)
#遗传多样性，通过共亲系数来体现
keep.genetic.diversity <- round(1 - mean(keep.kin[keep, keep]),2)
cat(paste0("遗传多样性为：",keep.genetic.diversity,", "))
g <- keep.completeness[keep, "equiGen"]
#g <- keep.completeness[Indiv %in% keep, "equiGen"]
N <- length(g)
n <-
  (matrix(g, N, N, byrow = TRUE) + matrix(g, N, N, byrow = FALSE)) /
  2
# 使用 Cervantes (2011) 原版公式（去除 Gutiérrez 的 -1 修正）
delta.c_matrix <- 1 - (1 - keep.kin[keep, keep]) ^ (1 / n)
delta.c_real <- delta.c_matrix[upper.tri(delta.c_matrix, diag = FALSE)]
effective.size_re <- 1 / (2 * mean(delta.c_real, na.rm = TRUE))
cat(paste0("有效群体大小为：",round(effective.size_re,2)))
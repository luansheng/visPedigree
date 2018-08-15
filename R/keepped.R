#' @import data.table
keepped <- function(ped, candidate=NULL, gen=NULL, num=FALSE){
  ped_check <- checkped(ped)
  #pruning the pedigree by candidate
  if (!is.null(candidate)) {
    ped_num <- numped(ped_check)
    i <- 1
    keep_ind <- match(candidate,ped_num$Ind)
    ind_n <- length(keep_ind) + 1

    while (length(keep_ind) != ind_n) {
      ind_n <- length(keep_ind)
      keep_ind <- unique(c(unlist(ped_num[keep_ind,c("SireNum","DamNum")]),keep_ind))
      keep_ind <- keep_ind[which(keep_ind > 0)]
      if (!is.null(gen)){
        if (i == gen) {
          break
        }
      }
      i <- i+1
    }
    ped_new <- ped_num[sort(keep_ind)]
    ped_new[,":="(IndNum=NULL,SireNum=NULL,DamNum=NULL)]
    if (length(keep_ind) != ind_n) {
      ped_new <- checkped(ped_new)
    }
  } else {
    ped_new <- ped_check
  }

  #converting individual, sire, and dam IDs to sequential integer
  if (num) {
    ped_new <- numped(ped_new)
  }
  return(ped_new)
}

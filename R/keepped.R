#' @import data.table
keepped <- function(ped, candidate=NULL, gen=NULL, num=FALSE){
  ped_check <- checkped(ped)
  #pruning the pedigree by candidate
  if (!is.null(candidate)) {
    ped_num <- numped(ped_check)
    i <- 1
    #backward traverse
    keep_ind_backward <- match(candidate,ped_num$Ind)
    ind_n <- length(keep_ind_backward) + 1
    while (length(keep_ind_backward) != ind_n) {
      ind_n <- length(keep_ind_backward)
      keep_ind_backward <- unique(c(unlist(ped_num[keep_ind_backward,c("SireNum","DamNum")]),keep_ind_backward))
      keep_ind_backward <- keep_ind_backward[which(keep_ind_backward > 0)]
      if (!is.null(gen)){
        if (i == gen) {
          break
        }
      }
      i <- i+1
    }

    i <- 1
    #forward traverse
    keep_ind_foreward <- match(candidate,ped_num$Ind)
    ind_n <- length(keep_ind_foreward) + 1
    while (length(keep_ind_foreward) != ind_n) {
      ind_n <- length(keep_ind_foreward)
      keep_ind_foreward <- unique(c(ped_num[which(SireNum %in% keep_ind_foreward),IndNum],
                                    ped_num[which(DamNum %in% keep_ind_foreward),IndNum],keep_ind_foreward))
      keep_ind_foreward <- keep_ind_foreward[which(keep_ind_foreward > 0)]
      if (!is.null(gen)){
        if (i == gen) {
          break
        }
      }
      i <- i+1
    }

    keep_ind <- unique(c(keep_ind_backward,keep_ind_foreward))
    ped_new <- ped_num[sort(keep_ind)]
    ped_new[,":="(IndNum=NULL,SireNum=NULL,DamNum=NULL)]

    #insure the pruned pedigree with the missing parents.
    ped_new <- checkped(ped_new)

  } else {
    ped_new <- ped_check
  }

  #converting individual, sire, and dam IDs to sequential integer
  if (num) {
    ped_new <- numped(ped_new)
  }
  return(ped_new)
}

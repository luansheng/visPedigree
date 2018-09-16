#' Tidy and prepare a pedigree
#'
#' \code{tidyped} function checks, sorts and returns a tidy pedigree.
#'
#' This function takes a pedigree, checks duplicated, bisexual individuals, detects pedigree loop, adds missing founders, and sorts the pedigree. If the parameter \emph{cand} contains individuals's IDs, then only these individuals and their ancestors will be kept in the pedigree. The tracing direction and tracing generation number can be provided when the parameter \emphs{cand} is not NULL. Individual virtual generation will be infered and added when the parameter \emph{addgen} is TRUE. Individual numeric pedigree will be added when the parameter \emph{addnum} is TRUE. Individual sex will be added if there is no "Sex" column in the pedigee. If the pedigree includes the column \emph{Sex}, then individuals's sexes need to be recoded as "male" and "female". Missing sexes will be identified from the pedigree structure if possible.
#'
#' The first three columns in the pedigree must be individual, sire and dam identification. Names of the three columns are assigend by users, but their orders must be not changded in the pedigree. Individual ID should not be coded as "", " ", "0", "*", and "NA", these individuals will be deleted from the pedigree. Missing parents should be denoted by either "NA", "0", "*". " " and "" should be recoded as missing parents, but not be recommended.
#'
#'
#' @import data.table
#' @param ped Data table from data.table or data frame including the pedigree, which have the first three columns individual, sire and dam IDs. More columns, such as sex, generation can be included in the pedigree file.
#' @param cand Vector of charactor including individual IDs, or NULL. Only the candidates and their ancestors will be kept in the pedigree if this parameter is not NULL.
#' @param trace A variate of charactor about how to trace pedigree, . "up" means tracing candidates's pedigree to ancestors, "down" means tracing candidates's pedigree to progenies, "all" means tracing candidaes's pedigree to ancestors and progenies simutaneously. This paramters can be used only if the cand parameter is not NULL. The defaulted value is "up".
#' @param tracegen An integer means the number of tracing generation.
#' @param addgen A logical parameter indicates whether individual generation number will be added in the retrued data.table. The default values is TRUE, then individual generation number will be added.
#' @param addnum A logical parameter indicates wheter numeric pedigree will be added in the returned data.table. The defaulted value is TRUE, then three new columns of IndNum, SireNum and DamNum will be added.
#'
#' @export
tidyped <-
  function(ped,
           cand = NULL,
           trace = "up",
           tracegen = NULL,
           addgen = TRUE,
           addnum = TRUE) {
    ped_colnames <- colnames(ped)
    if (c("Cand") %in% ped_colnames) {
      ped[,Cand:=NULL]
    }
    ped_check <- checkped(ped, addgen)
    #pruning the pedigree by candidate
    if (!is.null(cand)) {
      if (!all(class(cand) %in% c("character"))) {
        cand <- as.character(cand)
      }

      if (!any(cand %in% ped_check$Ind)) {
        stop("Some candidates are not in the pedigree!")
      } else {
        ped_check[Ind %chin% cand,Cand:=TRUE]
        ped_check[!(Ind %chin% cand),Cand:=FALSE]
      }
      ped_num <- numped(ped_check)
      if (trace %in% c("up","all")) {
        # Trace from candidate to ancestry
        i <- 1
        keep_ind_backward <- match(cand, ped_num$Ind)
        ind_n <- length(keep_ind_backward) + 1
        while (length(keep_ind_backward) != ind_n) {
          ind_n <- length(keep_ind_backward)
          keep_ind_backward <-
            unique(c(unlist(ped_num[keep_ind_backward, c("SireNum", "DamNum")]), keep_ind_backward))
          keep_ind_backward <-
            keep_ind_backward[which(keep_ind_backward > 0)]
          if (!is.null(tracegen)) {
            if (i == tracegen) {
              break
            }
          }
          i <- i + 1
        }
        keep_ind_backward <- unique(keep_ind_backward)
      }

      if (trace %in% c("down","all")) {
        #Trace from candidate to descendant
        i <- 1
        keep_ind_foreward <- match(cand, ped_num$Ind)
        ind_n <- length(keep_ind_foreward) + 1
        while (length(keep_ind_foreward) != ind_n) {
          ind_n <- length(keep_ind_foreward)
          keep_ind_foreward <-
            unique(c(ped_num[which(SireNum %in% keep_ind_foreward), IndNum],
                     ped_num[which(DamNum %in% keep_ind_foreward), IndNum], keep_ind_foreward))
          keep_ind_foreward <-
            keep_ind_foreward[which(keep_ind_foreward > 0)]
          if (!is.null(tracegen)) {
            if (i == tracegen) {
              break
            }
          }
          i <- i + 1
        }
        keep_ind_foreward <- unique(keep_ind_foreward)
      }

      if (trace %in% c("up","all")) {
        keep_ind <- keep_ind_backward
        ped_up <- ped_num[sort(keep_ind)]
        ped_up <- ped_up[((SireNum %in% keep_ind) | (DamNum %in% keep_ind))]
        ped_new <- ped_up
      }
      if (trace %in% c("down","all")) {
        keep_ind <- keep_ind_foreward
        ped_down <- ped_num[sort(keep_ind)]
        ped_down <- ped_down[((SireNum %in% keep_ind) | (DamNum %in% keep_ind))]
        ped_down[!(SireNum %in% keep_ind),Sire:=NA]
        ped_down[!(DamNum %in% keep_ind),Dam:=NA]
        ped_new <- ped_down
      }
      if (trace %in% c("all")) {
        ped_new <- unique(rbind(ped_up,ped_down))
      }

      ped_new[, ":="(IndNum = NULL,
                     SireNum = NULL,
                     DamNum = NULL)]

      #insure the pruned pedigree with the missing parents.
      ped_new <- checkped(ped_new, addgen)

    } else {
      ped_new <- ped_check
    }

    #converting individual, sire, and dam IDs to sequential integer
    if (addnum) {
      ped_new <- numped(ped_new)
    }
    return(ped_new)
  }

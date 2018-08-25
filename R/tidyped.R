#' @import data.table
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
      if (!(class(cand) %in% "character")) {
        cand <- as.character(cand)
      }

      if (!any(cand %in% ped_check$Ind)) {
        stop("not find candidate in the pedigree!")
      } else {
        ped_check[Ind %chin% cand,Cand:=TRUE]
        ped_check[!(Ind %chin% cand),Cand:=FALSE]
      }
      ped_num <- numped(ped_check)
      if (trace %in% c("all")) {
        i <- 1
        # Trace from candidate to ancestry
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

        i <- 1
        #Trace from candidate to descendant
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
        keep_ind <- unique(c(keep_ind_backward, keep_ind_foreward))
      }

      if (trace %in% c("up")) {
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
        keep_ind <- unique(keep_ind_backward)
      }

      if (trace %in% c("down")) {
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
        keep_ind <- unique(keep_ind_foreward)
      }


      ped_new <- ped_num[sort(keep_ind)]
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

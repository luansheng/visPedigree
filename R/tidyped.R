#' Tidy and prepare a pedigree
#'
#' \code{tidyped} function checks, sorts, traces, and returns a trimmed pedigree.
#'
#' This function takes a pedigree, checks duplicated, bisexual individuals, detects pedigree loop, adds missing founders, sorts the pedigree, and traces the pedigree of the candidates. If the parameter \emph{cand} contains individuals's IDs, then only these individuals and their ancestors or descendants will be kept in the pedigree. The tracing direction and tracing generation number can be provided when the parameters \emph{trace} and \emph{tracegen} are not NULL. Individual virtual generation will be infered and assigned when the parameter \emph{addgen} is TRUE. Numeric pedigree will be generated when the parameter \emph{addnum} is TRUE. All individuals' sex will be inferred if there is not sexual information in the pedigee. If the pedigree includes the column \strong{Sex}, then individuals' sexes need to be recoded as "male", "female", or NA (unknown sex). Missing sexes will be identified from the pedigree structure and be added if possible.
#'
#' @param ped A data.table or data frame including the pedigree, which including the first three columns: \strong{individual}, \strong{sire} and \strong{dam} IDs. More columns, such as sex, generation can be included in the pedigree file. Names of the three columns can be assigend as you would like, but their orders must be not changded in the pedigree. Individual ID should not be coded as "", " ", "0", asterisk, and "NA", otherwise these individuals will be deleted from the pedigree. Missing parents should be denoted by either "NA", "0", asterisk. Space and "" will also be recoded as missing parents, but not be recommended.
#' @param cand A character vector including individual IDs, or NULL. Only the candidates and their ancestors and offspring will be kept in the pedigree if this parameter is not NULL.
#' @param trace A character value implying how to trace the pedigree of the candidates, is equal to "\strong{up}", "\strong{down}", or "\strong{all}". "up" means tracing candidates's pedigree to ancestors; "down" means tracing candidates's pedigree to descendants, "all" means tracing candidaes' pedigree to ancestors and descendants simutaneously. This parameter can be used only if the \emph{cand} parameter is not NULL. The defaulted value is "up".
#' @param tracegen An integer means the number of tracing generation. This paramter can only be used when the \emph{trace} parameter is not NULL. All generations of the candidates will be traced when the parameter tracegen is NULL.
#' @param addgen A logical value indicates whether individual generation number will be generated. The default values is TRUE, then a new column named \strong{Gen} will be added in the returned data.table.
#' @param addnum A logical value indicates whether numeric pedigree will be generated. The defaulted value is TRUE, then three new columns of \strong{IndNum}, \strong{SireNum} and \strong{DamNum} will be added in the returned data.table.
#'
#' @return A data.table including the tidy pedigree is returned. Individual, sire and dam ID columns are renamed as \strong{Ind}, \strong{Sire} and \strong{Dam}. Missing parent is replacted with the default missing value \strong{NA}. The column \strong{Sex} includes individuals' sex (male or female, NA for unknown sex). The column \strong{Cand} will be included when the parameter \emph{cand} is not NULL. The column \strong{Gen} will be included when the parameter \emph{addgen} is TRUE. The columns \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} will be included when the parameter \emph{addnum} is TRUE. Ind, Sire, Dam and Sex columns are character; The column Cand is logical; The Gen, IndNum, SireNum and DamNum are integer.
#'
#' @examples
#' require(visPedigree)
#' require(data.table)
#' simple_ped
#' tidy_ped <- tidyped(simple_ped)
#' tidy_ped
#' # The pedigree of individual J5X804 to ancestors is pruned,
#' # and the column Cand is added and returned
#' tidy_ped_J5X804 <- tidyped(simple_ped,cand="J5X804")
#' tidy_ped_J5X804
#' # The pedigree of individual J5X804 to parents and grandparents is pruned
#' tidy_ped_J5X804_up_2 <- tidyped(simple_ped,cand="J5X804",trace="up",tracegen=2)
#' tidy_ped_J5X804_up_2
#' # The pedigree of individual J5X804 to offspring is pruned
#' tidy_ped_J0Z990_down <- tidyped(simple_ped,cand="J0Z990",trace="down")
#' tidy_ped_J0Z990_down
#' # The pedigree of individual J0Z990 to child and grandchild is pruned
#' tidy_ped_J0Z990_down_2 <- tidyped(simple_ped,cand="J0Z990",trace="down",tracegen=2)
#' tidy_ped_J0Z990_down_2
#'
#' @import data.table
#' @export
tidyped <-
  function(ped,
           cand = NULL,
           trace = "up",
           tracegen = NULL,
           addgen = TRUE,
           addnum = TRUE) {
    ped_is_DT <- "data.table" %in% class(ped)
    if (!ped_is_DT) {
      ped_inter <- as.data.table(ped)
    } else {
      ped_inter <- copy(ped)
    }
    attr(ped_inter,"tidyped") <- FALSE

    ped_colnames <- colnames(ped_inter)

    # Delete Cand column
    if (!is.null(cand)) {
      if (c("Cand") %in% ped_colnames) {
        ped_inter[,Cand:=NULL]
        warning("The column Cand of the original pedigree is deleted.")
      }
    }

    # The Gen column will be deleted
    if (c("Gen") %in% ped_colnames) {
        warning("The column Gen of the original pedigree is deleted.")
    }

    # IndNum SireNum or DamNum columns will be deleted
    three_num_columns <- c("IndNum","SireNum", "DamNum")
    if (any(three_num_columns %in% ped_colnames)) {
      exist_columns <- three_num_columns[three_num_columns %in% ped_colnames]
      warning("The columns ", paste(exist_columns,collapse = ","),
              " of the original pedigree are deleted.")
    }

    ped_check <- checkped(ped_inter, addgen)
    #pruning the pedigree by candidate
    if (!is.null(cand)) {
      if (!all(class(cand) %in% c("character"))) {
        cand <- as.character(cand)
      }

      if (!any(cand %in% ped_check$Ind)) {
        stop("No candidates are not in the pedigree!")
      } else {
        ped_check[Ind %chin% cand,Cand:=TRUE]
        ped_check[!(Ind %chin% cand),Cand:=FALSE]
      }
      ped_num <- numped(ped_check)
      cand_num <- match(cand, ped_num$Ind, nomatch = 0)
      if (trace %in% c("up","all")) {
        # Trace from candidate to ancestry
        i <- 1
        keep_ind_backward <- cand_num
        keep_ind_backward <-
          keep_ind_backward[which(keep_ind_backward > 0)]
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
        keep_ind_foreward <- cand_num
        keep_ind_foreward <-
          keep_ind_foreward[which(keep_ind_foreward > 0)]
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

    #Converting individual, sire, and dam IDs to sequential integer
    if (addnum) {
      ped_new <- numped(ped_new)
    }

    # add a new column Cand
    if (!is.null(cand)) {
      Cand <- NULL
      ped_new[,Cand := Ind %in% cand]
    }
    attr(ped_new,"tidyped") <- TRUE
    return(ped_new)
  }

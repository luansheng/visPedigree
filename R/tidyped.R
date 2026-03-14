#' Tidy and prepare a pedigree using graph theory
#'
#' This function takes a pedigree, checks for duplicate and bisexual individuals, detects pedigree loops using graph theory, adds missing founders, assigns generation numbers, sorts the pedigree, and traces the pedigree of specified candidates. If the \code{cand} parameter contains individual IDs, only those individuals and their ancestors or descendants are retained. Tracing direction and the number of generations can be specified using the \code{trace} and \code{tracegen} parameters.
#' 
#' Compared to the legacy version, this function handles cyclic pedigrees more robustly by detecting and reporting loops, and it is generally faster for large pedigrees due to the use of sparse graph algorithms from the \code{igraph} package. Generation assignment can be performed using either a top-down approach (default, aligning founders at the top) or a bottom-up approach (aligning terminal nodes at the bottom).
#'
#' @param ped A data.table or data frame containing the pedigree. The first three columns must be \strong{individual}, \strong{sire}, and \strong{dam} IDs. Additional columns, such as sex or generation, can be included. Column names can be customized, but their order must remain unchanged. Individual IDs should not be coded as "", " ", "0", "*", or "NA"; otherwise, they will be removed. Missing parents should be denoted by "NA", "0", or "*". Spaces and empty strings ("") are also treated as missing parents but are not recommended.
#' @param cand A character vector of individual IDs, or NULL. If provided, only the candidates and their ancestors/descendants are retained.
#' @param trace A character value specifying the tracing direction: "\strong{up}", "\strong{down}", or "\strong{all}". "up" traces ancestors; "down" traces descendants; "all" traces the union of ancestors and descendants. This parameter is only used if \code{cand} is not NULL. Default is "up".
#' @param tracegen An integer specifying the number of generations to trace. This parameter is only used if \code{trace} is not NULL. If NULL or 0, all available generations are traced.
#' @param addgen A logical value indicating whether to generate generation numbers. Default is TRUE, which adds a \strong{Gen} column to the output.
#' @param addnum A logical value indicating whether to generate a numeric pedigree. Default is TRUE, which adds \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} columns to the output.
#' @param inbreed A logical value indicating whether to calculate inbreeding coefficients. Default is FALSE. If TRUE, an \strong{f} column is added to the output. This uses the same optimized engine as \code{pedmat(..., method = "f")}.
#' @param selfing A logical value indicating whether to allow the same individual to appear as both sire and dam. This is common in plant breeding (monoecious species) where the same plant can serve as both male and female parent. If TRUE, individuals appearing in both the Sire and Dam columns will be assigned Sex = "monoecious" instead of triggering an error. Default is FALSE.
#' @param genmethod A character value specifying the generation assignment method: "\strong{top}" or "\strong{bottom}". "top" (top-aligned) assigns generations from parents to offspring, starting founders at Gen 1. "bottom" (bottom-aligned) assigns generations from offspring to parents, aligning terminal nodes at the bottom. Default is "top".
#' @param ... Additional arguments passed to \code{\link{inbreed}}.
#' 
#' @return A \code{tidyped} object (which inherits from \code{data.table}). Individual, sire, and dam ID columns are renamed to \strong{Ind}, \strong{Sire}, and \strong{Dam}. Missing parents are replaced with \strong{NA}. The \strong{Sex} column contains "male", "female", or NA. The \strong{Cand} column is included if \code{cand} is not NULL. The \strong{Gen} column is included if \code{addgen} is TRUE. The \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} columns are included if \code{addnum} is TRUE. The \strong{Family} and \strong{FamilySize} columns identify full-sibling families (e.g., "A x B" for offspring of sire A and dam B). The \strong{f} column is included if \code{inbreed} is TRUE.
#' 
#' @seealso 
#' \code{\link{summary.tidyped}} for summarizing tidyped objects
#' \code{\link{visped}} for visualizing pedigree structure
#' \code{\link{pedmat}} for computing relationship matrices
#' \code{\link{vismat}} for visualizing relationship matrices
#' \code{\link{splitped}} for splitting pedigree into connected components
#' \code{\link{inbreed}} for calculating inbreeding coefficients
#'
#' @examples
#' library(visPedigree)
#' library(data.table)
#'
#' # Tidy a simple pedigree
#' tidy_ped <- tidyped(simple_ped)
#' head(tidy_ped)
#'
#' # Trace ancestors of a specific individual within 2 generations
#' tidy_ped_tracegen <- tidyped(simple_ped, cand = "J5X804", trace = "up", tracegen = 2)
#' head(tidy_ped_tracegen)
#' 
#' # Trace both ancestors and descendants for multiple candidates
#' # This is highly optimized and works quickly even on 100k+ individuals
#' cand_list <- c("J5X804", "J3Y620")
#' tidy_ped_all <- tidyped(simple_ped, cand = cand_list, trace = "all")
#' 
#' # Check for loops (will error if loops exist)
#' try(tidyped(loop_ped))
#' 
#' # Example with a large pedigree: extract 2 generations of ancestors for 2007 candidates
#' cand_2007 <- big_family_size_ped[Year == 2007, Ind]
#' \donttest{
#' tidy_big <- tidyped(big_family_size_ped, cand = cand_2007, trace = "up", tracegen = 2)
#' summary(tidy_big)
#' }
#'
#' @import data.table
#' @importFrom igraph graph_from_data_frame is_dag topo_sort components subcomponent distances ego degree as_adj_list vcount add_vertices neighbors shortest_paths which_loop ends graph_from_edgelist neighborhood
#' @importFrom stats setNames
#' @importFrom utils head
#' @export
tidyped <- function(ped,
                          cand = NULL,
                          trace = "up",
                          tracegen = NULL,
                          addgen = TRUE,
                          addnum = TRUE,
                          inbreed = FALSE,
                          selfing = FALSE,
                          genmethod = "top",
                          ...) {
  # 1. Parameter Validation
  if (missing(ped) || is.null(ped)) {
    stop("The ped parameter cannot be NULL or missing")
  }
  if ((is.data.frame(ped) || is.matrix(ped)) && NROW(ped) == 0) {
    stop("The ped parameter cannot be empty")
  }
  
  # Boolean checks
  for (arg in c("addgen", "addnum", "inbreed", "selfing")) {
    val <- get(arg)
    if (!is.logical(val) || length(val) != 1 || is.na(val)) {
      stop(sprintf("The %s parameter only is assigned using TRUE or FALSE", arg))
    }
  }

  valid_trace <- c("up", "down", "all")
  if (length(trace) != 1 || !trace %in% valid_trace) {
    stop(sprintf("The trace parameter must be one of: %s", paste(valid_trace, collapse = ", ")))
  }

  valid_genmethod <- c("top", "bottom")
  if (length(genmethod) != 1 || !genmethod %in% valid_genmethod) {
    stop(sprintf("The genmethod parameter must be one of: %s", paste(valid_genmethod, collapse = ", ")))
  }
  
  if (!is.null(tracegen)) {
    if (!is.numeric(tracegen) || length(tracegen) != 1 || is.na(tracegen) || !is.finite(tracegen)) {
      stop("The tracegen parameter must be a single numeric value")
    }
    tracegen <- as.integer(tracegen)
  }

  # 2. Data Preparation
  res_prep <- validate_and_prepare_ped(ped, selfing = selfing)
  ped_dt <- res_prep$ped_dt
  bisexual_parents <- res_prep$bisexual_parents
  
  # 3. Build Initial Graph
  # We need the graph for cycle detection and tracing
  graph_res <- build_ped_graph(ped_dt)
  g <- graph_res$g
  
  # 4. Cycle Detection
  check_ped_loops(g)
  
  # 5. Trace Candidates (if needed)
  if (!is.null(cand)) {
    ped_dt <- trace_ped_candidates(g, ped_dt, cand, trace, tracegen)
    
    # Rebuild graph for the subsetted pedigree
    # This is necessary because generation assignment depends on the structure of the subset
    graph_res <- build_ped_graph(ped_dt)
    g <- graph_res$g
  }
  
  # 6. Generation Assignment
  # We need topo_order for sorting even if addgen=FALSE
  topo_order <- as.integer(topo_sort(g))
  
  # 6.5. Pre-calculate Family (moved from step 8.5) to assist generation assignment
  # Create family identifier based on parent combination
  ped_dt[, Family := ifelse(
    !is.na(Sire) & !is.na(Dam),
    paste0(Sire, "x", Dam),
    NA_character_
  )]
  
  # Calculate family sizes
  ped_dt[, FamilySize := .N, by = Family]
  # Individuals without both parents have FamilySize = 1
  ped_dt[is.na(Family), FamilySize := 1L]
  
  if (addgen) {
    ped_dt <- assign_ped_generations(g, ped_dt, topo_order, genmethod)
    
    # Optimization for bottom-up method:
    # 1. Align full siblings to the same generation (the highest/minimum Gen among them)
    # This prevents leaves (individuals with no progeny) from dropping to the bottom
    # when they have siblings who are parents of deep lineages.
    if (genmethod == "bottom") {
      ped_dt[!is.na(Family), Gen := min(Gen), by = Family]
      
      # 2. Align mates to the same generation
      # If Sire and Dam have different generations, align them to the earlier (min) generation.
      # This ensures pedigree plots look cleaner with horizontal mating lines.
      # We iterate this until stability or max depth 
      # (usually one pass is sufficient for most cases, but overlapping generations might need more).
      # But for speed, we do a single efficient pass that should catch the 'I' vs 'H' case.
      
      # We need to lookup mate generations. Join self to self on Sire/Dam.
      # Create a small lookup table of Ind -> Gen
      # Using match() is fast.
      
      # Iterate a few times to propagate changes.
      # For standard pedigrees, 1 or 2 passes cover the "mate alignment".
      # Using a while loop with a safe maximum to handle deep recursive alignments
      iter <- 0
      max_iter <- 10
      repeat {
        iter <- iter + 1
        # Get generations of parents for every individual
        sire_gen <- ped_dt$Gen[match(ped_dt$Sire, ped_dt$Ind)]
        dam_gen <- ped_dt$Gen[match(ped_dt$Dam, ped_dt$Ind)]
        
        # Identify families where SireGen != DamGen (and both exist)
        valid_mask <- !is.na(sire_gen) & !is.na(dam_gen) & (sire_gen != dam_gen)
        
        if (!any(valid_mask) || iter > max_iter) break
        
        # Get the target maximum generation for each mismatched couple
        # We push the "older/higher" parent DOWN to the "younger/lower" parent's generation.
        # This is safe from parent constraints (moving away from parents), but risky for child constraints.
        target_gen <- pmax(sire_gen[valid_mask], dam_gen[valid_mask])
        
        # Collect updates needed: (Ind, NewGen)
        sires_to_upd <- ped_dt$Sire[valid_mask]
        dams_to_upd <- ped_dt$Dam[valid_mask]
        
        upd_dt <- data.table(
          Ind = c(sires_to_upd, dams_to_upd),
          NewGen = c(target_gen, target_gen)
        )
        
        # Resolve multiple updates (take max to go deeper)
        upd_dt <- upd_dt[, .(NewGen = max(NewGen)), by = Ind]
        
        # Constraint check: Individuals cannot move down if they hit their children.
        # NewGen must be < Min(ChildGen).
        inds_vec <- upd_dt$Ind
        relevant_kids <- ped_dt[Sire %in% inds_vec | Dam %in% inds_vec, 
                                .(Sire, Dam, Gen)]
        
        limits <- rbind(
          relevant_kids[Sire %in% inds_vec, .(Ind=Sire, Limit=Gen)],
          relevant_kids[Dam %in% inds_vec, .(Ind=Dam, Limit=Gen)]
        )[, .(Limit = min(Limit)), by = Ind]
        
        limit_vals <- limits$Limit[match(upd_dt$Ind, limits$Ind)]
        
        curr_gen <- ped_dt$Gen[match(upd_dt$Ind, ped_dt$Ind)]
        
        # Valid only if pushing down (New > Curr) AND valid against children
        to_change <- (upd_dt$NewGen > curr_gen) & (is.na(limit_vals) | upd_dt$NewGen < limit_vals)
        
        if(any(to_change)) {
           update_inds <- upd_dt$Ind[to_change]
           update_gens <- upd_dt$NewGen[to_change]
           ped_dt[Ind %in% update_inds, Gen := update_gens[match(Ind, update_inds)]]
        } else {
           break
        }
      }
      
      # Final Sibling Alignment check:
      # Mate alignment might have pushed one sibling down (away from others) to match a mate.
      # But Sibling Consistency (P1) > Mate Consistency (P2). 
      # So we enforce sibling alignment again.
      ped_dt[!is.na(Family), Gen := min(Gen), by = Family]
    }
  }
  
  # 7. Sex Inference and Check
  ped_dt <- infer_and_check_sex(ped_dt, selfing = selfing)
  
  # 8. Sorting and Numeric IDs
  if (addgen) {
    setorder(ped_dt, Gen, Ind)
  } else {
    # If no Gen, preserve topo sort order
    sorted_inds <- V(g)$name[topo_order]
    ped_dt <- ped_dt[match(sorted_inds, Ind)]
  }
  
  if (addnum) {
    ped_dt[, IndNum := .I]
    ped_dt[, SireNum := match(Sire, Ind, nomatch = 0)]
    ped_dt[, DamNum := match(Dam, Ind, nomatch = 0)]
  }
  
  if (!is.null(cand)) {
    ped_dt[, Cand := Ind %in% cand]
  }

  # 9. Final S3 and Inbreeding
  ped_dt <- new_tidyped(ped_dt)
  attr(ped_dt, "bisexual_parents") <- bisexual_parents
  attr(ped_dt, "selfing") <- selfing
  
  if (inbreed) {
    ped_dt <- inbreed(ped_dt, ...)
  }
  
  return(ped_dt)
}

# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Validate and Prepare Pedigree Data
#' @noRd
validate_and_prepare_ped <- function(ped, selfing = FALSE) {
  ped_dt <- if (is.data.table(ped)) copy(ped) else as.data.table(ped)
  
  # Standardize primary column names
  setnames(ped_dt, 1:3, c("Ind", "Sire", "Dam"))
  
  # Ensure character types for IDs
  ped_dt[, c("Ind", "Sire", "Dam") := lapply(.SD, as.character), .SDcols = c("Ind", "Sire", "Dam")]
  
  # Remove records with missing Ind
  if (any(ped_dt$Ind %in% c("", " ", "0", "*", "NA", NA))) {
    warning("Missing values in 'Ind' column; records discarded.")
    ped_dt <- ped_dt[!Ind %in% c("", " ", "0", "*", "NA", NA)]
  }
  
  # Normalize missing parents as NA
  ped_dt[Sire %in% c("", " ", "0", "*", "NA"), Sire := NA_character_]
  ped_dt[Dam %in% c("", " ", "0", "*", "NA"), Dam := NA_character_]
  
  if (all(is.na(ped_dt$Sire)) && all(is.na(ped_dt$Dam))) {
    stop("All parents are missing; no pedigree can be built.")
  }

  # Duplicate checks
  if (anyDuplicated(ped_dt, by = c("Ind", "Sire", "Dam")) > 0) {
    n_duped <- sum(duplicated(ped_dt, by = c("Ind", "Sire", "Dam")))
    warning(sprintf("Removed %d records with duplicate Ind/Sire/Dam IDs.", n_duped))
    ped_dt <- unique(ped_dt, by = c("Ind", "Sire", "Dam"))
  }
  
  if (anyDuplicated(ped_dt, by = "Ind") > 0) {
    stop("Fatal error: Duplicate Ind IDs with different parents found!")
  }

  # Sex conflict check: same individual appears as both sire and dam
  sires <- unique(ped_dt[!is.na(Sire), Sire])
  dams <- unique(ped_dt[!is.na(Dam), Dam])
  bisexual_parents <- sort(intersect(sires, dams))
  
  if (length(bisexual_parents) > 0 && !selfing) {
    stop(sprintf(
      paste0("Sex conflict detected: The following individual(s) appear as both Sire and Dam: %s. ",
             "This is biologically impossible. Please check and correct the pedigree data. ",
             "If this is a plant pedigree with monoecious species, set selfing = TRUE."),
      paste(bisexual_parents, collapse = ", ")
    ), call. = FALSE)
  }
  
  if (length(bisexual_parents) > 0 && selfing) {
    message(sprintf(
      "Selfing mode: %d individual(s) appear as both Sire and Dam: %s. These will be assigned Sex = 'monoecious'.",
      length(bisexual_parents),
      paste(head(bisexual_parents, 10), collapse = ", ")
    ))
  }

  # Add missing founders
  all_parents <- unique(c(sires, dams))
  missing_founders <- setdiff(all_parents, ped_dt$Ind)
  
  if (length(missing_founders) > 0) {
    founder_dt <- data.table(Ind = missing_founders, Sire = NA_character_, Dam = NA_character_)
    ped_dt <- rbind(founder_dt, ped_dt, fill = TRUE)
  }
  
  list(ped_dt = ped_dt, bisexual_parents = bisexual_parents)
}

#' Build igraph Object from Pedigree
#' @noRd
build_ped_graph <- function(ped_dt) {
  node_names <- ped_dt$Ind
  
  # Use match() instead of named vector lookup for speed
  s_from <- match(ped_dt$Sire[!is.na(ped_dt$Sire)], node_names)
  s_to <- which(!is.na(ped_dt$Sire))
  
  d_from <- match(ped_dt$Dam[!is.na(ped_dt$Dam)], node_names)
  d_to <- which(!is.na(ped_dt$Dam))
  
  edge_mat <- matrix(c(s_from, d_from, s_to, d_to), ncol = 2)
  
  # Using graph_from_edgelist is faster usually
  g <- graph_from_edgelist(edge_mat, directed = TRUE)
  
  # Ensure all vertices exist (even isolated ones)
  num_needed <- length(node_names)
  if (vcount(g) < num_needed) {
    g <- add_vertices(g, num_needed - vcount(g))
  }
  
  V(g)$name <- node_names
  
  list(g = g)
}

#' Check for Cycles/Loops in Pedigree
#' @noRd
check_ped_loops <- function(g) {
  if (!is_dag(g)) {
    # 1. Check for self-loops (A -> A)
    # which_loop returns TRUE for edges that are self-loops
    loop_edges <- which(which_loop(g))
    self_loops <- character()
    if (length(loop_edges) > 0) {
      # ends() returns the vertex IDs of edge endpoints
      self_loops <- unique(V(g)[ends(g, loop_edges)[, 1]]$name)
    }
    
    # 2. Check for larger cycles via SCC
    scc <- components(g, mode = "strong")
    cyclic_nodes_idx <- which(scc$csize > 1)
    
    cycle_reports <- character()
    
    if (length(self_loops) > 0) {
      cycle_reports <- c(cycle_reports, paste(self_loops, "->", self_loops, "(self-loop)"))
    }
    
    if (length(cyclic_nodes_idx) > 0) {
      scc_reports <- sapply(cyclic_nodes_idx, function(i) {
        nodes_in_scc <- names(which(scc$membership == i))
        v_start <- nodes_in_scc[1]
        successors <- intersect(names(neighbors(g, v_start, mode = "out")), nodes_in_scc)
        
        if (length(successors) > 0) {
          v_next <- successors[1]
          p <- shortest_paths(g, from = v_next, to = v_start, mode = "out")$vpath[[1]]
          path_names <- c(v_start, names(p))
          paste(path_names, collapse = " -> ")
        } else {
          paste(nodes_in_scc, collapse = ", ")
        }
      })
      cycle_reports <- c(cycle_reports, scc_reports)
    }
    
    # Fallback
    if (length(cycle_reports) == 0) {
      cycle_reports <- "Complex cycle structure detected (try checking strong components)."
    }
    
    stop(paste("Pedigree error! Pedigree loops detected:\n", 
               paste(cycle_reports, collapse = "\n")), call. = FALSE)
  }
}

#' Trace Candidates and Subset Pedigree
#' @noRd
trace_ped_candidates <- function(g, ped_dt, cand, trace, tracegen) {
  # Drop NAs and empty strings from input
  cand_clean <- as.character(cand[!is.na(cand) & cand != "" & cand != " "])
  
  if (length(cand_clean) == 0) {
    stop("The cand parameter contains no valid individual IDs.")
  }

  valid_cand <- cand_clean[cand_clean %in% ped_dt$Ind]
  if (length(valid_cand) == 0) {
    stop("None of the specified candidates were found in the pedigree.")
  }

  if (length(valid_cand) < length(cand_clean)) {
    missing_cands <- setdiff(cand_clean, valid_cand)
    warning(sprintf("The following %d candidates were not found in the pedigree: %s", 
                    length(missing_cands), 
                    paste(head(missing_cands, 5), collapse = ", ")))
  }
  
  # Determine recursion depth
  # order must be numeric and >= 0 for neighborhood()
  order <- if (is.null(tracegen) || tracegen <= 0) vcount(g) else as.integer(tracegen)
  
  # For mode "all", we want the union of ancestors (in) and descendants (out).
  # We should NOT use igraph's mode="all" which is undirected search (connected component).
  if (trace == "all") {
    rel_in <- neighborhood(g, order = order, nodes = valid_cand, mode = "in")
    rel_out <- neighborhood(g, order = order, nodes = valid_cand, mode = "out")
    # Efficiently flatten and get unique IDs from both directions
    keep_indices <- unique(c(
      unlist(lapply(rel_in, as.integer)),
      unlist(lapply(rel_out, as.integer))
    ))
  } else {
    m <- if (trace == "up") "in" else "out"
    reaching_nodes_list <- neighborhood(g, order = order, nodes = valid_cand, mode = m)
    # Efficiently flatten and get unique IDs
    # lapply(list, as.integer) converts vertex sequences to integer IDs
    keep_indices <- unique(unlist(lapply(reaching_nodes_list, as.integer)))
  }
  
  keep_inds <- V(g)$name[keep_indices]
  
  ped_dt <- ped_dt[Ind %in% keep_inds]
  # Update parents not in set after tracing
  ped_dt[!(Sire %in% keep_inds), Sire := NA_character_]
  ped_dt[!(Dam %in% keep_inds), Dam := NA_character_]
  
  return(ped_dt[])
}

#' Assigns individual generation numbers based on topological sorting and parentage.
#'
#' @param g An igraph object representing the pedigree.
#' @param ped_dt A data.table containing the pedigree.
#' @param topo_order Integer vector specifying the topological order of vertices.
#' @param genmethod Character, "top" or "bottom".
#' @return A data.table with a \strong{Gen} column added.
#' @noRd
assign_ped_generations <- function(g, ped_dt, topo_order, genmethod) {
  # Get numeric parents for C++ (1-based, 0 for missing)
  # Pedigree is already complete (missing founders added)
  inds <- ped_dt$Ind
  sire_num <- match(ped_dt$Sire, inds, nomatch = 0)
  dam_num <- match(ped_dt$Dam, inds, nomatch = 0)
  
  if (genmethod == "top") {
    gen_vec <- cpp_assign_generations_top(sire_num, dam_num, as.integer(topo_order))
  } else {
    gen_vec <- cpp_assign_generations_bottom(sire_num, dam_num, as.integer(topo_order))
  }
  
  # Isolated nodes check
  in_deg <- degree(g, mode = "in")
  out_deg <- degree(g, mode = "out")
  iso_mask <- (in_deg == 0) & (out_deg == 0)
  if (any(iso_mask)) {
    # Match iso_mask (which is in graph order) to ped_dt order
    # topo_order is a permutation of 1:vcount(g)
    # The C++ gen_vec is in the same order as sire_num/dam_num, which is ped_dt order
    # BUT wait, the graph g was built from ped_dt. So V(g)$name is ped_dt$Ind.
    # So vertex i in the graph corresponds to row i in ped_dt.
    gen_vec[iso_mask] <- 0L
  }
  
  ped_dt[, Gen := gen_vec]
  return(ped_dt[])
}

#' Infer and Check Sex of Individuals
#' @noRd
infer_and_check_sex <- function(ped_dt, selfing = FALSE) {
  if (!("Sex" %in% names(ped_dt))) {
    ped_dt[, Sex := NA_character_]
  } else {
    ped_dt[, Sex := tolower(as.character(Sex))]
    ped_dt[Sex %in% c("", " ", "na"), Sex := NA_character_]
  }
  
  sires <- unique(ped_dt[!is.na(Sire), Sire])
  dams <- unique(ped_dt[!is.na(Dam), Dam])
  
  # Identify monoecious individuals (appear as both Sire and Dam)
  monoecious_ids <- intersect(sires, dams)
  
  if (selfing && length(monoecious_ids) > 0) {
    # In selfing mode: monoecious individuals can be both Sire and Dam
    # Check for conflicts with explicit sex annotation (only non-monoecious annotations)
    sex_conflicts <- ped_dt[
      !(Ind %in% monoecious_ids) & 
      ((Ind %in% sires & Sex == "female") | (Ind %in% dams & Sex == "male")), Ind]
    if (length(sex_conflicts) > 0) {
      stop(sprintf(
        paste0("Sex annotation conflicts detected for: %s. ",
               "These individuals have explicit sex that contradicts their role as Sire/Dam. ",
               "Please check and correct the pedigree data."),
        paste(sex_conflicts, collapse = ", ")
      ), call. = FALSE)
    }
    
    # Check that monoecious individuals are not explicitly annotated as male or female
    mono_sex_conflicts <- ped_dt[
      Ind %in% monoecious_ids & Sex %in% c("male", "female"), Ind]
    if (length(mono_sex_conflicts) > 0) {
      stop(sprintf(
        paste0("Sex annotation conflicts for monoecious individuals: %s. ",
               "These individuals appear as both Sire and Dam but have explicit male/female annotation. ",
               "Remove or change their Sex annotation to 'monoecious' or NA."),
        paste(mono_sex_conflicts, collapse = ", ")
      ), call. = FALSE)
    }
    
    # Assign monoecious sex
    ped_dt[Ind %in% monoecious_ids & (is.na(Sex) | Sex != "monoecious"), Sex := "monoecious"]
    
    # Infer sex for remaining individuals (excluding monoecious)
    sires_only <- setdiff(sires, monoecious_ids)
    dams_only <- setdiff(dams, monoecious_ids)
    ped_dt[is.na(Sex) & (Ind %in% sires_only), Sex := "male"]
    ped_dt[is.na(Sex) & (Ind %in% dams_only), Sex := "female"]
  } else {
    # Standard mode: no selfing allowed
    # Sex conflicts with explicit sex annotation
    sex_conflicts <- ped_dt[(Ind %in% sires & Sex == "female") | (Ind %in% dams & Sex == "male"), Ind]
    if (length(sex_conflicts) > 0) {
      stop(sprintf(
        paste0("Sex annotation conflicts detected for: %s. ",
               "These individuals have explicit sex that contradicts their role as Sire/Dam. ",
               "Please check and correct the pedigree data."),
        paste(sex_conflicts, collapse = ", ")
      ), call. = FALSE)
    }
    
    # Infer sex from roles
    ped_dt[is.na(Sex) & (Ind %in% sires), Sex := "male"]
    ped_dt[is.na(Sex) & (Ind %in% dams), Sex := "female"]
  }
  
  return(ped_dt[])
}

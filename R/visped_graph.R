#' Convert pedigree to igraph structure
#' @import data.table
#' @keywords internal
ped2igraph <- function(ped, compact = FALSE, highlight = NULL, trace = FALSE, showf = FALSE) {
  if (nrow(ped) == 0) {
    return(list(
      node = data.table(id = integer(), nodetype = character(), gen = integer(), layer = numeric(), label = character()),
      edge = data.table(from = integer(), to = integer())
    ))
  }
  
  # 1. Inject missing parents (for subsetted pedigrees)
  ped_new <- inject_missing_parents(ped)
  
  # 2. Resolve highlight IDs and relatives
  highlight_info <- get_highlight_ids(ped_new, highlight, trace)
  h_ids <- highlight_info$all_ids
  
  # 3. Prepare initial node table
  ped_node <- prepare_initial_nodes(ped_new)
  
  # 4. Compact pedigree (if requested)
  ped_node <- compact_pedigree(ped_node, compact, h_ids)
  
  # 5. Generate edges and virtual nodes
  graph_struct <- generate_graph_structure(ped_node, h_ids)
  ped_node <- graph_struct$node
  ped_edge <- graph_struct$edge
  
  # 6. Apply styles (colors, shapes, highlighting)
  ped_node <- apply_node_styles(ped_node, highlight_info)
  
  # 7. Finalize graph (reindex IDs, set edge colors)
  result <- finalize_graph(ped_node, ped_edge, highlight_info, trace, showf)
  
  return(result)
}

#' Inject missing parents for subsetted pedigrees
#' @param ped A data.table containing pedigree info.
#' @keywords internal
inject_missing_parents <- function(ped) {
  ped_new <- copy(ped)
  # Ensure keys for fast matching
  setkey(ped_new, IndNum)
  all_ind_nums <- ped_new$IndNum
  
  # Find missing Sires
  s_nums <- ped_new$SireNum
  missing_sire_nums <- unique(s_nums[s_nums != 0 & !is.na(s_nums) & !(s_nums %in% all_ind_nums)])
  
  # Find missing Dams (compute early to detect monoecious)
  d_nums <- ped_new$DamNum
  missing_dam_nums <- unique(d_nums[d_nums != 0 & !is.na(d_nums) & !(d_nums %in% all_ind_nums)])
  
  # Detect monoecious: individuals that are both missing sires and missing dams
  monoecious_nums <- intersect(missing_sire_nums, missing_dam_nums)
  
  if (length(missing_sire_nums) > 0) {
    sire_info <- unique(ped_new[SireNum %in% missing_sire_nums, .(SireNum, Sire)])
    # Determine sex: monoecious if also a dam, otherwise male
    sire_sex <- ifelse(sire_info$SireNum %in% monoecious_nums, "monoecious", "male")
    new_sires <- data.table(
      Ind = sire_info$Sire,
      Sire = NA_character_,
      Dam = NA_character_,
      Sex = sire_sex,
      Gen = min(ped_new$Gen, na.rm=TRUE) - 1,
      IndNum = sire_info$SireNum,
      SireNum = 0L,
      DamNum = 0L
    )
    extra_cols <- setdiff(names(ped_new), names(new_sires))
    if (length(extra_cols) > 0) new_sires[, (extra_cols) := NA]
    ped_new <- rbind(ped_new, new_sires, fill = TRUE)
    setkey(ped_new, IndNum)
    all_ind_nums <- ped_new$IndNum
  }

  # Find missing Dams (exclude monoecious already added as sires)
  missing_dam_nums_only <- setdiff(missing_dam_nums, monoecious_nums)
  
  if (length(missing_dam_nums_only) > 0) {
    dam_info <- unique(ped_new[DamNum %in% missing_dam_nums_only, .(DamNum, Dam)])
    new_dams <- data.table(
      Ind = dam_info$Dam,
      Sire = NA_character_,
      Dam = NA_character_,
      Sex = "female",
      Gen = min(ped_new$Gen, na.rm=TRUE) - 1,
      IndNum = dam_info$DamNum,
      SireNum = 0L,
      DamNum = 0L
    )
    extra_cols <- setdiff(names(ped_new), names(new_dams))
    if (length(extra_cols) > 0) new_dams[, (extra_cols) := NA]
    ped_new <- rbind(ped_new, new_dams, fill = TRUE)
  }
  
  return(ped_new[])
}

#' Prepare initial node table for igraph conversion
#' @param ped A data.table containing pedigree info.
#' @keywords internal
prepare_initial_nodes <- function(ped) {
  cols <- c("IndNum", "Ind", "SireNum", "DamNum", "Sire", "Dam", "Sex", "Gen")
  if ("Cand" %in% colnames(ped)) cols <- c(cols, "Cand")
  if ("f" %in% colnames(ped)) cols <- c(cols, "f")
  
  ped_node <- ped[, ..cols]
  setnames(ped_node, 
           old = c("IndNum", "Ind", "SireNum", "DamNum", "Sire", "Dam", "Sex", "Gen"),
           new = c("id", "label", "sirenum", "damnum", "sirelabel", "damlabel", "sex", "gen"))
  
  if ("Cand" %in% colnames(ped)) setnames(ped_node, "Cand", "cand")
  
  max_id <- max(ped_node$id, na.rm = TRUE)
  familylabel = NULL
  ped_node[!(is.na(sirelabel) & is.na(damlabel)),
           familylabel := paste(sirelabel, damlabel, sep = "x")]
  
  ped_node[!is.na(familylabel), 
           familynum := .GRP + max_id, 
           by = familylabel]
  ped_node[is.na(familynum), familynum := 0]
  
  nodetype = NULL
  ped_node[, nodetype := "real"]
  
  return(ped_node[])
}

#' Compact pedigree by merging full siblings
#' @param ped_node A data.table of nodes.
#' @param compact Logical, whether to compact.
#' @param h_ids Highlighted IDs to exempt from compaction.
#' @keywords internal
compact_pedigree <- function(ped_node, compact, h_ids) {
  if (!compact) return(ped_node)
  
  # Identify parents (individuals who appear in sirelabel or damlabel)
  # Efficiently find unique parents
  parents <- unique(c(ped_node$sirelabel, ped_node$damlabel))
  parents <- parents[!is.na(parents)]
  
  # Candidates for compaction: NOT parents, HAS a family, and NOT highlighted
  # Using %in% on unique labels is usually fast in data.table
  ped_node_1 <- ped_node[!(label %in% parents) & !is.na(familylabel)]
  
  if (!is.null(h_ids) && length(h_ids) > 0) {
    ped_node_1 <- ped_node_1[!(label %in% h_ids)]
  }
  
  if (nrow(ped_node_1) > 0) {
    familysize <- NULL
    # Efficiently count family sizes for eligible nodes
    ped_node_1[, familysize := .N, by = .(familylabel)]
    fullsib_id_DT <- ped_node_1[familysize >= 2]
    
    if (nrow(fullsib_id_DT) > 0) {
      # Take one representative from each family
      compact_family <- unique(fullsib_id_DT, by = c("familylabel"))
      compact_family[, `:=`(
        label = as.character(familysize),
        nodetype = "compact",
        sex = NA_character_
      )]
      
      next_id <- max(ped_node$id, na.rm = TRUE)
      compact_family[, id := seq.int(next_id + 1, length.out = .N)]
      
      # Mapping from old individual IDs to new compact family node IDs
      map_dt <- fullsib_id_DT[, .(from_id = id, familylabel)]
      map_dt <- compact_family[map_dt, on = .(familylabel), .(from_id, to_id = id, to_label = label)]
      
      # Remove individuals that were compacted
      ped_node <- ped_node[!(id %in% fullsib_id_DT$id)]
      
      # Update children of these individuals to point to the compact node instead
      if (nrow(map_dt) > 0) {
        # Note: In practice, these individuals are NOT parents (vetted above),
        # but we keep this for robustness if the logic ever changes.
        ped_node[map_dt, on = .(sirenum = from_id), `:=`(sirenum = to_id, sirelabel = to_label)]
        ped_node[map_dt, on = .(damnum = from_id), `:=`(damnum = to_id, damlabel = to_label)]
      }
      
      ped_node <- rbind(ped_node, compact_family, fill = TRUE)
      
      # Re-calculate family numbers to include new compact nodes
      max_id <- max(ped_node$id, na.rm = TRUE)
      ped_node[!is.na(familylabel), 
               familynum := .GRP + max_id, 
               by = familylabel]
      ped_node[is.na(familynum), familynum := 0]
    }
  }
  return(ped_node[])
}

#' Generate edges and virtual family nodes
#' @param ped_node A data.table of nodes.
#' @param h_ids Highlighted IDs.
#' @keywords internal
generate_graph_structure <- function(ped_node, h_ids) {
  # Individual → family edges
  edge_ind <- ped_node[familynum > 0, .(from = id, to = familynum, role = "ind")]

  # Family → sire edges (dedup siblings in same family)
  edge_sire <- unique(
    ped_node[sirenum > 0 & familynum > 0, .(from = familynum, to = sirenum, role = "sire")]
  )

  # Family → dam edges (dedup siblings in same family)
  edge_dam <- unique(
    ped_node[damnum > 0 & familynum > 0, .(from = familynum, to = damnum, role = "dam")]
  )

  # Detect selfing: same (from, to) appearing in both sire and dam roles
  edge_parents <- rbind(edge_sire, edge_dam)
  dup_ft <- edge_parents[, .N, by = .(from, to)]
  selfing_ft <- dup_ft[N > 1]
  if (nrow(selfing_ft) > 0) {
    edge_parents[selfing_ft, on = .(from, to), role := "selfing"]
    edge_parents <- unique(edge_parents, by = c("from", "to"))
  }

  ped_edge <- rbind(edge_ind, edge_parents)
  ped_edge <- ped_edge[order(from, to)]

  width = arrow.size = arrow.width = color = curved = NULL
  edge_default_color <- if (length(h_ids) > 0) "#3333334D" else "#333333"
  ped_edge[, ":="(width = 1, arrow.size = 1, arrow.width = 1, arrow.mode = 2, 
                  color = edge_default_color, curved = 0.10)]
  
  virtual_nodes <- unique(ped_node[familynum > 0, .(
    id = familynum,
    familylabel,
    label = familylabel,
    sirenum,
    damnum,
    sirelabel,
    damlabel,
    gen,
    familynum
  )])
  
  ped_node <- rbind(ped_node, virtual_nodes, fill = TRUE)
  ped_node[is.na(nodetype), nodetype := "virtual"]
  
  layer = NULL
  ped_node[nodetype %in% c("real", "compact"), layer := 2 * gen - 1]
  ped_node[nodetype %in% c("virtual"), layer := 2 * (gen - 1)]
  
  list(node = ped_node, edge = ped_edge)
}

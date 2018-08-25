#' @import visNetwork
#' @import data.table

visped <- function(ped,
                   fullsib = TRUE,
                   layout = NULL) {
  nodes <- ped2node(ped)
  ind_num_gen <- nodes$node[, .N, by = .(gen)]
  ind_num_to <- nodes$edge[, .N, by= .(to)]
  nodes$node <- merge(nodes$node,ind_num_to,by.x=c("id"),by.y=c("to"),all.x = TRUE)
  edge_to_num <- merge(nodes$edge,nodes$node[,.(id,N)],by.x="to",by.y="id",all.x = TRUE)
  edge_to_num_by_from <- edge_to_num[,.(meanN = mean(N)),by=.(from)]
  nodes$node <- merge(nodes$node,edge_to_num_by_from,by.x="id",by.y = "from",all.x = TRUE)
  nodes$node <- nodes$node[order(gen, -N, -meanN)]

  width <- max(ind_num_gen$N) * 100
  height <- max(ind_num_gen$gen) * 150
  if (width > 1920*3) {width <- 1920*3}
  if (height > 1080*3) {height <- 1080*3}

  if (height < width) {
    height <- round(width / 1.777778, 8)
  }

  if (height >= width) {
    width <- round(height * 1.7, 8)
  }
  gen_max <- max(ind_num_gen$gen)

  level_space = round(height / gen_max,8)*0.99
  if (gen_max %% 2 == 0) {
    gen_node_y_up <-
      seq(
        from = round(level_space / 2, 8),
        to = round(height / 2,8),
        by = level_space
      )
    gen_node_y_down <-
      seq(
        from = (-1) * round(level_space / 2, 8),
        to = (-1) * round(height / 2,8),
        by = (-1) * level_space
      )
    gen_node_y <-
      sort(c(gen_node_y_down[1:(gen_max / 2)],
             gen_node_y_up[1:(gen_max / 2)]))
  } else {
    gen_node_y_up <- seq(from = 0,
                         to = round(height / 2,8),
                         by = level_space)
    gen_node_y_down <-
      seq(
        from = 0,
        to = (-1) * round(height / 2,8),
        by = (-1) * level_space
      )
    gen_node_y <-
      sort(unique(c(gen_node_y_down[2:((gen_max + 1) / 2)],
                    0,
                    gen_node_y_up[2:((gen_max + 1) / 2)])))
  }
  for (g in 1:gen_max) {
    node_num_per_gen <- nrow(nodes$node[gen == g])
    node_space = round(width / node_num_per_gen, 8)*0.99
    if (node_num_per_gen %% 2 == 0) {
      node_x_right <-
        seq(
          from = round(node_space / 2, 8),
          to = round(width / 2,8),
          by = node_space
        )
      node_x_left <-
        seq(
          from = (-1) * round(node_space / 2, 8),
          to = (-1) * round(width / 2, 8),
          by = (-1) * node_space
        )
      node_x <- sort(c(node_x_left[1:(node_num_per_gen / 2)],
                       node_x_right[1:(node_num_per_gen / 2)]))
    } else {
      node_x_right <- seq(from = 0,
                          to = round(width / 2, 8),
                          by = node_space)
      node_x_left <-
        seq(
          from = 0,
          to = (-1) * round(width / 2,8),
          by = (-1) * node_space
        )

      node_x <-
        sort(unique(c(node_x_left[2:((node_num_per_gen + 1) / 2)], 0,
                      node_x_right[2:((node_num_per_gen + 1) / 2)])))
    }
    node_y <- rep(gen_node_y[g], node_num_per_gen)
    nodes$node[gen == g, ":="(x = node_x, y = node_y,
                              size=node_space*0.5)]
  }
  nodes$node[size>25,size:=25]
  if (sum(ind_num_gen$N,na.rm = TRUE) <= 200) {
    visNetwork(
      nodes$node,
      nodes$edge,
      width = paste(width, "px", sep = ""),
      height = paste(height, "px", sep = "")
    ) %>%
      visNodes(physics = FALSE,scaling = list(label=list(enabled=TRUE))) %>%
      visEdges(arrows = "to", smooth = TRUE) %>%
      visPhysics(stabilization = FALSE) %>%
      visInteraction(dragNodes = TRUE,
                     dragView = TRUE,
                     zoomView = TRUE)
  } else {
      visNetwork(
        nodes$node,
        nodes$edge,
        width = paste(width, "px", sep = ""),
        height = paste(height, "px", sep = "")
      ) %>%
        visNodes(physics = FALSE,scaling = list(label=list(enabled=TRUE))) %>%
        visEdges(physics = FALSE, arrows = "to", smooth = FALSE) %>%
        visPhysics(stabilization = FALSE) %>%
        visInteraction(dragNodes = FALSE,
                       dragView = TRUE,
                       zoomView = TRUE)
  }
}

ped2node <- function(ped) {
  ped_col_names <- names(ped)
  if (!all(c("IndNum","SireNum","DamNum") %in% ped_col_names)) {
    stop("The pedigree need to be firstly trimmed by the tidyped() function!")
  }
  ped_new <- copy(ped)
  if (c("Cand") %in% ped_col_names) {
    ped_node <- ped_new[,.(id=IndNum,label=Ind,cand=Cand,sex=Sex,gen=Gen)]
  } else {
    ped_node <- ped_new[,.(id=IndNum,label=Ind,sex=Sex,gen=Gen)]
  }

  ped_node[,shape:="diamond"]
  ped_node[,size:=25]
  ped_node[sex %in% c("male"),shape:="square"]
  ped_node[sex %in% c("female"),shape:="dot"]

  ped_node[sex %in% c("male"),color:="#119ecc"]
  ped_node[sex %in% c("female"),color:="#f4b131"]

  if (c("Cand") %in% ped_col_names) {
    ped_node[cand == TRUE,color:="red"]
  }

  ped_edge <- rbind(ped_new[,.(from=IndNum,to=SireNum)],ped_new[,.(from=IndNum,to=DamNum)])
  ped_edge <- ped_edge[!(to==0)]
  return(list(node=ped_node,edge=ped_edge))
}


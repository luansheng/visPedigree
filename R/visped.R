#' @import data.table
#' @import igraph
#' @export

`:=` = function(...) NULL

visped <- function(ped,
                   compact = FALSE, cex=NULL, show=TRUE, file="output.pdf") {
  ped_new <- copy(ped)
  # Convertting pedigree to node and edge two data.table
  ped_igraph <- ped2igraph(ped_new, compact)
  real_node <- ped_igraph$node[nodetype %in% c("real", "compact")]
  gen_node_num <- real_node[, .N, by = gen]
  gen_max_size <-  max(gen_node_num$N, na.rm = TRUE)


  #=== Obtaining the width of a node's label =====================================================
  # The maximum width of a pdf file is 200 inch
  pdf_max_width = 200
  cexs <-  seq(from = 0.1, to = 1, by = 0.05)
  best_cex <- 0
  max_strwidth_label <- real_node[which.max(strwidth(real_node$label, cex = 1, units = "inches")), label]
  for (i in length(cexs):1) {
    # Obtaining the maximum width of a node's label: inch
    label_max_width <- max(strwidth(max_strwidth_label, cex = cexs[i], units = "inches"),
          na.rm = TRUE)
    # Fixing the width of the node when the number of nodes (individuals) in one generation is small
    # The unit of 0.8 is inch, about 2cm for the width of one node
    if (gen_max_size <= 16 & label_max_width < 0.8) {
      label_max_width = 0.8
    }
    if ((label_max_width * gen_max_size) < pdf_max_width) {
      best_cex <- cexs[i] * 0.7
      break
    }
  }
  if (best_cex == 0) {
    stop(
      "Too many individuals (>=",
      gen_max_size,
      paste(") in one genration!!! please remove full-sib individuals using the parameter compact = TRUE,",
            "rerun visped() function",sep=" ")
    )
  }

  #=== Generating the hierarchy layout of all nodes using the sugiyama algorithm =======
  hgap <- round(1 / gen_max_size, 8)
  gen_num <- max(real_node$gen, na.rm = TRUE)
  max_layer <- max(ped_igraph$node$layer, na.rm = TRUE)
  # The mean space between two near nodes in each generation are obtained by the
  # normalized axis x from the layout_with_sugiyama
  layers <-  vector(mode = "list", length = max_layer)
  k <- 1
  for (i in max_layer:1) {
    layers[[k]] <- ped_igraph$node[layer == i, id]
    k <- k + 1
  }
  g <- graph_from_data_frame(ped_igraph$edge, directed = TRUE, ped_igraph$node)
  l <- layout_with_sugiyama(g,
        layers = apply(sapply(layers, function(x)
          V(g)$name %in% x)
          , 1
          , which),
        hgap = hgap,
        maxiter = 120,
        attributes = "all")$layout
  l <- norm_coords(l,
        xmin = 0,
        xmax = 1,
        ymin = 0,
        ymax = 1)
  ped_igraph$node <- cbind(ped_igraph$node, x = l[, 1], y = l[, 2])

  #=== Adjusting space between two nodes (individuals) in x axis for each generation ==
  real_node <- ped_igraph$node[nodetype %in% c("real", "compact")]
  if (gen_max_size >= 2) {
    x_stats_gen <-
      real_node[, .(.N, range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
                by =gen]
    meanspace = NULL # due to NSE notes in R CMD check
    x_stats_gen[range > 0 & N > 1, ":="(meanspace = range / (N - 1))]
    l_x_range <- range(l[, 1], na.rm = TRUE)
    l_x_distance <- diff(l_x_range)
    max_gen_mean_space <- max(x_stats_gen$meanspace, na.rm = TRUE)
    for (i in 1:gen_num) {
      v_rank <- rank(real_node[gen == i, x], na.last = TRUE)
      node_num <- length(v_rank)
      if (node_num >= 2) {
        x_distance_1 <- diff(range(real_node[gen == i, x], na.rm = TRUE))
        x_distance_2 <- (node_num - 1) * max_gen_mean_space
        if (x_distance_2 > l_x_distance) {
          x_distance_2 <- l_x_distance
        }
        mean_space <- round(x_distance_2 / (node_num - 1), 8)
        a <- x_distance_2 - x_distance_1
        x_min <- min(real_node[gen == i, x], na.rm = TRUE)
        b <- x_min - a
        if (b > 0) {
          x_min <- b
        } else {
          x_min <- 0
        }
        x_new <- x_min + seq(from = 0, to = node_num - 1) * mean_space
        real_node[gen == i, x := x_new[v_rank]]
      }
    }
    ped_igraph$node[nodetype %in% c("real", "compact")] <-
      real_node
  }

  #=== Matching a virtual node's x pos to the samllest position of the full-sib =======
  # A virutal node is a tie between two parents and their progenies
  virtual_node <- ped_igraph$node[nodetype %in% c("virtual")]
  for (i in 2:gen_num) {
    real_family_min_x <- real_node[gen == i, .(minx = min(x,na.rm=TRUE)),by=c("familylabel")]
    virtual_family_label <- virtual_node[gen == i, familylabel]
    min_x <- real_family_min_x[match(virtual_family_label, familylabel),minx]
    virtual_node[gen == i, x := min_x]
  }
  ped_igraph$node[nodetype %in% c("virtual")] <- virtual_node
  # l[,1] <- ped_igraph$node[match(V(g)$name,as.character(id)),x]
  l[, 1] <- ped_igraph$node[, x]


  #=== Rescale canvas' size, node's size and edge's size ==============================
  # calculate the width of each node: inch
  node_width <- label_max_width
  if (gen_max_size >= 2) {
    x_stats_gen <-
      real_node[, .(.N, range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
                by = gen]
    x_stats_gen[range > 0 & N > 1, ":="(meanspace = range / (N - 1))]
    min_node_space <-
      min(x_stats_gen[meanspace > 0, meanspace], na.rm = TRUE)
    f <- 1
    if (max(x_stats_gen$N) <= 16) {
      # increase large space between two nodes when the node number is small
      f <- 3 * round(node_width / min_node_space, 8)
    }
    else {
      # keep samll space between two nodes when the node number are big
      f <- round(node_width / min_node_space, 8)
    }
  } else {
    f <- 1
  }
  x_f <- f * (l[, 1])
  # Setting the width of the canvas
  # Adding extra 6 node width to the canvas's width to decerase node size
  # when only have one node in one layer
  # because node size is equal to the percentage of node width to the canvas width
  width <- max(x_f, na.rm = TRUE) - min(x_f, na.rm = TRUE) + 6 * node_width

  # The maximum width or height for a pdf file is 200 inch
  pdf_maximum_width <- pdf_maximum_height <- pdf_max_width
  if (width > pdf_maximum_width) {
    width <- pdf_maximum_width
  }

  height <- width * 0.618

  # inch
  gen_height <- 0.618
  if (height < gen_num * (node_width) + 3 * node_width) {
    height <- gen_num * (node_width) + 3 * node_width
  }
  if (height > pdf_maximum_height) {
    height <- pdf_maximum_height
  }

  # vertes_size is a percentage of the width of node to graph
  node_size <- round(node_width * 100 / width, 8)
  edge_size <- node_size * 0.001
  edge_arrow_size <- node_size * 0.002
  edge_arrow_width <- node_size * 0.006
  V(g)$size[V(g)$nodetype %in% c("real", "compact")] = node_size
  if (is.null(cex) & (best_cex > 0)) {
    V(g)$label.cex[V(g)$nodetype %in% c("real", "compact")] = best_cex
  } else {
    V(g)$label.cex[V(g)$nodetype %in% c("real", "compact")] = cex
  }
  E(g)$size = edge_size
  E(g)$arrow.size = edge_arrow_size
  E(g)$arrow.width = edge_arrow_width
  if (show) {
    plot.igraph(
      g,
      rescale = FALSE,
      xlim = c(0 - node_size / 100, 1 + node_size / 100),
      ylim = c(1 + node_size / 100, 0 - node_size / 100),
      layout = l,
      asp = 0
    )
  }
  pdf(file = file,
      width = width,
      height = height)
  plot.igraph(
    g,
    rescale = FALSE,
    xlim = c(0 - node_size / 100, 1 + node_size / 100),
    ylim = c(1 + node_size / 100, 0 - node_size / 100),
    layout = l,
    asp = 0
  )
  dev.off()
  message(paste("The HD graph of pedigree is saved in the ",
          getwd(),
          "/",file," file",sep=""))
  if (is.null(cex)) {
    message(paste("The cex for individual label is ", best_cex, ".", sep = ""))
  } else {
    message(paste("The cex for individual label is ", cex, ".", sep = ""))
  }
  message(
    paste("Please decease or increase the value of cex  paremter and rerun visped() function ",
          "when the label's width is longer or shorter than that of the circle in the ",file," file",
          sep=""))

}

ped2igraph <- function(ped,compact=TRUE) {
  # IndNum, SireNum, and DamNum columns are used as ids to node and edges
  ped_col_names <- names(ped)
  if (!all(c("IndNum", "SireNum", "DamNum") %in% ped_col_names)) {
    stop("The pedigree need to be firstly trimmed by the tidyped() function!")
  }
  ped_new <- copy(ped)
  # There is the Cand column in the pedigree if it is traced by the tidyped function
  if (c("Cand") %in% ped_col_names) {
    ped_node <-
      ped_new[, .(
        id = IndNum,
        label = Ind,
        sirenum = SireNum,
        damnum = DamNum,
        sirelabel = Sire,
        damlabel = Dam,
        cand = Cand,
        sex = Sex,
        gen = Gen
      )]
  } else {
    ped_node <-
      ped_new[, .(
        id = IndNum,
        label = Ind,
        sirenum = SireNum,
        damnum = DamNum,
        sirelabel = Sire,
        damlabel = Dam,
        sex = Sex,
        gen = Gen
      )]
  }

  max_id <- max(ped_node$id,na.rm = TRUE)

  # Adding two new columns family label (column name: familylabel) and it's numeric id
  # (column name: familynum) in the ped_node
  familylabel = NULL # due to NSE notes in R CMD check
  ped_node[!(is.na(sirelabel) &
               is.na(damlabel)),
           familylabel := paste(sirelabel, damlabel, sep = "")]
  family_label <- unique(ped_node$familylabel)
  family_label <- family_label[!is.na(family_label)]
  family_num <-
    setDT(list(
      familynum = seq(
        from = max(ped_node$id,na.rm=TRUE) + 1,
        to = max(ped_node$id,na.rm=TRUE) + length(family_label)
      ),
      familylabel = family_label
    ))
  ped_node <-
    merge(ped_node,
          family_num,
          by = c("familylabel"),
          all.x = TRUE)
  ped_node[is.na(familynum), familynum := 0]

  # There will be three node types in the ped_note, including real, compact, and virtual.
  # Real nodes are all individuals in the pedigree.
  # Compact nodes are full-sib individuals with parents, but without progeny,
  # they exist only when the "compact" paramete is TRUE
  nodetype = NULL # due to NSE notes in R CMD check
  ped_node[,nodetype:="real"]

  #=== Compact the pedigree============================================================
  # Full-sib individuals with parents but without progeny will be deleted from ped_note.
  # count individuals by family and sex as a number of node  replace full-sib individuals
  if (compact) {
    # Finding the individuals with parents, but without progeny
    sire_dam_label <- unique(c(ped_node$sirelabel,ped_node$damlabel))
    sire_dam_label <- sire_dam_label[!is.na(sire_dam_label)]
    ped_node_1 <- ped_node[!(label %in% sire_dam_label)]

    # Moreover, finding full-sib individuals
    ped_node_1[,N:=.N,by=.(familylabel,sex)]
    setnames(ped_node_1,c("N"),c("familysize"))
    if (max(ped_node_1$familysize,na.rm=TRUE)>=2) {
      # The full-sib individuals in a family will be compacted if the family size >= 2
      fullsib_id_DT <- ped_node_1[familysize >=2]
      fullsib_ids <- fullsib_id_DT$id
      familylabelsex = NULL # due to NSE notes in R CMD check
      fullsib_id_DT[,familylabelsex:=paste(familylabel,sex,sep="")]
      # Generating a compact family dataset, only including maximum three individuals for
      # each family: male, female and unknown sex individuals
      fullsib_family_label_sex <- unique(fullsib_id_DT$familylabelsex)
      compact_family <- fullsib_id_DT[match(fullsib_family_label_sex,familylabelsex)]
      # The compact families' id are the number of individuals by family and sex.
      compact_family[,":="(label=familysize,nodetype="compact")]
      # Deleting full-sib individuals from families with 2 and more full-sib individuals
      ped_node <- ped_node[!(id %in% fullsib_ids)]
      ped_node <- rbind(ped_node,compact_family,fill=TRUE)
    }
  }

  #=== Add virtual nodes between parents and progrenies================================
  # Add id to familynum and familynum to sirenum and damnum as new virtual edges
  ped_edge <-
    rbind(ped_node[, .(from = id, to = familynum)],
          ped_node[, .(from = familynum, to = sirenum)],
          ped_node[, .(from = familynum, to = damnum)])
  ped_edge <- ped_edge[!(to == 0)]
  # Delete duplicated edges from familynum to parents
  ped_edge <- unique(ped_edge)
  ped_edge <- ped_edge[order(from, to)]
  size = arrow.size = arrow.width = color = curved = NULL # due to NSE notes in R CMD check
  # grey
  #ped_edge[,":="(size=1,arrow.size=1,arrow.width=1,color="#9d96ad",curved=0.15)]
  #ped_edge[,":="(size=1,arrow.size=1,arrow.width=1,color="#a69f89",curved=0.15)]
  ped_edge[,":="(size=1,arrow.size=1,arrow.width=1,color="#afa8be",curved=0.10)]
  # Add familynum as new virtual nodes
  ped_node <-
    rbind(ped_node, unique(ped_node[familynum > 0, .(
      id = familynum,
      familylabel,
      label = familylabel,
      sirenum,
      damnum,
      sirelabel,
      damlabel,
      gen,
      familynum
    )]), fill = TRUE)
  ped_node[is.na(nodetype),nodetype:="virtual"]
  layer = NULL # due to NSE notes in R CMD check
  ped_node[nodetype %in% c("real","compact"),layer:=2*gen-1]
  ped_node[nodetype %in% c("virtual"),layer:=2*(gen-1)]


  #=== Set default shape, size and color for male and female===========================
  # Setting the default attributes of nodes
  # Notes: size = 15 means the width of a circle node account for 15% of the whole width
  # of the graph
  #ped_node[, ":="(shape = "circle", frame.color="#8495e8", color="#9daaea",size = 15)]
  #ped_node[, ":="(shape = "circle", frame.color="black", color="#aaa16c",size = 15)]
  shape = frame.color = color = size = label.color = NULL
  ped_node[, ":="(shape = "circle", frame.color="#7fae59", color="#9cb383",size = 15, label.color="#0d0312")]
  #ped_node[, ":="(shape = "circle", frame.color=NA, color="#9cb383",size = 15)]
  ped_node[nodetype %in% c("compact"), ":="(shape="square")]
  # Setting virtual size of nodes to 0.0001
  ped_node[id > max_id,":="(shape="none",label="",size=0)]
  # Setting male and female's color
  ped_node[sex %in% c("male"), ":="(frame.color="#0e8dbb", color = "#119ecc")]
  ped_node[sex %in% c("female"), ":="(frame.color="#e6a11f", color = "#f4b131")]

  # The edge color is same with the color of the it's "to" node.
  min_familynum <- min(family_num$familynum)
  ped_edge <- merge(ped_edge, ped_node[,.(id,tonodecolor=color)],by.x="to", by.y="id",all.x=TRUE)
  ped_edge[from >= min_familynum,":="(color=tonodecolor)]
  ped_edge[from < min_familynum,":="(curved=0)]


  # Sorting the from and to columns as the first two columns in the ped_edge
  old_names <- colnames(ped_edge)
  new_names <- c(c("from","to"),old_names[!(old_names %in% c("from","to"))])
  ped_edge <- ped_edge[, ..new_names]
  ped_edge <- ped_edge[order(from,to)]


  # Sorting the id column as the first column in the ped_node
  old_names <- colnames(ped_node)
  new_names <- c("id",old_names[!(old_names %in% c("id"))])
  ped_node <- ped_node[, ..new_names]
  ped_node <- ped_node[order(layer,id)]

  return(list(node = ped_node, edge = ped_edge))

}



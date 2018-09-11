#' @import data.table
#' @import igraph

visped <- function(ped,
                   compact = FALSE, cex=NULL) {
  ped_new <- copy(ped)
  # Convertting pedigree to node and edge two data.table
  ped_igraph <- ped2igraph(ped_new,compact)
  #ped_igraph$node <- sortnode(ped_igraph$node)
  real_node <- ped_igraph$node[nodetype %in% c("real","compact")]
  gen_node_num <- real_node[nodetype %in% c("real","compact")][,.N,by=gen]
  max_size_gen <-  max(gen_node_num$N,na.rm=TRUE)


  # The maximum width of a pdf file is 200 inch
  pdf_max_width = 200
  cexs <-  seq(from=0.1, to =1, by=0.05)
  best_cex <- 0
  max_strwidth_label <- real_node[which.max(strwidth(real_node$label,cex = 1, units="inches")),label]
  for (i in length(cexs):1) {
    # Obtaining the maximum width of a node's label: inch
    label_max_width <- max(strwidth(max_strwidth_label,cex=cexs[i], units="inches"),na.rm = TRUE)
    if (max_size_gen <= 16 & label_max_width < 0.8) {label_max_width = 0.8}
    if ((label_max_width * max_size_gen) < pdf_max_width ) {
      best_cex <- cexs[i]*0.7
      break
    }
  }
  if (best_cex == 0) {
    stop("Too many individuals (>=", max_size_gen,
         ") in one genrations!!! please visualize the pedigree by setting the compact parameter to TRUE.")
  }

  hgap <- round(1/max_size_gen,8)
  max_gen_num <- max(real_node$gen,na.rm=TRUE)
  ped_igraph$node[nodetype %in% c("real","compact"),layer:=2*gen]
  ped_igraph$node[nodetype %in% c("virtual"),layer:=2*gen-1]
  max_layer <- max(ped_igraph$node$layer, na.rm=TRUE)
  # The mean space between two near nodes in each generation are obtained by the
  # normalized axis x from the layout_with_sugiyama
  layers <-  vector(mode = "list",length = max_layer)
  k <- 1
  for (i in max_layer:1) {
    layers[[k]] <- ped_igraph$node[layer==i,id]
    k <- k + 1
  }
  g <- graph_from_data_frame(ped_igraph$edge,directed = TRUE,ped_igraph$node)
  l <- layout_with_sugiyama(g, layers=apply(sapply(layers,function(x) V(g)$name %in% x)
                                            ,1
                                            ,which),
                            hgap=hgap)$layout
  l <- norm_coords(l,xmin=0 ,xmax = 1,ymin=0,ymax = 1)
  ped_igraph$node <- cbind(ped_igraph$node,x=l[,1],y=l[,2])
  real_node <- ped_igraph$node[nodetype %in% c("real","compact")]
  x_stats_gen <- real_node[,.(.N,range= max(x,na.rm=TRUE)-min(x,na.rm=TRUE)),by=gen]
  x_stats_gen[range>0 & N>1,":="(meanspace=range/(N-1))]

  # adjust the space between two nodes to mean space in each generation
  real_node <- real_node[order(gen,id)]
  for (i in 1:max_gen_num) {
    x1 <- sort(real_node[gen == i, x])
    node_num <- length(x1)
    if (node_num > 1) {
      x2 <- x1[1] + seq(from=0,to=node_num-1)*x_stats_gen[gen==i,meanspace]
      real_node[gen == i,x:=x2]
    }
  }
  ped_igraph$node[nodetype %in% c("real","compact")] <-  real_node

  # adjust the virtual nodes's x axis pos to that of the first small real node in the same family
  virtual_node <- ped_igraph$node[nodetype %in% c("virtual")]
  virtual_node <- virtual_node[order(gen,id)]
  for (i in 2:max_gen_num) {
    label_1 <- real_node[gen==i,familylabel]
    label_2 <- virtual_node[gen==i,familylabel]

    x1 <- real_node[gen==i,x]
    x2 <- x1[match(label_2,label_1)]
    virtual_node[gen==i,x:= x2]
  }
  ped_igraph$node[nodetype %in% c("virtual")] <- virtual_node
  l[,1] <- ped_igraph$node[match(V(g)$name,as.character(id)),x]

  #=== Rescale canvas' size, node's size and edge's size ==============================
  # calculate the width of each node: inch
  node_width <- label_max_width
  min_node_space <- min(x_stats_gen[meanspace>0,meanspace],na.rm=TRUE)
  f <- 1
  if (!is.infinite(min_node_space)) {
    if (max(x_stats_gen$N) <= 16) {
      # increase large space between two nodes when the node number is small
      f <- 3*round(node_width / min_node_space, 8)
    }
    else {
      # keep samll space between two nodes when the node number are big
      f <- round(node_width / min_node_space, 8)
    }
  }
  x_f <- f*(l[,1])
  # Setting the width of the canvas
  # Adding extra 6 node width to the canvas's width to decerase node size when only have one node in one layer
  # because node size is equal to the percentage of node width to the canvas width
  width <- max(x_f,na.rm=TRUE) -min(x_f,na.rm = TRUE) + 6*node_width

  # The maximum width or height for a pdf file is 200 inch
  pdf_maximum_width <- pdf_maximum_height <- pdf_max_width
  if (width > pdf_maximum_width) {
    width <- pdf_maximum_width
  }

  height <- width * 0.618

  # inch
  gen_height <- 0.618
  if (height < max_gen_num * (gen_height+node_width)) {
    height <- max_gen_num * (gen_height+node_width)
  }
  if (height > pdf_maximum_height) {
    height <- pdf_maximum_height
  }

  # vertes_size is a percentage of the width of the graph
  vertex_size <- round(node_width*100/width,8)
  edge_size <- vertex_size*0.001
  edge_arrow_size <- vertex_size*0.002
  edge_arrow_width <- vertex_size*0.006
  V(g)$size[V(g)$nodetype %in% c("real","compact")] = vertex_size
  if (is.null(cex) & (best_cex >0) ) {
    V(g)$label.cex[V(g)$nodetype %in% c("real","compact")] = best_cex
  } else {
    V(g)$label.cex[V(g)$nodetype %in% c("real","compact")] = cex
  }
  E(g)$size = edge_size
  E(g)$arrow.size = edge_arrow_size
  E(g)$arrow.width = edge_arrow_width
  plot.igraph(g,rescale=FALSE,xlim = c(0-vertex_size/100,1+vertex_size/100),
              ylim = c(1+vertex_size/100,0-vertex_size/100),layout=l,asp=0)
  pdf(file="output.pdf",width = width, height = height)
  plot.igraph(g,rescale=FALSE,xlim = c(0-vertex_size/100,1+vertex_size/100),
               ylim = c(1+vertex_size/100,0-vertex_size/100),layout=l,asp=0)
  dev.off()
  message("The HD graph of pedigree is saved in the ", getwd(), "/output.pdf file")
  if (is.null(cex)) {
    message(paste("The cex for individual label is ",best_cex,".",sep=""))
  } else {
    message(paste("The cex for individual label is ",cex,".",sep=""))
  }
  message("Please decease or increase the value of cex  paremter and rerun visped() function when the label's width is longer or shorter than that of the circle in the output.pdf file")

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
    N_by_family <- ped_node_1[,.N,by=.(familylabel,sex)]
    if (max(N_by_family$N,na.rm=TRUE)>=2) {
      ped_node_1 <- merge(ped_node_1, N_by_family,by=c("familylabel","sex"),all.x=TRUE)
      setnames(ped_node_1,c("N"),c("familysize"))
      # The full-sib individuals in a family will be compacted if the family size >= 2
      fullsib_id_DT <- ped_node_1[familysize >=2]
      fullsib_ids <- fullsib_id_DT$id
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

  #=== Set default shape, size and color for male and female===========================
  # Setting the default attributes of nodes
  # Notes: size = 5 means the width of a circle node account for 5% of the whole width
  # of the graph
  #ped_node[, ":="(shape = "circle", frame.color="#8495e8", color="#9daaea",size = 15)]
  #ped_node[, ":="(shape = "circle", frame.color="black", color="#aaa16c",size = 15)]
  ped_node[, ":="(shape = "circle", frame.color="#7fae59", color="#9cb383",size = 15, label.color="#0d0312")]
  #ped_node[, ":="(shape = "circle", frame.color=NA, color="#9cb383",size = 15)]
  ped_node[nodetype %in% c("compact"), ":="(shape="square")]
  # Setting virtual size of nodes to 0.0001
  ped_node[id > max_id,":="(shape="none",label="",size=0.0001)]
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
  ped_edge <- ped_edge[,new_names,with=FALSE]
  ped_edge <- ped_edge[order(from,to)]


  # Sorting the id column as the first column in the ped_node
  old_names <- colnames(ped_node)
  new_names <- c("id",old_names[!(old_names %in% c("id"))])
  ped_node <- ped_node[,new_names,with=FALSE]
  ped_node <- ped_node[order(id)]

  return(list(node = ped_node, edge = ped_edge))

}

sortnode <- function(node) {
  if (!is.data.table(node)) {
    stop("nodes need to be data.table type.")
  }
  each_family_size <- node[, .N, by=.(sirelabel,damlabel)]
  setnames(each_family_size,c("N"),c("familysize"))
  new_node <- merge(node,each_family_size,by=c("sirelabel","damlabel"),all.x=TRUE)
  new_node[is.na(sirelabel) & is.na(damlabel), familysize:=0]
  old_names <- colnames(new_node)
  new_names <- c("id",old_names[!(old_names %in% c("id"))])
  new_node <- new_node[,new_names,with=FALSE]
  return(new_node[order(gen,familysize,sirelabel,damlabel)])
}


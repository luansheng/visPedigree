#' @import visNetwork
#' @import data.table
#' @import igraph

visped <- function(ped,
                   fullsib = FALSE,
                   layout = NULL) {
  if (fullsib) {
    ped_new <- cutped(ped)
  } else {ped_new <- copy(ped)}

    nodes <- ped2igraphnode(ped_new)
    nodes$node <- sortnode(nodes$node)

    # The mean space between two near nodes in each generation are obtained by the
    # normalized axis x from the layout_as_tree
    g <- graph_from_data_frame(nodes$edge,directed = TRUE,nodes$node)
    l <- layout_as_tree(g,flip.y = FALSE)
    l <- norm_coords(l,xmin=0 ,xmax = 1,ymin=0,ymax = 1)
    nodes$node <- cbind(nodes$node,x=l[,1],y=l[,2])
    real_node <- nodes$node[nodetype %in% c("real")]
    x_stats_gen <- real_node[,.(.N,range= max(x,na.rm=TRUE)-min(x,na.rm=TRUE)),by=gen]
    x_stats_gen[,":="(meanspace=range/(N-1))]

    # adjust the space between two nodes  to mean space in each generation
    real_node <- real_node[order(gen,id)]
    gen <- max(real_node$gen,na.rm=TRUE)
    for (i in 1:gen) {
      x1 <- sort(real_node[gen == i, x])
      node_num <- length(x1)
      if (node_num > 1) {
        x2 <- x1[1] + seq(from=0,to=node_num-1)*x_stats_gen[gen==i,meanspace]
        real_node[gen == i,x:=x2]
      }
    }
    nodes$node[nodetype %in% c("real")] <-  real_node

    # adjust the virtual nodes's axis x pos to the first small real node in the same family
    virtual_node <- nodes$node[nodetype %in% c("virtual")]
    virtual_node <- virtual_node[order(gen,id)]
    for (i in 2:gen) {
      label_1 <- real_node[gen==i,familylabel]
      label_2 <- virtual_node[gen==i,familylabel]

      x1 <- real_node[gen==i,x]
      x2 <- x1[match(label_2,label_1)]
      virtual_node[gen==i,x:= x2]
    }
    nodes$node[nodetype %in% c("virtual")] <- virtual_node
    l[,1] <- nodes$node[match(V(g)$name,as.character(id)),x]

    #=== Rescale canvas' size, node's size and edge's size ==============================
    # The number of pixel in each inch
    pdf_dpi <- 300
    # Obtain the maximum width of a node's label: inch
    label_max_width <- max(strwidth(real_node$label, units="inches"),na.rm = TRUE)
    # calculate the width of each node: px
    node_width <- label_max_width * pdf_dpi
    min_node_space <- min(x_stats_gen[meanspace>0,meanspace],na.rm=TRUE)
    if (max(x_stats_gen$N) <= 16) {
      # There are large spaces between two nodes when the node number is small
      f <- round(2*node_width / min_node_space, 8)
    }
    else {
      # There are small spaces between two nodes when the node number are big
      f <- round(node_width / min_node_space, 8)
    }
    x_f <- f*(l[,1])
    # get the width of the canvas
    width <- max(x_f,na.rm=TRUE) -min(x_f,na.rm = TRUE) + node_width

    # The maximum width or height for a pdf file is 200 inch
    pdf_maximum_width <- pdf_maximum_height <- 200 * pdf_dpi
    if (width > pdf_maximum_width) {
      width <- pdf_maximum_width
    }

    height <- width * 0.618
    gen_height <- 300
    if (height < gen * (gen_height+node_width)) {
      height <- gen * (gen_height+node_width)
    }
    if (height > pdf_maximum_height) {
      height <- pdf_maximum_height
    }

    # vertes_size is a percentage compared to the width of the graph
    vertex_size <- round(node_width*100/width,8)
    edge_size <- vertex_size*0.001
    edge_arrow_size <- vertex_size*0.02
    edge_arrow_width <- vertex_size*0.06
    V(g)$size[V(g)$nodetype %in% c("real")] = vertex_size
    E(g)$size = edge_size
    E(g)$arrow.size = edge_arrow_size
    E(g)$arrow.width = edge_arrow_width
    pdf(file="output.pdf",width = width / pdf_dpi, height = height / pdf_dpi)
    plot.igraph(g,rescale=FALSE,xlim = c(0-vertex_size/100,1+vertex_size/100),
                ylim = c(0-vertex_size/100,1+vertex_size/100),layout=l,asp=0)
    dev.off()

# todo


    #
    #   if (sum(ind_num_gen$N,na.rm = TRUE) <= 10) {
    #     gen_max <- max(ind_num_gen$gen)
    #     level_space = round(height / gen_max,8)*0.99
    #     visNetwork(
    #       nodes$node,
    #       nodes$edge,
    #       width = paste(width, "px", sep = ""),
    #       height = paste(height, "px", sep = "")
    #     ) %>%
    #       visNodes(physics = TRUE,scaling = list(label=list(enabled=TRUE))) %>%
    #       visEdges(arrows = "to", color=list(inherit="to"),width = 0.1,smooth = FALSE) %>%
    #       visHierarchicalLayout(direction = "DU",sortMethod = "directed",
    #                             levelSeparation = level_space,
    #                             edgeMinimization = TRUE)  %>%
    #       visPhysics(stabilization = FALSE) %>%
    #       visInteraction(dragNodes = TRUE,
    #                      dragView = TRUE,
    #                      zoomView = TRUE)
    #   } else {
}

ped2igraphnode <- function(ped) {
  ped_col_names <- names(ped)
  if (!all(c("IndNum", "SireNum", "DamNum") %in% ped_col_names)) {
    stop("The pedigree need to be firstly trimmed by the tnameyped() function!")
  }
  ped_new <- copy(ped)
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

  #=== add virtual nodes between parents and progrenies================================
  # add a new column familylabel and it's numeraic id in the ped_node
  ped_node[!(is.na(sirelabel) &
               is.na(damlabel)),
           familylabel := paste(sirelabel, damlabel, sep = "")]
  family_label <- unique(ped_node$familylabel)
  family_label <- family_label[!is.na(family_label)]
  family_DT <-
    setDT(list(
      familynum = seq(
        from = max(ped_node$id,na.rm=TRUE) + 1,
        to = max(ped_node$id,na.rm=TRUE) + length(family_label)
      ),
      familylabel = family_label
    ))
  ped_node <-
    merge(ped_node,
          family_DT,
          by = c("familylabel"),
          all.x = TRUE)
  ped_node[is.na(familynum), familynum := 0]

  # add id to familynum and familynum to sirenum and damnum as new virtual edges
  ped_edge <-
    rbind(ped_node[, .(from = id, to = familynum)],
          ped_node[, .(from = familynum, to = sirenum)],
          ped_node[, .(from = familynum, to = damnum)])
  ped_edge <- ped_edge[!(to == 0)]
  # delete duplicated edges from familynum to parents
  ped_edge <- unique(ped_edge)
  ped_edge <- ped_edge[order(from, to)]
  ped_edge[,":="(size=1,arrow.size=1,arrow.wnameth=1,color="darkgrey",curved=0.15)]

  # add familynum as new virtual nodes
  ped_node <- cbind(ped_node,nodetype=rep("real",nrow(ped_node)))
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

  #=== set default shape, size and color for male and female===========================
  # the shape size 10 may be percentage in igraph
  ped_node[, ":="(shape = "none", frame.color=NA, size = 5,label.degree=0)]
  # set virtual size of nodes to 0.0001
  ped_node[id > nrow(ped_new),":="(shape="none",label="",size=0.0001)]
  # set male and female's shape
  ped_node[sex %in% c("male"), shape := "square"]
  ped_node[sex %in% c("female"), shape := "circle"]
  # set male and female's color
  ped_node[sex %in% c("male"), color := "#119ecc"]
  ped_node[sex %in% c("female"), color := "#f4b131"]

  # edge color is same with that of the to node.
  min_familynum <- min(family_DT$familynum)
  ped_edge <- merge(ped_edge, ped_node[,.(id,tonodecolor=color)],by.x="to", by.y="id",all.x=TRUE)
  ped_edge[from >= min_familynum,":="(color=tonodecolor)]
  ped_edge[from < min_familynum,":="(curved=0)]

  # make sure from and to  are the first two columns in the ped_edge
  old_names <- colnames(ped_edge)
  new_names <- c(c("from","to"),old_names[!(old_names %in% c("from","to"))])
  ped_edge <- ped_edge[,new_names,with=FALSE]
  ped_edge <- ped_edge[order(from,to)]


  # make sure id is the first column in the ped_node
  old_names <- colnames(ped_node)
  new_names <- c("id",old_names[!(old_names %in% c("id"))])
  ped_node <- ped_node[,new_names,with=FALSE]
  ped_node <- ped_node[order(id)]

  return(list(node = ped_node, edge = ped_edge))
}

ped2visnetworknode <- function(ped) {
  ped_col_names <- names(ped)
  if (!all(c("IndNum", "SireNum", "DamNum") %in% ped_col_names)) {
    stop("The pedigree need to be firstly trimmed by the tidyped() function!")
  }
  ped_new <- copy(ped)
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

  ped_node[!(is.na(sirelabel) &
               is.na(damlabel)),
           familylabel := paste(sirelabel, damlabel, sep = "")]
  family_label <- unique(ped_node$familylabel)
  family_label <- family_label[!is.na(family_label)]
  family_DT <-
    setDT(list(
      familynum = seq(
        from = max(ped_node$id,na.rm=TRUE) + 1,
        to = max(ped_node$id,na.rm=TRUE) + length(family_label)
      ),
      familylabel = family_label
    ))
  ped_node <-
    merge(ped_node,
          family_DT,
          by = c("familylabel"),
          all.x = TRUE)
  ped_node[is.na(familynum), familynum := 0]

  ped_edge <-
    rbind(ped_node[, .(from = id, to = familynum)],
          ped_node[, .(from = familynum, to = sirenum)],
          ped_node[, .(from =familynum, to = damnum)])
  ped_edge <- ped_edge[!(to == 0)]
  ped_edge <- unique(ped_edge)
  ped_edge <- ped_edge[order(from, to)]

  ped_node <-
    rbind(ped_node, unique(ped_node[familynum > 0, .(
      familylabel,
      id = familynum,
      label = familylabel,
      sirenum,
      damnum,
      sirelabel,
      damlabel,
      gen,
      familynum
    )]), fill = TRUE)
  ped_node[, ":="(shape = "ellipse",size=15)]
  ped_node[id>nrow(ped_new),":="(shape="square",label="",size=0.1)]
  ped_node[sex %in% c("male"), shape := "square"]
  ped_node[sex %in% c("female"), shape := "dot"]

  ped_node[sex %in% c("male"), color := "#119ecc"]
  ped_node[sex %in% c("female"), color := "#f4b131"]

  if (c("Cand") %in% ped_col_names) {
    ped_node[cand == TRUE, color := "red"]
  }
  old_names <- colnames(ped_node)
  new_names <- c("id",old_names[!(old_names %in% c("id"))])
  ped_node <- ped_node[,new_names,with=FALSE]
  ped_node <- ped_node[order(id)]

  return(list(node = ped_node, edge = ped_edge))
}

cutped <- function(ped) {
  ped_col_names <- colnames(ped)
  if (!all(c("IndNum","SireNum","DamNum") %in% ped_col_names)) {
    stop("The pedigree need to be firstly trimmed by the tidyped(addnum = TRUE) function!")
  }
  ped_new <- copy(ped)
  ped_new[,":="(IndNum=NULL,SireNum=NULL,DamNum=NULL)]
  ped_new_col_names <- colnames(ped_new)

  # delete candidate individuals
  if (c("Cand") %in% ped_col_names) {
    ped_new <- ped_new[!(cand == TRUE)]
  }

  #去除双亲未知的个体
  ped_new <- ped_new[!(is.na(Sire) | is.na(Dam))]
  sires_dams <- unique(c(ped_new$Sire,ped_new$Dam))
  sires_dams <- sires_dams[!is.na(sires_dams)]

  #找那些没有后代，只作为后代出现的个体，寻找全同胞个体#
  ped_new <- ped_new[!(Ind %in% sires_dams)]

  #找出全同胞个体
  ped_new[,FamilyID:=paste(Sire,Dam,sep = "")]
  N_by_family <- ped_new[,.N,by=.(FamilyID)]
  if (max(N_by_family$N,na.rm=TRUE)>=2) {
    ped_new <- merge(ped_new,N_by_family,by="FamilyID",all.x=TRUE)
    setnames(ped_new,c("N"),c("FamilySize"))
    fullsib_ind_DT <- ped_new[FamilySize >=2]
    fullsib_ind_name <- fullsib_ind_DT$Ind
    fullsib_family_ID <- unique(fullsib_ind_DT$FamilyID)
    fullsib_family_ID_DT <- setDT(list(FamilyID=fullsib_family_ID,SeqNum=1:length(fullsib_family_ID)))
    fullsib_ind_DT <- merge(fullsib_ind_DT,fullsib_family_ID_DT,by="FamilyID",all.x=TRUE)

    fullsib_ind_DT <- fullsib_ind_DT[, ":="(Ind=paste("F",SeqNum,"size",FamilySize,sep=""), Sex=NA)]
    fullsib_ind_DT <- unique(fullsib_ind_DT[,ped_new_col_names, with=FALSE])
    ped_cut <- rbind(ped[!(Ind %in% fullsib_ind_name),ped_new_col_names, with=FALSE],fullsib_ind_DT)
    ped_cut <- tidyped(ped_cut)
    return(ped_cut)
  } else {
    return(ped)
  }
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


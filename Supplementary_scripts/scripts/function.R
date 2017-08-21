clusterCells_Density_Peak <- clusterCells

#' #' Scale pseudotime to be in the range from 0 to 100 (it works both for situations involving only one state and complex states)
#' #' @param cds the CellDataSet upon which to perform this operation
#' #' @param ordering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' #' @return an updated CellDataSet object which an
#' #' @export
#' scale_pseudotime <- function(cds, verbose = F) {
#'   pd <- pData(cds)
#'   pd$Cell_name <- row.names(pd)
#'   range_df <- plyr::ddply(pd, .(State), function(x) {
#'     min_max <- range(x$Pseudotime)
#'     min_cell <- subset(x, Pseudotime %in% min_max[1])
#'     max_cell <- subset(x, Pseudotime %in% min_max[2])
#'     min_cell$fate_type <- 'Start'
#'     max_cell$fate_type <- 'End'
#'     rbind(min_cell, max_cell)
#'   }) #pseudotime range for each state
#'   
#'   #1. construct a tree of the selected cells
#'   adj_list <- data.frame(Source = subset(range_df, length(Parent) > 0 & !is.na(Parent))[, 'Parent'], Target = subset(range_df, length(Parent) > 0 & !is.na(Parent))[, 'Cell_name'])
#'   #convert to cell fate adjancency list:
#'   adj_list$Source <- pd[as.character(adj_list$Source), 'State']
#'   adj_list$Target <- pd[as.character(adj_list$Target), 'State']
#'   
#'   uniq_cell_list <- unique(c(as.character(adj_list$Source), as.character(adj_list$Target)))
#'   adj_mat <- matrix(rep(0, length(uniq_cell_list)^2), nrow = length(uniq_cell_list), ncol = length(uniq_cell_list), dimnames = list(uniq_cell_list, uniq_cell_list))
#'   adj_mat[as.matrix(adj_list)] <- 1
#'   net <- graph.adjacency(as.matrix(adj_mat), mode = 'directed', weighted = NULL, diag = FALSE)
#'   
#'   # plot(net, layout = layout.fruchterman.reingold,
#'   #        vertex.size = 25,
#'   #        vertex.color="red",
#'   #        vertex.frame.color= "white",
#'   #        vertex.label.color = "white",
#'   #        vertex.label.cex = .5,
#'   #        vertex.label.family = "sans",
#'   #        edge.width=2,
#'   #        edge.color="black")
#'   
#'   #dfs_net <- graph.dfs(net, root = 1, order.out = T) #DFS search for the cell fate tree
#'   #net_order_out <- as.vector(dfs_net$order.out)
#'   net_leaves <- which(degree(net, v = V(net), mode = "out")==0, useNames = T)
#'   
#'   pd$scale_pseudotime <- NA
#'   
#'   #2. scale the psedutime:
#'   for(i in net_leaves) {
#'     path_vertex <- as.vector(get.all.shortest.paths(net, from = 1, to = i, mode="out")$res[[1]])
#'     pd_subset <- subset(pd, State %in% path_vertex & is.na(scale_pseudotime))
#'     
#'     #scale the pseudotime between the parent cell of the first cell to the last cell on the remaing branch
#'     # print(path_vertex)
#'     min_cell_name <- row.names(subset(pd_subset, Pseudotime == min(Pseudotime)))
#'     # print(min_cell_name)
#'     
#'     if(!is.na(pd[min_cell_name, 'Parent'])) {
#'       parent_min_cell <- as.character(pd[min_cell_name, 'Parent'])
#'       subset_min_pseudo <- pd[parent_min_cell, 'Pseudotime']
#'       scale_pseudotime_ini <- pd[parent_min_cell, 'scale_pseudotime']
#'       scaling_factor <- (100 - pd[parent_min_cell, 'scale_pseudotime']) / c(max(pd_subset$Pseudotime) - subset_min_pseudo)
#'     }
#'     else {
#'       subset_min_pseudo <- min(pd[, 'Pseudotime'])
#'       scale_pseudotime_ini <- 0
#'       scaling_factor <- 100 / c(max(pd_subset$Pseudotime) - min(pd_subset$Pseudotime))
#'     }
#'     
#'     pseudotime_scaled <- (pd_subset$Pseudotime - subset_min_pseudo) * scaling_factor + scale_pseudotime_ini
#'     
#'     if(verbose)
#'       message(i, '\t', range(pseudotime_scaled)[1],'\t', range(pseudotime_scaled)[2])
#'     
#'     pd[row.names(pd_subset), 'ori_pseudotime'] <- pd[row.names(pd_subset), 'Pseudotime']
#'     pd[row.names(pd_subset), 'Pseudotime'] <- pseudotime_scaled
#'   }	
#'   scale_pseudotime <- (pd_subset$Pseudotime - subset_min_pseudo) * scaling_factor + scale_pseudotime_ini
#'   message(i, '\t', range(scale_pseudotime)[1],'\t', range(scale_pseudotime)[2])
#'   pd[row.names(pd_subset), 'scale_pseudotime'] <- scale_pseudotime
#'   
#'   pData(cds) <- pd
#'   
#'   return(cds)
#' }
#' 
#' # Methods for PQ-tree based ordering
#' 
#' get_next_node_id <- function()
#' {
#'   next_node <<- next_node + 1
#'   return (next_node) 
#' }
#' 
#' #' Recursively builds and returns a PQ tree for the MST
#' #' @import igraph
#' pq_helper<-function(mst, use_weights=TRUE, root_node=NULL)
#' {
#'   new_subtree <- graph.empty()
#'   
#'   root_node_id <- paste("Q_", get_next_node_id(), sep="")
#'   
#'   new_subtree <- new_subtree + vertex(root_node_id, type="Q", color="black")
#'   
#'   if (is.null(root_node) == FALSE){
#'     sp <- get.all.shortest.paths(mst, from=V(mst)[root_node])
#'     #print(sp)
#'     sp_lengths <- sapply(sp$res, length)
#'     target_node_idx <- which(sp_lengths == max(sp_lengths))[1]
#'     #print(unlist(sp$res[target_node_idx]))
#'     diam <- V(mst)[unlist(sp$res[target_node_idx])]
#'     #print(diam)
#'   }else{
#'     if (use_weights){
#'       diam <- V(mst)[get.diameter(mst)]
#'     }else{
#'       diam <- V(mst)[get.diameter(mst, weights=NA)]
#'     }
#'   }
#'   
#'   
#'   #print (diam)
#'   
#'   V(new_subtree)[root_node_id]$diam_path_len = length(diam)
#'   
#'   diam_decisiveness <- igraph::degree(mst, v=diam) > 2
#'   ind_nodes <- diam_decisiveness[diam_decisiveness == TRUE]
#'   
#'   first_diam_path_node_idx <- head(as.vector(diam), n=1)
#'   last_diam_path_node_idx <- tail(as.vector(diam), n=1)
#'   if (sum(ind_nodes) == 0 || 
#'       (igraph::degree(mst, first_diam_path_node_idx) == 1 && 
#'        igraph::degree(mst, last_diam_path_node_idx) == 1))
#'   {
#'     ind_backbone <- diam
#'   }
#'   else 
#'   {
#'     last_bb_point <- names(tail(ind_nodes, n=1))[[1]]
#'     first_bb_point <- names(head(ind_nodes, n=1))[[1]]  
#'     #diam_path_vertex_names <- as.vector()
#'     #print (last_bb_point)
#'     #print (first_bb_point)
#'     diam_path_names <- V(mst)[as.vector(diam)]$name
#'     last_bb_point_idx <- which(diam_path_names == last_bb_point)[1]
#'     first_bb_point_idx <- which(diam_path_names == first_bb_point)[1]
#'     ind_backbone_idxs <- as.vector(diam)[first_bb_point_idx:last_bb_point_idx]
#'     #print (ind_backbone_idxs)
#'     ind_backbone <- V(mst)[ind_backbone_idxs]
#'     
#'     #ind_backbone <- diam[first_bb_point:last_bb_point]
#'   }
#'   
#'   
#'   
#'   mst_no_backbone <- mst - ind_backbone
#'   #print (V(mst_no_backbone)$name)
#'   
#'   for (backbone_n in ind_backbone)
#'   {
#'     #print (n)
#'     #backbone_n <- ind_backbone[[i]]
#'     
#'     if (igraph::degree(mst, v=backbone_n) > 2)
#'     {
#'       new_p_id <- paste("P_", get_next_node_id(), sep="")
#'       #print(new_p_id)
#'       new_subtree <- new_subtree + vertex(new_p_id, type="P", color="grey")
#'       new_subtree <- new_subtree + vertex(V(mst)[backbone_n]$name, type="leaf", color="white")
#'       new_subtree <- new_subtree + edge(new_p_id, V(mst)[backbone_n]$name)
#'       new_subtree <- new_subtree + edge(root_node_id, new_p_id)
#'       
#'       nb <- graph.neighborhood(mst, 1, nodes=backbone_n)[[1]]
#'       
#'       #print (E(nb))
#'       #print (V(nb))
#'       for (n_i in V(nb))
#'       {
#'         n <- V(nb)[n_i]$name			
#'         if (n %in% V(mst_no_backbone)$name)
#'         {	
#'           #print (n)
#'           
#'           sc <- subcomponent(mst_no_backbone, n)
#'           
#'           sg <- induced.subgraph(mst_no_backbone, sc, impl="copy_and_delete")
#'           
#'           
#'           if (ecount(sg) > 0)
#'           {
#'             #print (E(sg))	
#'             sub_pq <- pq_helper(sg, use_weights)
#'             
#'             
#'             # Works, but slow:
#'             for (v in V(sub_pq$subtree))
#'             {
#'               new_subtree <- new_subtree + vertex(V(sub_pq$subtree)[v]$name, type=V(sub_pq$subtree)[v]$type, color=V(sub_pq$subtree)[v]$color, diam_path_len=V(sub_pq$subtree)[v]$diam_path_len)
#'             }
#'             
#'             edge_list <- get.edgelist(sub_pq$subtree)
#'             for (i in 1:nrow(edge_list))
#'             {
#'               new_subtree <- new_subtree + edge(V(sub_pq$subtree)[edge_list[i, 1]]$name, V(sub_pq$subtree)[edge_list[i, 2]]$name)
#'             }   					
#'             #plot (new_subtree)
#'             
#'             new_subtree <- new_subtree + edge(new_p_id, V(sub_pq$subtree)[sub_pq$root]$name)  
#'           }
#'           else
#'           {
#'             new_subtree <- new_subtree + vertex(n, type="leaf", color="white")
#'             new_subtree <- new_subtree + edge(new_p_id, n)
#'           }
#'         }
#'         
#'       }
#'       #print ("##########################")
#'     }
#'     else
#'     {
#'       new_subtree <- new_subtree + vertex(V(mst)[backbone_n]$name, type="leaf", color="white")
#'       new_subtree <- new_subtree + edge(root_node_id, V(mst)[backbone_n]$name)
#'     }
#'   }
#'   # else
#'   # {
#'   #     for (backbone_n in diam)
#'   #     {
#'   #           new_subtree <- new_subtree + vertex(backbone_n, type="leaf")
#'   #           new_subtree <- new_subtree + edge(root_node_id, backbone_n)
#'   #     }
#'   # }
#'   
#'   return (list(root=root_node_id, subtree=new_subtree))
#' }
#' 
#' make_canonical <-function(pq_tree)
#' {
#'   nei <- NULL
#'   
#'   canonical_pq <- pq_tree
#'   
#'   V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out") == 2]$color="black"
#'   V(canonical_pq)[type == "P" &  igraph::degree(canonical_pq, mode="out")== 2]$type="Q"
#'   
#'   single_child_p <- V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")== 1]
#'   V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")==1]$color="blue"
#'   
#'   for (p_node in single_child_p)
#'   {
#'     child_of_p_node <- V(canonical_pq) [ suppressWarnings(nei(p_node, mode="out")) ]
#'     parent_of_p_node <- V(canonical_pq) [ suppressWarnings(nei(p_node, mode="in")) ]
#'     
#'     for (child_of_p in child_of_p_node)
#'     {
#'       canonical_pq[parent_of_p_node, child_of_p] <- TRUE
#'       # print (p_node)
#'       # print (child_of_p)
#'       # print (parent_of_p_node)
#'       # print ("*********")
#'     }
#'   }
#'   
#'   canonical_pq <- delete.vertices(canonical_pq, V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")==1])
#'   #print (V(canonical_pq)[type == "Q" & igraph::degree(canonical_pq, mode="in")==0])
#'   return (canonical_pq)
#' }
#' 
#' count_leaf_descendents <- function(pq_tree, curr_node, children_counts)
#' {
#'   nei <- NULL
#'   
#'   if (V(pq_tree)[curr_node]$type == "leaf")
#'   {
#'     #V(pq_tree)[curr_node]$num_desc = 0
#'     children_counts[curr_node] = 0
#'     return(children_counts)
#'   } else {
#'     children_count = 0
#'     for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
#'     {
#'       children_counts <- count_leaf_descendents(pq_tree, child, children_counts)
#'       if (V(pq_tree)[child]$type == "leaf")
#'       {
#'         children_count <- children_count + 1
#'       }
#'       else
#'       {
#'         children_count <- children_count + children_counts[child]
#'       }
#'     }
#'     #print (curr_node)
#'     children_counts[curr_node] = children_count
#'     return(children_counts)
#'   }
#' }
#' 
#' 
#' #' Return an ordering for a P node in the PQ tree
#' #' @importFrom combinat permn
#' order_p_node <- function(q_level_list, dist_matrix)
#' { 
#'   q_order_res <- combinat::permn(q_level_list, fun=order_q_node, dist_matrix)
#'   #print (q_order_res)
#'   all_perms <- lapply(q_order_res, function(x) { x$ql } )
#'   #print ("perm ql:")
#'   #print(all_perms)
#'   all_perms_weights <- unlist(lapply(q_order_res, function(x) { x$wt }))
#'   #print ("perm weights:")
#'   #print (all_perms_weights)
#'   
#'   opt_perm_idx <- head((which(all_perms_weights == min(all_perms_weights))), 1)
#'   opt_perm <- all_perms[[opt_perm_idx]]
#'   
#'   #print ("opt_path:")
#'   #print (opt_perm)
#'   #print ("opt_all_weight:")
#'   #print (min(all_perms_weights))
#'   #print ("weights:")
#'   #print (all_perms_weights)
#'   # print ("q_level_list:")
#'   # print (q_level_list)
#'   stopifnot (length(opt_perm) == length(q_level_list))
#'   
#'   return(opt_perm)
#' }
#' 
#' order_q_node <- function(q_level_list, dist_matrix)
#' {
#'   new_subtree <- graph.empty()
#'   
#'   if (length(q_level_list) == 1)
#'   {
#'     return (list(ql=q_level_list, wt=0))
#'   }
#'   for (i in 1:length(q_level_list))
#'   {
#'     new_subtree <- new_subtree + vertex(paste(i,"F"), type="forward")
#'     new_subtree <- new_subtree + vertex(paste(i,"R"), type="reverse")
#'   }
#'   
#'   for (i in (1:(length(q_level_list)-1)))
#'   {
#'     cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], q_level_list[[i+1]][1]]
#'     new_subtree <- new_subtree + edge(paste(i,"F"), paste(i+1,"F"), weight=cost)
#'     
#'     cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], q_level_list[[i+1]][length(q_level_list[[i+1]])]]
#'     new_subtree <- new_subtree + edge(paste(i,"F"), paste(i+1,"R"), weight=cost)
#'     
#'     cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i+1]][1]]
#'     new_subtree <- new_subtree + edge(paste(i,"R"), paste(i+1,"F"), weight=cost)
#'     
#'     cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i+1]][length(q_level_list[[i+1]])]]
#'     new_subtree <- new_subtree + edge(paste(i,"R"), paste(i+1,"R"), weight=cost)
#'   }
#'   
#'   first_fwd = V(new_subtree)[paste(1,"F")]
#'   first_rev = V(new_subtree)[paste(1,"R")]
#'   last_fwd = V(new_subtree)[paste(length(q_level_list),"F")]
#'   last_rev = V(new_subtree)[paste(length(q_level_list),"R")]
#'   
#'   FF_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_fwd), to=as.vector(last_fwd), mode="out", output="vpath")$vpath)
#'   FR_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_fwd), to=as.vector(last_rev), mode="out", output="vpath")$vpath)
#'   RF_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_rev), to=as.vector(last_fwd), mode="out", output="vpath")$vpath)
#'   RR_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_rev), to=as.vector(last_rev), mode="out", output="vpath")$vpath)
#'   
#'   # print (FF_path)
#'   # print (FR_path)
#'   # print (RF_path)
#'   # print (RR_path)
#'   
#'   FF_weight <- sum(E(new_subtree, path=FF_path)$weight)
#'   FR_weight <- sum(E(new_subtree, path=FR_path)$weight)
#'   RF_weight <- sum(E(new_subtree, path=RF_path)$weight)
#'   RR_weight <- sum(E(new_subtree, path=RR_path)$weight)
#'   
#'   # print (FF_weight)
#'   # print (FR_weight)
#'   # print (RF_weight)
#'   # print (RR_weight)
#'   
#'   paths <- list(FF_path, FR_path, RF_path, RR_path)
#'   path_weights <- c(FF_weight, FR_weight, RF_weight, RR_weight)
#'   opt_path_idx <- head((which(path_weights == min(path_weights))), 1)
#'   opt_path <- paths[[opt_path_idx]]
#'   
#'   # print ("opt_path:")
#'   # print (opt_path)
#'   # print ("q_level_list:")
#'   # print (q_level_list)
#'   stopifnot (length(opt_path) == length(q_level_list))
#'   
#'   directions <- V(new_subtree)[opt_path]$type
#'   #print (directions)
#'   q_levels <- list()
#'   for (i in 1:length(directions))
#'   {
#'     if (directions[[i]] == "forward"){
#'       q_levels[[length(q_levels)+1]] <- q_level_list[[i]]
#'     }else{
#'       q_levels[[length(q_levels)+1]] <- rev(q_level_list[[i]])
#'     }
#'   }
#'   
#'   return(list(ql=q_levels, wt=min(path_weights)))
#' }
#' 
#' measure_diameter_path <- function(pq_tree, curr_node, path_lengths)
#' {
#'   nei <- NULL
#'   
#'   if (V(pq_tree)[curr_node]$type != "Q")
#'   {
#'     #V(pq_tree)[curr_node]$num_desc = 0
#'     path_lengths[curr_node] = 0
#'     return(path_lengths)
#'   } else {
#'     
#'     children_count = 0
#'     for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
#'     {
#'       children_counts <- count_leaf_descendents(pq_tree, child, children_counts)
#'       if (V(pq_tree)[child]$type == "leaf")
#'       {
#'         children_count <- children_count + 1
#'       }
#'       else
#'       {
#'         children_count <- children_count + children_counts[child]
#'       }
#'     }
#'     
#'     
#'     path_lengths[curr_node] = children_count
#'     return(children_counts)
#'   }
#' }
#' 
#' # Assign leaf nodes reachable in pq_tree from curr_node to assigned_state
#' assign_cell_lineage <- function(pq_tree, curr_node, assigned_state, node_states)
#' {
#'   nei <- NULL
#'   
#'   if (V(pq_tree)[curr_node]$type == "leaf")
#'   {
#'     #V(pq_tree)[curr_node]$num_desc = 0
#'     #print (curr_node)
#'     node_states[V(pq_tree)[curr_node]$name] = assigned_state
#'     return(node_states)
#'   } else {
#'     for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
#'     {
#'       node_states <- assign_cell_lineage(pq_tree, child, assigned_state, node_states)
#'     }
#'     return(node_states)
#'   }
#' }
#' 
#' 
#' extract_good_ordering <- function(pq_tree, curr_node, dist_matrix)
#' {
#'   nei <- NULL
#'   
#'   if (V(pq_tree)[curr_node]$type == "leaf")
#'   {
#'     #print ("ordering leaf node")
#'     return (V(pq_tree)[curr_node]$name)
#'   }else if (V(pq_tree)[curr_node]$type == "P"){
#'     #print ("ordering P node")
#'     p_level <- list()
#'     for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
#'     {
#'       p_level[[length(p_level)+1]] <- extract_good_ordering(pq_tree, child, dist_matrix)
#'     }
#'     p_level <- order_p_node(p_level, dist_matrix)
#'     p_level <- unlist(p_level)
#'     #print (p_level)
#'     return (p_level)
#'   }else if(V(pq_tree)[curr_node]$type == "Q"){
#'     #print ("ordering Q node")
#'     q_level <- list()
#'     for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
#'     {
#'       q_level[[length(q_level)+1]] <- extract_good_ordering(pq_tree, child, dist_matrix)
#'     }
#'     q_level <- order_q_node(q_level, dist_matrix)
#'     q_level <- q_level$ql
#'     q_level <- unlist(q_level)
#'     #print (q_level)
#'     return (q_level)
#'   }
#' }
#' 
#' 
#' 
#' #' Extract a linear ordering of cells from a PQ tree
#' #' @importFrom plyr arrange
#' extract_good_branched_ordering <- function(orig_pq_tree, curr_node, dist_matrix, num_branches, reverse_main_path=FALSE)
#' {
#'   nei <- NULL
#'   type <- NULL
#'   pseudo_time <- NULL
#'   
#'   pq_tree <- orig_pq_tree
#'   
#'   # children_counts <- rep(0, length(as.vector(V(pq_tree))))
#'   #     names(children_counts) <- V(pq_tree)$name
#'   # children_counts <- count_leaf_descendents(pq_tree, curr_node, children_counts)
#'   # 
#'   # branch_node_counts <- children_counts[V(res$subtree)[type == "P"]]
#'   # branch_node_counts <- sort(branch_node_counts, decreasing=TRUE)
#'   # print (branch_node_counts)
#'   
#'   
#'   branch_node_counts <- V(pq_tree)[type == "Q"]$diam_path_len
#'   names(branch_node_counts) <- V(pq_tree)[type == "Q"]$name
#'   if(length(names(branch_node_counts)) < num_branches)
#'     stop('Number of branches attempted is larger than the branches constructed from pq_tree algorithm')
#'   
#'   branch_node_counts <- sort(branch_node_counts, decreasing=TRUE)
#'   #print (branch_node_counts)
#'   
#'   
#'   cell_states <- rep(NA, length(as.vector(V(pq_tree)[type=="leaf"])))
#'   names(cell_states) <- V(pq_tree)[type=="leaf"]$name
#'   
#'   cell_states <- assign_cell_lineage(pq_tree, curr_node, 1, cell_states)
#'   
#'   branch_point_roots <- list()
#'   
#'   # Start building the ordering tree. Each pseudo-time segment will be a node.
#'   branch_tree <- graph.empty()
#'   #root_branch_id <- "Q_1"
#'   #branch_tree <- branch_tree + vertex(root_branch_id)
#'   
#'   for (i in 1:num_branches)
#'   {
#'     #cell_states <- assign_cell_lineage(pq_tree, names(branch_node_counts)[i], i+1, cell_states)
#'     #print (head(cell_states))
#'     #print(names(branch_node_counts)[i])
#'     
#'     branch_point_roots[[length(branch_point_roots) + 1]] <- names(branch_node_counts)[i]
#'     branch_id <- names(branch_node_counts)[i]
#'     #print (branch_id)
#'     branch_tree <- branch_tree + vertex(branch_id)
#'     parents <- V(pq_tree)[suppressWarnings(nei(names(branch_node_counts)[i], mode="in"))]
#'     if (length(parents) > 0 && parents$type == "P")
#'     {
#'       p_node_parent <- V(pq_tree)[suppressWarnings(nei(names(branch_node_counts)[i], mode="in"))]
#'       parent_branch_id <- V(pq_tree)[suppressWarnings(nei(p_node_parent, mode="in"))]$name
#'       #print (parent_branch_id)
#'       #print (branch_id)
#'       branch_tree <- branch_tree + edge(parent_branch_id, branch_id)
#'     }
#'     pq_tree[V(pq_tree) [ suppressWarnings(nei(names(branch_node_counts)[i], mode="in")) ], names(branch_node_counts)[i] ] <- FALSE
#'   }
#'   
#'   #branch_point_roots[[length(branch_point_roots) + 1]] <- curr_node
#'   #branch_point_roots <- rev(branch_point_roots)
#'   branch_pseudotimes <- list()
#'   
#'   for (i in 1:length(branch_point_roots))
#'   {
#'     branch_ordering <- extract_good_ordering(pq_tree, branch_point_roots[[i]], dist_matrix)
#'     branch_ordering_time <- weight_of_ordering(branch_ordering, dist_matrix)
#'     names(branch_ordering_time) <- branch_ordering
#'     branch_pseudotimes[[length(branch_pseudotimes) + 1]] = branch_ordering_time
#'     names(branch_pseudotimes)[length(branch_pseudotimes)] = branch_point_roots[[i]]
#'   }
#'   
#'   cell_ordering_tree <- graph.empty()
#'   curr_branch <- "Q_1"
#'   
#'   extract_branched_ordering_helper <- function(branch_tree, curr_branch, cell_ordering_tree, branch_pseudotimes, dist_matrix, reverse_ordering=FALSE)
#'   {
#'     nei <- NULL
#'     
#'     curr_branch_pseudotimes <- branch_pseudotimes[[curr_branch]]
#'     #print (curr_branch_pseudotimes)
#'     curr_branch_root_cell <- NA
#'     for (i in 1:length(curr_branch_pseudotimes))
#'     {
#'       cell_ordering_tree <- cell_ordering_tree + vertex(names(curr_branch_pseudotimes)[i])
#'       if (i > 1)
#'       {
#'         if (reverse_ordering == FALSE){
#'           cell_ordering_tree <- cell_ordering_tree + edge(names(curr_branch_pseudotimes)[i-1], names(curr_branch_pseudotimes)[i])
#'         }else{
#'           cell_ordering_tree <- cell_ordering_tree + edge(names(curr_branch_pseudotimes)[i], names(curr_branch_pseudotimes)[i-1])
#'         }
#'       }
#'     }
#'     
#'     if (reverse_ordering == FALSE)
#'     {
#'       curr_branch_root_cell <- names(curr_branch_pseudotimes)[1]
#'     }else{
#'       curr_branch_root_cell <- names(curr_branch_pseudotimes)[length(curr_branch_pseudotimes)]
#'     }
#'     
#'     for (child in V(branch_tree) [ suppressWarnings(nei(curr_branch, mode="out")) ])
#'     {
#'       child_cell_ordering_subtree <- graph.empty()
#'       
#'       child_head <- names(branch_pseudotimes[[child]])[1]
#'       child_tail <- names(branch_pseudotimes[[child]])[length(branch_pseudotimes[[child]])]
#'       
#'       # find the closest cell in the parent branch for each of the head and the tail
#'       
#'       curr_branch_cell_names <- names(branch_pseudotimes[[curr_branch]])
#'       head_dist_to_curr <- dist_matrix[child_head, curr_branch_cell_names]
#'       closest_to_head <- names(head_dist_to_curr)[which(head_dist_to_curr == min(head_dist_to_curr))]
#'       
#'       head_dist_to_anchored_branch = NA
#'       branch_index_for_head <- NA
#'       
#'       head_dist_to_anchored_branch <- dist_matrix[closest_to_head, child_head]
#'       
#'       tail_dist_to_curr <- dist_matrix[child_tail, curr_branch_cell_names]
#'       closest_to_tail <- names(tail_dist_to_curr)[which(tail_dist_to_curr == min(tail_dist_to_curr))]
#'       
#'       tail_dist_to_anchored_branch = NA
#'       branch_index_for_tail <- NA
#'       
#'       tail_dist_to_anchored_branch <- dist_matrix[closest_to_tail, child_tail]
#'       
#'       if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch)
#'       {
#'         reverse_child <- TRUE
#'       }else{
#'         reverse_child <- FALSE
#'       }
#'       
#'       res <- extract_branched_ordering_helper(branch_tree, child, child_cell_ordering_subtree, branch_pseudotimes, dist_matrix, reverse_child)
#'       child_cell_ordering_subtree <- res$subtree
#'       child_subtree_root <- res$root
#'       
#'       # Works, but slow:
#'       for (v in V(child_cell_ordering_subtree))
#'       {
#'         cell_ordering_tree <- cell_ordering_tree + vertex(V(child_cell_ordering_subtree)[v]$name)
#'       }
#'       
#'       edge_list <- get.edgelist(child_cell_ordering_subtree)
#'       for (i in 1:nrow(edge_list))
#'       {
#'         cell_ordering_tree <- cell_ordering_tree + edge(V(cell_ordering_tree)[edge_list[i, 1]]$name, V(cell_ordering_tree)[edge_list[i, 2]]$name)
#'       }   					
#'       
#'       if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch)
#'       {
#'         cell_ordering_tree <- cell_ordering_tree + edge(closest_to_tail, child_subtree_root)
#'       }else{
#'         cell_ordering_tree <- cell_ordering_tree + edge(closest_to_head, child_subtree_root)
#'       }
#'       
#'     }
#'     
#'     return (list(subtree=cell_ordering_tree, root=curr_branch_root_cell, last_cell_state=1, last_cell_pseudotime=0.0))
#'   }
#'   
#'   res <- extract_branched_ordering_helper(branch_tree, curr_branch, cell_ordering_tree, branch_pseudotimes, dist_matrix, reverse_main_path)
#'   cell_ordering_tree <- res$subtree
#'   
#'   curr_state <- 1
#'   
#'   assign_cell_state_helper <- function(ordering_tree_res, curr_cell)
#'   {
#'     nei <- NULL
#'     
#'     cell_tree <- ordering_tree_res$subtree
#'     V(cell_tree)[curr_cell]$cell_state = curr_state
#'     
#'     children <- V(cell_tree) [ suppressWarnings(nei(curr_cell, mode="out")) ]
#'     ordering_tree_res$subtree <- cell_tree
#'     
#'     if (length(children) == 1){
#'       ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[children]$name)
#'     }else{
#'       for (child in children)	{
#'         curr_state <<- curr_state + 1
#'         ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[child]$name)
#'       }
#'     }
#'     return (ordering_tree_res)
#'   }
#'   
#'   res <- assign_cell_state_helper(res, res$root)
#'   
#'   assign_pseudotime_helper <- function(ordering_tree_res, dist_matrix, last_pseudotime, curr_cell)
#'   {
#'     nei <- NULL
#'     
#'     cell_tree <- ordering_tree_res$subtree
#'     curr_cell_pseudotime <- last_pseudotime
#'     V(cell_tree)[curr_cell]$pseudotime = curr_cell_pseudotime
#'     V(cell_tree)[curr_cell]$parent =  V(cell_tree)[ suppressWarnings(nei(curr_cell, mode="in")) ]$name
#'     #print (curr_cell_pseudotime)
#'     
#'     ordering_tree_res$subtree <- cell_tree
#'     children <- V(cell_tree) [ suppressWarnings(nei(curr_cell, mode="out")) ]
#'     
#'     for (child in children)	{
#'       next_node <- V(cell_tree)[child]$name
#'       delta_pseudotime <- dist_matrix[curr_cell, next_node]
#'       ordering_tree_res <- assign_pseudotime_helper(ordering_tree_res, dist_matrix, last_pseudotime + delta_pseudotime, next_node)
#'     }
#'     
#'     return (ordering_tree_res)
#'   }
#'   
#'   res <- assign_pseudotime_helper(res, dist_matrix, 0.0, res$root)
#'   
#'   cell_names <- V(res$subtree)$name
#'   cell_states <- V(res$subtree)$cell_state
#'   cell_pseudotime <- V(res$subtree)$pseudotime
#'   cell_parents <- V(res$subtree)$parent
#'   # print (cell_names)
#'   # print (cell_states)
#'   # print (cell_pseudotime)
#'   ordering_df <- data.frame(sample_name = cell_names,
#'                             cell_state = factor(cell_states),
#'                             pseudo_time = cell_pseudotime,
#'                             parent = cell_parents)
#'   
#'   ordering_df <- plyr::arrange(ordering_df, pseudo_time)
#'   return(list("ordering_df"=ordering_df, "cell_ordering_tree"=cell_ordering_tree))
#' }
#' 
#' reverse_ordering <- function(pseudo_time_ordering)
#' {
#'   pt <- pseudo_time_ordering$pseudo_time
#'   names(pt) <- pseudo_time_ordering$sample_name
#'   rev_pt <- -((pt - max(pt)))
#'   rev_df <- pseudo_time_ordering
#'   rev_df$pseudo_time <- rev_pt
#'   return(rev_df)
#' }
#' 
#' 
#' weight_of_ordering <- function(ordering, dist_matrix)
#' {
#'   time_delta <- c(0)
#'   curr_weight <- 0
#'   ep <- 0.01
#'   for (i in 2:length(ordering))
#'   {
#'     d <- dist_matrix[ordering[[i]], ordering[[i-1]]]
#'     curr_weight <- curr_weight + d + ep
#'     time_delta <- c(time_delta, curr_weight)
#'   }
#'   
#'   return(time_delta)
#' }
#' 
#' 
#' #' Sets the features (e.g. genes) to be used for ordering cells in pseudotime.
#' #' @param cds the CellDataSet upon which to perform this operation
#' #' @param ordering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' #' @return an updated CellDataSet object
#' #' @export
#' setOrderingFilter <- function(cds, ordering_genes){
#'   fData(cds)$use_for_ordering <- row.names(fData(cds)) %in% ordering_genes
#'   cds
#' }
#' 
#' #' Run the fastICA algorithm on a numeric matrix.
#' #' @importFrom fastICA  ica.R.def ica.R.par
#' #' @importFrom irlba irlba
#' #' 
#' ica_helper <- function(X, n.comp, alg.typ = c("parallel", "deflation"), fun = c("logcosh", "exp"), alpha = 1, 
#'                        row.norm = TRUE, maxit = 200, tol = 1e-4, verbose = FALSE, w.init = NULL, use_irlba=TRUE){
#'   dd <- dim(X) 
#'   #FIXME: This will internally convert to a dense matrix
#'   d <- dd[dd != 1L]
#'   if (length(d) != 2L) 
#'     stop("data must be matrix-conformal")
#'   X <- if (length(d) != length(dd)) 
#'     matrix(X, d[1L], d[2L])
#'   else as.matrix(X)
#'   if (alpha < 1 || alpha > 2) 
#'     stop("alpha must be in range [1,2]")
#'   alg.typ <- match.arg(alg.typ)
#'   fun <- match.arg(fun)
#'   n <- nrow(X)
#'   p <- ncol(X)
#'   if (n.comp > min(n, p)) {
#'     message("'n.comp' is too large: reset to ", min(n, p))
#'     n.comp <- min(n, p)
#'   }
#'   if (is.null(w.init)) 
#'     w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
#'   else {
#'     if (!is.matrix(w.init) || length(w.init) != (n.comp^2)) 
#'       stop("w.init is not a matrix or is the wrong size")
#'   }
#'   
#'   if (verbose) 
#'     message("Centering")
#'   X <- scale(X, scale = FALSE)
#'   X <- if (row.norm) 
#'     t(scale(X, scale = row.norm))
#'   else t(X)
#'   if (verbose) 
#'     message("Whitening")
#'   V <- X %*% t(X)/n
#'   
#'   if (verbose) 
#'     message("Finding SVD")
#'   
#'   s <- irlba::irlba(V, n.comp, n.comp)  
#'   svs <- s$d  
#'   
#'   D <- diag(c(1/sqrt(s$d)))
#'   K <- D %*% t(s$u)
#'   K <- matrix(K[1:n.comp, ], n.comp, p)
#'   X1 <- K %*% X
#'   
#'   if (verbose) 
#'     message("Running ICA")
#'   if (alg.typ == "deflation") {
#'     a <- fastICA::ica.R.def(X1, n.comp, tol = tol, fun = fun, 
#'                             alpha = alpha, maxit = maxit, verbose = verbose, 
#'                             w.init = w.init)
#'   }
#'   else if (alg.typ == "parallel") {
#'     a <- fastICA::ica.R.par(X1, n.comp, tol = tol, fun = fun, 
#'                             alpha = alpha, maxit = maxit, verbose = verbose, 
#'                             w.init = w.init)
#'   }
#'   w <- a %*% K
#'   S <- w %*% X
#'   A <- t(w) %*% solve(w %*% t(w))
#'   return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S), svs=svs))
#' }
#' 
#' #
#' # extract_ddrtree_ordering <- function(cds, root_cell, verbose=T)
#' # {
#' #   nei <- NULL
#' #   type <- NULL
#' #   pseudo_time <- NULL
#' #   
#' #   dp_mst <- minSpanningTree(cds) 
#' #   
#' #   Pseudotime <- shortest.paths(dp_mst, v=root_cell, to=V(dp_mst))
#' #   
#' #   curr_state <- 1
#' #   #' a function to assign pseudotime for the MST
#' #   assign_cell_state_helper <- function(ordering_tree_res, curr_cell, visited_node = curr_cell) {
#' #     nei <- NULL
#' #     
#' #     cell_tree <- ordering_tree_res$subtree
#' #     V(cell_tree)[curr_cell]$cell_state = curr_state
#' #     
#' #     children <- V(cell_tree) [ suppressWarnings(nei(curr_cell, mode="all")) ]$name 
#' #     children <- setdiff(children, V(cell_tree)[visited_node]$name)
#' #     
#' #     ordering_tree_res$subtree <- cell_tree
#' #     #V(ordering_tree_res$subtree)[curr_cell]$parent = visited_node
#' #     
#' #     if (length(children) == 0){
#' #       return (ordering_tree_res)
#' #     }else if (length(children) == 1){
#' #       #visited_node <- union(children, visited_node)
#' #       V(ordering_tree_res$subtree)[children]$parent = V(cell_tree)[curr_cell]$name
#' #       ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[children]$name, curr_cell)
#' #     }else{
#' #       for (child in children) {
#' #         #visited_node <- union(child, visited_node)
#' #         V(ordering_tree_res$subtree)[children]$parent = rep(V(cell_tree)[curr_cell]$name, length(children))
#' #         curr_state <<- curr_state + 1
#' #         ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[child]$name, curr_cell)
#' #       }
#' #     }
#' #     return (ordering_tree_res)
#' #   }
#' #   
#' #   res <- list(subtree = dp_mst, root = root_cell)
#' #   #V(res$subtree)$parent <- rep(NA, nrow(pData(cds)))
#' #   res <- assign_cell_state_helper(res, res$root)
#' #   
#' #   states <- V(res$subtree)$cell_state
#' #   
#' #   cell_names <-  colnames(Pseudotime)
#' #   cell_states <- states
#' #   cell_pseudotime <- Pseudotime
#' #   cell_parents <- V(res$subtree)$parent
#' #   
#' #   ordering_df <- data.frame(sample_name = cell_names,
#' #                             cell_state = factor(cell_states),
#' #                             pseudo_time = as.vector(cell_pseudotime),
#' #                             parent = cell_parents)
#' #   row.names(ordering_df) <- ordering_df$sample_name
#' #   # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
#' #   return(ordering_df)
#' # }
#' 
#' extract_ddrtree_ordering <- function(cds, root_cell, verbose=T)
#' {
#'   nei <- NULL
#'   type <- NULL
#'   pseudo_time <- NULL
#'   
#'   dp <- cellPairwiseDistances(cds) 
#'   dp_mst <- minSpanningTree(cds) 
#'   
#'   curr_state <- 1
#'   
#'   res <- list(subtree = dp_mst, root = root_cell)
#'   
#'   states = rep(1, ncol(dp))
#'   names(states) <- V(dp_mst)$name
#'   
#'   pseudotimes = rep(0, ncol(dp))
#'   names(pseudotimes) <- V(dp_mst)$name
#'   
#'   parents = rep(NA, ncol(dp))
#'   names(parents) <- V(dp_mst)$name
#'   
#'   mst_traversal <- graph.dfs(dp_mst, 
#'                              root=root_cell, 
#'                              neimode = "all", 
#'                              unreachable=FALSE, 
#'                              father=TRUE)
#'   mst_traversal$father <- as.numeric(mst_traversal$father)
#'   curr_state <- 1
#'   
#'   for (i in 1:length(mst_traversal$order)){
#'     curr_node = mst_traversal$order[i]
#'     curr_node_name = V(dp_mst)[curr_node]$name
#'     
#'     if (is.na(mst_traversal$father[curr_node]) == FALSE){
#'       parent_node = mst_traversal$father[curr_node]
#'       parent_node_name = V(dp_mst)[parent_node]$name
#'       parent_node_pseudotime = pseudotimes[parent_node_name]
#'       parent_node_state = states[parent_node_name]
#'       curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
#'       if (degree(dp_mst, v=parent_node_name) > 2){
#'         curr_state <- curr_state + 1
#'       }
#'     }else{
#'       parent_node = NA
#'       parent_node_name = NA
#'       curr_node_pseudotime = 0
#'     }
#'     
#'     curr_node_state = curr_state
#'     pseudotimes[curr_node_name] <- curr_node_pseudotime
#'     states[curr_node_name] <- curr_node_state
#'     parents[curr_node_name] <- parent_node_name
#'   }
#'   
#'   ordering_df <- data.frame(sample_name = names(states),
#'                             cell_state = factor(states),
#'                             pseudo_time = as.vector(pseudotimes),
#'                             parent = parents)
#'   row.names(ordering_df) <- ordering_df$sample_name
#'   # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
#'   return(ordering_df)
#' }
#' 
#' 
select_root_cell <- function(cds, root_state=NULL, reverse=FALSE){
  if (is.null(root_state) == FALSE) {
    if (is.null(pData(cds)$State)){
      stop("Error: State has not yet been set. Please call orderCells() without specifying root_state, then try this call again.")
    }
    # FIXME: Need to gaurd against the case when the supplied root state isn't actually a terminal state in the tree.
    root_cell_candidates <- subset(pData(cds), State == root_state)
    if (nrow(root_cell_candidates) == 0){
      stop(paste("Error: no cells for State =", root_state))
    }

    # build a local MST to find a good root cell for this state
    dp <- as.matrix(dist(t(reducedDimS(cds)[,row.names(root_cell_candidates)])))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)

    # Make sure to use the real MST here
    tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
    #root_cell_candidates <- root_cell_candidates[row.names(root_cell_candidates) %in% tip_leaves,]
    #sg <- make_ego_graph(dp_mst, nodes=row.names(root_cell_candidates))[[1]]

    diameter <- get.diameter(dp_mst)

    if (length(diameter) == 0){
      stop(paste("Error: no valid root cells for State =", root_state))
    }

    #root_cell = names(diameter)[tip_leaves %in% names(diameter)]
    root_cell_candidates <- root_cell_candidates[names(diameter),]
    if (root_state == 1){
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == min(root_cell_candidates$Pseudotime))]
    }else{
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == max(root_cell_candidates$Pseudotime))]
    }
    if (length(root_cell) > 1)
      root_cell <- root_cell[1]

    # If we used DDRTree, we need to go from this actual cell to the nearst
    # point on the principal graph
    if (cds@dim_reduce_type == "DDRTree"){
      #root_cell_idx <- which(V(minSpanningTree(cds))$name == root_cell, arr.ind=T)
      graph_point_for_root_cell <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[root_cell,]
      root_cell = V(minSpanningTree(cds))[graph_point_for_root_cell]$name
    }

  }else{
    if (is.null(minSpanningTree(cds))){
      stop("Error: no spanning tree found for CellDataSet object. Please call reduceDimension before calling orderCells()")
    }
    diameter <- get.diameter(minSpanningTree(cds))
    if (is.null(reverse) == FALSE && reverse == TRUE){
      root_cell = names(diameter[length(diameter)])
    } else {
      root_cell = names(diameter[1])
    }
  }
  return(root_cell)
}

#' 
#' normalize_expr_data <- function(cds, 
#'                                 norm_method = c("vstExprs", "log", "none"), 
#'                                 pseudo_expr = NULL){
#'   FM <- exprs(cds)
#'   
#'   # If the user has selected a subset of genes for use in ordering the cells
#'   # via setOrderingFilter(), subset the expression matrix.
#'   if (is.null(fData(cds)$use_for_ordering) == FALSE && 
#'       nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
#'     FM <- FM[fData(cds)$use_for_ordering, ]
#'   }
#'   
#'   norm_method <- match.arg(norm_method)
#'   if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
#'     
#'     # If we're going to be using log, and the user hasn't given us a pseudocount
#'     # set it to 1 by default.
#'     if (is.null(pseudo_expr)){
#'       if(norm_method == "log")
#'         pseudo_expr = 1
#'       else
#'         pseudo_expr = 0
#'     }
#'     
#'     if (norm_method == "vstExprs") {
#'       checkSizeFactors(cds)
#'       
#'       if (is.null(fData(cds)$use_for_ordering) == FALSE && 
#'           nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
#'         VST_FM <- vstExprs(cds[fData(cds)$use_for_ordering,], round_vals = FALSE)
#'       }else{
#'         VST_FM <- vstExprs(cds, round_vals = FALSE)
#'       }
#'       
#'       if (is.null(VST_FM) == FALSE) {
#'         FM <- VST_FM
#'       }
#'       else {
#'         stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
#'       }
#'     }else if (norm_method == "log") {
#'       # If we are using log, normalize by size factor before log-transforming
#'       FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
#'       FM <- FM + pseudo_expr
#'       FM <- log2(FM)
#'     }else if (norm_method == "none"){
#'       # If we are using log, normalize by size factor before log-transforming
#'       FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
#'       FM <- FM + pseudo_expr
#'     }
#'   }else if (cds@expressionFamily@vfamily == "binomialff") {
#'     if (norm_method == "none"){
#'       #If this is binomial data, transform expression values into TF-IDF scores.
#'       ncounts <- FM > 0
#'       ncounts[ncounts != 0] <- 1
#'       FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
#'     }else{
#'       stop("Error: the only normalization method supported with binomial data is 'none'")
#'     }
#'   }else if (cds@expressionFamily@vfamily == "Tobit") {
#'     FM <- FM + pseudo_expr
#'     if (norm_method == "none"){
#'       
#'     }else if (norm_method == "log"){
#'       FM <- log2(FM)
#'     }else{
#'       stop("Error: the only normalization methods supported with Tobit-distributed (e.g. FPKM/TPM) data are 'log' (recommended) or 'none'")
#'     }
#'   }else if (cds@expressionFamily@vfamily == "gaussianff") {
#'     if (norm_method == "none"){
#'       FM <- FM + pseudo_expr
#'     }else{
#'       stop("Error: the only normalization method supported with gaussian data is 'none'")
#'     }
#'   }
#'   return (FM)
#' }

#this helper function is used to convert the data into the format for reducing dimension reduction in R
convert2DRData <- function(cds,
         max_components=2,
         reduction_method=c("DDRTree", "ICA"),
         norm_method = c("log", "vstExprs", "none"),
         residualModelFormulaStr=NULL,
         pseudo_expr=NULL,
         verbose=FALSE,
         scaling = TRUE, 
         ...){

  FM <- monocle:::normalize_expr_data(cds, norm_method, pseudo_expr)

  #FM <- FM[unlist(sparseApply(FM, 1, sd, convert_to_dense=TRUE)) > 0, ]
  xm <- Matrix::rowMeans(FM)
  xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
  FM <- FM[xsd > 0,]

  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose)
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
                                       data = pData(cds), drop.unused.levels = TRUE)

    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
  }else{
    X.model_mat <- NULL
  }

  if(scaling)
   FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))

  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  return(FM)
}

custom_ordering <- function(DDRTree_res, root_cell = NULL, branch_num = NULL){
  #library(devtools)
  #load_all

  dp <- as.matrix(dist(t(DDRTree_res$Y)))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)

  if(is.null(root_cell)){
    root_cell <- which(degree(dp_mst) == 1)[1]
  }
  else{
    root_cell_candidates <- which(degree(dp_mst) == 1)[1]
    Pseudotime <- shortest.paths(dp_mst, v=root_cell, to=V(dp_mst)[root_cell_candidates])
    root_cell <- root_cell_candidates[Pseudotime == min(Pseudotime)]
  }
  Pseudotime <- shortest.paths(dp_mst, v=root_cell, to=V(dp_mst))
  names(Pseudotime) <- 1:ncol(DDRTree_res$Y)

  cc_ordering <- extract_ddrtree_ordering_xj(dp_mst = dp_mst, dp = dp, root_cell = root_cell)
  #assign the branches as well as the branch time point:
  # next_node <<- 0
  # 
  #   res <- pq_helper(dp_mst, use_weights=FALSE, root_node=root_cell)
  # 
  #   if(is.null(branch_num))
  #   branch_num <- sum(degree(dp_mst) > 2) + 1
  # 
  #   order_list <- extract_good_branched_ordering(res$subtree, res$root, dp, branch_num, FALSE)
  # 
  #   cc_ordering <- order_list$ordering_df
  #   row.names(cc_ordering) <- cc_ordering$sample_name

  mst_branch_nodes <- V(dp_mst)[which(degree(dp_mst) > 2)]$name #calculate the branch point

  return(list(Pseudotime = Pseudotime, cc_ordering = cc_ordering, mst_branch_nodes = mst_branch_nodes, dp_mst = dp_mst))
}

find_root_cells <- function(ordering_res, root_state) {
  root_cell_candidates <- subset(ordering_res$cc_ordering, cell_state == root_state)
  if (root_state == 1){
    root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$pseudo_time == min(root_cell_candidates$pseudo_time))]
  }else{
    root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$pseudo_time == max(root_cell_candidates$pseudo_time))]
  }
  return(root_cell)
}

calSparsity <- function(matrix){
  1 - sum(matrix != 0) / length(matrix)
}

cal_cross_map <- function(ordered_exprs_mat, lib_colum, target_column) {
  lib_xmap_target <- ccm(ordered_exprs_mat, E = 3, random_libs = TRUE, lib_column = lib_colum, #ENSG00000122180.4
                         target_column = target_column, lib_sizes = seq(10, 75, by = 5), num_samples = 300)
  target_xmap_lib <- ccm(ordered_exprs_mat, E = 3, random_libs = TRUE, lib_column = target_column, #ENSG00000122180.4
                         target_column = lib_colum, lib_sizes = seq(10, 75, by = 5), num_samples = 300)

  lib_xmap_target_means <- ccm_means(lib_xmap_target)
  target_xmap_lib_means <- ccm_means(target_xmap_lib)

  return(list(lib_xmap_target = lib_xmap_target, target_xmap_lib = target_xmap_lib,
              lib_xmap_target_means = lib_xmap_target_means, target_xmap_lib_means = target_xmap_lib_means))
}

#parallel the CCM algorithm:
parallelCCM <- function(ordered_exprs_mat, cores = detectCores()) {
  ordered_exprs_mat <- ordered_exprs_mat
  combn_mat <- combn(1:ncol(ordered_exprs_mat), 2)

  combn_mat_split <- split(t(combn_mat), 1:ncol(combn_mat))
  CCM_res <- mclapply(combn_mat_split, function(x, ordered_exprs_mat){
    col_names <- colnames(ordered_exprs_mat)[x]
    cal_cross_map(ordered_exprs_mat[, col_names], col_names[1], col_names[2])
  }, ordered_exprs_mat = ordered_exprs_mat, mc.cores = cores)

  return(list(CCM_res = CCM_res, combn_mat_split = combn_mat_split))
}

################################################################################################################
################################################################################################################

#following are functions for different dimension reduction methods:
PCA <- function(data, ...) {
  res <- prcomp(t(data), center = F, scale = F)
  res$x
}
ICA <- function(data, max_components = 3, ...) {
  data <- data[, !(duplicated(t(data)))]
  monocle:::ica_helper(t(data), n.comp = max_components,  use_irlba = F)$S #dimension need to be less than variable
}

LLE2 <- function(data, m = 2, k = 5){
  lle_reduced_data <- lle(t(data), m = m, k = k)
  res <- t(lle_reduced_data$Y)
  colnames(res) <- colnames(data)
  return(t(res))
}

LLE <- function (data, max_components = 3, num_neigh = NULL, reg = 2, ss = FALSE,
                 id = TRUE, v = 0.9, iLLE = FALSE) {
  if(is.null(num_neigh)) {
    # num_neigh_list <- calc_k(t(data), m = max_components, kmin = 1,
    #                          kmax = 20, plotres = TRUE, parallel = T, cpus = detectCores(), iLLE = iLLE)
    # num_neigh <- num_neigh_list$k[which(num_neigh_list$rho ==
    #                                       min(num_neigh_list$rho))]
    
    #use slicer's approach: 
    k = SLICER::select_k(t(data), kmin=5)
    message('k is ', k)
  }
  lle_reduced_data <- lle(t(data), m = max_components, k = k,
                          reg = reg, ss = ss, id = id, v = v, iLLE = iLLE)

  res <- t(lle_reduced_data$Y)
  colnames(res) <- colnames(data)
  return(t(res))
}
ISOMAP <- function(data, max_components = 3){
  tmp <- isomap(dist(t(data)), ndim = max_components, k = 3, fragmentedOK = T)
  res <- tmp$points
  row.names(res) <- colnames(X)
}

#run the default tSNE: 
tSNE <- function (data, max_components = 3, max_iter = 500, min_cost = 0,
                  whiten = TRUE, epoch = 100) {
  tsne_data <- tsne(t(data), k = max_components, epoch_callback = NULL,
                    max_iter = max_iter, min_cost = min_cost, whiten = whiten,
                    epoch = epoch)
  row.names(tsne_data) <- colnames(data)
  return(tsne_data)
}

#run Van der Maaten's Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding
R_tSNE <- function (data, max_components = 3, initial_dims = 50, perplexity = 30, 
                    theta = 0.5, check_duplicates = TRUE, pca = TRUE, max_iter = 1000,
                    verbose = FALSE, is_distance = FALSE, Y_init = NULL, ... ) {
  data <- data[, !duplicated(as.matrix(t(data)))]
  tsne_res <- Rtsne(as.matrix(t(data)), dims = max_components, initial_dims = initial_dims, perplexity = perplexity,
                                  theta = theta, check_duplicates = check_duplicates, pca = pca, max_iter = max_iter,
                                  verbose = verbose, is_distance = is_distance, Y_init = Y_init, ...)
  tsne_data <- tsne_res$Y[, 1:max_components]
  row.names(tsne_data) <- colnames(tsne_data)
  return(tsne_data)
}

#determine how many clusters we want: 
determine_cluster_num <- function(data, gaussian = TRUE, rho=2, delta=2){
  dataDist <- dist(t(data))
  dataClust <- densityClust(dataDist, gaussian=gaussian)
  plot(dataClust) # Inspect clustering attributes to define thresholds
  
  #automatically pick up the rho and delta values: 
  
  dataClust <- findClusters(dataClust)
  #find the number of clusters: 
  cluster_num <- length(unique(dataClust$clusters))
  
  return(list(dataClust = dataClust, cluster_num = cluster_num))
}

DM <- function(data, eps.val = NULL, neigen = NULL, t = 0, max_components = 3, delta=10^-5){
  D <- dist(t(data))
  if(is.null(eps.val))
    eps.val = epsilonCompute(D)
  tmp <- diffuse(D, eps.val = eps.val, neigen = neigen, t = t, maxdim = max_components, delta=delta)
  res <- tmp$X
  row.names(res) <- colnames(data)
  return(res)
}

diffusion_maps <- function (data, bw_ini = 0, pseudo_cnt = 1, neighbours = 0.2,
                            log2_data = F, max_components = 3) {
  if (log2_data)
    data <- t(log2(data + pseudo_cnt))
  else data <- t(data)
  data_dm_bw_res <- diffusion_maps_bw(data, pseudo_cnt = pseudo_cnt,
                                      neighbours = neighbours, bw_ini = 0, iter = 100, step = 0.02,
                                      log2_data = log2_data)
  bw = 10^(0.2 * (which(data_dm_bw_res$av_d_sigma == max(data_dm_bw_res$av_d_sigma)) +
                    1))
  nn <- ceiling(nrow(data) * neighbours)
  d2 <- as.matrix(dist(data)^2)
  sigma <- bw^2
  W <- exp(-d2/(2 * sigma))
  R <- apply(d2, 2, function(x) sort(x)[nn])
  R <- matrix(rep(R, ncol(d2)), ncol = ncol(d2))
  W <- (d2 < R) * W
  W <- W + t(W)
  D <- colSums(W, na.rm = T)
  q <- D %*% t(D)
  diag(W) <- 0
  H <- W/q
  colS <- colSums(H)
  Hp <- t(t(H)/colS)
  E <- eigen(Hp)
  eigOrd <- order(Re(E$values), decreasing = TRUE)
  E$values <- E$values[eigOrd][-1]
  E$vectors <- E$vectors[, eigOrd][, -1]
  rownames(E$vectors) <- rownames(data)
  colnames(E$vectors) <- 1:ncol(E$vectors)
  diffMap <- t(E$vectors)
  return(t(diffMap[1:max_components, ]))
}

#use destiny to perform the diffusion map dimension reduction
destiny_diffusionMaps <- function(data, max_components = 3, sigma = NULL, k = find.dm.k(nrow(data) - 1L),
                                  n.eigs = min(20L, nrow(data) - 2L), density.norm = TRUE,
                                  ..., distance = c("euclidean", "cosine", "rankcor"), censor.val = NULL,
                                  censor.range = NULL, missing.range = NULL, vars = NULL, verbose = !is.null(censor.range)){
  data <- t(data)
  tmp <- destiny::DiffusionMap(data, sigma = sigma, k = k,
                        n.eigs = n.eigs, density.norm = density.norm, distance = distance, censor.val = censor.val,
                        censor.range = censor.range, missing.range = missing.range, vars = vars, verbose = verbose)
  res <- tmp@eigenvectors
  row.names(res) <- row.names(data)
  
  return(res[, 1:max_components])
}

reduceDimension_custom <- function (cds, max_components = 2, method = c("ICA", "DDRTree"),
                                    pseudo_expr = NULL, residualModelFormulaStr = NULL, use_vst = NULL,
                                    verbose = FALSE, use_irlba = NULL, ...)
{
  FM <- exprs(cds)
  if (is.null(use_irlba)) {
    message("Warning: argument 'use_irlba' is deprecated and will be removed in a future release")
  }
  if (is.null(use_vst) && cds@expressionFamily@vfamily == "negbinomial") {
    use_vst = TRUE
    pseudo_expr = 0
  }
  if (cds@expressionFamily@vfamily == "negbinomial") {
    if (is.null(use_vst))
      use_vst = TRUE
    if (is.null(pseudo_expr))
      pseudo_expr = 0
  }
  else {
    if (is.null(use_vst))
      use_vst = FALSE
    if (is.null(pseudo_expr))
      pseudo_expr = 1
  }
  if (use_vst == FALSE && cds@expressionFamily@vfamily == "negbinomial") {
    checkSizeFactors(cds)
    size_factors <- sizeFactors(cds)
    FM <- Matrix::t(Matrix::t(FM)/size_factors)
  }
  if (is.null(fData(cds)$use_for_ordering) == FALSE && nrow(subset(fData(cds),
                                                                   use_for_ordering == TRUE)) > 0)
    FM <- FM[fData(cds)$use_for_ordering, ]
  if (cds@expressionFamily@vfamily == "binomialff") {
    ncounts <- FM > 0
    ncounts[ncounts != 0] <- 1
    FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
  }
  if (cds@expressionFamily@vfamily != "binomialff") {
    FM <- FM + pseudo_expr
  }
  FM <- FM[apply(FM, 1, sd) > 0, ]
  if (cds@expressionFamily@vfamily != "binomialff") {
    if (use_vst) {
      VST_FM <- vstExprs(cds, expr_matrix = FM, round_vals = FALSE)
      if (is.null(VST_FM) == FALSE) {
        FM <- VST_FM
      }
      else {
        stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
      }
    }
    else {
      FM <- log2(FM)
    }
  }
  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose)
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
                                       data = pData(cds), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
    if (cds@expressionFamily@vfamily != "binomialff") {
      if (use_vst == FALSE) {
        FM <- 2^FM
      }
    }
  }
  if (verbose)
    message("Reducing to independent components")
  FM <- Matrix::t(scale(Matrix::t(FM)))
  FM <- FM[apply(FM, 1, sd) > 0, ]
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  if (is.function(method)) {
    reducedDim <- method(FM, ...)
    return(reducedDim)
  }
}

################################################################################################################
###big version of PCA, ICA, LLE, ISOMAP, tSNE, Diffusion map, GPLVM
#‘bigpca’
# install.packages(c('reader', 'NCmisc', 'bigmemory', 'biganalytics', 'bigmemory.sri', 'bigpca'))

MAPPER <- function(mapper_func = mapper2D,
                   distance_matrix = dist(data.frame(x = 2 * cos(0.5 * (1:100)), y = sin(1:100))),
                   filter_values = sin(1:100), #2 * cos(0.5 * (1:100)),
                   num_intervals = 10,
                   percent_overlap = 50,
                   num_bins_when_clustering = 10){
  res <- mapper_func(distance_matrix = distance_matrix,
                     filter_values = filter_values, #2 * cos(0.5 * (1:100)),
                     num_intervals = num_intervals,
                     percent_overlap = percent_overlap,
                     num_bins_when_clustering = num_bins_when_clustering)

  return(res)
}

mapperVertices <- function(m, pt_labels) {

  # Hovering over vertices gives the point labels:
  # convert the list of vectors of point indices to a list of vectors of labels
  labels_in_vertex <- lapply( m$points_in_vertex, FUN=function(v){ pt_labels[v] } )
  nodename <- sapply( sapply(labels_in_vertex, as.character), paste0, collapse=", ")
  nodename <- paste0("V", 1:m$num_vertices, ": ", nodename )

  # Hovering over vertices gives the point indices:
  # list the points in each vertex
  # nodename <- sapply( sapply(m$points_in_vertex, as.character), paste0, collapse=", ")
  # concatenate the vertex number with the labels for the points in each vertex
  #nodename <- paste0("V", 1:m$num_vertices, ": ", nodename )

  nodegroup <- m$level_of_vertex
  nodesize <- sapply(m$points_in_vertex, length)

  return(data.frame( Nodename=nodename,
                     Nodegroup=nodegroup,
                     Nodesize=nodesize ))

}

mapperEdges <- function(m) {
  linksource <- c()
  linktarget <- c()
  linkvalue <- c()
  k <- 1
  for (i in 2:m$num_vertices) {
    for (j in 1:(i-1)) {
      if (m$adjacency[i,j] == 1) {
        linksource[k] <- i-1
        linktarget[k] <- j-1
        linkvalue[k] <- 2
        k <- k+1
      }
    }
  }
  return( data.frame( Linksource=linksource,
                      Linktarget=linktarget,
                      Linkvalue=linkvalue ) )

}

#use the ordering from DDRTree:
#calculate the ordering of cells from DDRTree:

#DDRTree_res$Z for all the downstream analysis:
DDRTree_analysis <- function(DDRTree_res, root_cell = NULL) {
  dp <- as.matrix(dist(t(DDRTree_res$Y)))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)

  qplot(DDRTree_res$Y[1, ], DDRTree_res$Y[2, ]) #
  if(is.null(root_cell)) {
    root_cell <- which(degree(dp_mst, mode = "total") == 1)[1] #root cell
  }

  Pseudotime <- shortest.paths(dp_mst, v=root_cell, to=V(dp_mst))

  #assign the branches as well as the branch time point:
  next_node <<- 0
  res <- pq_helper(dp_mst, use_weights=FALSE, root_node=root_cell)

  order_list <- extract_good_branched_ordering(res$subtree, res$root, dp, 2, FALSE)

  #sample_name cell_state pseudo_time parent
  #1          32          1   0.0000000   <NA>

  cc_ordering <- order_list$ordering_df
  row.names(cc_ordering) <- cc_ordering$sample_name

  mst_branch_nodes <- V(dp_mst)[which(degree(dp_mst) > 2)]$name #calculate the branch point

  order_list$ordering_df$index <- as.numeric(as.character(order_list$ordering_df$sample_name)) #convert factor to numeric
  ori_order_list <- order_list$ordering_df[sort(order_list$ordering_df$index, index.return = T)$ix, ]
}

order_cell_tree <- function(order_genes, cds, initial_method = PCA, norm_method = 'log',order_by="Time") {
  cds <- setOrderingFilter(cds, order_genes)
  plot_ordering_genes(cds)
  
  cds <- reduceDimension(cds, max_components = 2, norm_method = norm_method, initial_method = initial_method, verbose = T, auto_param_selection = F) #, initial_method = DM , initial_method = DM
  cds <- orderCells(cds, num_paths=1, root_state = NULL)
  PD <- pData(cds)
  
  if(!is.null(order_by)){
    avg_pseudotime <- ddply(PD, .(Time), function(x) mean(x$Pseudotime))
    starting_cell_states <- apply(table(PD[, c("Time", "State")]), 1, function(x) which(x == max(x)))[1]
    cds <- orderCells(cds, num_paths=1, root_state = starting_cell_states)
  }
  # return(pData(HSMM_myo))
  return(cds)
}

#calculate smoothness of the trajectory: standard deviation of the derviative divided by the mean of the derivative
cal_smoothness <- function(cds, order_by = 'Pseudotime', States = unique(pData(cds)$State)){
  cds <- cds[, pData(cds)$State %in% States]
  sorted_cds <- cds[, order(pData(cds)[, order_by])]
  esApply(sorted_cds, 1, function(x) {
    x <- log(x + 1)
    sd(diff(x))/abs(mean(diff(x))) #scale free statistics; coefficient of variance 
  })
}

#calculate the smoothness based on the average expression of a group of genes: 
cal_smoothness_gene_clusters <- function(cds, clusters, order_by = 'Pseudotime', States = unique(pData(cds)$State)){
  cds <- cds[, pData(cds)$State %in% States]
  ordered_cds <- cds[names(clusters), order(pData(cds)[, order_by])]
  smoothness_list <- lapply(unique(clusters), function(x){
    avg_cluster_genes <- esApply(ordered_cds[clusters == x, ], 2, mean)
    sd(diff(avg_cluster_genes))/abs(mean(diff(avg_cluster_genes)))
  })
  return(unlist(smoothness_list))
}

#calculate tightness of the trajectory: tightness is defined as square of residual error divided by variance of residual error
cal_tightness_var <- function(cds, order_by = 'Pseudotime', States = unique(pData(cds)$State)){
  cds <- cds[, pData(cds)$State %in% States]
  sum_residual_error <- esApply(cds, 1, function(x, Pseudotime) {
    tryCatch({
      x <- log(x + 1)
      plx <- predict(loess(Pseudotime ~ x), se=T)
      sum_error <- sum((plx$fit - x)^2)
      vars <- var(x)
      vars <- sum_error / vars
      percentage <- mean((abs(plx$fit - x) / x)[x != 0]) #average percentage while avoiding the zero values 
      res <- data.frame(vars = vars, percentage = percentage)
      return(res)
    }, error = function(e){
      print(e)
      res <- data.frame(vars = NA, percentage = NA)
      return(res)
    })
    
  }, Pseudotime = pData(cds)[, order_by])
  
}

#calculate the lag 1 auto-correlation: 
cal_autocorrelation_lag <- function(cds, ignore_zero = F, order_by = 'Pseudotime', States = unique(pData(cds)$State)){
  cds <- cds[, pData(cds)$State %in% States]
  ordered_cds <- cds[, order(pData(cds)[, order_by])]
  autocorrelation <- esApply(ordered_cds, 1, function(x){
    if(ignore_zero)
      x <- x[x > 0]
    cor(x[-length(x)], x[-1])
  })
  return(autocorrelation)
}

#calculate the autocorrelation based on the average expression of a group of genes: 
cal_autocorrelation_gene_clusters <- function(cds, clusters, order_by = 'Pseudotime', States = unique(pData(cds)$State)){
  cds <- cds[, pData(cds)$State %in% States]
  ordered_cds <- cds[names(clusters), order(pData(cds)[, order_by])]
  autocorrelation_list <- lapply(unique(clusters), function(x){
    avg_cluster_genes <- esApply(ordered_cds[clusters == x, ], 2, mean)
    cor(avg_cluster_genes[-length(avg_cluster_genes)], avg_cluster_genes[-1])
  })
  return(unlist(autocorrelation_list))
}

#function to run slicer algorithm: 
run_slicer <- function(cds, select_genes = F, start = NULL, min_branch_len = 10){
  traj <- t(convert2DRData(cds, norm_method = 'log'))
  if(select_genes)
    genes = select_genes(traj)
  else
    genes = 1:ncol(traj)
  k = select_k(traj[,genes], kmin=5)
  traj_lle = lle(traj[,genes], m=2, k)$Y
  traj_graph = conn_knn_graph(traj_lle,5)
  ends = find_extreme_cells(traj_graph, traj_lle)
  
  if(is.null(start))
    start <- ends[1]
    
  cells_ordered = cell_order(traj_graph, start)
  branches = assign_branches(traj_graph,start, min_branch_len = min_branch_len)
  
  return(list(traj_lle = traj_lle, ends = ends, order_df = data.frame(cells_ordered = cells_ordered, branches = branches)))
}

#function to run slicer algorithm: 
run_dpt <- function(cds, branching = T, normalize = T, root = NULL, color_by = 'Hours'){
  if(normalize)
    data <- t(convert2DRData(cds[, ], norm_method = 'log'))
  else
    data <- t(as.matrix(log(exprs(cds)[fData(cds)$use_for_ordering, ] + 1)))
  ts <- Transitions(data)
  M <- dpt:::propagation_matrix(ts)

  if(is.null(root))
    pt <- dpt(ts, branching = branching)
  else
    pt <- dpt(ts, branching = branching, root = root)
  
  ev <- eigen(as.matrix(ts@transitions), TRUE)$vectors
  dm <- as.data.frame(ev[, -1])
  colnames(dm) <- paste0('DC', seq_len(ncol(dm)))
  
  if(!('Hours' %in% colnames(pData(cds))))
    pData(cds)$Hours <- pData(cds)[, color_by]
  p1 <- qplot(DC1, DC2, data = dm, colour = pData(cds)$Hours)
  #p1 <-plot_dpt_update(ts, pt, 1:2) # DPT and average path , 1:2
  dm <- DiffusionMap(data)
  DPT_res <- DPT(dm)
  dp_res <- list(dm = dm, pt = pt, ts = ts, M = M, ev = ev, p1 = p1, DPT = DPT_res)
  
  return(dp_res)
}

run_new_dpt <- function(cds, branching = T, normalize = T, root = NULL, color_by = 'Hours'){
  message('root should be the id to the cell not the cell name ....')
  
  if(benchmark_type == 'na_sim_data') {
    norm_data <- t(as.matrix(exprs(cds)[fData(cds)$use_for_ordering, ]))
  } 
  else {
    if(normalize)
      norm_data <- t(convert2DRData(cds[, ], norm_method = 'log'))
    else
      norm_data <- t(as.matrix(log(exprs(cds)[fData(cds)$use_for_ordering, ] + 1)))
  }
  
  dm <- DiffusionMap(norm_data)
  dpt <- DPT(dm)
  
  ts <- dm@transitions
  M <- destiny:::accumulated_transitions(dm)
  
  if(is.null(root)){
    
    # plot(dpt, root = 2, paths_to = c(1,3), col_by = 'branch', pch = 20)
    # plot(dpt, col_by = 'branch', divide = 3, dcs = c(-1,3,-2), pch = 20)
  }
  else{
    dm <- DiffusionMap(norm_data)
    dpt <- DPT(dm, tips = root)
    # plot(dpt, root = root, paths_to = c(1,3), col_by = 'branch', pch = 20)
    # plot(dpt, col_by = 'branch', divide = 3, dcs = c(-1,3,-2), pch = 20)
  }
  
  if('Hours' %in% colnames(pData(cds)))
    pData(cds)$Hours <- pData(cds)[, color_by]
  p1 <- qplot(DM$DC1, DM$DC2, colour = pData(cds)$Hours)
  
  branch <- dpt@branch
  row.names(branch) <- row.names(norm_data[!duplicated(norm_data), ])
  
  if(is.null(root))
    root <- which.min(pData(cds)$Pseudotime)
  pt <- dpt[root, ]
  dp_res <- list(dm = dm, pt = pt, ts = ts, M = M, ev = dm@eigenvectors, p1 = p1, branch = branch)
  
  return(dp_res)
}

plot_dpt_update <- function (ts, pt, branches, x = 1L, y = 2L, w_width = 0.1, path_col = "red", size = 1, 
          ..., dcs = eig_decomp(ts@transitions, max(x, y))$vectors[, 
                                                                   -1]) 
{
  stopifnot(all(c("Branch", "DPT", "DPT.1", "DPT.2") %in% names(pt)))
  stopifnot(is(ts, "Transitions"))
  stopifnot(length(x) == 1L, length(y) == 1L)
  stopifnot(is.integer(branches))
  idx <- as.integer(pt$Branch) %in% c(1L, 2L, branches + 2L)
  if (is.null(colnames(dcs))) {
    colnames(dcs) <- paste0("DC", seq_len(ncol(dcs)))
  }
  evs <- as.data.frame(as.matrix(dcs))[, c(x, y)]
  nms <- names(evs)
  DPT <- switch(branches[[1L]], pt$DPT, pt$DPT.1, pt$DPT.2)
  path <- average_path(DPT[idx], evs[idx, ], w_width)
  (ggplot(cbind(evs, DPT = DPT), aes_string(nms[[1L]], nms[[2L]], 
                                            colour = "DPT")) + geom_point(size = size) + scale_size(range = c(0.1, size)) + geom_path(data = path, 
                                                                                        colour = path_col) + monocle:::monocle_theme_opts() + theme(axis.text.x = element_text(angle = 30, hjust = 1)))
}

#run density clustering algorithm: 
#this script is an implementation of the density cluster method recently published in Science 
densityClustering <- function(cds) {
  dist_mat <- as.matrix(dist(t(exprs(absolute_cds)))) 
  locDen_vec <- locDen(t(exprs(absolute_cds)), dist_mat)
  sigma_vec <- dist2HighDensity(dist_mat, locDen_vec)
  DecisionGraph(locDen_vec, sigma_vec)
  cluster_vec <- assignCluster(dist_mat, locDen_vec, sigma_vec)
  
  p1 <- 
    qplot(x = 1:199, pData(absolute_cds)$State, color = c('red', 'blue', 'green')[cluster_vec]) +
    xlab('sample') + ylab('State calculated from monocle') + 
    scale_color_discrete(name = 'clusters based on densityClustering', label = sort(unique(cluster_vec)))
  
  p2 <- 
    qplot(x = 1:199, pData(absolute_cds)$State, color = c('red', 'blue', 'green')[cluster_vec], shape = pData(absolute_cds)$Time) +
    xlab('sample') + ylab('State calculated from monocle') + 
    scale_color_discrete(name = 'clusters based on densityClustering', label = sort(unique(cluster_vec))) + 
    scale_shape_discrete(name = 'Time')
}

calDensity <- function(dist_vec, dist_c) {
  sum(dist_vec < dist_c)
}

#calculat the local density 
locDen <- function(exprs_mat, dist_mat, method = 'euclidean', dist_c_method = mean) {
  if(ncol(exprs_mat) < nrow(exprs_mat)) exprs_mat <- t(exprs_mat) #rows are sample, columns are features 
  
  locDen <- apply(dist_mat, 1, function(x) calDensity(x, dist_c_method(dist_mat)))
  
  return(locDen)
}

#calculate the dist for each point to any other points with higher density
dist2HighDensity <- function(dist_mat, locDen_vec) {
  sampleName <- rownames(dist_mat)
  if(!identical(sampleName, colnames(dist_mat))) dist_mat <- dist_mat[, sampleName]
  if(!identical(sampleName, names(locDen_vec))) locDen <- locDen_vec[sampleName]
  
  nSample <- length(locDen_vec)
  sigma_vec <- rep(0, nSample)
  for(i in 1:nSample) {
    HighDenIndex <- which(locDen_vec > locDen_vec[i]) #any other points with higher density
    
    if(length(HighDenIndex))
      sigma_vec[i] <- min(dist_mat[i, HighDenIndex])
    else #point with highest density (global center)
      sigma_vec[i] <- max(dist_mat[i, ])
  }
  
  names(sigma_vec) <- sampleName
  return(sigma_vec)
}

#plot the decision graph to find the centers of clusters
DecisionGraph <- function(locDen_vec, sigma_vec) {
  sampleName <- names(locDen_vec)
  if(identical(sampleName, names(sigma_vec))) sigma_vec <- sigma_vec[sampleName]
  # qplot(locDen_vec, sigma_vec, label = sampleName, size = 2, geom = c("point", "text")) + geom_text(size = 0.2)
  
  plot_df <- data.frame(locDensity = locDen_vec, sigma = sigma_vec)
  ggplot(plot_df, aes(x = locDen_vec, y = sigma_vec, label = sampleName)) + geom_text(size = .5) + geom_point()
}

#assign each sample to different clusters
assignCluster <- function(dist_mat, locDen_vec, sigma_vec, clustCent_criteria = NULL) {
  sampleName <- rownames(dist_mat)
  if(!identical(sampleName, colnames(dist_mat))) dist_mat <- dist_mat[, sampleName]
  if(!identical(sampleName, names(locDen_vec))) locDen_vec <- locDen_vec[sampleName]
  if(!identical(sampleName, names(sigma_vec))) sigma_vec <- sigma_vec[sampleName]
  
  if(!is.null(clustCent_criteria))
    clusterCenter <- which(locDen_vec > clustCent_criteria[1] & sigma_vec > clustCent_criteria[2])
  else 
    clusterCenter <- which(locDen_vec > 1/2 * max(locDen_vec) & sigma_vec > 1/2 * max(sigma_vec))
  
  cluster_vec <- rep(0, length(clusterCenter))
  # cluster <- list(center = names(clusterCenter), clusterNum = length(clusterCenter), clusters = list())
  # length(cluster$clusters) <- cluster$clusterNum
  nSample <- nrow(dist_mat)
  for(i in 1:nSample) {
    cluster_id <- which.min(dist_mat[i, clusterCenter]) #find clustering center with smallest dist
    # cluster$clusters[[cluster_id]] <- c(cluster$clusters[[cluster_id]], sampleName[i])
    
    cluster_vec[i] <- cluster_id
  }
  
  names(cluster_vec) <- sampleName
  return(cluster_vec)
}

#clustCent_criteria <- c(40, 10000)
densityCluster <- function(cds) {
  data_matrix <- t(convert2DRData(cds))
  Dist <- dist(data_matrix)
  Clust <- densityClust(Dist, gaussian=TRUE)
  plot(Clust) # Inspect clustering attributes to define thresholds
  
  Clust <- findClusters(Clust, rho=2, delta=2)
  # split(iris[,5], irisClust$clusters)
}
#run densityClustering for the tSNE result and compare with the DDRTree clustering? 


#calculate the Dunn index: 

#calculate the correct assignment of the trunks: 
cal_percentage_correct_trunk_assignment <- function(assignment_to_compare, correct_assignment){
  sum(assignment_to_compare ==  correct_assignment) / length(correct_assignment)
}

#calculate the number of cells need to traverse on the graph to reach the true branch time point: 
cal_distance_branch_time_point <- function(mst, correct_branch_point){
  from <- degree(mst) == 3
  shortest_path <- shortest_paths(mst, from, to = correct_branch_point, weights = NULL, predecessors = FALSE, inbound.edges = FALSE)
  return(length(shortest_path))
}

#run andrew's pipeline code to prepare a cds based on files generated from cellranger: 
tenx_to_cds = function(pipeline_dirs, genome="hg19") {
  # Takes a list of 10X pipeline output directories and generates a cellDataSet containing all cells in these experiments
  #
  # Args:
  #	pipeline_dirs: Directory name or list of directory names of the top level 10X output directory for an experiment(s)
  # 	genome: String with genome name specified for 10X run (such as hg19)
  #
  # Returns:
  # 	A cellDataSet object containing data from all experiments.
  
  # Outer scope variables for collecting data
  expression_matrices = list()
  metadata_dfs = list()
  gene_table = NULL # will keep first gene table so can match ordering for each dataset
  
  lapply(pipeline_dirs, function(pipeline_dir) {
    # Check initial user input
    if( ! file.exists(pipeline_dir) ) { stop(paste("Specified 10X output directory does not exist:", pipeline_dir)) }
    
    # Construct paths for pipeline output files and check that they exist
    base_path = file.path(pipeline_dir, "outs", "filtered_gene_bc_matrices", genome)
    
    if( ! file.exists(base_path) ) { stop(paste("Specified genome does not appear in 10X output:", base_path)) }
    
    matrix_path = file.path(base_path, "matrix.mtx")
    genes_path = file.path(base_path, "genes.tsv")
    barcodes_path = file.path(base_path, "barcodes.tsv")
    analysis_path = file.path(pipeline_dir, "outs", "analysis")
    
    if( ! file.exists(matrix_path) ) { stop(paste("Expression matrix not found in 10X output:", matrix_path)) }
    if( ! file.exists(genes_path) ) { stop(paste("Genes file not found in 10X output:", genes_path)) }
    if( ! file.exists(barcodes_path) ) { stop(paste("Barcodes file not found in 10X output:", barcodes_path)) }
    if( ! file.exists(analysis_path) ) { stop(paste("Analysis  path not found in 10X output:", analysis_path)) }
    
    # All files exist, read them in
    matrix = Matrix::readMM(matrix_path)
    barcodes = read.table(barcodes_path, header=F, as.is=T)[,1]
    current_gene_table = read.table(genes_path, header=F, as.is=T) ## saves for later
    tsne = read.delim(file.path(analysis_path, "tsne", "projection.csv"), sep=',')
    
    if ( is.null(gene_table) ) { gene_table <<- current_gene_table } ## store the first gene table so can match ordering between all experiments
    
    genes = current_gene_table[, 1]
    
    # Add gene and sample names to expression matrix (adding dataset post-fix in case barcodes appear in multiple samples)
    row.names(matrix) = genes
    colnames(matrix) = paste(barcodes, "_", length(expression_matrices), sep="") ## adds dataset post-fix
    matrix = matrix[gene_table[, 1], ] ## ensures order of genes matches between experiments
    
    # Construct metadata table that includes directory samples come from and other stats
    total_umis = colSums(matrix)
    sample = basename(pipeline_dir)
    
    metadata_df = data.frame(
      cell = colnames(matrix),
      total_umis = total_umis,
      sample = sample,
      TSNE.1 = tsne$TSNE.1,
      TSNE.2 = tsne$TSNE.2
    )
    
    # Add both matrices to the running list
    expression_matrices[[length(expression_matrices) + 1]] <<- matrix
    metadata_dfs[[length(metadata_dfs) + 1]] <<- metadata_df
  })
  
  # Now combine all the dataframes into one and make CDS
  combined_expression_matrix <- do.call(cBind, expression_matrices)
  row.names(combined_expression_matrix) <- gene_table[, 1]
  
  combined_metadata_df <- do.call(rbind, metadata_dfs)
  row.names(combined_metadata_df) = combined_metadata_df$cell
  
  colnames(gene_table) = c("id", "gene_short_name")
  row.names(gene_table) = gene_table$id
  
  pd = new("AnnotatedDataFrame", data = combined_metadata_df)
  fd = new("AnnotatedDataFrame", data = gene_table)
  cds = newCellDataSet(combined_expression_matrix, 
                       phenoData=pd, 
                       featureData=fd,
                       expressionFamily=negbinomial(), 
                       lowerDetectionLimit=1)
  return(cds)
}

#create CDS based on the directory we provide: 
createCellRangeCDS <- function(matrix_path, barcodes_path, genes_path){
  mtx_matrix = Matrix::readMM(matrix_path)
  barcodes = read.table(barcodes_path, header=F, as.is=T)[,1]
  gene_table = read.table(genes_path, header=F, as.is=T) ## saves for later
  row.names(gene_table) <- gene_table$V1
  genes = gene_table[, 1]
  
  # Add gene and sample names to expression matrix (adding dataset post-fix in case barcodes appear in multiple samples)
  row.names(mtx_matrix) = genes
  colnames(mtx_matrix) = paste(barcodes, "_", length(mtx_matrix), sep="") ## adds dataset post-fix
  mtx_matrix = mtx_matrix[gene_table[, 1], ] ## ensures order of genes matches between experiments
  
  metadata_df = data.frame(
    cell = colnames(mtx_matrix),
    total_umis = apply(mtx_matrix, 2, sum),
    sample = 'sample'
  )
  
  pd = new("AnnotatedDataFrame", data = metadata_df)
  fd = new("AnnotatedDataFrame", data = gene_table)
  cds = newCellDataSet(as.matrix(mtx_matrix), 
                            phenoData=pd, 
                            featureData=fd,
                            expressionFamily=negbinomial(), 
                            lowerDetectionLimit=1)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
}

#functions to calculate the statistics for matching of two clusters: 
#require 

calClusteringMetrics <- function(cl1, cl2){
  tab_1 <- table(cl1, cl2)
  randi <- flexclust::randIndex(tab_1)
  
  vi <- vi.dist(cl1, cl2, parts = FALSE, base = 2)
  
  #arandi
  arandi <- arandi(cl1, cl2, adjust = TRUE)
  
  #add all result into a data frame
  randIndex_df <- data.frame(randIndex = c(randi, vi, arandi), 
                             Type = c("rand index", "variation of information", "adjusted rand index"))
  
  return(randIndex_df)
}

#run the function: 
benchmark_cluster_peformance <- function(c1, c2){
  tab_1 <- table(c1, c2)
  graceCluster_cellType <- flexclust::randIndex(tab_1)
  
  vi.dist(pData(subset_pbmc_cds)$grace_cluster, pData(subset_pbmc_cds)$cell_type, parts = FALSE, base = 2)
  
  arandi(pData(subset_pbmc_cds)$grace_cluster, pData(subset_pbmc_cds)$cell_type, adjust = TRUE)
}

#function to find the distance, path segment as well as branch point 
# - give two states or two cells, obtain the short path between these two cells states
# - calculate the distance between them in terms of graph
# - return the branch point traversed by the short path

traverseTree <- function(g, initial_vertex, terminal_vertex){ 
  distance <- shortest.paths(g, v=initial_vertex, to=terminal_vertex)
  branchPoints <- which(degree(g) == 3)
  path <- shortest_paths(g, from = initial_vertex, terminal_vertex)
  
  return(list(shortest_path = path$vpath, distance = distance, branch_points = intersect(branchPoints, unlist(path$vpath))))
}

#
# buildBranchCellDataSet <- function(cds,
#                                    progenitor_method = c('sequential_split', 'duplicate'), 
#                                    branch_states = NULL, 
#                                    branch_point = 1,
#                                    branch_labels = NULL, 
#                                    stretch = TRUE)
# {
#   # TODO: check that branches are on the same paths
#   
#   if(is.null(pData(cds)$State) | is.null(pData(cds)$Pseudotime)) 
#     stop('Please first order the cells in pseudotime using orderCells()')
#   if(is.null(branch_point) & is.null(branch_states)) 
#     stop('Please either specify the branch_point or branch_states to select subset of cells')
#   #if(ncol(cds@reducedDimS) != ncol(cds))
#   #  stop('You probably used clusterCells function which should be used together with buildBranchCellDataSet, try re-run reduceDimension without clustering cells again')
#   
#   if (!is.null(branch_labels) & !is.null(branch_states)) {
#     if(length(branch_labels) != length(branch_states))
#       stop("length of branch_labels doesn't match with that of branch_states")
#     branch_map <- setNames(branch_labels, as.character(branch_states))
#   }
#   
#   if(cds@dim_reduce_type == "DDRTree") {
#     pr_graph_cell_proj_mst <- minSpanningTree(cds)
#   }
#   else {
#     pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
#   }
#   
#   root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
#   root_state <- pData(cds)[root_cell,]$State
#   #root_state <- V(pr_graph_cell_proj_mst)[root_cell,]$State
#   
#   pr_graph_root <- subset(pData(cds), State == root_state)
#   
#   if (cds@dim_reduce_type == "DDRTree"){
#     closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
#     root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root),]
#   }else{
#     root_cell_point_in_Y <- row.names(pr_graph_root)
#   }
#   
#   root_cell <- names(which(degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y, mode = "all")==1, useNames = T))[1]
#   
#   paths_to_root <- list()
#   if (is.null(branch_states) == FALSE){
#     
#     # If the user didn't specify a branch point,
#     # let's walk back from the branch states
#     for (leaf_state in branch_states){
#       
#       curr_cell <- subset(pData(cds), State == leaf_state)
#       #Get all the nearest cells in Y for curr_cells:
#       
#       if (cds@dim_reduce_type == "DDRTree"){
#         closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
#         curr_cell_point_in_Y <- closest_vertex[row.names(curr_cell),] 
#       }else{
#         curr_cell_point_in_Y <- row.names(curr_cell)
#       }
#       
#       # Narrow down to a single tip cell in Y:
#       curr_cell <- names(which(degree(pr_graph_cell_proj_mst, v = curr_cell_point_in_Y, mode = "all")==1, useNames = T))[1]
#       
#       path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst,curr_cell, root_cell)
#       path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
#       
#       if (cds@dim_reduce_type == "DDRTree"){
#         closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
#         ancestor_cells_for_branch <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_ancestor)]
#       }else if (cds@dim_reduce_type == "ICA"){
#         ancestor_cells_for_branch <- path_to_ancestor
#       }
#       ancestor_cells_for_branch <- intersect(ancestor_cells_for_branch, colnames(cds))
#       paths_to_root[[as.character(leaf_state)]] <- ancestor_cells_for_branch
#     }
#   }else{
#     if(cds@dim_reduce_type == "DDRTree")
#       pr_graph_cell_proj_mst <- minSpanningTree(cds)
#     else
#       pr_graph_cell_proj_mst <- cds@auxOrderingData$ICA$cell_ordering_tree
#     
#     mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
#     branch_cell <- mst_branch_nodes[branch_point]
#     mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_cell]
#     
#     path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_cell, root_cell)
#     path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
#     
#     #post_branch_cells <- c()
#     for (backbone_nei in V(pr_graph_cell_proj_mst)[suppressWarnings(nei(branch_cell))]$name){
#       descendents <- bfs(mst_no_branch_point, V(mst_no_branch_point)[backbone_nei], unreachable=FALSE)
#       descendents <- descendents$order[!is.na(descendents$order)]
#       descendents <- V(mst_no_branch_point)[descendents]$name
#       if (root_cell %in% descendents == FALSE){
#         path_to_root <- unique(c(path_to_ancestor, branch_cell, descendents))
#         
#         if (cds@dim_reduce_type == "DDRTree"){
#           closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
#           path_to_root <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_root)]
#         }else{
#           path_to_root <- path_to_root
#         }
#         
#         closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
#         #branch_state <- unique(pData(cds)[backbone_nei, ]$State)[1]
#         
#         path_to_root <- intersect(path_to_root, colnames(cds))
#         paths_to_root[[backbone_nei]] <- path_to_root
#         #post_branch_cells <- c(post_branch_cells, backbone_nei)
#       }
#     }
#   }
#   all_cells_in_subset <- c()
#   
#   if (is.null(branch_labels) == FALSE){
#     if (length(branch_labels) != 2)
#       stop("Error: branch_labels must have exactly two entries")
#     names(paths_to_root) <- branch_labels
#   }
#   
#   for (path_to_ancestor in paths_to_root){
#     if (length(path_to_ancestor) == 0){
#       stop("Error: common ancestors between selected State values on path to root State")
#     }
#     all_cells_in_subset <- c(all_cells_in_subset, path_to_ancestor)
#   }
#   all_cells_in_subset <- unique(all_cells_in_subset)
#   
#   common_ancestor_cells <- intersect(paths_to_root[[1]], paths_to_root[[2]])
#   # if (length(paths_to_root) > 2){
#   #   for (i in seq(3,length(paths_to_root))){
#   #     common_ancestor_cells <- intersect(intersect(paths_to_root[i], paths_to_root[i-1]), common_ancestor_cells)
#   #   }
#   # }
#   
#   #when n-center used, this creates problems
#   cds <- cds[, row.names(pData(cds[,all_cells_in_subset]))] #or just union(ancestor_cells, branch_cells)
#   
#   #State <- pData(cds)$State 
#   Pseudotime <- pData(cds)$Pseudotime 
#   
#   pData <- pData(cds)
#   
#   if(stretch) {
#     max_pseudotime <- -1
#     for (path_to_ancestor in paths_to_root){
#       max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime)  
#       if (max_pseudotime < max_pseudotime_on_path){
#         max_pseudotime <- max_pseudotime_on_path
#       }
#     }
#     
#     branch_pseudotime <- max(pData[common_ancestor_cells,]$Pseudotime)
#     #ancestor_scaling_factor <- branch_pseudotime / max_pseudotime
#     
#     for (path_to_ancestor in paths_to_root){
#       max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime) 
#       path_scaling_factor <-(max_pseudotime - branch_pseudotime) / (max_pseudotime_on_path - branch_pseudotime)
#       if (is.finite(path_scaling_factor)){
#         branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
#         pData[branch_cells,]$Pseudotime <- ((pData[branch_cells,]$Pseudotime - branch_pseudotime) * path_scaling_factor + branch_pseudotime)
#       }
#     }
#     #pData[common_ancestor_cells,]$Pseudotime <- pData[common_ancestor_cells,]$Pseudotime / max_pseudotime
#     
#     pData$Pseudotime <- 100 * pData$Pseudotime / max_pseudotime
#   }
#   pData$original_cell_id <- row.names(pData)
#   
#   pData$original_cell_id <- row.names(pData)
#   
#   if(length(paths_to_root) != 2)
#     stop('more than 2 branch states are used!')
#   
#   pData[common_ancestor_cells, "Branch"] <- names(paths_to_root)[1] #set progenitors to the branch 1
#   
#   progenitor_pseudotime_order <- order(pData[common_ancestor_cells, 'Pseudotime'])
#   
#   if (progenitor_method == 'duplicate') {
#     ancestor_exprs <- exprs(cds)[,common_ancestor_cells]
#     expr_blocks <- list()
#     
#     # Duplicate the expression data
#     for (i in 1:length(paths_to_root)) { #duplicate progenitors for multiple branches
#       if (nrow(ancestor_exprs) == 1)
#         exprs_data <- t(as.matrix(ancestor_exprs))
#       else exprs_data <- ancestor_exprs
#       
#       colnames(exprs_data) <- paste('duplicate', i, 1:length(common_ancestor_cells), sep = '_')
#       expr_lineage_data <- exprs(cds)[,setdiff(paths_to_root[[i]], common_ancestor_cells)]
#       exprs_data <- cbind(exprs_data, expr_lineage_data)
#       expr_blocks[[i]] <- exprs_data
#     }
#     
#     # Make a bunch of copies of the pData entries from the common ancestors
#     ancestor_pData_block <- pData[common_ancestor_cells,]
#     
#     pData_blocks <- list()
#     
#     weight_vec <- c()
#     for (i in 1:length(paths_to_root)) {
#       weight_vec <- c(weight_vec, rep(1, length(common_ancestor_cells)))
#       weight_vec_block <- rep(1, length(common_ancestor_cells))
#       
#       #pData <- rbind(pData, pData[common_ancestor_cells, ])
#       new_pData_block <- ancestor_pData_block
#       # new_pData_block$Lineage <- lineage_states[i]
#       # new_pData_block$State <- lineage_states[i]
#       
#       row.names(new_pData_block) <- paste('duplicate', i, 1:length(common_ancestor_cells), sep = '_')
#       
#       pData_lineage_cells <- pData[setdiff(paths_to_root[[i]], common_ancestor_cells),]
#       # pData_lineage_cells$Lineage <- lineage_states[i]
#       # pData_lineage_cells$State <- lineage_states[i]
#       
#       weight_vec_block <- c(weight_vec_block, rep(1, nrow(pData_lineage_cells)))
#       
#       weight_vec <- c(weight_vec, weight_vec_block)
#       
#       new_pData_block <- rbind(new_pData_block, pData_lineage_cells)
#       pData_blocks[[i]] <- new_pData_block
#     }
#     pData <- do.call(rbind, pData_blocks)
#     exprs_data <- do.call(cbind, expr_blocks)
#   }
#   else if(progenitor_method == 'sequential_split') {
#     pData$Branch <- names(paths_to_root)[1]
#     
#     branchA <- progenitor_pseudotime_order[seq(1, length(common_ancestor_cells), by = 2)]
#     pData[common_ancestor_cells[branchA], 'Branch'] <- names(paths_to_root)[1]
#     branchB <- progenitor_pseudotime_order[seq(2, length(common_ancestor_cells), by = 2)]
#     pData[common_ancestor_cells[branchB], 'Branch'] <- names(paths_to_root)[2]   
#     
#     # Duplicate the root cell to make sure both regression start at pseudotime zero:
#     zero_pseudotime_root_cell <- common_ancestor_cells[progenitor_pseudotime_order[1]]
#     exprs_data <- cBind(exprs(cds), 'duplicate_root' = exprs(cds)[, zero_pseudotime_root_cell])
#     pData <- rbind(pData, pData[zero_pseudotime_root_cell, ])
#     row.names(pData)[nrow(pData)] <- 'duplicate_root'
#     pData[nrow(pData), 'Branch'] <- names(paths_to_root)[2]
#     
#     weight_vec <- rep(1, nrow(pData))
#     
#     for (i in 1:length(paths_to_root)){
#       path_to_ancestor <- paths_to_root[[i]]
#       branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
#       pData[branch_cells,]$Branch <- names(paths_to_root)[i]
#     }
#   }
#   
#   pData$Branch <- as.factor(pData$Branch)
#   
#   pData$State <- factor(pData$State)
#   Size_Factor <- pData$Size_Factor
#   
#   fData <- fData(cds)
#   
#   colnames(exprs_data) <- row.names(pData) #check this 
#   cds_subset <- newCellDataSet(as.matrix(exprs_data),
#                                phenoData = new("AnnotatedDataFrame", data = pData),
#                                featureData = new("AnnotatedDataFrame", data = fData),
#                                expressionFamily=cds@expressionFamily,
#                                lowerDetectionLimit=cds@lowerDetectionLimit)
#   pData(cds_subset)$State <- as.factor(pData(cds_subset)$State)
#   pData(cds_subset)$Size_Factor <- Size_Factor
#   
#   cds_subset@dispFitInfo <- cds@dispFitInfo
#   
#   return (cds_subset)
# }

selectOrderingGenes2 <- function(cds, dim = c(2, 3), top_gene_num = 1000, method = c('PCA', 'FSTree'), initial_method = PCA, verbose = F, ...) {
  cds <- detectGenes(cds)
  cds <- estimateSizeFactors(cds)
  exprs_filtered <- t(t(exprs(cds)/pData(cds)$Size_Factor))
  nz_genes <- which(as.matrix(exprs_filtered) != 0)
  exprs_filtered[nz_genes] <- log(exprs_filtered[nz_genes] + 1)
  
  # # Calculate the variance across genes without converting to a dense
  # # matrix:
  # expression_means <- Matrix::rowMeans(exprs_filtered)
  # expression_vars <- Matrix::rowMeans((exprs_filtered - expression_means)^2)
  # 
  # # Filter out genes that are constant across all cells:
  # genes_to_keep <- expression_vars > 0
  # exprs_filtered <- exprs_filtered[genes_to_keep,]
  # expression_means <- expression_means[genes_to_keep]
  # expression_vars <- expression_vars[genes_to_keep]
  # 
  # # Here's how to take the top PCA loading genes, but using
  # # sparseMatrix operations the whole time, using irlba.
  # irlba_pca_res <- irlba(t(exprs_filtered), center=expression_means, scale=sqrt(expression_vars), right_only=TRUE)$v # nu=0,
  # row.names(irlba_pca_res) <- row.names(exprs_filtered)
  # 
  # # Select the top_gene_num from each selected pca dimension for performing the DEG test: 
  # 
  # # Add function to perform GO enrichment analysis to determine the best components for running fstree 
  fData(cds)$use_for_ordering[fData(cds)$num_cells_expressed > 5] <- T
  data <- convert2DRData(cds[fData(cds)$num_cells_expressed > 5, ], norm_method = 'log')
  
  pca_res <- prcomp(t(data), center = T, scale = T) #data: row is the gene, column is the cell 
  # print(plot(pca_res)) #plot the pca results
  # print(biplot(pca_res))
  sort(abs(pca_res$rotation[, 1]), decreasing = T)
  
  top_gene_num <- min(top_gene_num, nrow(data))
  gene_vec <- c()
  for(i in dim) {
    tmp <- names(sort(abs(pca_res$rotation[, i]), decreasing = T))[1:top_gene_num]
    gene_vec <- c(gene_vec, tmp)
  }
  # gene_vec <- c()
  # for(i in dim) {
  #   tmp <-  names(sort(abs(irlba_pca_res[, i]), decreasing = T))[1:top_gene_num]
  #   gene_vec <- c(gene_vec, tmp)
  # }
  
  # fstree_res <- fstree(log2(exprs_filtered[gene_vec, ] + 1), initial_method = PCA, maxIter = maxIter, eps = 1e-5, dim = d, cell_time = NULL, 
  #                              C_constant = 0.1,  lambda = 1,  gamma = 1, sigma = 1 , num_landmarks = NULL, verbose = T)
  
  if(ncol(exprs_filtered) > 10000)
    num_landmarks <- 1000
  else
    num_landmarks <- NULL
  
  if(method == 'FSTree'){
    res <- fstree(as.matrix(exprs_filtered[gene_vec, ]), d = 2, num_landmarks = num_landmarks, maxIter = 20, verbose = verbose, initial_method = initial_method)
    res$gene_vec <- gene_vec
  }
  else{
    res <- list(gene_vec = gene_vec, irlba_pca_res = NULL) #irlba_pca_res
  }
  
  return(res)
}

selectOrderingGenes <- function(cds, dim = c(2, 3), top_gene_num = 1000, method = c('PCA', 'FSTree'), verbose = F, ...) {
  exprs_filtered <- t(t(exprs(cds)/pData(cds)$Size_Factor))
  nz_genes <- which(exprs_filtered != 0)
  exprs_filtered[nz_genes] <- log(exprs_filtered[nz_genes] + 1)
  
  # Calculate the variance across genes without converting to a dense
  # matrix:
  expression_means <- Matrix::rowMeans(exprs_filtered)
  expression_vars <- Matrix::rowMeans((exprs_filtered - expression_means)^2)
  
  # Filter out genes that are constant across all cells:
  genes_to_keep <- expression_vars > 0
  exprs_filtered <- exprs_filtered[genes_to_keep,]
  expression_means <- expression_means[genes_to_keep]
  expression_vars <- expression_vars[genes_to_keep]
  
  # Here's how to take the top PCA loading genes, but using
  # sparseMatrix operations the whole time, using irlba.
  irlba_pca_res <- irlba(t(exprs_filtered), nu=0, center=expression_means, scale=sqrt(expression_vars), right_only=TRUE)$v
  row.names(irlba_pca_res) <- row.names(exprs_filtered)
  
  # Select the top_gene_num from each selected pca dimension for performing the DEG test: 
  
  # Add function to perform GO enrichment analysis to determine the best components for running fstree 
  
  gene_vec <- c()
  for(i in dim) {
    tmp <-  names(sort(abs(irlba_pca_res[, i]), decreasing = T))[1:top_gene_num]
    gene_vec <- c(gene_vec, tmp)
  }
  
  d = 2;
  lambda = 1;
  C = 0.1;
  maxIter <- 100
  
  # fstree_res <- fstree(log2(exprs_filtered[gene_vec, ] + 1), initial_method = PCA, maxIter = maxIter, eps = 1e-5, dim = d, cell_time = NULL, 
  #                              C_constant = 0.1,  lambda = 1,  gamma = 1, sigma = 1 , num_landmarks = NULL, verbose = T)
  
  if(ncol(exprs_filtered) > 10000)
    num_landmarks <- 1000
  else
    num_landmarks <- NULL
  
  if(method == 'FSTree'){
    res <- fstree(log2(exprs_filtered[gene_vec, ] + 1), num_landmarks = num_landmarks, verbose = verbose)
  }
  else{
    res <- list(gene_vec = gene_vec, irlba_pca_res = irlba_pca_res) #
  }
  
  return(res)
}

#create a new cds based on two states or cells: 
SubSet_cds <- function(cds, cells){
  cells <- unique(cells)
  if(ncol(reducedDimK(cds)) != ncol(cds))
    stop("SubSet_cds doesn't support cds with ncenter run for now. You can try to subset the data and do the construction of trajectory on the subset cds")
    
  exprs_mat <- as(as.matrix(exprs(cds)[, cells]), "sparseMatrix")
  cds_subset <- newCellDataSet(exprs_mat, 
                                 phenoData = new("AnnotatedDataFrame", data = pData(cds)[colnames(exprs_mat), ]), 
                                 featureData = new("AnnotatedDataFrame", data = fData(cds)), 
                                 expressionFamily=negbinomial.size(), 
                                 lowerDetectionLimit=1)
  sizeFactors(cds_subset) <- sizeFactors(cds[, cells])
  cds_subset@dispFitInfo <- cds@dispFitInfo
  # if(nrow(cds_subset@reducedDimS) > 0 ){
  #   cds_subset@reducedDimW <- cds_subset@reducedDimW[, cells]      
  # }
  cds_subset@reducedDimW <- cds@reducedDimW   
  cds_subset@reducedDimS <- cds@reducedDimS[, cells]
  cds_subset@reducedDimK <- cds@reducedDimK[, cells]
  # cds_subset@minSpanningTree <- subgraph(cds@minSpanningTree, cells)
  cds_subset@cellPairwiseDistances <- cds@cellPairwiseDistances[cells, cells]
  
  adjusted_K <- Matrix::t(reducedDimK(cds_subset))
  dp <- as.matrix(dist(adjusted_K))
  cellPairwiseDistances(cds_subset) <- dp
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  minSpanningTree(cds_subset) <- dp_mst
  cds_subset@dim_reduce_type <- "DDRTree"
  cds_subset <- findNearestPointOnMST(cds_subset)
  
  cds_subset <- orderCells(cds_subset)
  #determine the root state: 
  root_cell <- row.names(subset(pData(cds[, cells]), Pseudotime == min(Pseudotime)))
  root_state <- as.numeric(pData(cds_subset[, root_cell])$State)
  cds_subset <- orderCells(cds_subset, root_state = root_state)

  return(cds_subset)
} 
  
# if(ncol(ddrtree_res$Y) == ncol(cds))
#   colnames(ddrtree_res$Y) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
# else
#   colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")

# colnames(ddrtree_res$Z) <- colnames(FM)
# reducedDimW(cds) <- ddrtree_res$W       
# reducedDimS(cds) <- ddrtree_res$Z
# reducedDimK(cds) <- ddrtree_res$Y
# cds@auxOrderingData[["DDRTree"]]$objective_vals <- ddrtree_res$objective_vals

# adjusted_K <- Matrix::t(reducedDimK(cds))
# dp <- as.matrix(dist(adjusted_K))
# cellPairwiseDistances(cds) <- dp
# gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
# dp_mst <- minimum.spanning.tree(gp)
# minSpanningTree(cds) <- dp_mst
# cds@dim_reduce_type <- "DDRTree"
# cds <- findNearestPointOnMST(cds)

traverseTreeCDS <- function(cds, starting_cell, end_cells){
  subset_cell <- c()
  dp_mst <- cds@minSpanningTree
  
  for(end_cell in end_cells) {
    traverse_res <- traverseTree(dp_mst, starting_cell, end_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])
    
    subset_cell <- c(subset_cell, path_cells)
  }
  
  # cds_subset <- subset_cds[, subset_cell]
  cds_subset <- SubSet_cds(cds, subset_cell)
# 
#   # subset_dp_mst <- subgraph(cds_subset@minSpanningTree, subset_cell)
#   Pseudotime <- shortest.paths(dp_mst, v=starting_cell, to=subset_cell)
#   pData(cds_subset)$Pseudotime <- as.numeric(Pseudotime)
  
  root_state <- pData(cds_subset[, starting_cell])[, 'State']
  cds_subset <- orderCells(cds_subset, root_state = as.numeric(root_state))
  
  return(cds_subset)
}

#extended BEAM test: 
multiple_BEAM <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch", 
         reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
         branch_states = NULL,
         branch_points=c(1),
         relative_expr = TRUE, 
         # branch_labels = NULL, 
         verbose = FALSE,
         cores = 1, 
         ...) {
  #performing the tree substraction and do BEAM on each branch 
  branchTest_res_list <- list()
  for(branch_point_ind in 1:length(branch_points)) {
    branchTest_res <- branchTest(cds, fullModelFormulaStr = fullModelFormulaStr,
                                 reducedModelFormulaStr = reducedModelFormulaStr, 
                                 branch_states = branch_states, 
                                 branch_point=branch_points[branch_point_ind],
                                 relative_expr = relative_expr,
                                 cores = cores, 
                                 branch_labels = NULL, 
                                 verbose=verbose, 
                                 ...)
    branchTest_res_list[[i]] <- branchTest_res
  }
  
  names(branchTest_res_list) <- branch_points
  
  return(branchTest_res_list)
  # cmbn_df <- branchTest_res[, 1:4] 
  # 
  # #make a newCellDataSet object with the smoothed data? 
  # if(verbose)
  #   message('pass branchTest')
  # 
  # fd <- fData(cds)[row.names(cmbn_df),]
  # 
  # #combined dataframe: 
  # cmbn_df <- cbind(cmbn_df, fd)
  # 
  # if(verbose)
  #   message('return results')
  # 
  # return(cmbn_df)
}

#add the calABC after running calILRs for branchTimePoint detection: 
##already did: 

selectFeatureByDP <- function(cds, num_cells_expressed = 5, num_dim = 5, rho_threshold = NULL, delta_threshold = NULL, qval_threshold = 0.01){
  #1. determine how many pca dimension you want:
  cds <- detectGenes(cds)
  fData(cds)$use_for_ordering[fData(cds)$num_cells_expressed > num_cells_expressed] <- T

  if(is.null(num_dim)){
    lung_pc_variance <- plot_pc_variance_explained(cds, return_all = T)
    ratio_to_first_diff <- diff(lung_pc_variance$variance_explained) / diff(lung_pc_variance$variance_explained)[1]
    num_dim <- which(ratio_to_first_diff < 0.1) + 1
  }
  
  #2. run reduceDimension with tSNE as the reduction_method
  # absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
  cds <- reduceDimension(cds, return_all = F, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = num_dim,  verbose = T)

  #3. initial run of clusterCells_Density_Peak
  cds <- clusterCells_Density_Peak(cds, rho_threshold = rho_threshold, delta_threshold = delta_threshold, verbose = T)

  #perform DEG test across clusters: 
  cds@expressionFamily <- negbinomial.size()
  pData(cds)$Cluster <- factor(pData(cds)$Cluster)
  clustering_DEG_genes <- differentialGeneTest(cds, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
  clustering_DEG_genes_subset <- lung_clustering_DEG_genes[fData(cds)$num_cells_expressed > num_cells_expressed, ]
  
  #use all DEG gene from the clusters
  clustering_DEG_genes_subset <- clustering_DEG_genes_subset[order(clustering_DEG_genes_subset$qval), ]
  ordering_genes <- row.names(subset(clustering_DEG_genes, qval < qval_threshold))
  
  cds <- setOrderingFilter(cds, ordering_genes = lung_ordering_genes)
  cds <- reduceDimension(cds, norm_method = 'log', verbose = T)
  plot_cell_trajectory(cds, color_by = 'Time')
  
  return(list(new_cds = cds, ordering_genes = ordering_genes))
}

#compare with the algorithm from SLICER: 
slicer_feature_genes <- function(cds){ 
  data <- log(as.matrix(exprs(cds) + 1))
  # lle_res <- LLE(data)
  res <- select_genes(t(data))
  return(res)
}

#recreate a cds from the original cds: 
recreate_cds <- function(cds) {
  pd <- new("AnnotatedDataFrame", data = pData(cds))
  fd <- new("AnnotatedDataFrame", data = fData(cds))
  
  new_cds <- newCellDataSet(as(as.matrix(cds), "sparseMatrix"), 
                                 phenoData = pd, 
                                 featureData = fd, 
                                 expressionFamily=cds@expressionFamily, 
                                 lowerDetectionLimit=1)
  
  new_cds <- estimateSizeFactors(new_cds)
  new_cds <- estimateDispersions(new_cds)
  new_cds
}

select_loading_genes_in_gates <- function(cds, dim = c(1, 2), gene_num = 1000) {
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  fData(cds)$use_for_ordering <- T
  data <- convert2DRData(cds, norm_method = 'log')
  
  pca_res <- prcomp(t(data), center = T, scale = T) #data: row is the gene, column is the cell 
  print(plot(pca_res)) #plot the pca results
  # print(biplot(pca_res))
  sort(abs(pca_res$rotation[, 1]), decreasing = T)
  
  gene_vec <- c()
  for(i in dim) {
    tmp <- names(sort(abs(pca_res$rotation[, i]), decreasing = T))[1:gene_num]
    gene_vec <- c(gene_vec, tmp)
  }
  return(gene_vec)
}

#http://stats.stackexchange.com/questions/113485/weighted-principal-components-analysis
# Let XX be the data matrix with variables in columns and nn observations xixi in rows. If each observation has an associated weight wiwi, then it is indeed straightforward to incorporate these weights into PCA.
# 
# First, one needs to compute the weighted mean μ=1n∑wixiμ=1n∑wixi and subtract it from the data in order to center it.
# 
# Then we compute the weighted covariance matrix X⊤WX/(n−1)X⊤WX/(n−1), where W=diag(wi)W=diaga(wi) is the diagonal matrix of weights, and apply standard PCA to analyze it.
# 
# This is equivalent to multiplying each row of the appropriately centered data matrix by the corresponding wi‾‾√wi and proceeding with the standard PCA, because X⊤WX/(n−1)X⊤WX/(n−1) is the covariance matrix of W1/2XW1/2X.
# 
# Note that this is conceptually related to rescaling the variables (e.g. standardizing them), when one multiplies each column of the data matrix by a certain value.

#
normalized_PCA <- function(data, max_components = 3) {
  gene_weight_number_cells_expressed <- 1 / matrix(rep(apply(data, 1, function(x) sum(x > 1)), each = ncol(data)), ncol = ncol(data), byrow = T)
  normalized_Y_ori <- t(scale(t(data * gene_weight_number_cells_expressed), center = T, scale = F))
  W <- diag(gene_weight_number_cells_expressed[, 1])
  n <- nrow(data)
  weighted_cor_mat <- t(normalized_Y_ori) %*% W %*% normalized_Y_ori / (n - 1)
  
  res <- prcomp(t(weighted_cor_mat), center = T, scale = T)
  res$x[, 1:max_components]
}

order_by_original_states_custom_ordering_function <- function(data, root_state, cells_state_2, cells_state_3, initial_method = PCA, root_cell = NULL) {
  
  DDRTree_res <- DDRTree(data, dimensions = 2, verbose = F, initial_method = initial_method) #, maxIter = 5, sigma = 1e-2, lambda = 1, ncenter = 3, param.gamma = 10, tol = 1e-2
  
  #set the name for the ordering dataframe: 
  if(!is.null(root_cell)) {
    dpt_ordering <- custom_ordering(DDRTree_res, root_cell == which(colnames(data) ==  root_cell))
    cc_ordering <- dpt_ordering$cc_ordering
  }
  else {
    cc_ordering <- dpt_ordering$cc_ordering
    row.names(cc_ordering) <- colnames(data)[as.numeric(as.character(cc_ordering$sample_name))]
    
    #determine the mapping between original state 1/2/3 and new state 1/2/3:  
    overlap_state_1 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 1, ]), root_state))
    overlap_state_2 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 2, ]), root_state))
    overlap_state_3 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 3, ]), root_state))
    
    #find the state corresponding to the original root state
    overlap_vec <- c(overlap_state_1, overlap_state_2, overlap_state_3)
    max_ind <- which(overlap_vec == max(overlap_vec))[1]
    
    dp <- as.matrix(dist(t(DDRTree_res$Y)))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    root_cell <- intersect(which(degree(dp_mst) == 1), cc_ordering[cc_ordering$cell_state == max_ind, "sample_name"]) #find the tip cell of the root state cells
    dpt_ordering <- custom_ordering(DDRTree_res, root_cell = root_cell)
    
    #set the name for the ordering dataframe: 
    cc_ordering <- dpt_ordering$cc_ordering
    row.names(cc_ordering) <- colnames(data)[as.numeric(as.character(cc_ordering$sample_name))]
    
    overlap_state_2 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 2, ]), cells_state_2))
    overlap_state_3 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 3, ]), cells_state_2))
    
    #find the new state corresponding to the original state 2: 
    overlap_state_2_vec <- c(overlap_state_2, overlap_state_3)
    state_2_max_ind <- which(overlap_state_2_vec == max(overlap_state_2_vec)) #root state
    
    if(state_2_max_ind != 1) {
      State_vec <- cc_ordering$cell_state
      cc_ordering$cell_state[State_vec == 2] <- 3
      cc_ordering$cell_state[State_vec == 3] <- 2
    }
  }
  cc_ordering
}

extract_ddrtree_ordering_xj <- function(dp_mst, root_cell, verbose=T)
{
  
  curr_state <- 1
  
  res <- list(subtree = dp_mst, root = root_cell)
  
  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name
  
  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name
  
  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name
  
  mst_traversal <- graph.dfs(dp_mst, 
                             root=root_cell, 
                             neimode = "all", 
                             unreachable=FALSE, 
                             father=TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  
  for (i in 1:length(mst_traversal$order)){
    curr_node = mst_traversal$order[i]
    curr_node_name = V(dp_mst)[curr_node]$name
    
    if (is.na(mst_traversal$father[curr_node]) == FALSE){
      parent_node = mst_traversal$father[curr_node]
      parent_node_name = V(dp_mst)[parent_node]$name
      parent_node_pseudotime = pseudotimes[parent_node_name]
      parent_node_state = states[parent_node_name]
      curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
      if (degree(dp_mst, v=parent_node_name) > 2){
        curr_state <- curr_state + 1
      }
    }else{
      parent_node = NA
      parent_node_name = NA
      curr_node_pseudotime = 0
    }
    
    curr_node_state = curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  
  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = factor(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents)
  row.names(ordering_df) <- ordering_df$sample_name
  # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  return(ordering_df)
}

extract_ddrtree_ordering_xj <- function(dp_mst, dp = dp, root_cell, verbose=T) # dp,
{
  nei <- NULL
  type <- NULL
  pseudo_time <- NULL
  
  curr_state <- 1
  
  res <- list(subtree = dp_mst, root = root_cell)
  
  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name
  
  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name
  
  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name
  
  mst_traversal <- graph.dfs(dp_mst, 
                             root=root_cell, 
                             neimode = "all", 
                             unreachable=FALSE, 
                             father=TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  
  for (i in 1:length(mst_traversal$order)){
    curr_node = mst_traversal$order[i]
    curr_node_name = V(dp_mst)[curr_node]$name
    
    if (is.na(mst_traversal$father[curr_node]) == FALSE){
      parent_node = mst_traversal$father[curr_node]
      parent_node_name = V(dp_mst)[parent_node]$name
      parent_node_pseudotime = pseudotimes[parent_node_name]
      parent_node_state = states[parent_node_name]
      curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
      if (degree(dp_mst, v=parent_node_name) > 2){
        curr_state <- curr_state + 1
      }
    }else{
      parent_node = NA
      parent_node_name = NA
      curr_node_pseudotime = 0
    }
    
    curr_node_state = curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  
  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = factor(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents)
  row.names(ordering_df) <- ordering_df$sample_name
  # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  return(ordering_df)
}

select_top_graph_test_genes <- function(cds, graphTest_res, cor_graph_res = NULL, cth, top_num = 10){ #how to consider the cth when selecting genes? 
  #get the branchpoints as well as tip cells
  g <- cth@classificationTree
  branchpoints <- which(degree(g) == 3)
  tips <- which(degree(g) == 1)[-1]
  
  #calculate correlated genes: 
  cell_type_num <- length(V(cth@classificationTree)) - 1
  cell_type_vec <- V(cth@classificationTree)$name[-1]
  if(is.null(cor_graph_res)) {
    cor_graph_res <-  matrix(list(), nrow = cell_type_num, ncol = cell_type_num, dimnames = list(cell_type_vec, cell_type_vec )) #data.frame of data.frame
    for(child1 in V(cth@classificationTree)$name[-1]) {
      markers <- V(cth@classificationTree) [ child1 ]$markers[[1]]
      marker_ids <- row.names(subset(fData(cds), gene_short_name %in% markers))
      for(child2 in V(cth@classificationTree)$name[-1]) {
        cell_ids <- row.names(subset(pData(cds), Cell_Type == child2))
        
        mat <- exprs(cds)[, cell_ids]
        # Center each variable
        mat <- mat - rowMeans(mat);
        # Standardize each variable
        mat <- mat / sqrt(rowSums(mat^2));   
        # Calculate correlations
        cr <- tcrossprod(mat);
        
        cor_graph_res[[child1, child2]] <- cr[marker_ids, ]
        rm(cr) #remove cr
      }
    }
  }
  
  #get the top genes from the graphTest result
  valid_ids <- which(!sapply(graphTest_res$graph_deg_res, is.null))
  valid_deg_list <- graphTest_res$graph_deg_res[valid_ids]
  
  top_gene_lists <- lapply(valid_ids, function(x) {
    deg_res <- graphTest_res$graph_deg_res[[x]]
    branch_genes <- row.names(deg_res[order(deg_res$qval), ])#[1:(4 * max(top_num))]
    
    #select genes based on the cth 
    row_col_ind <- c((x %% 5), round(x / 5) + 1) #convert to index for the matrix
    #1. diag 
    if(row_col_ind[1] == row_col_ind[2]) {
      cur_node <- row.names(graphTest_res$graph_deg_res)[row_col_ind[1]]
      out_nodes <- neighborhood(g, nodes = cur_node, order = 1, mode = 'out')[[1]]$name[-1]
      
      corA <- cor_graph_res[[out_nodes[1], out_nodes[2]]]; corB <- cor_graph_res[[out_nodes[2], out_nodes[1]]]; 
      mutual_exclusive_genes <- apply(corA[, branch_genes], 2, function(x) all(x > 0, na.rm = T)) & apply(corB[, branch_genes], 2, function(x) all(x < 0, na.rm = T))
      mutual_exclusive_names <- names(mutual_exclusive_genes[mutual_exclusive_genes])

      mutual_exclusive_names[1:(4 * max(top_num))]
    }
    #2. off diag
    else {
      start_node <- row.names(graphTest_res$graph_deg_res)[row_col_ind[1]]
      end_node <- row.names(graphTest_res$graph_deg_res)[row_col_ind[2]]
      transition_path <- shortest_paths(g, start_node, end_node)$vpath[[1]]$name
      transition_cells <- transition_path[2:(length(transition_path) - 1)]

      transition_cor_list <- lapply(transition_cells, function(cell) {
        cor <- cor_graph_res[[cell, cell]]
        apply(cor[, branch_genes], 2, function(x) all(x > 0, na.rm = T))
      })
      
      transition_cor_df <- do.call(rbind, transition_cor_list)
      transition_genes <- apply(transition_cor_df, 2, function(x) all(x > 0, na.rm = T))
      transition_genes_names <- names(transition_genes[transition_genes])
      
      transition_genes_names[1:(4 * max(top_num))]
    }
    
    })
  
  #obtain unique top 250 deg genes from each test
  if(length(top_num) != length(valid_ids))
    top_num <- rep(max(top_num), length(top_gene_lists))
  
  top_gene_vec <- c()
  for(i in 1:length(top_gene_lists)) {
    top_gene_vec <-  c(top_gene_vec, setdiff(top_gene_lists[[i]], unlist(top_gene_lists[-i]))[1:top_num[i]])
  }
  return(top_gene_vec)
}

#Considerations: 
#1. directly provide major branch points 
#2. create a directed graph based on the root cell from the cds 
#3. combine group test with BEAM or pseudotime test in the same framework 
#4. engineer the function to nicely return the data 

# a function to perform DEG test based on graph topology for cds: 
# fullModelFormulaStr = "~Cell_Type", reducedModelFormulaStr = "~1", top_num = 4, 
graphTest <- function(cds, cth = NULL, major_branch_points = NULL, branch_length = round(ncol(cds) / 25), cores = 1, verbose = F, ...) {
  if(is.null(cth)) {
    if(verbose)
      message('Running graph based test ...')
    ####################################################################################################################################################################################
    #perform DEG test based on the principal graph 
    ####################################################################################################################################################################################
    #retrieve all the branch points and tips 
    root_cell <- row.names(subset(pData(cds), Pseudotime == 0))
    mst <- minSpanningTree(cds)
    all_bps <-  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points #V(mst)[which(degree(mst) > 2)]$name
    all_tips <- V(mst)[which(degree(mst) == 1)]$name

    #convert to directed graph:  
    #which(colnames(cds) %in% c('H7hESC.p7c12r7', 'APS.p1c6r3'))
    mst_adj <- get.adjacency(mst, type = 'lower')
    directed_mst <- graph.adjacency(mst_adj, mode=c("directed"), weighted=T)
    shortest_path_df <- shortest.paths(directed_mst, all_bps, all_tips, weights = NULL)

    #identify major branch points: 
    is_major_bp <- rep(rep(TRUE, length(all_bps)))
    names(is_major_bp) <- all_bps

    #create a graph for the major bps and the the corresponding direct tips: 
    bp_tip_adj_mat <- matrix(rep(0, (length(all_bps) +  length(all_tips))^2 ), 
      nrow = length(all_bps) +  length(all_tips), 
      dimnames = list(c(all_bps, all_tips), 
        c(all_bps, all_tips) ))
    
    #detect the major branch point: 
    direct_tip <- c()
    for(branch_point_ind in 1:length(all_bps)) {
      branch_point <- all_bps[branch_point_ind]
      shortest_paths_list <- all_shortest_paths(mst, branch_point, all_tips, weights = NULL)$res
      overlap_list <- lapply(shortest_paths_list, function(x) length(intersect(all_bps, x$name)))
      direct_tip_tmp <- all_tips[which(unlist(overlap_list) == 1)]
      direct_tip[branch_point_ind] <- direct_tip_tmp
      is_major_bp[branch_point_ind] <- shortest_path_df[branch_point, direct_tip[branch_point_ind]] > branch_length

      print(is_major_bp)
      if(is_major_bp[branch_point_ind])
        bp_tip_adj_mat[branch_point, direct_tip_tmp] <- 1
    }

    #connect the major branch point: 
    major_bps <- names(is_major_bp[is_major_bp])
    
    #add the link between root cell to the first major branch point: 
    first_major_branch_point <- major_bps[which(pData(cds[, major_bps])$Pseudotime == min(pData(cds[, major_bps])$Pseudotime))]
    bp_tip_adj_mat[root_cell, first_major_branch_point] <- 1
    
    #store the regular BEAM results of the major bps: 
    regular_major_BEAM_res <- list()

    for(major_bp_ind in 1:length(major_bps)){
      major_bp <- major_bps[major_bp_ind]
      shortest_paths_list <- all_shortest_paths(mst, major_bp, major_bps)$res
      overlap_list <- lapply(shortest_paths_list, function(x) length(intersect(major_bps, x$name)))
      direct_major_bp <- major_bps[which(unlist(overlap_list) == 2)]
      bp_tip_adj_mat[major_bp, direct_major_bp] <- 1

      #regular BEAM on the major branch points: 
      branchpoint <- which(cds@auxOrderingData[['DDRTree']]$branch_points %in% major_bp)
      regular_major_BEAM_res[[major_bp_ind]] <- BEAM(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
                                                    reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
                                                    branch_states = NULL,
                                                    branch_point=branchpoint,
                                                    cores = cores, ...)
    }

    #create a major bps and tips graph: 
    diag(bp_tip_adj_mat) <- 0

    major_bp_tip_graph <- graph_from_adjacency_matrix(bp_tip_adj_mat, mode = 'undirected', diag = F)

    #retrieve only the connected graph 
    graph_decompose <- decompose.graph(major_bp_tip_graph)
    conn_comp_id <- which(unlist(lapply(graph_decompose, function(x) length(V(x)$name))) > 1)
    major_bp_tip_graph <- graph_decompose[[conn_comp_id]]
    plot(major_bp_tip_graph)

    #perform branch test / two graph test on the branchpoint: 
    graph_res <-  matrix(list(), nrow = length(V(major_bp_tip_graph)), ncol = length(V(major_bp_tip_graph)), dimnames = list(V(major_bp_tip_graph)$name, V(major_bp_tip_graph)$name)) #data.frame of data.frame

    #BEAM only on the neighbor trunk segment to the next major branch point
    for(i in V(major_bp_tip_graph)$name) {
      print(i)
      out_nodes <- neighborhood(major_bp_tip_graph, nodes = i, order = 1)[[1]]$name
      
      if(length(out_nodes) > 2) {
        #branch test:     
        nb_ct <- out_nodes[-1]
        cells <- c()
        
        for(nb_ct_tmp in out_nodes) {
          cells_1 <- shortest_paths(mst, out_nodes[1], nb_ct_tmp)$vpath[[1]]$name
          cells <- unique(c(cells, cells_1)) 
        }
        
        cds_subset <- SubSet_cds(cds, cells)
        
        print(cds_subset)
        #perform BEAM on the subset cds: 
        bp_ind <- which(cds_subset@auxOrderingData[['DDRTree']]$branch_points %in% out_nodes)
        #find the branch point 
        bp_cell <- cds_subset@auxOrderingData[['DDRTree']]$branch_points[[bp_ind]]
        BEAM_test_res <- BEAM(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch", 
                              reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                              branch_states = NULL,
                              branch_point = bp_ind, 
                              cores = cores, ...)
        graph_res[[bp_cell, bp_cell]] <- BEAM_test_res   
      }
      #pseduotime test on the transition path between major branch points
      else {
        cells <- shortest_paths(mst, i, out_nodes[2])$vpath[[1]]$name
        cds_subset <- cds[, cells]
        
        #checked the pval when subtract the smallest pseudotime: the same     
        pseudotime_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr = "~1", cores =cores)
        graph_res[[i, out_nodes[2]]] <- pseudotime_test_res   
        # top_gene_lists[[list_ind]] <- row.names(branchpoint_test_res[order(branchpoint_test_res$qval), ])[1:(4 * top_num)]; list_ind <- list_ind + 1
      }
    }
    res <- list(GraphTest_g = major_bp_tip_graph, GraphTest_res = graph_res,  RegularTest_res = regular_major_BEAM_res)
  }
  else {
  ####################################################################################################################################################################################
  #perform DEG test based on the cth
  ####################################################################################################################################################################################
    if(verbose)
      message('Running regular test ...')
    g <- cth@classificationTree
    branchpoints <- which(degree(g) == 3)
    tips <- which(degree(g) == 1)[-1]
    
    adj_mat <- matrix(rep(0, (length(branchpoints) +  length(tips))^2 ), nrow = length(branchpoints) +  length(tips), dimnames = list(c(names(branchpoints), names(tips)), c(names(branchpoints), names(tips)) ))
    
    adj_g <- get.adjacency(g)
    
    #store the regular BEAM results of the major bps: 
    regular_group_res <- list()
    #create the graph for the deg tests based on the cth: 
    for(index in 1:length(branchpoints)) {
      i <- names(branchpoints)[index]
      print(i)
      s_paths <- shortest_paths(g, i, c(tips, branchpoints))
      adj_mat[i, names(c(tips, branchpoints))] <- as.numeric(unlist(lapply(s_paths$vpath, length)) > 0 & 
              unlist(lapply(s_paths$vpath, function(x) if(length(x) > 1) { #not consider paths pass through another branch point: 
                !any(setdiff(names(branchpoints), i) %in% x$name[1:(length(x)- 1)])
                } else TRUE )))
      
      #regular group test on the branch points (intermediate states): 
      #separate cells into multiple groups based on CellTypeHiearchy:  
      first_branch_cell <- which(degree(g) > 2)[1]
      states <- rep(1, length(g) - first_branch_cell)
      names(states) <- V(g)$name[-(1:first_branch_cell)]
      root_cell <- V(g)[first_branch_cell]$name
      root_cells_nei <- neighbors(g, root_cell)$name
      cell_type_traversal <- names(states)

      mst_traversal <- graph.dfs(g, 
                             root=root_cell, 
                             neimode = "out", 
                             unreachable=FALSE, 
                             father=TRUE)
      curr_state <- 0
      # avoid this
      #      names(mst_traversal$order)
      # [1] "H7hESC"          "APS"             "DLL1PXM"         "D2_25somitomere"
      # [5] "Earlysomite"     "cDM"             "Sclerotome"      "MPS3"
      # [9] "LatM"            NA
      valid_ids <- 2:sum((!is.na(names(mst_traversal$order))))
      for (valid_ids_ind in 1:length(valid_ids)){
        i <- valid_ids[valid_ids_ind]
        curr_node = mst_traversal$order[i]
        curr_node_name = V(dp_mst)[curr_node]$name
        
        ####already avoid above? 
        if (is.na(mst_traversal$father[curr_node]) == FALSE){ 
          if (curr_node$name %in% root_cells_nei){
            curr_state <- curr_state + 1
          }
        }else{}
        
        states[valid_ids_ind] <- curr_state
      }

      cds_subset <- cds[, pData(cds)$CellType %in% cell_type_traversal]
      pData(cds_subset)$cth_state <- "0"
      for(i in names(states)){
          pData(cds_subset)[pData(cds_subset)$CellType == i, 'cth_state'] <- as.character(states[i])
      }
      regular_group_res[[index]] <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~cth_state", reducedModelFormulaStr = "~1", cores =cores)
    }
    diag(adj_mat) <- 0
    graph_test_g <- graph_from_adjacency_matrix(adj_mat, mode = 'directed', diag = F)
    plot(graph_test_g)
    
    #assign the deg test to the edge of the graph: 
    tip_branch_names <- c(names(branchpoints), names(tips))
    tip_branch_len <- length(tip_branch_names)
    graph_res <-  matrix(list(), nrow = tip_branch_len, ncol = tip_branch_len, dimnames = list(tip_branch_names, tip_branch_names )) #data.frame of data.frame
    # top_gene_lists <- list(); list_ind <- 1
    
    #perform branch test / two graph test on the branchpoint: 
    for(i in V(graph_test_g)$name) {
      out_nodes <- neighborhood(graph_test_g, nodes = i, order = 1, mode = 'out')[[1]]$name
      in_nodes <- neighborhood(graph_test_g, nodes = i, order = 1, mode = 'in')[[1]]$name

      if(length(out_nodes) > 2) {
        #branch test:     
        nb_ct_tmp <- neighborhood(g, nodes = i, order = 1, mode = 'out')[[1]]$name
        nb_ct <- nb_ct_tmp[-1]
        if(all(nb_ct %in% pData(cds)$CellType)){
          branchpoint_test_res <- differentialGeneTest(cds[, pData(cds)$CellType %in% nb_ct], fullModelFormulaStr = "~CellType", reducedModelFormulaStr = "~1", cores =cores)
        }
        else branchpoint_test_res <- data.frame()
        graph_res[[i, i]] <- branchpoint_test_res
        # top_gene_lists[[list_ind]] <- row.names(branchpoint_test_res[order(branchpoint_test_res$qval), ])[1:(4 * top_num)]; list_ind <- list_ind + 1
        #perform trunk test
        for(out_nodes_names in out_nodes[-1]) {
          nb_ct <- shortest_paths(g, i, out_nodes_names)$vpath[[1]]$name
          
          if(length(nb_ct) > 2) {
            branchpoint_test_res <- differentialGeneTest(cds[, pData(cds)$CellType %in% nb_ct], fullModelFormulaStr = "~CellType", reducedModelFormulaStr = "~1", cores =cores)
            graph_res[[i, out_nodes_names]] <- branchpoint_test_res   
            # top_gene_lists[[list_ind]] <- row.names(branchpoint_test_res[order(branchpoint_test_res$qval), ])[1:(4 * top_num)]; list_ind <- list_ind + 1
          }
        }
      }
    }
    
    # #obtain unique top 250 deg genes from each test
    # top_gene_vec <- c()
    # for(i in 1:length(top_gene_lists)) {
    #   top_gene_vec <-  c(top_gene_vec, setdiff(top_gene_lists[[i]], unlist(top_gene_lists[-i]))[1:top_num])
    # }
    # 
    res <- list(GraphTest_g = graph_test_g, GraphTest_res = graph_res, RegularTest_res = regular_group_res)
  }
  #return the results as a graph: 
  return(res)
}

#function to project the points into originall projected pints: 

ProjectNewData <- function(cds, new_mat){
  if(row.names(new_mat) != row.names(row.names)){
    stop("genes names from new data doesn't match with that from cell dataset")
  }
    
  W <- URMM_all_fig1b@reducedDimW
  Z <- URMM_all_fig1b@reducedDimS
  Y <- URMM_all_fig1b@reducedDimK
  
  new_lowdim_Data <- MASS::ginv(W)[, which(is.finite(new_mat[, 1]))] %*% new_mat[is.finite(new_mat[, 1]), ] # log(new_mat + 1)
  URMM_all_fig1b@reducedDimS <- cbind(Z, new_lowdim_Data)
  return(cds)
  
}

# URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, initial_method = PCA, norm_method = 'log', max_components = 2, verbose = T, auto_param_selection = F, scaling = F) #, initial_method = DM , initial_method = DM
# URMM_all_fig1b <- orderCells(URMM_all_fig1b, num_paths=1, root_state = NULL)
# 
# knockout_cds <- URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout', 'Gfi1_knockout', 'Irf8_knockout')]
# fData(knockout_cds)$use_for_ordering <- fData(URMM_all_fig1b)$use_for_ordering
# knockout_cds <- estimateSizeFactors(knockout_cds)
# 
# #test the strategy: 
# new_mat <- monocle::normalize_expr_data(URMM_all_fig1b, 'log')
# new_mat <- as.matrix(Matrix::t(scale(Matrix::t(new_mat))))
# new_mat <- new_mat[!is.na(row.names(new_mat)), ]
# 
# new_mat <- normalize_expr_data(knockout_cds, 'log')
# new_mat <- as.matrix(Matrix::t(scale(Matrix::t(new_mat))))
# new_mat <- new_mat[!is.na(row.names(new_mat)), ]
# 
# 
# qplot(Z[1, ], Z[2, ]) + geom_point(aes(new_lowdim_Data[1, ], new_lowdim_Data[2, ], color = pData(knockout_cds)$Type))
# 
# trimTree <- function(cds, num_paths = 2, min_branch_thrsld = 0.1){ #
#   mst <- minSpanningTree(cds)
#   branch_len_thrsld <- min_branch_thrsld * length(V(mst))
#   
#   mst_traversal <- graph.dfs(dp_mst, 
#                              root=root_cell, 
#                              neimode = "all", 
#                              unreachable=FALSE, 
#                              father=TRUE)
#   mst_traversal$father <- as.numeric(mst_traversal$father)
#   curr_state <- 1
#   
#   states = rep(1, length(V(mst)))
#   names(states) <- V(mst)$name
#   parents = rep(NA, length(V(mst)))
#   names(parents) <- V(mst)$name
#   
#   for (i in 1:length(mst_traversal$order)){
#     curr_node = mst_traversal$order[i]
#     curr_node_name = V(mst)[curr_node]$name
#     
#     if (is.na(mst_traversal$father[curr_node]) == FALSE){
#       parent_node = mst_traversal$father[curr_node]
#       parent_node_name = V(mst)[parent_node]$name
#       parent_node_state = states[parent_node_name]
#       
#       
#       if (degree(mst, v=parent_node_name) > 2){ 
#         #the two tip cells has length larger than a certain ratio 
#         branch_tip_inds <- which(degree(mst, mst_traversal$order$name[-c(1:i)]) == 1)[1:2]
#         
#         #update state only when the cells have more than a certain threshold 
#         if(branch_tip_inds[1] >  branch_len_thrsld & diff(branch_tip_inds) > branch_len_thrsld){
#           curr_state <- curr_state + 1
#         }
#       }
#     }else{
#       parent_node = NA
#       parent_node_name = NA
#     }
#     
#     curr_node_state = curr_state
#     states[curr_node_name] <- curr_node_state
#     parents[curr_node_name] <- parent_node_name
#   }
#   
#   pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex
#   pr_graph_cell_proj_closest_vertex[, 1] <- paste("Y_", pr_graph_cell_proj_closest_vertex[, 1], sep = '')
#   pData(cds)$State <- states[pr_graph_cell_proj_closest_vertex[, 1]]
#   
#   return(cds)
# }

#compare graph: 
calGraphSimilarity <- function(g1, g2) {
  # source('./delta_con.R', echo = T)
  graph1 <- get.edgelist(g1)
  graph2 <- get.edgelist(g2)
  
  graph1 <- matrix(as.numeric(get.edgelist(g1)), ncol = 2)
  graph2 <- matrix(as.numeric(get.edgelist(g2)), ncol = 2)
  nnodes <- max(rbind(graph1, graph2))
  
  delta_con(data.frame(graph1), data.frame(graph2), nnodes = nnodes, percent = 1)
}

# #proven version of trimTree: (too slow and will be replaced) 
# trimTree <- function(cds, min_branch_thrsld = 0.1){ #
#   tmp <- colnames(cds)[pData(cds)$Pseudotime == 0]
#   pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex
#   
#   root_cell <- paste('Y', pr_graph_cell_proj_closest_vertex[tmp, ], sep = '_')
#   
#   adjusted_S <- t(cds@reducedDimK)
#   
#   dp <- as.matrix(dist(adjusted_S))
#   
#   # Build an MST of the cells in ICA space.
#   gp <- graph.adjacency(dp, mode="undirected", weighted=TRUE)
#   dp_mst <- minimum.spanning.tree(gp)
#   
#   # Build the PQ tree
#   next_node <<- 0
#   res <- pq_helper(dp_mst, use_weights=FALSE, root_node=root_cell)
#   
#   tip_cells <- names(which(degree(dp_mst) == 1))
# 
#   cc_ordering_df <- extract_ddrtree_ordering_xj(dp_mst, dp, root_cell)
#   tip_states <- as.numeric(cc_ordering_df[tip_cells, 'cell_state'])
#   state_tbl <- table(as.numeric(subset(cc_ordering_df, cell_state %in% tip_states)$cell_state))
#   tiny_branches_num <- sum(state_tbl < min_branch_thrsld * nrow(cc_ordering_df))
#     
#   num_paths <- length(tip_states)  - tiny_branches_num - 1
# 
#   order_list <- monocle::extract_good_branched_ordering(res$subtree, res$root, dp, num_paths, FALSE)
#   cc_ordering <- order_list$ordering_df
#   row.names(cc_ordering) <- cc_ordering$sample_name
#   
#   tmp <- pr_graph_cell_proj_closest_vertex[row.names(pData(cds)), ]
#   projection_points_names <- paste('Y_', tmp, sep = '')
#   
#   pData(cds)$State <- cc_ordering[projection_points_names, ]$cell_state
#   return(cds)
# }

# trimTree1 <- function(cds, num_paths = 2, min_branch_thrsld = 0.1){ #
#   mst <- minSpanningTree(cds)
#   branch_len_thrsld <- min_branch_thrsld * length(V(mst))
#   root_cell <- colnames(cds)[which(pData(cds)$Pseudotime == 0)]
#   
#   mst_traversal <- graph.dfs(mst, 
#                              root=root_cell, 
#                              neimode = "all", 
#                              unreachable=FALSE, 
#                              father=TRUE)
#   mst_traversal$father <- as.numeric(mst_traversal$father)
#   curr_state <- 1
#   
#   states = rep(1, length(V(mst)))
#   names(states) <- V(mst)$name
#   parents = rep(NA, length(V(mst)))
#   names(parents) <- V(mst)$name
#   
#   for (i in 1:length(mst_traversal$order)){
#     curr_node = mst_traversal$order[i]
#     curr_node_name = V(mst)[curr_node]$name
#     
#     if (is.na(mst_traversal$father[curr_node]) == FALSE){
#       parent_node = mst_traversal$father[curr_node]
#       parent_node_name = V(mst)[parent_node]$name
#       parent_node_state = states[parent_node_name]
#       
#       
#       if (degree(mst, v=parent_node_name) > 2){ 
#         #the two tip cells has length larger than a certain ratio 
#         branch_tip_inds <- which(degree(mst, mst_traversal$order$name[-c(1:i)]) == 1)[1:2]
#         
#         #update state only when the cells have more than a certain threshold 
#         if(branch_tip_inds[1] >  branch_len_thrsld & diff(branch_tip_inds) > branch_len_thrsld){
#           curr_state <- curr_state + 1
#         }
#       }
#     }else{
#       parent_node = NA
#       parent_node_name = NA
#     }
#     
#     curr_node_state = curr_state
#     states[curr_node_name] <- curr_node_state
#     parents[curr_node_name] <- parent_node_name
#   }
#   
#   pData(cds)$State <- states
#   return(cds)
# }
# 
# #prune tree based on the branch length
# trimTree <- function(cds, num_paths = 2, min_branch_thrsld = 0.1){ #
#   mst <- minSpanningTree(cds)
#   branch_len_thrsld <- min_branch_thrsld * length(V(mst))
#   root_state <- as.numeric(as.character(pData(cds)[colnames(cds)[which(pData(cds)$Pseudotime == 0)], 'State']))
#   root_cell <- select_root_cell(cds, root_state, reverse)
#   
#   #find the tip cells as well as the cells on the branch point
#   tip_cells <- V(mst)[degree(mst) == 1]$name
#   branch_cells <- V(mst)[degree(mst) > 2]$name
#   
#   ignore_tip_branch_cells_df <- do.call(rbind.data.frame, lapply(tip_cells, function(x) {
#     vpaths <- shortest_paths(mst, x, branch_cells)$vpath
#     paths_len <- unlist(lapply(vpaths, length))
#     if(any(paths_len < branch_len_thrsld)){
#       return(data.frame(tip_cell = x, branch_cell = branch_cells[which.min(paths_len)]))
#     }
#   }))
#   
#   ignore_branch_cells <- as.character(ignore_tip_branch_cells_df$branch_cell)
#   ignore_tip_cells <- as.character(ignore_tip_branch_cells_df$tip_cell)
#   
#   mst_traversal <- graph.dfs(mst, 
#                              root=root_cell, 
#                              neimode = "all", 
#                              unreachable=FALSE, 
#                              father=TRUE)
#   mst_traversal$father <- as.numeric(mst_traversal$father)
#   curr_state <- 1
#   
#   states = rep(1, length(V(mst)))
#   names(states) <- V(mst)$name
#   parents = rep(NA, length(V(mst)))
#   names(parents) <- V(mst)$name
#   
#   for (i in 1:length(mst_traversal$order)){
#     curr_node = mst_traversal$order[i]
#     curr_node_name = V(mst)[curr_node]$name
#     
#     if (is.na(mst_traversal$father[curr_node]) == FALSE){
#       parent_node = mst_traversal$father[curr_node]
#       parent_node_name = V(mst)[parent_node]$name
#       parent_node_state = states[parent_node_name]
#       
#       #update state only when the cells have more than a certain threshold 
#       #the two tip cells has length larger than a certain ratio 
#       if (degree(mst, v=parent_node_name) > 2 & !(parent_node_name %in% ignore_branch_cells)){ 
#           curr_state <- curr_state + 1
#       }
#     }else{
#       parent_node = NA
#       parent_node_name = NA
#     }
#     
#     curr_node_state = curr_state
#     states[curr_node_name] <- curr_node_state
#     parents[curr_node_name] <- parent_node_name
#   }
#   
#   #post-process those off branches with miss assigned states
#   for(i in 1:nrow(ignore_tip_branch_cells_df)) {
#     x <- ignore_tip_branch_cells_df[i, ]
#     off_branch_path <- names(shortest_paths(mst, as.character(x[[1]]), as.character(x[[2]]))$vpath[[1]])
#     states[off_branch_path] <- states[as.character(x[[2]])]
#   }
#     
#   pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex
#   
#   match_states <- paste('Y', pr_graph_cell_proj_closest_vertex, sep = '_')
#   names(match_states) <- row.names(pr_graph_cell_proj_closest_vertex)
#   
#   pData(cds)[names(match_states), 'State'] <- as.character(pData(cds)[names(match_states), 'State']) 
#   pData(cds)[names(match_states), 'State'] <- factor(states[match_states])
#   return(cds)
# }

# min_branch_thrsld is not used for now
# trimTree <- function(cds, num_paths = 2, min_branch_thrsld = 0.1){ #
#   mst <- minSpanningTree(cds)
#   branch_len_thrsld <- min_branch_thrsld * length(V(mst))
#   root_state <- as.numeric(as.character(subset(pData(cds), Pseudotime == 0)[, 'State']))
#   root_cell <- monocle:::select_root_cell(cds, root_state)
#   
#   #find the tip cells as well as the cells on the branch point
#   pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex
#   
#   if(length(V(minSpanningTree(cds))) < ncol(cds))
#     samples_pc_pairs <- paste('Y', pr_graph_cell_proj_closest_vertex, sep = '_')
#   else
#     samples_pc_pairs <- colnames(cds)[pr_graph_cell_proj_closest_vertex]
#   
#   names(samples_pc_pairs) <- row.names(pr_graph_cell_proj_closest_vertex)
#   
#   tip_cells <- V(mst)[degree(mst) == 1]$name
#   branch_cells <- V(mst)[degree(mst) > 2]$name
#   
#   #find the length of any tip branches, and only select the longest branches
#   # cells_project_tips <- names(samples_pc_pairs)[which(samples_pc_pairs %in% tip_cells)]
#   
#   #map the tip principal point to a state
#   state_of_tip_cells <- lapply(tip_cells, function(x){
#     cells_project_tips <- names(samples_pc_pairs)[which(samples_pc_pairs %in% x)]
#     names(sort(table(as.character(pData(cds)[cells_project_tips, 'State'])), decreasing = T))#[1] #find only the state with largest number of cells
#   })
#   state_of_tip_cells <- unlist(state_of_tip_cells)
#   names(state_of_tip_cells) <- tip_cells
#   
#   #number for each state
#   tip_states_order <- sort(table(as.character(subset(pData(cds), State %in% state_of_tip_cells)[, 'State'])), decreasing = T)
#   top_states <- names(tip_states_order)[1:c(num_paths + 1)] #only select the first num_paths branches with largest number of cells
#   tips_for_branch <- names(state_of_tip_cells)[state_of_tip_cells %in% top_states]
#   
#   #find the tip cells and their nearest branch points, then record the tip cells/branch points which are not included in the longest pathes
#   ignore_tip_branch_cells_df <- do.call(rbind.data.frame, lapply(tip_cells, function(x) {
#     vpaths <- shortest_paths(mst, x, branch_cells)$vpath
#     paths_len <- unlist(lapply(vpaths, length))
#     if(!(x %in% tips_for_branch)){ #any(paths_len < branch_len_thrsld) |
#       return(data.frame(tip_cell = x, branch_cell = branch_cells[which.min(paths_len)]))
#     }
#   }))
#   
#   ignore_branch_cells <- as.character(ignore_tip_branch_cells_df$branch_cell)
#   ignore_tip_cells <- as.character(ignore_tip_branch_cells_df$tip_cell)
#   
#   mst_traversal <- graph.dfs(mst, 
#                              root=root_cell, 
#                              neimode = "all", 
#                              unreachable=FALSE, 
#                              father=TRUE)
#   mst_traversal$father <- as.numeric(mst_traversal$father)
#   curr_state <- 1
#   
#   states = rep(1, length(V(mst)))
#   names(states) <- V(mst)$name
#   parents = rep(NA, length(V(mst)))
#   names(parents) <- V(mst)$name
#   
#   for (i in 1:length(mst_traversal$order)){
#     curr_node = mst_traversal$order[i]
#     curr_node_name = V(mst)[curr_node]$name
#     
#     if (is.na(mst_traversal$father[curr_node]) == FALSE){
#       parent_node = mst_traversal$father[curr_node]
#       parent_node_name = V(mst)[parent_node]$name
#       parent_node_state = states[parent_node_name]
#       
#       #update state after branch point only when the cells located at the branches with num_path length
#       if (degree(mst, v=parent_node_name) > 2 & !(parent_node_name %in% ignore_branch_cells)){ 
#         curr_state <- curr_state + 1
#       }
#     }else{
#       parent_node = NA
#       parent_node_name = NA
#     }
#     
#     curr_node_state = curr_state
#     states[curr_node_name] <- curr_node_state
#     parents[curr_node_name] <- parent_node_name
#   }
#   
#   #post-process those off branches with miss assigned states
#   for(i in 1:nrow(ignore_tip_branch_cells_df)) {
#     x <- ignore_tip_branch_cells_df[i, ]
#     off_branch_path <- names(shortest_paths(mst, as.character(x[[1]]), as.character(x[[2]]))$vpath[[1]])
#     states[off_branch_path] <- states[as.character(x[[2]])]
#   }
# 
#   pData(cds)[names(samples_pc_pairs), 'State'] <- as.character(pData(cds)[names(samples_pc_pairs), 'State']) 
#   pData(cds)[names(samples_pc_pairs), 'State'] <- factor(states[samples_pc_pairs])
#   return(cds)
# }
# the real one is here!! 
trimTree <- function(cds, num_paths = 2, min_branch_thrsld = 0.1) {
  dp_mst <- cds@minSpanningTree
  pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex
  
  root_sample <- row.names(subset(pData(cds), Pseudotime == 0))
  root_cell_id <- cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex[root_sample, ]
  
  if(length(V(dp_mst)) < ncol(cds)) {
    root_cell  <- paste('Y', root_cell_id, sep = '_')
    pr_graph_cell_proj_closest_vertex[, 1]  <- paste('Y', pr_graph_cell_proj_closest_vertex, sep = '_')
    cell_pairwise_dist <- cellPairwiseDistances(cds)
  }
  else {
    root_cell <- colnames(cds)[root_cell_id]
    cell_pairwise_dist <- as.matrix(dist(t(cds@reducedDimK)))
    pr_graph_cell_proj_closest_vertex[, 1]  <- row.names(pr_graph_cell_proj_closest_vertex)
  }
  # Build the PQ tree
  next_node <<- 0
  res <- monocle:::pq_helper(dp_mst, use_weights=FALSE, root_node=root_cell)
  
  # cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
  order_list <- monocle:::extract_good_branched_ordering(res$subtree, res$root, cell_pairwise_dist, num_paths, FALSE)
  cc_ordering <- order_list$ordering_df
  row.names(cc_ordering) <- cc_ordering$sample_name
  
  pData(cds)$State <- cc_ordering[pr_graph_cell_proj_closest_vertex[, 1], 'cell_state']
  
  return(cds)
}

make_cds <- function (exprs_matrix, pd, fd, expressionFamily) {
  cds <- newCellDataSet(exprs_matrix, 
                        phenoData = new("AnnotatedDataFrame", data = pd), 
                        featureData = new("AnnotatedDataFrame", data = fd), 
                        expressionFamily = expressionFamily, 
                        lowerDetectionLimit = 0.1)
  
  if(identical(expressionFamily, tobit())) {
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
  }
  return(cds)
}
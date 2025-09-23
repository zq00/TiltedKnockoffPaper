# Compute a node ordering for the metropolized knockoff method
# Input: A graph defined in the class in the class `graphNEL` (from the package `graph`)

order_nodes <- function(g){
  junction_tree <- get_junction_tree(g)
  # edges from the junction tree
  
  ordered_nodes <- get_node_order(junction_tree$junction_tree, junction_tree$cliques, p)
  
  return(ordered_nodes)
}

# compute the junction tree from a graph
get_junction_tree <- function(g){
  # triangulated graph
  triangulated_gR <- gRbase::triangulate(g)
  edges_triangulated <- graph::edges(triangulated_gR)
  # clique graph from triangulated graph
  edge_matrix <- NULL
  nodes <- graph::nodes(triangulated_gR)
  for(i in 1:length(nodes)){
    new_edges <- cbind(nodes[i], edges_triangulated[[i]])
    edge_matrix <- rbind(edge_matrix, new_edges)
  }
  triangulated_graph <- graph_from_edgelist(edge_matrix, directed = F)
  cliques_triangulated <- max_cliques(triangulated_graph)
  # construct edges in the clique graph and compute weights for each edge
  ncliques <- length(cliques_triangulated)
  edge_weight <- NULL
  for(i in 1:length(nodes)){
    clique_containing_i <- NULL
    cat(i, ":")
    for(j in 1:ncliques){
      if(nodes[i] %in% names(cliques_triangulated[[j]])){
        clique_containing_i <- c(clique_containing_i, j)
      }
    }
    # add edges between cliques if they have not showed up yet
    new_edge_weight <- NULL
    for(k in 1:length(clique_containing_i)){
      for(l in 1:length(clique_containing_i)){
        if(k != l){
          # Compute the weights 
          ind1 <- names(cliques_triangulated[[clique_containing_i[k]]])
          ind2 <- names(cliques_triangulated[[clique_containing_i[l]]])
          new_weight <- length(intersect(ind1, ind2)) 
          new_edge <- sort(c(clique_containing_i[l], clique_containing_i[k]))
          # order the two nodes in increasing order
          new_edge_weight <- rbind(new_edge_weight, c(new_edge, new_weight))
        }
      }
    }
    new_edge_weight <- base::unique(new_edge_weight)
    edge_weight <- rbind(edge_weight, new_edge_weight)
    cat("\n")
  }
  edge_weight <- base::unique(edge_weight)
  # construct a clique graph
  clique_graph <- igraph::graph_from_edgelist(edge_weight[,1:2], directed = F)
  junction_tree <- mst(clique_graph, weights = -edge_weight[,3]) # the junction tree
  
  return(list(junction_tree = junction_tree,
              cliques = cliques_triangulated))
}

# compute a node ordering from the junction tree and the cliques (the junction tree is a clique graph)
get_node_order <- function(junction_tree, cliques, p){
  ncliques <- length(cliques)
  tree_edges <- as_edgelist(junction_tree)
  # list of ordered vertices
  ordered_vertices <- NULL
  neighbors <- list() # the neighbors of the i-th vertex, not the one in the ordered list, e.g. if the first one in order is 89, then find its neighbor by calling neighbors[[89]]
  
  i <- 0
  V <- as.numeric(names(which(table(as.vector(tree_edges)) == 1)[1]))
  while(i < p){
    which_edge <- which(rowSums(tree_edges == V) == 1)
    if(sum(which_edge) == 0){
      # if no V', add all of the vertices in V to the list (we will remove those that are already in the tree later)
      new_vertices <- names(cliques[[V]])
    }else{
      VPrime <- setdiff(tree_edges[which_edge,], V)
      # what are the vertices in V but not in VPrime?
      new_vertices <- setdiff(names(cliques[[V]]), names(cliques[[VPrime]]))
    }
    
    if(length(new_vertices) == 0){
      # remove V from the list and continue
      stop("There's no vertex in V that is not in VPrime!")
    }else{
      # add the new vertex to the list 
      for(j in 1:length(new_vertices)){
        i <- i + 1
        ordered_vertices[i] <- as.numeric(new_vertices[j])
        new_nbd <- list(vertex = as.numeric(new_vertices[j]), nbd = as.numeric(setdiff(names(cliques[[V]]), new_vertices[1:j])))
        neighbors[[as.numeric(new_vertices[j])]]  <- new_nbd
      }
      
      # remove V from the edge list 
      if(sum(which_edge) >  0){
        tree_edges <- tree_edges[-which_edge, , drop = F]
      }
    }
    if(sum(tree_edges == VPrime) == 1 | nrow(tree_edges) == 0){
      V = VPrime
    }else{
      V <- as.numeric(names(which(table(as.vector(tree_edges)) == 1)[1]))
    }
  }
  
  return(list(order = ordered_vertices, 
              nbd = neighbors))
}


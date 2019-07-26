# automatically recode amalgamated traits 
# using the matrices returned by amalgamate_deps, recode this data into td and return td
recode_traits <- function(td, M){
  # iterate through each new amalgamated trait and recode the matrix for that trait into td
  for(i in seq_along(M$new_traits)){
    trait_name <- M$new_traits[[i]]
    
    # if trait is not connected (ie has no dependencies)
    if(strcmp(M$traits[[i]], M$new_traits[[i]]) == TRUE){
      states <- setNames(c("0", "1"), 1:2)
      td <- recode_td(td, trait_name, states)
    }
    else{ # if it has dependencies, recode the combined matrix for that subgraph into td
      # set gtraits 
      gtraits <- M$traits[[i]]
      # set states
      states <- colnames(M$M[[trait_name]])
      
      td <-recode_td(td, gtraits, states) 
    }
  }
  return(td)
}


# given a a dependency matrix, generate a graph for the traits and combined matrices for each set of dependencies on that graph
# return the matrices, our graph, the graph with organized subgraphs, and lists of the old and new trait names
amalgamate_deps <- function(dep_mat) {
  M <- list() # store all of our combined matrices
  traits <- list() 
  new_traits <- list() # store new trait names
  
  # remove redundant dependencies
  dep_mat <- remove_indirect(dep_mat)
  
  # generate our graph
  g <- graph_from_adjacency_matrix(t(as.matrix(dep_mat)))
  con.comp <- components(g, "weak")
  
  # break up our graph into subgraphs, and combine matrices for the traits in those graphs
  for(i in 1:con.comp$no){
    g_sub <- induced_subgraph(g, which(con.comp$membership == i))
    bin_mats <- init_binary_matrix(g_sub) 
    gtraits <- names(bin_mats) # get names of traits
    newTraitName <- paste(gtraits, collapse="+") # name of new, amalgamated trait
    
    traits[[i]] <- gtraits # store names of traits before they were amalgamated
    new_traits[[i]] <- newTraitName # store name of new trait
    M[[newTraitName]] <- combine_subgraph(g_sub) # store new matrix
  }
  return(list(M=M, graph=g, con.comp=con.comp, traits=traits, new_traits=new_traits))
}


# combine binary matrices for a subgraph on a dependency graph
# return the final combined matrix of all traits
combine_subgraph <- function(subgraph){
  # sort graph so that the "root" comes first
  sorted_g <- topo_sort(subgraph, mode = "out")
  
  bin_mats <- init_binary_matrix(subgraph) 
  
  # setting up:
  M <- list() 
  node_list <- names(sorted_g) 
  # trait_order <- list()  # keeps track of our traits in order of their dependencies
  ancestor_v <- node_list[[1]] # set to root
  M[[1]] <- bin_mats[[ancestor_v]]
  # trait_order[[1]] <- ancestor_v
  
  # now combine matrices based on dependencies
  # for each node in subgraph, combine trait matrix with its ancestor matrix
  for (i in seq_along(node_list)){
    if (i != 1) { # if we are on the root node don't combine
      node_name <- node_list[[i]]
    
      if(i == 2) { # if we are currently on the first non-root vertex, we want to combine our current binary matrix with the binary matrix of the root vertex
        # set dependent state
        dep_state <- get_dep_state(subgraph, V(subgraph)[name == node_name], bin_mats[[ancestor_v]], node_list) #trait_order)
        M[[i]] <- comb2matrices(bin_mats[[ancestor_v]], bin_mats[[node_name]], dependent.state = dep_state)
      }
      else { # otherwise, combine with the resulting matrix from the combination we just did
        # set dependent state
        dep_state <- get_dep_state(subgraph, V(subgraph)[name == node_name], M[[i - 1]], node_list) #trait_order)
        M[[i]] <- comb2matrices(M[[i - 1]], bin_mats[[node_name]], dependent.state = dep_state)
      }
      # for debugging:
      # print(M[[i]]) 
      # print(dep_state)
    }
  }
  M[[newTraitName]] <- M[[i]]
  return(M[[newTraitName]])
}


# given a graph and a vertex v on that graph, the trait matrix for the ancestor node of v,
# and a list containing the traits in order of their combination into the matrix, 
# return the columns in the ancestor trait matrix which v depends on
get_dep_state <- function(subgraph, v, ancestor_mat, trait_order) {
  # get colnames of ancestor matrix
  mat_cols <- colnames(ancestor_mat)
  
  ## find which position the dependent trait(s) occupies in the colnames of our trait matrix
  adj_v <- adjacent_vertices(subgraph, v, mode="in") 
  ancestor_v <- names(adj_v[[v$name]]) # store the list of ancestor vertices of v
  
  # create a list containing positions in the ancestor_mat character where we have dependencies
  trait_pos <- list() 
  for(i in seq_along(ancestor_v)){
    pos <- which(trait_order == ancestor_v[[i]]) 
    trait_pos[[i]] <- pos
  }
  
  # for each position in the ancestor_mat character that v depends on,
  # store the corresponding column and add to our list
  depstate_list <- list() 
  for(i in seq_along(trait_pos)){
    depstate <- substr(mat_cols, trait_pos[[i]], trait_pos[[i]]) # extract values at that position for each column 
    depstate_list[[i]] <- c(which(depstate == 1)) # return the columns where the value for our dependent trait is 1
  }
  depstate_list <- unlist(depstate_list)
  return(unique(depstate_list))
}


# function for removing redundant dependency edges from our graph
remove_indirect <- function(dependency_matrix){
  z <- dependency_matrix
  diag(z) <- NA
  tmp <-  z & (!logical_mult(z,z))
  tmp[is.na(tmp)] <- 1
  res <- z*tmp
  return(res)
}
logical_mult = function(x,y)
{
  I = dim(x)[1]
  K = dim(x)[2]
  J = dim(y)[2]
  
  z = rep(FALSE,I*J)
  dim(z) = c(I,J)
  for(i in 1:I)
  {
    for(j in 1:J)
    {
      for(k in 1:K)
      {
        z[i,j] = z[i,j] || (x[i,k] && y[k,j])
      }
    }
  }
  z
}



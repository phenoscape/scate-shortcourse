#' Function to display the number of matching species and average value across them from an Ontotrace matrix
#' @export
print_coverage <- function(x){
  coverage <- apply(x, 2, function(x) sum(!is.na(x)))
  average <- sapply(x, function(x) mean(as.numeric(na.omit(x[x %in% c("0", "1")])), na.rm=TRUE))
  cover <- cbind(coverage, average)
  tmp <- filter(data.frame(traits=rownames(cover), cover), coverage > 0, average < 1, average > 0) %>% arrange(., desc(coverage))
  print(tmp)
}

#' Utility function for stripping off url from IRI
#' @export
strip_IRI <- function(x){
  x <- gsub("http://purl.obolibrary.org/obo/", "", x)
  x <- gsub("_", ":", x)
  return(x)
}

#' Utility function for cleaning up character data table after amalgamating characters
#' @export
dropDependentTraits <- function(char_info, dep.mat, td){
  char_info_comb <- char_info[which(apply(dep.mat, 1, sum, na.rm=TRUE)==0), c(1,5)]
  new.traits <- colnames(td$dat)
  old.traits <- sapply(new.traits, function(x) strsplit(x, "+", fixed=TRUE)[[1]][1])
  trait.trans <- setNames(new.traits, old.traits)
  char_info_comb$ID <- unname(trait.trans[as.character(char_info_comb$ID)])
  return(char_info_comb)
}

#' Utility function for filtering based on missing traits and taxa
#' @export
filter_coverage <- function(td, traits=0, taxa=0){
  taxa_coverage <- apply(td$dat, 1, function(x) mean(as.numeric(!is.na(x))))
  trait_coverage <- apply(td$dat, 2, function(x) mean(as.numeric(!is.na(x))))
  td <- dplyr::filter(td, taxa_coverage > taxa)
  td <- dplyr::select(td, which(trait_coverage > traits))
  return(td)
}

#' Data object for Rev scripts
#' @export
ratemat1 <- function() {
  '\nfor (i in 1:NUM_STATES) {\nfor (j in 1:NUM_STATES) {\nrates[i][j] <-0.0\n}\n}\n#rate prior\nr1 ~ dnExp(20)\nr2 ~ dnExp(20)\n\nmoves[++mvi] = mvScale(r1, lambda=1, tune=true, weight=2)\nmoves[++mvi] = mvScale(r2, lambda=1, tune=true, weight=2)
\n\n# place rate categories into matrix\nrates[2][1]:=r1\nrates[1][2]:=r2\n\n\nrate_matrix := fnFreeK(transition_rates=rates, rescaled=false, matrixExponentialMethod="eigen")\n\nroot_freq <- simplex(1, 1)\n\n'
  
}

#' recodes a treedata object based on amalgamated characters
#' @export
recode_td <- function(td, traits, states, hidden0=numeric(0)){
  tmp <- select(td, traits)
  obs_char <- which(!(1:length(traits) %in% hidden0))
  for(i in 1:ncol(tmp$dat)){
    recode0 <- ifelse(i %in% hidden0, "*", "0")
    missingObs <- ifelse(i %in% obs_char, "X", "*")
    tmp$dat[[i]] <- recode(as.character(tmp$dat[[i]]), "1"="1", "0"=recode0, "1 and 0"=missingObs, "0 and 1"=missingObs, .missing=missingObs)
  }
  new.char <- unite(tmp$dat, "new", sep="")
  new.char <- unname(sapply(new.char[[1]], function(x) paste(which(grepl(glob2rx(x), states))-1, collapse="&")))
  new.char[which(new.char=="")] <- "?"
  
  new.td <- select(td, -one_of(traits))
  new.td$dat[[paste(traits, collapse="+")]] <- new.char
  return(new.td)
}

#' Plots a heatmap along with a phylogeny and trait tree
#' @export
ontologyHeatMap <- function(td, njt, start=3, margs=c(0.2, 0.25), ...){
  #vals <- na.omit(unique(do.call(c, lapply(3:ncol(td$dat), function(x) unique(as.character(td$dat[[x]]))))))
  
  X <- do.call(cbind, lapply(start:ncol(td$dat), function(x) as.numeric(recode(td$dat[[x]], "0 and 1"=0.5, "1 and 0"=0.5, "1"=1, "0"=0, "2"=2, "3"=3))))
  colnames(X) <- colnames(td$dat)[start:ncol(td$dat)]
  
  .vals <- sort(na.omit(unique(as.vector(X))))
  dimx <- dim(X)
  
  tree1 <- force.ultrametric(td$phy,method = "extend")
  
  if(!is.null(njt)){
    X <- X[,njt$edge[njt$edge[,2] <= length(njt$tip.label),2]]
    tree2 <- chronopl(njt, 1)
    tree2$edge.length <- tree2$edge.length/(max(branching.times(tree2)))*margs[2]*dimx[1]
    png(tempfile())
    invisible(h2 <- plot(tree2, plot = FALSE, direction = "downwards", show.tip.label=FALSE))
    dev.off()
  } else{
    h2 <- list(x.lim=c(1,dimx[2]+1), y.lim=c(0,0.2*dimx[1]))
  }
  
  #alters the length of the tips of the trees
  tree1$edge.length <- tree1$edge.length/(max(branching.times(tree1)))*margs[1]*dimx[2]
  
  #Changes the direction of the top plot
  png(tempfile())
  invisible(h1 <- plot(tree1, plot = FALSE, cex=0.5))
  dev.off()
  
  # adjustible color palette for the plot and legend
  colors <- c("#ffeaa7","#fab1a0", "#e17055")
  
  
  # this is all setting boundaries for the different plots and combining them into one image
  par(mar = c(0,0,0,0))
  plot(0,0, type = 'n', xlim = c(0,h1$x.lim[2]+h2$x.lim[2]), ylim=c(0,h1$y.lim[2]+h2$y.lim[2]))
  
  image(seq(h1$x.lim[2]+1,h1$x.lim[2]+h2$x.lim[2], length.out=ncol(X)), seq(1, h1$y.lim[2], length.out=nrow(X)), t(X),xlim=c(1+h1$x.lim[2],h1$x.lim[2]+h2$x.lim[2]+1) ,ylim=c(0, h1$y.lim[2]-1), add=TRUE, col=colors)
  
  legend(0, (h1$y.lim[2]+h2$y.lim[2])*.99, legend=.vals ,pch=22, pt.bg=colors)
  
  par(new = TRUE)
  
  plot(tree1, x.lim=c(0,(1+margs[2])*(h2$x.lim[2]+h1$x.lim[1])), y.lim=c(0,h1$y.lim[2]+h2$y.lim[2]),...)
  
  if(!is.null(njt)){
    par(new = TRUE)
    plot(tree2, direction = "downwards", x.lim=c(-h1$x.lim[2],h2$x.lim[2]), y.lim=c((-h1$y.lim[2])-0.01*dimx[1],h2$y.lim[2]))
  }
  
  return ()
}

#' Makes a trait tree using semantic similarity
#' @export
makeTraitTree <- function (td, skip=1:2){
  traits <- colnames(td$dat)
  traits <- traits[-(1:2)] #delete otu data
  traits
  
  semanticSimilarityMatrix <- jaccard_similarity(terms = traits, .colnames = "label", .labels = traits)
  
  #rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- traits
  
  #Neighbor-joining tree of SS matrix
  njt <- nj(1-semanticSimilarityMatrix)
  
  return (njt)
}

#' Function for processing revbayes stochastic maps
#' @export
prepareMaps <- function(td, dirR, dirW, discretization_level=100, start_tree=1, end_tree=2) {
  characters <- colnames(td$dat)
  characters <- gsub(" ", "_", characters)
  #####################################
  # Read a sample of 2 maps from .stm files and save them in the proper format .stmR
  #####################################
  
  for (i in 1:length(characters))
  {
    .tree<-read_Simmap_Rev(paste0(dirR, characters[i], ".stm"),
                           start=start_tree, end=end_tree,
                           save = NULL) 
    tree  <- read.simmap(text=.tree, format="phylip")
    
    
    write.simmap(tree, file=paste0(dirW, characters[i], ".stmR"))
  }
  ##########
  
  #####################################
  # Read stmR, discretize maps, and save each map as a separate rds file; 
  #all rds filea for a chracter are stored in a zip archive
  #####################################
  
  for (i in 1:length(c))
  { 
    # read in undesritezed trees
    #print(paste0("Reading ", characters[i]))
    sim=read.simmap(file=paste0(dirW, characters[i], ".stmR"), format="phylip")
    
    # discretize trees by looping over sample and saving as rds
    
    for (j in 1:length(sim)){
      tryCatch({
        
        #print(paste0("Discretizing tree ", j))
        
        ## errors with na
        
        ##
        
        ##### make trees equal with template
        sim.d<-make_tree_eq(td$phy, sim[[j]], round=5)
        ###
        
        #sim.d<-discr_Simmap_all(sim[[j]], 1000)
        sim.d<-discr_Simmap_all(sim.d, discretization_level)
        
        saveRDS(sim.d, file =  paste0(dirW,characters[i], "_", j, ".rds") )
        
      }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        #errors<-rbind(errors, c(ii,jj))
      }  )
      
    } 
    
    # putting rds files into archive
    files<-paste0(dirW, characters[i], "_", c(1:length(sim)), ".rds")
    zip(paste0(dirW, characters[i], ".zip"), files=files)
    file.remove(files)
    
  }
  
  # close connections
  showConnections (all=T)
  closeAllConnections()
  
}

#' Function for aggregating characters under a specified set of terms (for example, body regions)
#' @export
RAC_query <- function(char_info, ONT, names){
  c <-  char_info$ID
  c <- gsub(" ", "_", c)
  
  annot <- as.list(as.character(char_info$IRI))
  names(annot) <- as.character(char_info$ID)
  ONT$terms_selected_id <- annot
  
  
  parts <- do.call(rbind,lapply(names,  pk_get_iri, as="uberon"))
  parts$IRI <- sapply(parts[,1], strip_IRI)
  levelA <- setNames(parts$IRI, names)     
  
  res <- lapply(levelA, function(x)
    get_descendants_chars(ONT, annotations="manual", terms=x)  )
  
  cat("\nAggregations by :\n")
  print(res)
  return(res)
}

#' Function to write Ontology CTMC models as a .Rev file for execution in Revbayes
#' @export
writeRevOntologyModels <- function(td, M, dir, dirW, dirR, dirD) {
  MT <- as.data.frame(td$dat)
  colnames(MT) <- gsub(" ", "_", colnames(MT))
  rownames(MT) <- td$phy$tip.label
  for (i in 1:ncol(MT)){
    C.rev<-MT[,i]
    C.rev<-gsub("&", " ", C.rev)
    o <- order(nchar(C.rev))
    
    out<-cbind(rownames(MT), C.rev)
    out <- out[o,]
    write.table(file=paste0(dirD, colnames(MT[i]), ".char"), out, quote=F, sep=" ", 
                row.names=F, col.names=F)
  }
  
  data(ParamoRevTemplate)
  
  
  for (i in 1:ncol(MT)){
    fl.in  <- ParamoRevTemplate
    fl.in  <- gsub(pattern = "Hymenoptera_br_resolved", replace = "fishtree",
                   x = fl.in)
    fl.in  <- gsub(pattern = "@analysis_name@", replace = paste0(colnames(MT[i])),
                   x = fl.in)
    
    fl.in <- gsub(pattern = "@chrs_2_read@", 
                  replace = paste0("data/", colnames(MT[i]), ".char"), x = fl.in)
    
    if(colnames(MT)[i] %in% gsub(" ", "_", names(M))){
      in.rev<-Mk_Rev(M[[gsub("_", " ", colnames(MT)[i])]])
      
      fl.in <- gsub(pattern = "@numstates@", 
                    replace = as.character(max(dim(M[[gsub("_", " ", colnames(MT)[i])]]))), x = fl.in)
      
      fl.in <- gsub(pattern = "@ratematrix@", 
                    replace = in.rev, x = fl.in)
    } else {
      
      fl.in <- gsub(pattern = "@numstates@", 
                    replace = "2", x = fl.in)
      
      fl.in <- gsub(pattern = "@ratematrix@", 
                    replace = ratemat1(), x = fl.in)
    }
    
    cat(file=paste0(dir, colnames(MT[i]), ".Rev"), sep="\n", fl.in)
  }
  
  
}



##############################################################################################
# FUNCTIONS
##############################################################################################

#' @title Combining two matrices
#' @description Combining two matrices. The parametric schem of matrice is defined by nattural
#' numbers; different numbers = different rate parameters
#' @param M1 matrix; if dependency true thenM1 controls M2
#' @param M2 matrix; if dependency true then: M2 depends on those states of M1 specified in dependent.state
#' @param dependent.state state(s) of M1 that switches on matrix M2 
#' @param name.sep separator for state names
#' @param diag.as hpopulate main diagonal with
#' @return Matrix
#' @examples
#' M1<-matrix(c(-1,1,  2,-2),2,2,byrow=TRUE)
#' rownames(M1)<-colnames(M1)<-c("0","1")
#' M2<-matrix(c(-3,3,  4,-4),2,2,byrow=TRUE)
#' rownames(M2)<-colnames(M2)<-c("0","1")
#' comb2matrices(M1, M2, dependent.state=NULL)
#' comb2matrices(M1, M2, dependent.state=2)
# if dependency true then: M2 depends on M1 states specified in dependent.state
comb2matrices<-function(M1,M2, dependent.state=NULL, name.sep="", diag.as=""){
  
  if (!is.null(dependent.state)){
    matrix.diag<-rep(0, ncol(M1))
    matrix.diag[dependent.state]<-1
    matrix.diag<-diag(matrix.diag)
  }
  
  if (is.null(dependent.state)){
    matrix.diag<-diag(nrow(M1))
  }
  
  M_kr=(M1%x%diag(nrow(M2)))+ (matrix.diag%x%M2)
  
  
  #getting colnames
  
  col=paste(colnames(kronecker(M1, diag(nrow(M2)), make.dimnames = T)),
            colnames(kronecker(diag(nrow(M1)), M2, make.dimnames = T)), sep="")
  col=gsub("::", name.sep, col, fixed=T)
  
  # merging two names
  rownames(M_kr)<-colnames(M_kr)<-col
  if (diag.as!="") diag(M_kr)<-diag.as
  
  return(M_kr)
}



#' @title Combining multiple matrices
#' @description Liely as independently evolving
#' @param list.matrices Lsit of matrices
#' @return Matrix
#' @examples
combNmatrices<-function(list.matrices,  ...){
  comb.matrix<-list.matrices[[1]]
  
  for (i in 1:(length(list.matrices)-1)){
    comb.matrix=comb2matrices(comb.matrix, list.matrices[[i+1]], dependent.state=NULL, diag.as = 0)
  }
  
  return(comb.matrix)
}

#' @title Initialize binary matrices given graph
#' @description Call matrices are populated with different parameters
#' @param graph igraph object
#' @return List of matrices
#' @examples
#' init_binary_matrix(g)
init_binary_matrix<-function(graph){
  matrix.list=list()
  n.matrix=vcount(graph)
  vertices=V(graph)$name
  
  param.scheme<-matrix(seq(1:(2* n.matrix )), ncol = 2, byrow=TRUE)
  for (i in 1:n.matrix){
    matrix.list[[vertices[i]]]<-matrix(c(-param.scheme[i,1],param.scheme[i,1], 
                                         param.scheme[i,2],-param.scheme[i,2]),2,2,byrow=TRUE)
    rownames(matrix.list[[vertices[i]]])<-colnames(matrix.list[[vertices[i]]])<-c("0","1")
  }
  
  return(matrix.list)
}


#' @title Get all dependency matrices given a dependecy graph
#' @description Construct dependency matrices and their correponding attributes
#' @param graph igraph object of ontology terms
#' @return List of matrices and their attributes
#' @examples
#' get_graph_matrix(g)

########## Structure of the List
#' $binary.matrices # intial binary matrices assigned to each node of graph
#' $comb.matrices$matrix # combined matrix for each node
#' $comb.matrices$state.string # vector of states [1] "00" "01" "10" "11"
#' $comb.matrices$state.ident # specifies the order of ontology terms in each state [1] "UBERON:0007829" "UBERON:2000663"
#' $comb.matrices$state.observable # ids and names of "observable" states. In red-blue tail notation refers to blue and red states
#' $comb.matrices$state.hidden # ids and names of "hidden" states. In red-blue tail notation refers to "blue absent"
# and "red absent"

#' $nodes.sorted # topologically sorted nodes
#' $vertex.hier # hierrachy of the nodes
###########

get_graph_matrix<-function(graph){
  
  g=graph  
  
  # dependent chars object
  complex.char<-list()
  complex.char$binary.matrices<-init_binary_matrix(g)
  complex.char$comb.matrices<-list()
  
  # traverse graph
  .topo <- topo_sort(g, mode = c("out"))
  topo= names(.topo)
  complex.char$nodes.sorted=topo
  vertex.hier=ego(g, order=1, nodes = topo, mode = c("in"), mindist = 1)
  names(vertex.hier)<-topo
  complex.char$vertex.hier=vertex.hier
  
  
  for (i in seq_along(topo)){
    focal.v=complex.char$vertex.hier[[i]]
    
    if (length(focal.v)==0){
      complex.char$comb.matrices[[topo[i]]]$matrix=complex.char$binary.matrices[[topo[i]]]
      complex.char$comb.matrices[[topo[i]]]$state.string=row.names(complex.char$binary.matrices[[topo[i]]])
      complex.char$comb.matrices[[topo[i]]]$state.ident=topo[i]
      #complex.char$comb.matrices[[topo[i]]]$dependency.true=2
      complex.char$comb.matrices[[topo[i]]]$state.observable=integer(0)
      complex.char$comb.matrices[[topo[i]]]$state.hidden=integer(0)
    }
    
    if (length(focal.v)>0){
      
      if (length(focal.v)==1){ # if length =1 the dependency is chain like
        MC=complex.char$comb.matrices[[names(focal.v)]]$matrix
        M=complex.char$binary.matrices[[topo[i]]]
        dps=ncol(MC)
        cmb=comb2matrices(MC, M, dependent.state=dps, diag.as = 0)
      }
      
      
      if (length(focal.v)>1){
        # sequentially combine multiple matrices as independently coevolving
        list.matrices=lapply(names(focal.v), function(x) complex.char$comb.matrices[[x]]$matrix)
        #names(list.matrices)=names(focal.v)
        #list.matrices=list.matrices[1:2]
        comb.mt=combNmatrices(list.matrices)
        
        # combine matrices  from above with the focal node matrix;
        # the state bearing dependency is where all entities=1, i.e. the last state
        M=complex.char$binary.matrices[[topo[i]]]
        dps=ncol(comb.mt)
        cmb=comb2matrices(comb.mt, M, dependent.state=dps, diag.as = 0)
      }
      
      # adding attributes
      complex.char$comb.matrices[[topo[i]]]$matrix=cmb
      r.name <- row.names(cmb)
      complex.char$comb.matrices[[topo[i]]]$state.string<-r.name #<-cmb %>% row.names()
      
      st.iden=unlist(lapply(names(focal.v), function(x) complex.char$comb.matrices[[x]]$state.ident))
      complex.char$comb.matrices[[topo[i]]]$state.ident<-st.iden<-c(st.iden, topo[i])
      
      # get observable and "hidden" states of the focal node
      ln=length(st.iden)-1 # observable a/p are only those which have all states from other chrs=1
      obs=which(substr(r.name, 1, ln)==paste(rep(1, ln), collapse=""))
      names(obs)=r.name[obs]
      complex.char$comb.matrices[[topo[i]]]$state.observable=obs
      
      hid=(1:length(r.name))[-obs]
      names(hid)<-r.name[hid]
      complex.char$comb.matrices[[topo[i]]]$state.hidden=hid
      
      
      
    } #end if (length(focal.v)>0)
  } #end all
  
  return(complex.char)
} # end function


#################
# END FUNCTIONS
################################################################################################################################



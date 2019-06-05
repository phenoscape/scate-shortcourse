
print_coverage <- function(x){
  coverage <- apply(x, 2, function(x) sum(!is.na(x)))
  average <- sapply(x, function(x) mean(as.numeric(na.omit(x[x %in% c("0", "1")])), na.rm=TRUE))
  cover <- cbind(coverage, average)
  kableExtra::kable(filter(data.frame(traits=rownames(cover), cover), coverage > 0, average < 1, average > 0) %>% arrange(., desc(coverage)), format = "markdown")
}


strip_IRI <- function(x){
  x <- gsub("http://purl.obolibrary.org/obo/", "", x)
  x <- gsub("_", ":", x)
  return(x)
}





filter_coverage <- function(td, traits=0, taxa=0){
  taxa_coverage <- apply(td$dat, 1, function(x) mean(as.numeric(!is.na(x))))
  trait_coverage <- apply(td$dat, 2, function(x) mean(as.numeric(!is.na(x))))
  td <- dplyr::filter(td, taxa_coverage > taxa)
  td <- dplyr::select(td, which(trait_coverage > traits))
  return(td)
}

plotData <- function(td, njt, start=3, margs=c(0.2, 0.25), ...){
  vals <- na.omit(unique(do.call(base::c, lapply(start:ncol(td$dat), function(x) unique(as.character(td$dat[[x]]))))))
  X <- do.call(cbind, lapply(start:ncol(td$dat), function(x) as.numeric(recode(td$dat[[x]], "0 and 1"=0, "1 and 0"=0, "1"=1, "0"=0, "2"=2, "3"=3))))
  colnames(X) <- colnames(td$dat)[start:ncol(td$dat)]
  X <- X[,njt$edge[njt$edge[,2] <= length(njt$tip.label),2]]
  .vals <- sort(na.omit(unique(as.vector(X))))
  dimx <- dim(X)
  
  tree <- njt
  tree1 <- chronopl(td$phy,1)
  tree2 <- chronopl(njt, 1)
  
  #alters the length of the tips of the trees
  tree1$edge.length <- tree1$edge.length/(max(branching.times(tree1)))*margs[1]*dimx[2]
  tree2$edge.length <- tree2$edge.length/(max(branching.times(tree2)))*margs[2]*dimx[1]
  
  #Changes the direction of the top plot
  h1 <- plot(tree1, plot = FALSE, cex=0.5, show.tip.label=FALSE)
  h2 <- plot(tree2, plot = FALSE, direction = "downwards", show.tip.label=FALSE)
  
  # this is all setting boundaries for the different plots and combining them into one image
  par(mar = c(0,0,0,0))
  plot(0,0, type = 'n', xlim = c(0,h1$x.lim[2]+h2$x.lim[2]), ylim=c(0,h1$y.lim[2]+h2$y.lim[2]))
  
  imgcol <- (brewer.pal(max(c(length(.vals),3)), "YlOrRd"))
  
  image(seq(h1$x.lim[2]+1,h1$x.lim[2]+h2$x.lim[2], length.out=ncol(X)), seq(1, h1$y.lim[2], length.out=nrow(X)), t(X),xlim=c(1+h1$x.lim[2],h1$x.lim[2]+h2$x.lim[2]+1) ,ylim=c(0, h1$y.lim[2]-1), add=TRUE, col=imgcol)
  
  legend(0, (h1$y.lim[2]+h2$y.lim[2])*.99, legend=.vals ,pch=22, pt.bg=imgcol)
  
  par(new = TRUE)
  
  
  plot(tree1, x.lim=c(0,(1+margs[2])*(h2$x.lim[2]+h1$x.lim[1])), y.lim=c(0,h1$y.lim[2]+h2$y.lim[2]), show.tip.label=FALSE)
  
  par(new = TRUE)
  plot(tree2, direction = "downwards", x.lim=c(-h1$x.lim[2],h2$x.lim[2]), y.lim=c((-h1$y.lim[2])-0.01*dimx[1],h2$y.lim[2]), ...)
  
  return ()
  
}

makeTree <- function (td, skip=1:2)
{
  
  
  traits <- colnames(td$dat)
  if(length(skip)>0){
    traits <- traits[-(skip)] #delete otu data
  }
  traits
  
  
  
  #Get IRI ids for each trait.
  
  
  traitDetails <- lapply(traits, function(x) pk_anatomical_detail(x, verbose=TRUE))
  
  
  
  traitDetails[1:5]
  traitIDs <- unname(do.call(base::c, sapply(traitDetails, function(x) x[,'@id'])))
  
  
  
  irisPhenotypes <- sapply(traitIDs, url_encode)
  
  
  
  filename <- paste0(tempdir(), "/a.txt")
  
  cat("iris=%5B%0A%20%20", file=filename)
  irisPhenotypes <- lapply(irisPhenotypes, function(x) gsub("/", "%2F", x, fixed=TRUE))
  irisPhenotypes <- lapply(irisPhenotypes, function(x) gsub(":", "%3A", x, fixed=TRUE))
  irisPhenotypes <- lapply(irisPhenotypes, function(x) gsub("=", "%3D%0A", x, fixed=TRUE))
  dum <- lapply(irisPhenotypes[1:(length(irisPhenotypes)-1)],function(x) cat(paste0('%22', x,'%22%2C', sep=""), file=filename, append=TRUE))
  cat(paste0('%22', irisPhenotypes[[length(irisPhenotypes)]],'%22',"%5D%0A", sep=""), file=filename, append=TRUE)
  
  
  
  
  
  api.semanticSimilarity_query <- paste0("curl -X POST -d @", filename," 'https://kb.phenoscape.org/api/similarity/jaccard'")
  semanticSimilarityAPIResults <- system(api.semanticSimilarity_query, intern=TRUE)
  
  
  
  
  
  results <- fromJSON(semanticSimilarityAPIResults)
  scores <- lapply(results$results, function(x) x$score)
  scores <- sapply(scores, function(x) if(is.null(x)) NA else(x))
  result_terms <- do.call(rbind, lapply(results$results, function(x) do.call(cbind, lapply(x$terms, curlUnescape))))
  semanticSimilarityMatrix <- matrix(NA, nrow=length(irisPhenotypes), ncol=length(irisPhenotypes))
  diag(semanticSimilarityMatrix) <- 1
  rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- curlUnescape(irisPhenotypes)
  
  for(i in 1:nrow(result_terms)){
    semanticSimilarityMatrix[result_terms[i,1], result_terms[i,2]] <- semanticSimilarityMatrix[result_terms[i,2], result_terms[i,1]] <- scores[i]
  }
  
  rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- traits
  
  write.csv(semanticSimilarityMatrix, file="siluriformesSemanticSimMatrix.csv")
  
  
  
  #Check to see if semantic similarity matrix makes sense.
  
  
  maxSS <- list()
  for(i in 1:ncol(semanticSimilarityMatrix)){
    ss <- semanticSimilarityMatrix[-i,]
    j <- which(ss[,i]==max(ss[,i]))[1]
    maxSS[[i]] <- cbind(colnames(semanticSimilarityMatrix)[i], rownames(ss)[j], round(ss[j, i],4))
  }
  maxSS <- do.call(rbind, maxSS)
  as.data.frame(maxSS)
  
  
  minSS <- list()
  for(i in 1:ncol(semanticSimilarityMatrix)){
    ss <- semanticSimilarityMatrix[-i,]
    j <- which(ss[,i]==min(ss[,i]))[1]
    minSS[[i]] <- cbind(colnames(semanticSimilarityMatrix)[i], rownames(ss)[j], round(ss[j, i],4))
  }
  minSS <- do.call(rbind, minSS)
  as.data.frame(minSS)
  
  
  
  #Neighbor-joining tree of SS matrix
  
  njt <- nj(1-semanticSimilarityMatrix)
  pdf("njTreeSiluriformesSemanticMatrix.pdf", height=30, width=30)
  plot(njt, type="unrooted", cex=0.35)
  dev.off()
  
  return (njt)
}
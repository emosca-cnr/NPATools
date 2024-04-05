#' Find topological communities
#' 
#' @description Find topological communities
#' @param g graph
#' @param eWeights edge weights
#' @param vWeights vertex weights
#' @param verbose TRUE/FALSE
#' @param methods one or more of "fastgreedy", "labprop", "walktrap", 
#' "eigen", "multilev", "infomap"
#' @return list od community objects
#' @importFrom igraph fastgreedy.community modularity label.propagation.community walktrap.community leading.eigenvector.community multilevel.community infomap.community
#' @export
find_communities <- function(g=NULL, eWeights=NULL, vWeights=NULL, 
                             verbose=TRUE, methods=c("fastgreedy", "multilev")){

    commList <- setNames(vector('list', length = length(methods)), methods)
    
    commInfo <- data.frame(algorithm=methods,
                            modularity=NA,
                            n=NA,
                            stringsAsFactors=FALSE)

    if("fastgreedy" %in% methods){
        if(verbose)
            print("fastgreedy")
        idx <- which(grepl("fastgreedy", methods))
        commList[[idx]] <- fastgreedy.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            modularity(g, commList[[idx]]$membership,
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
    if("labprop" %in% methods){
        if(verbose)
            print("labprop")
        idx <- which(grepl("labprop", methods))
        commList[[idx]] <- 
            label.propagation.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
    if("walktrap" %in% methods){
        if(verbose)
            print("walktrap")
        idx <- which(grepl("walktrap", methods))
        commList[[idx]] <- walktrap.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
    if("eigen" %in% methods){
        if(verbose)
            print("eigen")
        idx <- which(grepl("eigen", methods))
        commList[[idx]] <- 
            leading.eigenvector.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
    if("multilev" %in% methods){
        if(verbose)
            print('multilev')
        idx <- which(grepl("multilev", methods))
        commList[[idx]] <-
            multilevel.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
  if("infomap" %in% methods){
        if(verbose)
            print('infomap')
        idx <- which(grepl("infomap", methods))
        commList[[idx]] <- 
            infomap.community(g, e.weights = eWeights, v.weights = vWeights)
        commInfo$modularity[idx] <- 
            modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
  }

  return(list(comm=commList, info=commInfo))

}

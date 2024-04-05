#' Create a network of communities
#'
#' @param g igraph object, vertex attribute "comm_id" is mandatory
#' @param remove.multiple whether to remove multiple links among communities
#' @description Function for reduce the network as a community network
#' @return igraph object
#' @export
#' @importFrom igraph V E V<- E<- get.vertex.attribute remove.vertex.attribute contract.vertices as_ids get.edge.attribute simplify
#' 
comm_net <- function(g=NULL, remove.multiple=FALSE) {
    
    #set weights to count vertices and edges in the community
    V(g)$w <- 1
    E(g)$w <- 1
    
    comm_id <- V(g)$comm_id
    v_attr <- names(get.vertex.attribute(g))
    if(any(!v_attr %in% "w")){
        v_to_rm <-  v_attr[v_attr != "w"]
        for(i in 1:length(v_to_rm)){
            g <- remove.vertex.attribute(g, name = v_to_rm[i])
        }
    }
    commNet <- contract.vertices(graph = g, mapping = comm_id, vertex.attr.comb = "sum")
    
    V(commNet)$name <- V(commNet)$label <- as_ids(V(commNet))
    V(commNet)$size <- seq(10, 20, length.out=5)[cut(V(commNet)$w, 5)]
    
    e_attr <- names(get.edge.attribute(commNet))
    if(any(!e_attr %in% "w")){
        e_to_rm <-  e_attr[!e_attr %in% c("w")]
        for(i in 1:length(e_to_rm)){
            commNet <- remove.edge.attribute(commNet, name = e_to_rm[i])
        }
    }
    commNet <- simplify(graph = commNet, remove.multiple = remove.multiple, remove.loops = T, edge.attr.comb = "sum")
    if(remove.multiple){
        E(commNet)$width <- c(1:10)[cut(E(commNet)$w, 5)]
    }
    
    return(commNet)
    
}

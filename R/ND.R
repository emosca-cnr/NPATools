#' Network Diffusion 
#' 
#' @param X0 vector or matrix composed of column vectors with initial 
#' distribution of information
#' @param W symmetrically normalized adjacency matrix W = D^-1 A D^-1, 
#' see normalize_adj_mat function
#' @param alpha numeric, the smothing factor
#' @param nMax numeric, maximum number of iterations
#' @param eps numeric, the iteration will stop when the maximum difference 
#' between matrix Xs between two consecutive iteraction is 
#' smaller than \code{eps}
#' @param finalSmooth TRUE/FALSE, whether to do the final step of smoothing
#' @param fullOutput, TRUE/FALSE, whether to output all steps
#' @param verbose, TRUE/FALSE
#' @usage ND(X0, W, alpha=0.7, nMax=1e4, eps=1e-6, finalSmooth=FALSE, 
#' fullOutput=FALSE, verbose=FALSE)
#' @examples 
#' \dontrun{ND(X0, W, alpha=0.7, nMax=1e4, eps=1e-6, finalSmooth=FALSE, 
#' fullOutput=FALSE, verbose=FALSE)}
#' @return A matrix with steady state values. If fullOutput is TRUE, a list with:
#' \item{\code{Xs}}{ the smoothed matrix;}
#' \item{\code{eps}}{ see above;}
#' \item{\code{maxAbsDiff}}{ maximum absolute difference between F(t) and F(t+1);}
#' \item{\code{XsAll}}{ transient Xs matrices.}
#' @export
#' @importFrom methods is
#' @importFrom Matrix Matrix

ND <- function(X0=NULL, W=NULL, alpha=0.7, nMax=1e4, eps=1e-6, finalSmooth=FALSE, 
               fullOutput=FALSE, verbose=FALSE){

    if(!identical(rownames(X0), rownames(W))){
        cat("Row names of W and X0 are not identical, trying to match...\n")
        X0 <- X0[match(rownames(W), rownames(X0)), ]
        cat("nrows of X0 and W:", nrow(X0), "\n")
    }
    
    if(!is(X0, "dgCMatrix")){
        X0 <- Matrix(X0, sparse = TRUE)
    }
    if(!is(W, "dgCMatrix")){
        W <- Matrix(W, sparse = TRUE)
    }
    

    Xs <- X0
    Fprev <- X0

    if(fullOutput){
        XsAll <- list()
        XsAll[[1]] <- X0
    }

    #X0 and A multiplied by their weight
    X0a <- (1 - alpha) * X0
    Wa <- alpha * W

    for(i in 2:nMax){

        if(i %% 5 == 0 & verbose) {
            cat(i, " ")            
        }

        #current iteration
        Xs <- Wa %*%  Fprev + X0a

        if(fullOutput){
            XsAll[[i]] <- Xs
        }

        maxAbsDiff <- max(abs(Xs-Fprev))

        #update Fprev for next iteration
        Fprev <- Xs

        if(maxAbsDiff < eps){
            if(finalSmooth){
                Xs <- Wa %*% Fprev
            }
            break
        }
    }
    
    if(verbose)
        cat('\n')

    if(fullOutput){
        
        return(list(Xs=Xs, maxAbsDiff=maxAbsDiff, XsAll=XsAll))
        
    }else{
        
        return(Xs)
    }
}




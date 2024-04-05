#' assess enrichment of the top networks
#' @param G igraph object
#' @param topList ranked list of vertex names that will be used to define top networks
#' @param ranks ranks of topList that will be assessed
#' @param type gsea or ora
#' @param X0Vector named numeric vector that will be tested with GSEA or ORA. In case of GSEA it will be ranked by decreasing orderg, while in the case of ORA the names of the X0 values grater than 0 will be tested for enrichment, while all X0 names will be the universe.
#' @param k number of permutations
#' @param minComponentSize size of the smallest graph component
#' @param minNetSize minimum network size
#' @param minKNes minimum k for considering a NES "confident"
#' @param BPPARAMGsl BiocParallelParam instance to parallelize over gene sets. See `BiocParallel::bplapply()`
#' @param BPPARAMK BiocParallelParam instance to paralleliz over permutations. `BiocParallel::bplapply()`
#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom igraph V induced_subgraph
#' @importFrom stats setNames rnorm

assess_enrichment <-
  function(G = NULL,
           topList = NULL,
           ranks = NULL,
           X0Vector = NULL,
           type = c("gsea", "ora"),
           k = 99,
           minComponentSize = 2,
           minNetSize = 10,
           minKNes = 10,
           BPPARAMGsl = NULL,
           BPPARAMK = NULL) {
    
    type <- match.arg(type, c("gsea", "ora"))
    
    if (is.null(ranks)) {
      ranks <- seq(10, length(topList), by = 10)
    }
    en_res <- setNames(vector("list", length(ranks)), ranks)
    net_ccf <- en_res
    
    if (is.null(BPPARAMGsl)) {
      BPPARAMGsl <- SerialParam()
    }
    
    if (is.null(BPPARAMK)) {
      BPPARAMK <- SerialParam()
    }
    
    
    
    #defines the increasing top networks to be tested
    for (i in 1:length(ranks)) {
      net_ccf[[i]] <-
        get_nconn_comp(induced_subgraph(G, V(G)$name[V(G)$name %in% topList[1:ranks[i]]]), n =
                         minComponentSize)
    }
    
    gsl <- lapply(net_ccf, function(i_net)
      V(i_net)$name)
    gsl <- gsl[lengths(gsl) >= minNetSize]
    net_ccf <- net_ccf[names(net_ccf) %in% names(gsl)]
    
    if (length(gsl) > 0) {
      #check X0 to choose ORA or GSEA
      if (type == "ora") {
        cat("ORA\n")
        
        wb <- names(X0Vector)[X0Vector > 0]
        if (length(wb) > 0) {
          #bb <- names(X0Vector)[!names(X0Vector) %in% wb]
          en_res <- ora(
            wb = list(wb),
            universe = names(X0Vector),
            gsl = gsl,
            out_file_prefix = NULL
          )
          en_res <- en_res[[1]][, 1:11]
        } else{
          stop("Can' find any element of X0Vector > 0\n.")
        }
        
      } else{
        cat("GSEA\n")
        
        #if there are zeros... a noise is added to avoid ties in the ranked gene list
        if (sum(X0Vector == 0) > 1) {
          cat("Detected multiple null values in X0. Adding a noise to avoid ties.\n")
          idx_zeros <- which(X0Vector == 0)
          min_pos <- min(X0Vector[-idx_zeros])
          X0Vector[idx_zeros] <-
            abs(rnorm(
              n = length(idx_zeros),
              mean = 0,
              sd = min_pos * 10 ^ -6
            ))
          print(summary(X0Vector[-idx_zeros]))
          print(summary(X0Vector[idx_zeros]))
        }
        
        en_res <-
          gsea(
            rl = matrix(
              X0Vector,
              ncol = 1,
              dimnames = list(names(X0Vector), "X0")
            ),
            gsl = gsl,
            k = k,
            min.k = minKNes,
            BPPARAMGsl = BPPARAMGsl,
            BPPARAMK = BPPARAMK,
            out_file_prefix = NULL
          )
        #en_res <- en_res$X0
        en_res <- en_res[[1]][, c(1:9, 11:15)]
      }
      
      colnames(en_res) <- replace(colnames(en_res), colnames(en_res) == "id", "rank")
      
      return(list(en_summary = en_res, net_ccf = net_ccf))
      
    } else{
      cat("No networks with at least ", minNetSize, ".\n")
      return(NULL)
    }
  }

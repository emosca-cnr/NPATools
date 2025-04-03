#' Get KEGG Pathways using KEGGREST
#' @param org organism, following KEGG terminology
#' @param mapping data.frame with two columns: NCBI Entrez gene_id and symbol
#' @param out_file output file name
#' @export
#' @importFrom KEGGREST keggList keggLink keggConv keggInfo

get_kegg_db <- function(org = "hsa", mapping = NULL, out_file = "kegg_pathway.rds") {
  
  writeLines(keggInfo("pathway"))
  
  kl <- keggList("pathway", org)
  klink <- keggLink(org, "pathway")
  kconv <- keggConv(org, "ncbi-geneid")
  
  kegg_db <- data.frame(pathway_id = names(klink), kegg_gene = as.vector(klink), stringsAsFactors = F)
  kegg_db$pathway_id <- gsub("path:", "", kegg_db$pathway_id)
  
  kegg_desc <- data.frame(pathway_id = names(kl), desc = as.vector(kl), stringsAsFactors = F)
  #kegg_desc$desc <- gsub(" - Homo sapiens (human)", "", kegg_desc$desc, fixed = T)
  kegg_gene <- data.frame(kegg_gene = as.vector(kconv), entrez_id = names(kconv), stringsAsFactors = F)
  kegg_gene$entrez_id <- gsub("ncbi-geneid:", "", kegg_gene$entrez_id)
  
  kegg_db <- merge(kegg_db, kegg_desc)
  kegg_db <- merge(kegg_db, kegg_gene)
  
  eg2sym <- unique(eg2sym[, c(1, 2)])
  kegg_db <- merge(kegg_db, eg2sym, by.x = "entrez_id", by.y = "gene_id", all.x = T)
  
  if(!is.null(out_file)){
    saveRDS(kegg_db, out_file, compress = "bzip2")
  }
  
  return(kegg_db)
}

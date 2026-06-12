suppressPackageStartupMessages({library(tidyverse); library(fgsea); library(msigdbr)})
run_functional_enrichment <- function(){
 dir.create("results/enrichment",recursive=TRUE,showWarnings=FALSE); files <- list.files("results/deseq2", pattern="DEGs_.*csv$", full.names=TRUE); hallmark <- msigdbr(species="Homo sapiens", category="H") |> split(x=.$gene_symbol, f=.$gs_name)
 for(f in files){ tab <- read.csv(f); if(!all(c("gene_id","stat") %in% colnames(tab))) next; ranks <- tab$stat; names(ranks) <- tab$gene_id; ranks <- sort(na.omit(ranks), decreasing=TRUE); fg <- fgsea(pathways=hallmark, stats=ranks, minSize=10, maxSize=500); write.csv(fg, file.path("results/enrichment", paste0(tools::file_path_sans_ext(basename(f)),"_Hallmark_fgsea.csv")), row.names=FALSE) }
}

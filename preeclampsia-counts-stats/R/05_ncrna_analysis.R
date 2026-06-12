source("R/00_setup.R")
run_ncrna_analysis <- function(annotation_path=NULL){
 message("ncRNA module: provide annotation_path with gene_id and gene_biotype to activate biotype-specific outputs."); if(is.null(annotation_path) || !file.exists(annotation_path)) return(invisible(NULL))
 ann <- read.csv(annotation_path); dir.create("results/deseq2/ncRNA", showWarnings=FALSE); files <- list.files("results/deseq2", pattern="DEGs_.*csv$", full.names=TRUE)
 for(f in files){ tab <- read.csv(f); out <- left_join(tab,ann,by="gene_id") |> filter(grepl("lncRNA|miRNA|antisense|non_coding|pseudogene", gene_biotype, ignore.case=TRUE)); write.csv(out,file.path("results/deseq2/ncRNA",basename(f)),row.names=FALSE) }
}

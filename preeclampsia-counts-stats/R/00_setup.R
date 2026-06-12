suppressPackageStartupMessages({library(tidyverse); library(data.table); library(DESeq2); library(pheatmap); library(ggrepel); library(yaml)})
load_config <- function(path="config/analysis_config.yaml") yaml::read_yaml(path)
make_dirs <- function(){ dirs <- c("results/tables","results/deseq2","results/enrichment","figures/qc","figures/volcano","figures/heatmaps","figures/enrichment"); invisible(lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE)) }
read_counts <- function(path){ x <- data.table::fread(path) |> as.data.frame(); rownames(x) <- x[[1]]; x <- x[,-1,drop=FALSE]; round(as.matrix(x)) }
read_metadata <- function(path) read.csv(path, check.names=FALSE)
filter_low_counts <- function(counts, min_count=10, min_samples=3) counts[rowSums(counts >= min_count) >= min_samples,,drop=FALSE]

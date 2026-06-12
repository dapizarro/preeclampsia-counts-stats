#!/usr/bin/env Rscript
source("R/01_import_and_qc.R"); source("R/02_deseq2_by_placenta.R"); source("R/03_maternal_vs_fetal.R"); source("R/04_interaction_model.R"); source("R/05_ncrna_analysis.R"); source("R/06_isoform_analysis.R"); source("R/07_functional_enrichment.R"); source("R/08_gradient_analysis.R"); source("R/09_make_report.R")
run_import_qc(); run_deseq2_by_placenta(); run_maternal_vs_fetal(); run_interaction_model(); run_ncrna_analysis(); run_isoform_analysis(); run_functional_enrichment(); run_gradient_analysis(); run_make_report()

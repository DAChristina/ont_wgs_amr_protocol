# Data is updated per-batch
# Columns MUST NOT BE MERGED

library(tidyverse)
source_all <- function(path = "global"){
  files <- list.files(path,
                      pattern = "\\.R$",
                      full.names = TRUE
  )
  invisible(lapply(files, source))
}

source_all()

# Batch list:
# Batch 1 (22 April 2026)
# Batch 2 (xx May 2026)


current_batch <- 1
batch_prep <- prepare_batch(
  df1 = "raw_data/ACORN_HAI_2026_log_book_batch1_22APRIL2026.xlsx",
  df2 = "raw_data/demux_batch1.csv",
  batch = current_batch,
  batch_address = "/srv/nfs_share/2026_ACORNHAI/raw_data/ACORNHAI_3_22APRIL2026/ACORNHAI_3_22APRIL2026/20260422_1507_X4_FBD04135_d247130b/output_demux/"
  ) %>% 
  glimpse()

# prepare_batch automatically creates GHRU-mod samplesheet & compressed target list
# run GHRU-mod pipeline for getting final tables.

# compile QC workLab & workSeq
batch_qc <- compile_qc_batch(
  df1 = batch_prep,
  # TEMPORARY read.csv data on my local PC
  # df2 = "raw_data/ACORNHAI_3_22APRIL2026/ACORNHAI_3_22APRIL2026/20260422_1507_X4_FBD04135_d247130b/output_ghru-mod_lambda_filter_off/final_tables/combined.post.csv"
  df2 = "raw_data/output_ghru-mod_lambda_filter_off/final_tables/combined.post.csv",
  batch = current_batch
) %>% 
  glimpse()

# Compiled QC workLab & workSeq automatically saved in:
# "inputs/worklab_workSeq_compiled_batch<batch>.csv"

################################################################################
# compile batches by using bind_rows
compile_batches <- function(path = "inputs"){
  files <- list.files(path,
                      pattern = "worklab_workSeq_compiled_.*\\.csv$",
                      full.names = TRUE
  )
  
  return(files %>% 
           lapply(read.csv) %>% 
           dplyr::bind_rows()
  )
}

compiled_batches <- compile_batches() %>% 
  glimpse()


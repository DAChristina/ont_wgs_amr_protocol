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
# Batch 2 (06 May 2026)
# Batch 3 (20 May 2026)
# Batch 4 FAILED due to insufficient active pore
# Batch 5 (29 June 2026)
# Batch 6 (06 July 2026)


current_batch <- 6
ont_address <- "ACORNHAI_BATCH6_12_06072026/ACORNHAI_BATCH6_12_06072026/20260706_1621_X4_FBD01789_7551ff55/"

# use ls <folder>/* > raw_data/demux_batch<batch>.csv
batch_prep <- prepare_batch(
  df1 = "raw_data/ACORN_HAI_2026_log_book_batch6_06JULY2026.xlsx",
  df2 = "raw_data/demux_batch6.csv",
  batch = current_batch,
  batch_address = paste0("/srv/nfs_share/2026_ACORNHAI/raw_data/", ont_address, "output_demux/")
) %>% 
  glimpse()

# prepare_batch automatically creates GHRU-mod samplesheet & compressed target list
# run GHRU-mod pipeline for getting final tables.

# compile QC workLab & workSeq
batch_qc <- compile_qc_batch(
  df1 = batch_prep,
  # TEMPORARY read.csv data on my local PC
  df2 = paste0("raw_data/", ont_address, "output_ghru-mod_lambda_filter_on/final_tables/combined.post.csv"),
  # df2 = "raw_data/output_ghru-mod_lambda_filter_off/final_tables/combined.post.csv",
  batch = current_batch
) %>% 
  glimpse()

# Compiled QC workLab & workSeq automatically saved in:
# "inputs/worklab_workSeq_compiled_batch<batch>.csv"

################################################################################
# compile batches by using bind_rows
workLab_workSeq_all <- compile_csvfiles() %>% 
  # generate personal id
  dplyr::mutate(
    id = paste0(workLab_isolat_id, "_batch", workLab_batch, "_", workLab_native_barcode)
  ) %>% 
  glimpse()

write.csv(workLab_workSeq_all,
          "inputs/workLab_workSeq_compiled_all.csv",
          row.names = FALSE)

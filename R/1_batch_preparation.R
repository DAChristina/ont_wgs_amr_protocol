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
# Batch 7 (14 July 2026)


current_batch <- 7
ont_address <- "ACORNHAI_BATCH7_7_140726/ACORNHAI_BATCH7_7_140726/20260714_1518_X5_FBD01789_ae6a5c4d/"

# use ls <folder>/* > raw_data/demux_batch<batch>.csv
batch_prep <- prepare_batch(
  df1 = "raw_data/ACORN_HAI_2026_log_book_batch7_14JULY2026.xlsx",
  df2 = "raw_data/demux_batch7.csv",
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
# compile batches by using left_join & bind_rows
compiled_labSeq <- compile_csvfiles() %>% 
  glimpse()

workLab_workSeq_all <- dplyr::left_join(
  
  # received isolates from RS (duplication filtered)
  readxl::read_excel("raw_data/update acorn hai_isolasi DNA.xlsx",
                     sheet = "per June 26",
                     skip = 1) %>% 
    dplyr::select(2:9) %>% 
    janitor::clean_names() %>% 
    dplyr::rename_all(~ paste0("received_", .)) %>% 
    # fix species name
    dplyr::mutate(
      received_isolate = case_when(
        stringr::str_detect(received_isolate,
                            regex("KPN|K\\.\\s*pneumoniae|Klebsiella\\s+pneumoniae",
                                  ignore_case = TRUE)) ~ "Klebsiella pneumoniae",
        stringr::str_detect(received_isolate,
                            regex("E\\.\\s*coli|E\\.coli|Escherichia\\s+coli",
                                  ignore_case = TRUE)) ~ "Escherichia coli",
        stringr::str_detect(received_isolate,
                            regex("A\\.\\s*baumannii|A\\.baumannii|A\\.\\s*baumanii|A\\.baumanii|Acinetobacter\\s+baumannii",
                                  ignore_case = TRUE)) ~ "Acinetobacter baumannii",
        stringr::str_detect(received_isolate,
                            regex("P\\.\\s*aeruginosa|P\\.aeruginosa|P\\.\\s*aeroginosa|P\\.aeroginosa|Pseudomonas\\s+aeruginosa",
                                  ignore_case = TRUE)) ~ "Pseudomonas aeruginosa",
        TRUE ~ received_isolate
      )
    ) %>% 
    dplyr::mutate(dplyr::across(everything(),
                                as.character)) %>% 
    dplyr::filter(!duplicated(
      received_isolate_id,
      fromLast = TRUE))
  ,
  # worklab_workSeq_compiled_batch
  compiled_labSeq
  ,
  by = c("received_isolate_id" = "workLab_isolat_id")
) %>% 
  
  # adjust QC based on tolerable isolates
  dplyr::mutate(
    # adjust depth cutoff as 30x (yellow) 50x (green)
    workSeq_Depth.all_checks_passed_post = case_when(
      as.numeric(workSeq_Depth.Depth_post) < 30 ~ "False",
      as.numeric(workSeq_Depth.Depth_post) >= 30 & 
        as.numeric(workSeq_Depth.Depth_post) <= 50 ~ "Tolerated",
      as.numeric(workSeq_Depth.Depth_post) > 50 ~ "True",
      TRUE ~ workSeq_Depth.all_checks_passed_post
    ),
    # adjust completeness & contamination cutoff to 95% & 5% (up to 8% tolerance)
    # update: I don't think I can tolerate 8% contamination;
    # people use ~1-2% contamination cutoff instead of 5%
    workSeq_Checkm.Contamination.check_post = ifelse(
      as.numeric(workSeq_Checkm.Contamination_post) > 5, "False",
      workSeq_Checkm.Contamination.check_post
    ),
    workSeq_Checkm.all_checks_passed_post = ifelse(
      (as.numeric(workSeq_Checkm.Completeness_post) < 90 | 
         as.numeric(workSeq_Checkm.Contamination_post) > 5), "False",
      workSeq_Depth.all_checks_passed_post
    ),
    # set all checks
    workSeq_all_checks_passed_post = case_when(
      (as.numeric(workSeq_Sylph.top_adjusted_ani_post) >= 95 &
         workSeq_Checkm.all_checks_passed_post == "True" &
         workSeq_Depth.all_checks_passed_post == "True" &
         workSeq_Quast.all_checks_passed_post == "True" &
         workSeq_Speciator.all_checks_passed_post == "True"
      ) ~ "True",
      # adjust QC for tolerated isolates
      (as.numeric(workSeq_Sylph.top_adjusted_ani_post) >= 95 &
         workSeq_Checkm.all_checks_passed_post == "Tolerated" |
         workSeq_Depth.all_checks_passed_post == "Tolerated" |
         workSeq_Quast.all_checks_passed_post == "Tolerated" |
         workSeq_Speciator.all_checks_passed_post == "Tolerated"
      ) ~ "Tolerated",
      TRUE ~ "False"
    ),
  ) %>% 
  
  # generate personal id
  dplyr::mutate(
    id = paste0(received_isolate_id, "_batch", workLab_batch, "_", workLab_native_barcode)
  ) %>%
  glimpse()

write.csv(workLab_workSeq_all,
          "inputs/workLab_workSeq_compiled_all.csv",
          row.names = FALSE)

# Compile results from ABRicate
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
# Batch 3 (xx May 2026)
# Open new terminal and change owner to let the files be accessed with R

current_batch <- 2
patterns <- c(
  "ncbi" = ".*\\_ncbi.tab$",
  "resfinder" = ".*\\_resfinder.tab$",
  "plasmidfinder" = ".*\\_plasmidfinder.tab$",
  "vfdb" = ".*\\_vfdb.tab$"
)

for (db_name in names(patterns)) {
  compiled <- compile_tabfiles(
    path = paste0("inputs/output_abricate_batch", current_batch),
    name_file = patterns[[db_name]]
  )
  
  write.table(
    compiled,
    file = paste0("inputs/output_abricate_batch", current_batch,
                       "/abricate_", db_name,"_batch", current_batch, 
                       "_long.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


################################################################################
# correct ID following workLab workSeq; combine batches and import to inputs/
workLab_workSeq_all <- read.csv("inputs/workLab_workSeq_compiled_all.csv") %>% 
  dplyr::select(id,
                workDmx_filename,
                workSeq_Speciator.genusName_post,
                workSeq_Speciator.speciesName_post
                ) %>% 
  glimpse()

abricate_fol <- list.dirs("inputs",
                          recursive = FALSE,
                          full.names = TRUE)
abricate_fol <- abricate_fol[grepl("output_abricate_batch",
                                   abricate_fol)]

for (ft in c("long", "wide")){
  patterns <- c(
    "ncbi" = paste0(".*_ncbi_batch[0-9]+_", ft, "\\.t.*$"),
    "resfinder" = paste0(".*_resfinder_batch[0-9]+_", ft, "\\.t.*$"),
    "plasmidfinder" = paste0(".*_plasmidfinder_batch[0-9]+_", ft, "\\.t.*$"),
    "vfdb" = paste0(".*_vfdb_batch[0-9]+_", ft, "\\.t.*$")
  )
  
  for (db_name in names(patterns)){
    compiled <- compile_tabfiles(
      path = abricate_fol,
      name_file = patterns[[db_name]]
    ) %>% 
      # delete info in the last "/" & weird extention
      dplyr::mutate(
        filename = gsub(".*/|\\.fastq.gz.long.fasta$|*_abricate_.*$||\\.tab$", "",
                        X.FILE)
      ) %>% 
      dplyr::left_join(
        workLab_workSeq_all
        ,
        by = c("filename" = "workDmx_filename")
      ) %>% 
      dplyr::rename_with(~ paste0("workAbr_", db_name, "_", .))
    
    write.table(
      compiled,
      file = paste0("inputs/abricate_", db_name,"_compiled_", ft, "_all.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}

# I don't think full data acquired for "_wide" data
# (bind_rows with possibly different columns)
# How did Seeman compile the info to wide format again?
# nah it is perfectly combined ;)

test_wide <- read.table(
  "inputs/abricate_plasmidfinder_compiled_wide_all.tsv",
  header = TRUE,
  sep = "\t",
) %>% 
  glimpse()


################################################################################
# Analyse both acquired & point mutation-related AMR using abriTAMR
# This samplesheet can also be used for species-specific MLST scheme
# Prepare samplesheet specified for species from all pass QC samples
# Available MLST sheme option: mlst --info

workLab_workSeq_all <- read.csv("inputs/workLab_workSeq_compiled_all.csv") %>% 
  dplyr::transmute(#id = id,
                   #workDmx_filename = workDmx_filename,
                   workAMR_fasta = gsub("output_demux", "output_ghru-mod_lambda_filter_on/assemblies",
                                        workDmx_address),
                   workAMR_fasta = gsub(".fastq", ".fastq.gz.long.fasta", workAMR_fasta),
                   mlst_scheme = tolower(stringr::str_replace(workSeq_Speciator.speciesName_post,
                                                      "(\\w)\\w+ (\\w+)", "\\1\\2")),
                   abritamr_scheme = gsub(" ", "_", workSeq_Speciator.speciesName_post)
  ) %>%
  # specify both MLST & abriTAMR scheme
  dplyr::mutate(
    mlst_scheme = ifelse(mlst_scheme == "kpneumoniae",
                         "klebsiella",
                         mlst_scheme),
    abritamr_scheme = ifelse(abritamr_scheme == "Escherichia_coli",
                             "Escherichia",
                             abritamr_scheme)
  ) %>% 
  glimpse()

workLab_workSeq_all <- workLab_workSeq_all %>% 
  dplyr::filter(mlst_scheme == "ecoli") %>%
  dplyr::mutate(mlst_scheme = paste0(mlst_scheme, "_achtman_4")) %>%
  dplyr::bind_rows(workLab_workSeq_all, .) %>%
  dplyr::arrange(mlst_scheme) %>% 
  # view() %>% 
  glimpse()


write.csv(workLab_workSeq_all,
          "inputs/species_samplesheet.csv",
          row.names = FALSE,
          quote = FALSE
          )



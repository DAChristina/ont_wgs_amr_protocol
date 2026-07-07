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
# Batch 3 (20 May 2026)
# Batch 4 FAILED due to insufficient active pore
# Batch 5 (29 June 2026)
# Batch 6 (06 July 2026)

################################################################################
# Analyse both acquired & point mutation-related AMR using abriTAMR
# ALL data simultaneously (not per-batch)
# This samplesheet can also be used for species-specific MLST scheme
# Prepare samplesheet specified for species from all pass QC samples
# Available MLST scheme option: mlst --info

# COPY all fasta files to "/srv/nfs_share/2026_ACORNHAI/raw_data/0_compiled_fasta/" first!
workLab_workSeq_all <- read.csv("inputs/workLab_workSeq_compiled_all.csv") %>% 
  
  # filter just only for PASS QC
  # need to be adjusted following QC report
  # dplyr::filter() %>% 
  dplyr::transmute(#id = id,
    #workDmx_filename = workDmx_filename,
    workAMR_fasta = paste0("/srv/nfs_share/2026_ACORNHAI/raw_data/0_compiled_fasta/",
                           workDmx_filename, ".fastq.gz.long.fasta"),
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

# additional row for MLST scheme if required
workLab_workSeq_all <- workLab_workSeq_all %>% 
  dplyr::filter(mlst_scheme == "ecoli") %>%
  dplyr::mutate(mlst_scheme = paste0(mlst_scheme, "_achtman_4")) %>%
  dplyr::bind_rows(workLab_workSeq_all, .) %>%
  dplyr::arrange(mlst_scheme) %>% 
  
  dplyr::filter(mlst_scheme == "abaumannii") %>%
  dplyr::mutate(mlst_scheme = paste0(mlst_scheme, "_2")) %>%
  dplyr::bind_rows(workLab_workSeq_all, .) %>%
  dplyr::arrange(mlst_scheme) %>% 
  # view() %>% 
  glimpse()

write.csv(workLab_workSeq_all, "inputs/species_samplesheet.csv",
          row.names = FALSE, quote = FALSE)


# prepare A. baumannii-specialised for plasmid samplesheet
write.csv(workLab_workSeq_all %>% 
            dplyr::filter(abritamr_scheme == "Acinetobacter_baumannii") %>% 
            dplyr::distinct(workAMR_fasta, .keep_all = TRUE)
          ,
          "inputs/species_samplesheet_abaumannii_plasmid.csv",
          row.names = FALSE, quote = FALSE)


################################################################################
# compile abriTAMR results into one df
workLab_workSeq_all <- read.csv("inputs/workLab_workSeq_compiled_all.csv") %>% 
  dplyr::select(id,
                workDmx_filename,
                workSeq_Speciator.genusName_post,
                workSeq_Speciator.speciesName_post
  ) %>% 
  glimpse()

patterns <- c(
  "abritamr" = ".*\\_abritamr.txt$",
  # "amrfinder" = ".*\\_amrfinder.out$",
  # amrfinder is special case coz' isolates not included in their output table
  "summary_matches" = ".*\\_summary_matches.txt$",
  "summary_partials" = ".*\\_summary_partials.txt$",
  "summary_virulence" = ".*\\_summary_virulence.txt$"
)

for (db_name in names(patterns)) {
  compiled <- compile_tabfiles(
    path = paste0("inputs/output_abritamr/"),
    name_file = patterns[[db_name]]
  )
  
  write.table(
    compiled,
    file = paste0("inputs/abritamr_", db_name, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

# fixing amrfinder results
files <- list.files("inputs/output_abritamr",
                    pattern = ".*\\_amrfinder.out",
                    full.names = TRUE,
                    recursive = TRUE)

compiled <- lapply(files, function(f) {
  file_id <- basename(f) 
  read.delim(f, sep = "\t", header = TRUE) %>% 
    dplyr::select(-any_of(c("text_source", "filename"))
                  ) %>% 
    dplyr::mutate(
      dplyr::across(everything(), as.character),
      Class = stringr::str_to_title(Class),
      Subclass = stringr::str_to_title(Subclass),
      text_source = gsub("_amrfinder.out", "", file_id)
    )
}) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(
    filename = text_source,
    
  ) %>%
  dplyr::left_join(
    workLab_workSeq_all
    ,
    by = c("filename" = "workDmx_filename")
  ) %>% 
  dplyr::filter(!is.na(id)) %>% 
  dplyr::rename_with(~ paste0("workAbr_", .))

write.table(compiled, "inputs/abritamr_amrfinder.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# test AMR result
amr <- read.table("inputs/abritamr_amrfinder.tsv",
                  sep = "\t",
                  header = TRUE,
                  quote = "") %>% # use quote = "" to avoid tab misalignment
  glimpse()
length(unique(amr$workAbr_id))

################################################################################
# compile ABRicate results into one df (long format)
# Just use _long data for viz and further analyses; standardise the IDs
# correct ID following workLab workSeq; combine batches and import to inputs/
workLab_workSeq_all <- read.csv("inputs/workLab_workSeq_compiled_all.csv") %>% 
  dplyr::select(id,
                workDmx_filename,
                workSeq_Speciator.genusName_post,
                workSeq_Speciator.speciesName_post
  ) %>% 
  glimpse()

patterns <- c(
  "plasmidfinder" = ".*\\_plasmidfinder.tab$",
  "vfdb" = ".*\\_vfdb.tab$"
)

for (db_name in names(patterns)) {
  compiled <- compile_tabfiles(
    path = paste0("inputs/output_abricate"),
    name_file = patterns[[db_name]]
  ) %>% 
    dplyr::mutate(
      filename = gsub(".*/|\\.fastq.gz.long.fasta$|*_abricate_.*$||\\.tab$", "",
                      X.FILE),
      gene_class = str_extract(PRODUCT, "(?<= - )[^\\(]+(?= \\()"),  # " - " & " ("
      species_reference = str_extract(PRODUCT, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" & "]"
    ) %>%
    dplyr::left_join(
      workLab_workSeq_all
      ,
      by = c("filename" = "workDmx_filename")
    ) %>% 
    dplyr::rename_with(~ paste0("workAbr_", .))
  
  write.table(
    compiled,
    file = paste0("inputs/abricate_", db_name, "_long.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

# test plasmid result; check wide file to see detected plasmids
plasmid_long <- read.table("inputs/abricate_plasmidfinder_long.tsv",
                          sep = "\t",
                          header = TRUE,
                          quote = "") %>% # use quote = "" to avoid tab misalignment
  glimpse()
length(unique(plasmid_long$workAbr_id))

plasmid_wide <- read.delim("inputs/abricate_plasmidfinder_wide.tab",
                           sep = "\t",
                           header = TRUE,
                           # quote = ""
                           ) %>% # use quote = "" to avoid tab misalignment
  glimpse()
length(unique(plasmid_wide$X.FILE))


# test VFDB 
vfdb_long <- read.table("inputs/abricate_vfdb_long.tsv",
                           sep = "\t",
                           header = TRUE,
                           quote = "") %>% # use quote = "" to avoid tab misalignment
  glimpse()
length(unique(vfdb_long$workAbr_id))

vfdb_wide <- read.delim("inputs/abricate_vfdb_wide.tab",
                           sep = "\t",
                           header = TRUE,
                           # quote = ""
) %>% # use quote = "" to avoid tab misalignment
  glimpse()
length(unique(plasmid_wide$X.FILE))

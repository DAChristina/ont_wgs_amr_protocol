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


current_batch <- 2
patterns <- c(
  "card"         = ".*\\_card.tab$",
  "ncbi"         = ".*\\_ncbi.tab$",
  "resfinder"    = ".*\\_resfinder.tab$",
  "plasmidfinder"= ".*\\_plasmidfinder.tab$",
  "vfdb"         = ".*\\_vfdb.tab$"
)

for (db_name in names(patterns)) {
  compiled <- compile_tabfiles(
    path      = paste0("inputs/output_abricate_batch", current_batch),
    name_file = patterns[[db_name]]
  )
  
  write.table(
    compiled,
    file      = paste0("inputs/output_abricate_batch", current_batch,
                       "/abricate_", db_name,"_batch", current_batch, 
                       "_long.tsv"),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )
}


test_card <- read.delim("inputs/output_abricate_batch2/abricate_card_batch2_long.tsv",
                        sep = "\t",
                        header = T
                        ) %>% 
  glimpse()

test_ncbi <- read.delim("inputs/output_abricate_batch2/abricate_ncbi_batch2_long.tsv",
                        sep = "\t",
                        header = T
) %>% 
  glimpse()


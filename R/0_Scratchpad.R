# test AMR-Plasmid data compilation
library(tidyverse)

current_batch <- 1

# compiled abricate data
test1 <- read.delim(
  "inputs/abricate_plasmidfinder_batch1.tab",
  sep = "\t",
  header = T
) %>% 
  glimpse()



# real abricate data per-ID
test2 <- read.delim(
  "raw_data/ACORNHAI_3_22APRIL2026/ACORNHAI_3_22APRIL2026/20260422_1507_X4_FBD04135_d247130b/output_abricate/d247130b-acd0-4eb3-9f8e-f90fd6634eeb_SQK-NBD114-24_barcode09_abricate_plasmidfinder.tab",
  sep = "\t",
  header = T
) %>% 
  glimpse()


# I'm not sure we really want to use the compiled abricate data for every single batch (plasmids from PlasmidFinder might be different depending on the samples).
tests <- readxl::read_excel("raw_data/update acorn hai_isolasi DNA.xlsx",
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
    fromLast = TRUE)) %>% 
  glimpse()

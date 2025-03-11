library(tidyverse)

HRISPLEX_CSV <- "/home/lab/Downloads/hirisplexs.csv"
QUERY_PATH <- "/media/lab/Data1/longReads/R10_1KG/skin/hrisplex/vcfs/hrisplex_query.txt"
OUTPUT_PATH <- "/media/lab/Data1/longReads/R10_1KG/skin/hrisplex/vcfs/hrisplex.csv"


hrisplex_format <- read_csv(HRISPLEX_CSV)

hrisplex_tbl <-
  hrisplex_format %>% 
  add_row() %>% 
  pivot_longer(-sampleid,
               names_to = "markers",
               values_to = "id") %>% 
  select(markers)

hrisplex_tbl_tidy <-
  hrisplex_tbl %>% 
  separate(col = markers,into = c("ID","target"))


data <- read_delim(QUERY_PATH, delim = "\t", col_names = FALSE)

# Ajustar os nomes das colunas
colnames(data) <- c("SNP", paste0("Sample", 1:(ncol(data) - 1)))

data_long <- pivot_longer(data, cols = -SNP, names_to = "Sample", values_to = "Genotype")

data_long <-
  data_long %>% 
  separate(col = Genotype,sep = "\\=",into = c("sampleid","Genotype")) %>% 
  select(-Sample) %>% 
  rename(ID = SNP)

data_long %>% 
  right_join(hrisplex_tbl_tidy) %>% 
  mutate(cname = paste0(ID,"_",target)) %>% 
  mutate(Counts = str_count(Genotype,target)) %>% 
  select(sampleid, markers = cname,Counts) %>% 
  right_join(hrisplex_tbl,na_matches = "never") %>% 
  pivot_wider(names_from = markers, values_from = Counts) %>% 
  filter(!is.na(sampleid)) %>% 
  select(all_of(names(hrisplex_format))) %>% 
  write_csv(file = OUTPUT_PATH)

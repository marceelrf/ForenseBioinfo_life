library(tidyverse)
library(glue)
library(fs)

samples_fastq_path="/media/lab/Data1/longReads/R10_1KG/skin/PECAT"
fasta_path="/media/lab/Data1/longReads/R10_1KG/skin/assembly"
output_path="/media/lab/Data1/longReads/R10_1KG/skin/assembly/retro_BAMs"
threads = 6

# list.files(path = samples_fastq_path,recursive = T,pattern = ".fastq")

samples <- 
  fs::dir_ls(samples_fastq_path,type = "directory") %>% 
  str_remove_all(glue("{samples_fastq_path}/"))

fastq <- fs::dir_ls(path = glue("{samples_fastq_path}"),type = "file",regexp = ".fastq$",recurse = 1)

fasta <- fs::dir_ls(path = glue("{fasta_path}"),type = "file",regexp = ".fasta$")

fastq_tbl <- tibble(file_fastq = fastq) %>% 
  mutate(bas = str_remove_all(file_fastq,glue("{samples_fastq_path}/|_reads.fastq"))) %>% 
  separate(col = bas,into = c("Sample","Gene"))

fasta_tbl <- tibble(file_fasta = fasta) %>% 
  mutate(bas = str_remove_all(file_fasta,glue("{fasta_path}/|.fasta"))) %>% 
  mutate(Gene = str_split_i(bas,"_",2),
         Sample = str_split_i(bas,"_",1),
         Hap = str_split_i(bas,"_",3),
         ctg = str_split_i(bas,"_",4)) %>% 
  select(-bas)

all_tbl <- 
  full_join(fastq_tbl,fasta_tbl) %>% 
  drop_na() %>% 
  mutate(cmd = glue("minimap2 -ax map-ont {file_fasta} {file_fastq} -t {threads} > {output_path}/{Sample}_{Gene}_{Hap}_{ctg}.sam"))

tictoc::tic()
for (i in seq_along(all_tbl$cmd)) {
  print(paste0(i," of ",length(all_tbl$cmd), " done!"))
  system(all_tbl$cmd[i])
}
tictoc::toc()
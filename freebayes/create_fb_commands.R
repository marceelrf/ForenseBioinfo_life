library(tidyverse)
library(glue)


BAM_LIST_PATH <- "/media/lab/Data1/longReads/R10_1KG/skin/hrisplex/bam_list.txt"
BED_PATH <- "/media/lab/Data1/longReads/R10_1KG/skin/BEDs/bedsrel2/hirisplex_base1.bed"
OUTPUT_PATH <- "/media/lab/Data1/longReads/R10_1KG/skin/hrisplex/vcfs"
GENOME_PATH <- "/media/lab/Data1/genomes/1KG_ONT_VIENNA_hg38.fa"

bed <- read_tsv(file = BED_PATH,col_names = c("chr","start","end","id"))

if(!dir.exists(OUTPUT_PATH)){
  dir.create(OUTPUT_PATH)
}


cmd <- NULL
for (i in 1:nrow(bed)) {
  cmd <- append(cmd, glue("freebayes -f {GENOME_PATH} -L {BAM_LIST_PATH} -r {bed$chr[i]}:{bed$start[i]}-{bed$end[i]} > {OUTPUT_PATH}/{bed$id[i]}.vcf"))
}

write_lines(cmd,file = glue("{OUTPUT_PATH}/run_fb_commands.sh"),sep = "\n")
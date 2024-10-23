library(tidyverse)
library(seqinr)
library(glue)

"\\wsl.localhost/Ubuntu/home/marceelrf/Projs/Proteomes"

hsa <- read.fasta(file = "UP000005640_9606.fasta.gz",
                  seqtype = "AA",as.string = T)


nm <- hsa %>% 
  seqinr::getName.SeqFastaAA()

seq <- hsa %>% 
  seqinr::getSequence(as.string = T)

unlist(seq)


hsa_tib <- tibble(name = nm, seq = unlist(seq))

hsa_tib_BMP <- 
  hsa_tib %>% 
  filter(str_detect(name,"BMP|GDF|ACV|TGF|INHB|INHA")) %>% 
  mutate(fasta = glue(">{name}\n{seq}")) 


clipr::write_clip(hsa_tib_BMP$name)

hsa_tib_BMP %>% 
  pull(fasta) %>% 
  write_lines(file = "BMP_human.fasta",sep = "\n")

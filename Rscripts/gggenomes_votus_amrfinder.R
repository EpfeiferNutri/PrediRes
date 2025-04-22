# libraries
library(tidyverse)
library(gggenomes)
library(readxl)

### path to vOTUs table
    vOTU_tbl = read_excel("Table_S1_vOTUs_final.xlsx") 

    # AMR finderout vOTUs
amrfinder_out = read_tsv("non_public/amrfinder_out.tsv") %>%
  mutate(vOTU_ID = str_remove(pattern="_[0-9]{1,5}$", `Protein identifier`)) %>% 
  left_join(vOTU_tbl , by = "vOTU_ID") 
  

#### make a plot of the genomes using gggenomes

  # subset the genomad protein table (see DOI: 10.5281/zenodo.14758483) for the plot 
  genomad_annotation = read_tsv("vOTUs/vOTU_genomad_protein_tbl.tsv") %>%
     inner_join(amrfinder_out, by = c("vOTU_ID")) %>%
     mutate( strand = ifelse( strand=="1","+","-"), seq_id = vOTU_ID
            ,gene_col = case_when(vOTU_protein_id %in% amrfinder_out$`Protein identifier`~"tomato", !(is.na(marker)) ~ "skyblue", T ~ "grey" )) %>%
    select(seq_id,vOTU_protein_id, start, end, strand, width=length, size, gene_col)
  
  # vOTUs to plot
  votu_plot = genomad_annotation %>% select(seq_id, length=size) %>% distinct()

  ### plot blast comparison to make a simple plot
  gggenomes(genes = genomad_annotation, seqs = votu_plot) +  
    geom_seq() + geom_bin_label() + geom_gene(aes(fill = gene_col))+
    scale_fill_manual(values = c("grey","skyblue","red"))+
    scale_x_continuous(breaks= seq(0,100e4,10e3))+
    theme(legend.position = "none", text = element_text(size = 20))
  
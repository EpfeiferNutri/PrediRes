# libraries
library(tidyverse)
library(gggenomes)
library(readxl)

### path to vOTUs table
    vOTU_tbl = read_tsv("../Table_S1_vOTUs_final_v2.tsv") 

    # AMR finderout vOTUs
amrfinder_out = read_tsv("amrfinder_out.tsv") %>%
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

  ### To plot highly abundant phages:
  
  ## subset from abundance tbl  
  
  # read signal overview
  Read_signal_grouped = read_csv("../Table_S2_overview_readsignal.csv")  

  # take only the viral
  Viral_signal = Read_signal_grouped %>% filter(Sequence_set == "viral") %>% filter(mapped_read_count > 1e6, `fraction in %` > 10) %>% 
    select(donor, day, mapped_read_count, fraction=`fraction in %`) %>%
    mutate(Sampling_day = ifelse( grepl(day,pattern="^-"), sprintf("%+04d", day), sprintf("%03d", day))
      ,Donor = ifelse(nchar(donor)==1, sprintf("%02d", donor), donor)) %>% select(-donor, -day) %>%
    mutate( Donor = ifelse(Donor == "11", "No_Ref_11", Donor))

  ## load abundance table and select highly abundant phages
  abundance_vOTU_tbl = read_tsv("../Table_S4_abundance_tbl_final_pub_v2.tsv") %>%
      inner_join(Viral_signal, by = c("Donor","Sampling_day")) %>%
      mutate(vOTU_Donor = str_c(vOTU_ID, Donor, sep="_")) %>%
      group_by(Donor)  %>%
      mutate(Sample_ID = str_c(Donor, Sampling_day, sep = "_")
           # down size to the smallest count per donor
          ,min_read_set = min(mapped_read_count)
          ,rary_abundance = ab_abundance*(min_read_set/mapped_read_count)
           # normalise to the size
           ,norm_rare_abundance = rary_abundance/size) %>% 
            group_by(Sample_ID) %>% mutate(
           # compute relative abundance (sum = 100) per sample 
           rel_abund  = round(100*norm_rare_abundance/sum(norm_rare_abundance),9)) %>%
      select(vOTU_ID, Sample_ID,Sampling_day,Donor, rel_abund) %>%
      group_by(vOTU_ID,Donor) %>% 
      mutate(max_abdundance = max(rel_abund) ,abdundance_metric= ifelse(max_abdundance> 25, vOTU_ID, "Low_abundance"))
  
  ## make tbl with dominant vOTUs
  
  # vOTU_0984 is the only found in two donors = 11, 20 (at a cutoff of 25%)
  dominant_vOTUs = abundance_vOTU_tbl %>% filter(abdundance_metric != "Low_abundance") %>% 
    left_join(vOTU_tbl) %>% group_by(vOTU_ID,Donor) %>% arrange(desc(rel_abund)) %>% slice(1) %>% ungroup()
  
  
  
  ### plot tbl for gggenomes
  genomad_annotation_dominant = read_tsv("../vOTUs/vOTU_genomad_protein_tbl.tsv") %>%
     inner_join(dominant_vOTUs %>% distinct(vOTU_ID, size), by = c("vOTU_ID")) %>%
     mutate( strand = ifelse( strand=="1","+","-"), seq_id = vOTU_ID
            ,geNomad_marker_type = str_extract(marker, pattern ="[A-Z]{2}$")
             ,gene_col = case_when( !is.na(geNomad_marker_type)~ geNomad_marker_type, T ~ "NA" )
            ) %>%
    select(seq_id,vOTU_protein_id, start, end, strand, width=length, size, gene_col, marker, geNomad_marker_type)
  
  
  # dominant vOTUs to plot
    votu_plot_dominant = genomad_annotation_dominant %>% select(seq_id, length=size) %>% arrange(desc(length)) %>% distinct() %>%
      # choose how many to plot
      slice(1:10)
  
    ### plot blast comparison to make a simple plot
    gggenomes(genes = genomad_annotation_dominant, seqs = votu_plot_dominant) +  
      geom_seq() + geom_bin_label() + geom_gene(aes(fill = gene_col))+
      scale_fill_manual(values = c("grey","darkviolet","limegreen","skyblue"))+
      scale_x_continuous(breaks= seq(0,100e4,10e3))+
      theme(text = element_text(size = 20), legend.title = element_blank())
  
  
    votu_plot_dominant = genomad_annotation_dominant %>% select(seq_id, length=size) %>% arrange(desc(length)) %>% distinct() %>%
      # choose how many to plot
      slice(11:20)
  
    ### plot blast comparison to make a simple plot
    gggenomes(genes = genomad_annotation_dominant, seqs = votu_plot_dominant) +  
      geom_seq() + geom_bin_label() + geom_gene(aes(fill = gene_col))+
      scale_fill_manual(values = c("tomato","grey","darkviolet","limegreen","skyblue"))+
      scale_x_continuous(breaks= seq(0,100e4,10e3))+
      theme(text = element_text(size = 20), legend.title = element_blank())
  

  votu_plot_dominant = genomad_annotation_dominant %>% select(seq_id, length=size) %>% arrange(desc(length)) %>% distinct() %>%
      # choose how many to plot
      slice(21:30)
  
    ### plot blast comparison to make a simple plot
    gggenomes(genes = genomad_annotation_dominant, seqs = votu_plot_dominant) +  
      geom_seq() + geom_bin_label() + geom_gene(aes(fill = gene_col))+
     scale_fill_manual(values = c("grey","darkviolet","limegreen","skyblue"))+
      scale_x_continuous(breaks= seq(0,100e4,10e3))+
      theme(text = element_text(size = 20), legend.title = element_blank())
    
  
  
  
  
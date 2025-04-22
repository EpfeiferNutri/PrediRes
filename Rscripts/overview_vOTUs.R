### libraries

library(tidyverse)
library(readxl)
library(UpSetR)

# path to vOTUs table
vOTU_tbl = read_tsv("../Table_S1_vOTUs_final_v2.tsv") 

### Source plot 
   vOTU_tbl %>% dplyr::count(source, name = "counts") %>% mutate(Name="source") %>%
      ggplot(aes(x=Name, y=counts, fill=source)) + geom_col(col="black")+
      scale_fill_manual(values = alpha(c("grey40","white","royalblue"), alpha = 0.8))+
      theme_bw()+theme(text = element_text(size = 14))+xlab("")+ylab("#counts")+coord_flip()

### checkV_quality plot

vOTU_tbl %>% 
  # make counts and plot order
  dplyr::count(checkv_quality, name = "counts") %>% mutate(Name="checkV") %>%
      mutate(plot_order = case_when(checkv_quality=="Complete"~"01",checkv_quality=="High-quality"~"02"
                                    ,checkv_quality=="Medium-quality"~"03", T ~"04")) %>%
  arrange(plot_order) %>% filter(plot_order %in% c("01","02","03")) %>%
  mutate(checkv_quality = factor(checkv_quality, levels = checkv_quality)) %>%
  # plot    
  ggplot(aes(x=Name, y=counts, fill=checkv_quality)) + geom_col(col="black")+
      scale_fill_manual(values = alpha(c("skyblue","limegreen","orange","grey"), alpha = 0.8))+
      theme_bw()+theme(text = element_text(size = 14))+xlab("")+ylab("#counts")+coord_flip()


### upsetR plot of viral prediction tools and CheckV
df_to_plot = vOTU_tbl %>% select(vOTU_ID, geNomad, CheckV=`checkV (at least medium)`, VIBRANT, Viral_verify) %>% 
  column_to_rownames("vOTU_ID")

# plot of intersection 
upset(data=df_to_plot 
      #,nsets = 4
      ,mainbar.y.label = "Viral contigs"
      #,show.numbers = T
      ,order.by = "freq"
      ,sets.x.label = "
      " 
      ,decreasing = T, 
      main.bar.color = alpha("tomato", alpha = 0.8),
      sets.bar.color = "orange",
      text.scale = 1.8,
      point.size = 6)

# Sizes and qualities 

vOTU_tbl_sizes =  
  # make size categories
  vOTU_tbl %>% mutate( size_bins = cut(size, breaks = c(3e3,10e3,20e3,30e3,50e3,100e3,200e3,275e3))
                       ,size_bins_form = gsub(",", "-", gsub("\\(|\\]", "", as.character(size_bins)))) %>%
    dplyr::count(size_bins_form, checkv_quality) %>% 
  mutate(plot_order = case_when(checkv_quality=="Complete"~"01 Complete",checkv_quality=="High-quality"~"02 High-quality"
                                    ,checkv_quality=="Medium-quality"~"03 Medium-quality", T ~"04 Low-quality+ND")) %>% arrange(plot_order)
#plot 
vOTU_tbl_sizes %>% ggplot(aes(x=size_bins_form,y=n, fill = plot_order)) + geom_col()+
    scale_x_discrete(limits = c("3e+03-1e+04","1e+04-2e+04","2e+04-3e+04","3e+04-5e+04","5e+04-1e+05","1e+05-2e+05","2e+05-2.75e+05")
                     ,labels = c("3-10","10-20","20-30","30-50","50-100","100-200","200-275"))+
      scale_y_continuous(limits = c(0,2500), breaks = seq(0,2500,600))+
      scale_fill_manual(values = c("skyblue","limegreen","orange","grey80")) +
    theme_bw()+theme(text = element_text(size = 14))+xlab("Sizes (kb)")+ylab("#Counts")+coord_flip()
    
# Temperate phages 

vOTU_tbl %>% filter(checkv_quality %in% c("Complete","High-quality")) %>% 
  # count
    dplyr::count(viral_lifecycle) %>% mutate(Fraction = round(100*n/sum(n),1)) %>% 
  # plot  
  ggplot(aes(y=Fraction, x="BACPHLIP", fill = viral_lifecycle)) + 
    geom_bar(position = "stack", stat = "identity", width = 0.8, col = "grey30") +
    scale_fill_manual(values = alpha(c("orange","skyblue"), alpha = 1))+
    #scale_x_discrete(limits = c("Bacphlip\nn=1130","Vibrant\nn=1232","Combined\n(consistent)\nn=1118") )+
    theme_bw()+ylab("Fraction (%)")+
    scale_y_continuous(breaks = c( seq(0,100,20))) + xlab("Prediction Tool (complete+high quality)")+
    theme(text = element_text(size = 16), legend.position = "top",legend.title = element_blank())

### novel vOTUs, fraction quality

vOTU_tbl %>% filter(Novel_species %in% c("Yes")) %>% 
  # count
    dplyr::count(checkv_quality) %>% mutate(Fraction = round(100*n/sum(n),1)) %>% 
    mutate(plot_order = case_when(checkv_quality=="Complete"~"01",checkv_quality=="High-quality"~"02"
                                    ,checkv_quality=="Medium-quality"~"03", T ~"04")) %>%
  arrange(plot_order) %>%
  mutate(checkv_quality = factor(checkv_quality, levels = checkv_quality)) %>%
  # plot  
  ggplot(aes(x=Fraction, y="Quality", fill = checkv_quality)) + 
      geom_col(col="black")+
      scale_fill_manual(values = alpha(c("skyblue","limegreen","orange","grey","grey"), alpha = 0.8))+
      theme_bw()+theme(text = element_text(size = 14))+xlab("")+ylab("#counts")
  


#### Overview of read data set (Figure S3)

# read in read signal 

Read_signal_grouped = read_csv("Table_S2_overview_readsignal.csv")  

# adjust order of the bars, and label of x-axis 
  Read_signal_grouped_plot = Read_signal_grouped %>%
  #make plot order (bottom to up)
  mutate(plot_order = case_when( Sequence_set=="unmapped"~"01"
                                ,Sequence_set=="non_viral"~"02"
                                ,Sequence_set=="proviral"~"03"
                                ,Sequence_set=="UHGG"~"04"
                                ,Sequence_set=="UHGV"~"05"
                                ,Sequence_set=="viral"~"06")) %>%
  arrange(plot_order) %>% 
  mutate(Sequence_set = factor(Sequence_set, levels=unique(Sequence_set))
        ,Day_Reads = str_c(day,"\n",mill_reads_per_sample)) %>% arrange(day) %>%
  mutate(Day_Reads_factored = factor(Day_Reads, levels = unique(Day_Reads)))

# make the plot
  Read_signal_grouped_plot %>% 
    ggplot(aes(x=Day_Reads_factored, y = `fraction in %`, fill = Sequence_set )) +
    geom_bar(position = "stack", stat = "identity", width = 0.8)+
    scale_fill_manual(values = c("grey","black","orange","tomato","darkviolet","skyblue") )+
    theme_bw()+ylab("Mapped Total Reads Per Sample (%)")+
    scale_y_continuous(breaks = c( seq(0,98,33),100)) + xlab("Day // Total Reads (Mio)")+
    theme(text = element_text(size = 14)  ,axis.text.x = element_text(size = 10)
        , legend.position = "top",legend.title = element_blank()  , strip.background = element_blank()) + 
    facet_wrap(~donor, scales = "free_x", nrow = 3) + guides(fill=guide_legend(nrow = 1))

#####
  
  ### host classes plot
  
  # color code 
  host_class_colors_extended = data.frame(
    Class = c("Actinomycetia","Bacilli","Bacteroidia","Clostridia","Coriobacteriia","Gammaproteobacteria","Alphaproteobacteria","Methanobacteria","Negativicutes"
              ,"Peptococcia","Synergistia","Brachyspirae","Desulfovibrionia","Vampirovibrionia"
              ,"Elusimicrobia","Fusobacteriia","Lentisphaeria","Campylobacteria","Paceibacteria"
              ,"Saccharimonadia","Spirochaetia","Verrucomicrobiae")
    
    ,Colors = c("red","lightcoral","royalblue","orange","pink","lightgreen","#228B22","#A020F0"
                ,"#FF1493","#FFD700","#2E8B57","#90EE90","#8B4513","#B22222","#FF69B4","#BEBEBE","#40E0D0","#CDB79E","#B0E0E6","#BC8F8F","#7F7F7F","#7B68EE"))
  
  
  ### vOTUs
  vOTUs_host_classes = vOTU_tbl %>% mutate(Class = str_extract(host_taxonomy_lineage, "c__[^;]+") %>% str_remove("^c__") ) %>% 
    select(ID = vOTU_ID, Class) %>% mutate(Source = "vOTUs", Class = ifelse(Class == "Clostridia_A", "Clostridia", Class))
  
  vOTUs_host_classes_tbl = vOTUs_host_classes %>% dplyr::count(Class)
  
  
  vOTUs_host_species = vOTU_tbl %>% mutate(Species = str_extract(host_taxonomy_lineage, "s__[^;]+") %>% str_remove("^s__") ) %>% 
    select(ID = vOTU_ID, Species) %>% mutate(Source = "vOTUs")
  
  vOTUs_host_species_tbl = vOTUs_host_species %>% dplyr::count(Species)
  
  vOTUs_Caudo = vOTU_tbl %>% mutate(Caudo = str_extract(viral_taxonomy_lineage, "Caudoviricetes$") ) %>% 
    select(ID = vOTU_ID, Caudo,viral_taxonomy_lineage) %>% mutate(Source = "vOTUs")
  
  vOTUs_Caudo_tbl = vOTUs_Caudo %>% dplyr::count(Caudo)
  
  
  ### UHGV 
  UHGV_meta = read_tsv("metadata_UHGV.tsv") %>% mutate(host_class = str_extract(host_lineage, "c__[^;]+") %>% str_remove("^c__") ) %>%
    select(ID =uhgv_genome,Class=host_class) %>% mutate(Source = "UHGV", Class = ifelse(Class == "Clostridia_A", "Clostridia", Class))
  
  ### UHGG
  UHGG_meta = read_tsv("genomes-all_metadata_UHGG.tsv") %>% mutate(Class = str_extract(Lineage, "c__[^;]+") %>% str_remove("^c__") ) %>%
    select(ID = Genome, Class) %>% mutate(Source = "UHGG", Class = ifelse(Class == "Clostridia_A", "Clostridia", Class)) 
  
  ### bacterial species
  bacterial_spec = read_tsv("../Table_S5_microbial_species_overview.tsv") %>% mutate(Class = str_extract(gtdb_classification, "c__[^;]+") %>% str_remove("^c__") ) %>%
    select(ID = id_mgs, Class) %>% mutate(Source = "MSP", Class = ifelse(Class == "Clostridia_A", "Clostridia", Class)) 
  
  ### put together and make plot tbl
  
  class_plot_tbl = bind_rows(vOTUs_host_classes,UHGV_meta,bacterial_spec,UHGG_meta) %>% 
    # calculate the fractions 
    group_by(Source) %>% dplyr::count(Class) %>% mutate(Fraction = n/sum(n)) %>% ungroup() %>% 
    left_join(host_class_colors_extended, by = c("Class")) %>% 
     mutate( Colors = case_when( is.na(Class) ~ "grey40", is.na(Colors) ~ "grey40", T ~ Colors))
  
   # make a color vector
  category_colors <- setNames(class_plot_tbl$Colors, class_plot_tbl$Class)
  
  class_plot_tbl %>% ggplot(aes(x=Fraction, y=Source, fill = Class ))+
  geom_bar(stat = "identity",width = 0.8, col = NA)+xlab("Fraction in %")+ylab("")+
  scale_fill_manual(values = category_colors)+
  theme_bw() + 
    theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") 
  
  # look on separate studies
  UHGG_meta_studies = read_tsv("genomes-all_metadata_UHGG.tsv") %>% dplyr::count(Study_accession) %>% filter(n > 1160)
  
  ### keep only > 1160 species
  UHGG_1160 = read_tsv("genomes-all_metadata_UHGG.tsv") %>% semi_join(UHGG_meta_studies, by = "Study_accession") %>% mutate(Class = str_extract(Lineage, "c__[^;]+") %>% str_remove("^c__") ) %>%
    select(ID = Genome, Class,Study_accession) %>% mutate(Source = "UHGG", Class = ifelse(Class == "Clostridia_A", "Clostridia", Class)) %>%
    # add PrediRes project
    bind_rows(read_tsv("../Table_S5_microbial_species_overview.tsv") %>% mutate(Class = str_extract(gtdb_classification, "c__[^;]+") %>% str_remove("^c__") ) %>%
    select(ID = id_mgs, Class) %>% mutate(Study_accession = "PrediRes", Class = ifelse(Class == "Clostridia_A", "Clostridia", Class)) ) %>%
    # count
    group_by(Study_accession) %>% dplyr::count(Class) %>% mutate(Fraction = n/sum(n), sum_n = sum(n)) %>% ungroup() %>%
    left_join(host_class_colors_extended, by = c("Class")) %>% 
    mutate( Colors = case_when( is.na(Class) ~ "grey40", is.na(Colors) ~ "grey40", T ~ Colors)
            ,Study_accession_n = str_c( Study_accession,", n=",sum_n)) %>%
    drop_na()
  
  ## stats
  UHGG_1160_stats = UHGG_1160 %>% group_by(Class) %>% 
    summarise( mean_fraction = round(mean(Fraction),3), sd_fraction = round(sd(Fraction),3))
  
     # make a color vector
  category_colors <- setNames(UHGG_1160$Colors, UHGG_1160$Class)
  
  UHGG_1160 %>% ggplot(aes(x=Fraction, y=Study_accession_n, fill = Class ))+
  geom_bar(stat = "identity",width = 0.8, col = NA)+
  xlab("Fraction in %")+ylab("")+
  scale_fill_manual(values = category_colors)+
  theme_bw() + 
    theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") 
  
  # vOTU per donor
  
  vOTUs_host_classes_per_donor = vOTU_tbl %>% mutate(Class = str_extract(host_taxonomy_lineage, "c__[^;]+") %>% str_remove("^c__") ) %>% 
    select(ID = vOTU_ID, Class) %>% mutate(Source = "vOTUs", Class = ifelse(Class == "Clostridia_A", "Clostridia", Class))
  
  
  
  
  ### msp per donor
  
  ### bacterial species
  bacterial_spec_per_donor = read_tsv("../Table_S5_microbial_species_overview.tsv") %>% mutate(Class = str_extract(gtdb_classification, "c__[^;]+") %>% str_remove("^c__") ) %>%
    mutate(Class = ifelse(Class == "Clostridia_A", "Clostridia", Class)) %>%
    separate_rows(donorID_nsamples, sep = ";") %>% group_by(donorID_nsamples) %>%
    dplyr::count(Class) %>% mutate(Fraction = n/sum(n), sum_n = sum(n)) %>% ungroup() %>%
    left_join(host_class_colors_extended, by = c("Class")) %>% 
    mutate( Colors = case_when( is.na(Class) ~ "grey40", is.na(Colors) ~ "grey40", T ~ Colors)
            ,donorID_nsamples = str_c( donorID_nsamples,", n=",sum_n)) %>%
    drop_na()
  
  
      # make a color vector
  category_colors <- setNames(bacterial_spec_per_donor$Colors, bacterial_spec_per_donor$Class)
  
  bacterial_spec_per_donor %>% filter(sum_n > 100) %>%
  ggplot(aes(x=Fraction, y=donorID_nsamples, fill = Class ))+
  geom_bar(stat = "identity",width = 0.8, col = NA)+
  xlab("Fraction in %")+ylab("")+
  scale_fill_manual(values = category_colors)+
  theme_bw() + 
    theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") 
  
  
  
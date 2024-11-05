####
library(tidyverse)
library(readxl)
library(vegan)

####

### path to vOTUs table
vOTU_tbl = read_excel("path/to/Table_S1_vOTUs_final.xlsx") 

#### path to abundance table
abundance_tbl = read_tsv("path/to/Table_S4_abundance_tbl_final.tsv") %>%
  left_join(vOTU_tbl)

# Donors according to AB treatment
  donor_treatment = 
    data.frame(Donor_Treatment = c("01_ceftriaxone","02_cefotaxime","03_cefotaxime","06_ceftriaxone","08_ceftriaxone"
                                   ,"09_ceftriaxone","11_cefotaxime","12_cefotaxime","13_cefotaxime","14_cefotaxime",
                                   "16_ceftriaxone","17_cefotaxime","18_ceftriaxone","20_cefotaxime","22_ceftriaxone",
                                   "23_ceftriaxone","24_ceftriaxone","25_cefotaxime","27_cefotaxime","28_ceftriaxone","34_cefotaxime")) %>%
    mutate(Donor = str_remove(Donor_Treatment, pattern = "_.*"),Treatment = str_remove(Donor_Treatment, pattern = ".*_")) %>%
    select(-Donor_Treatment)

# read signal overview
Read_signal_grouped = read_csv("NutriPhage/Projects/predires_phageome_22_volunteers/data/Table_SX_overview_readsignal.csv")  

  # take only the viral
  Viral_signal = Read_signal_grouped %>% filter(sample_type == "viral") %>% filter(read_count > 1e6, fraction > 10) %>% 
    select(donor, day, read_count, fraction) %>%
    mutate(Sampling_day = ifelse( grepl(day,pattern="^-"), sprintf("%+04d", day), sprintf("%03d", day))
      ,Donor = ifelse(nchar(donor)==1, sprintf("%02d", donor), donor)) %>% select(-donor, -day)

  
## richness per sample
richness_per_sample = abundance_tbl %>% 
    select(vOTU_ID,Donor, Sampling_day, norm_rel_abund)  %>% 
    filter(norm_rel_abund > 0) %>%
    group_by(Donor, Sampling_day) %>% summarise(n_vOTUs_sample = n())

  #plot 
  richness_per_sample %>% 
    ggplot(aes(x=n_vOTUs_sample)) + geom_histogram(col="black", fill = "royalblue", bins = 9)+
    theme_bw()+theme(text = element_text(size =14))+
    xlab("vOTU richness")+ylab("Count")

# get the number of samples per donor, remove the one with n > 3
  n_samples_per_donor = richness_per_sample %>% 
    group_by(Sampling_day) %>% summarise(n=n()) %>% filter(n>3) %>% ungroup()
                                                  
# richness per day  
  richness_per_day =  abundance_tbl %>% 
  filter(norm_rel_abund > 0) %>% select(-Donor) %>% group_by(Sampling_day, vOTU_ID) %>% 
  summarise() %>% 
  #remove sample of day -7 (one sample; unclear if -14 or -1)
  filter(Sampling_day %in% n_samples_per_donor$Sampling_day ) %>%
  mutate(Sampling_day = case_when(Sampling_day %in% c("015","030") ~"15+30", 
                                  Sampling_day %in% c("004","007") ~ "4+7", 
         T ~ Sampling_day)) %>%
    group_by(Sampling_day) %>% summarise(vOTUs_per_day=n()) %>% ungroup()

  # plot
  richness_per_day %>% ggplot(aes(x=Sampling_day, y=vOTUs_per_day ))+
  geom_col(width = 0.8)+xlab("Day")+ylab("Richness vOTUs")+
  scale_x_discrete(limits = c("-014","-001","4+7","010","15+30","180")
                  ,labels = c("-14","-1","4+7","10","15+30","180"))+
  theme_bw() + 
    theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") +
  guides(color = guide_legend(ncol = 2)) 

###### compute shannon diversity per sample
shannon_per_sample = abundance_tbl %>% 
    semi_join(Viral_signal, by = c("Donor","Sampling_day")) %>%
    mutate(Sample_ID = str_c(Donor, Sampling_day, sep = "_")) %>% 
    select(vOTU_ID, Sample_ID, norm_rel_abund) %>%
    pivot_wider(names_from = Sample_ID, values_from = norm_rel_abund) %>% 
   column_to_rownames("vOTU_ID") %>% summarise_all(.funs = "diversity") %>% 
   t() %>% as.data.frame() %>% 
   rownames_to_column("Sample_ID") %>% rename(Shannon=V1) %>%
   mutate(Donor = str_remove(Sample_ID, pattern="_.*"), Sampling_day = str_remove(Sample_ID, pattern = ".*_")) %>%
   left_join(donor_treatment)
   
# normalize shannon to one day
shannon_of_min_day = shannon_per_sample %>% group_by(Donor) %>% filter(Sampling_day %in% c("-001")) %>% 
   select(Donor,min_shannon=Shannon) %>% ungroup()  
   
    #norm all shannon values
    Norm_shannon = shannon_per_sample %>% 
      inner_join(shannon_of_min_day) %>% 
      mutate(Norm_shannon = Shannon/min_shannon)
    
## plot norm shannon, Figure 2B
  Norm_shannon %>% ggplot(aes(x=Sampling_day, y=Norm_shannon))+
  geom_line(linewidth = 1, aes(group = Donor, color = Donor))  + 
  scale_color_manual(values=  colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set3"))(20))+
  geom_point(size = 2, color="grey70") +
  xlab("Day")+ylab("Shannon (norm.)\ndiversity")+
  scale_x_discrete(limits = c("-014","-007","-001","004","007","010","015","030","090","180")
                   ,labels = c("-14","-7","-1","4","7","10","15","30","90","180"))+
  theme_bw() + theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") +
  guides(color = guide_legend(ncol = 2)) +
  facet_wrap(~Treatment, ncol = 2)  

    
### compute the beta diversity between samples using bray-curtis distance, normalize all samples to one sample (and not per donor as for Shannon)
  
Abundance_tbl_bc = abundance_tbl %>%
  inner_join(Viral_signal, by = c("Donor","Sampling_day")) %>%  
  mutate(min_read_set = min(read_count)
         ,Donor = ifelse(Donor == "11", "No_Ref_11", Donor)
         ) %>% 
  # rarefy to all samples
  mutate(rary_abundance = ab_abundunance*(min_read_set/read_count),rel_rare_abundance = rary_abundance/size
         ,norm_rel_abund  = round(100*rel_rare_abundance/sum(rel_rare_abundance),9) ) %>% 
  ungroup() %>%
  select(vOTU_ID, Donor, Sampling_day, rel_rare_abundance) %>%
  mutate(Sample_ID = str_c(Donor, Sampling_day, sep = "_") ) %>% select(-Donor, -Sampling_day) %>%
  pivot_wider(names_from = Sample_ID, values_from = rel_rare_abundance)

  ## compute BC distances  
  bray_curt_distances = Abundance_tbl_bc %>% column_to_rownames("vOTU_ID") %>% 
  t(.) %>% vegdist(method = "bray", diag = T)
  
  # pcoa
  pcoa = cmdscale(bray_curt_distances, eig = T, k = 10)

  # # best proportion (just tot check)
  # prop_var_df = data.frame(var = round(100*pcoa$eig / sum(pcoa$eig),1)) %>% rownames_to_column("row_IDs")

  ## make a heatmap 
  
  ### heat map
  mat_bray_curt = as.data.frame(as.matrix(bray_curt_distances)) %>% 
    rownames_to_column("Sample_ID") %>%
    mutate(Donor = str_remove(Sample_ID, pattern="_.*"), Donor= ifelse(Donor == "11", "No_Ref_11", Donor)) %>%
    select(-Donor) %>% column_to_rownames("Sample_ID") %>% 
    t(.) %>% as.data.frame() %>%
    rownames_to_column("Sample_ID") %>%
    mutate(Donor = str_remove(Sample_ID, pattern="_.*"), Donor= ifelse(Donor == "11", "No_Ref_11", Donor)) %>%
    select(-Donor) %>% column_to_rownames("Sample_ID") %>%
    as.matrix(.)
  
  # anno col according to donor, day, AB  
  # donor
  anno_col = Viral_signal %>% select(Donor, Sampling_day) %>% 
     mutate(Sample_ID = str_c(Donor, Sampling_day, sep = "_") 
    ,Donor = ifelse(Donor == "11", "No_Ref_11", Donor), Donor = as.factor(Donor)) %>% 
    # remove Donor_29 (no sample before AB treatment)
    filter(Donor != "29") %>%
    column_to_rownames("Sample_ID")
  
  anno_colors = list(Donor= setNames( c(colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set3"))(21),"grey50"),levels(anno_col$Donor)))
  
  library(pheatmap)
  pheatmap( mat = mat_bray_curt
                     ,cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F
                     ,color = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 5,name = "Greys"))(100))   
                     , border_color = NA 
                     ,annotation_colors = anno_colors
                    ,annotation_col = anno_col, annotation_legend = T )  

  

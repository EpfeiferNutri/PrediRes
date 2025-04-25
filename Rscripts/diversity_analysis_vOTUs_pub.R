####
library(tidyverse)
library(readxl)
library(vegan)
library(pheatmap)

####

### path to vOTUs table
vOTU_tbl = read_tsv("Table_S1_vOTUs_final_v2.tsv") 

#### path to abundance table
abundance_tbl = read_tsv("Table_S5_abundance_tbl_final_pub_v2.tsv") %>%
  left_join(vOTU_tbl) 

# Donors according to AB treatment
  donor_treatment = 
    data.frame(Donor_Treatment = c("01_ceftriaxone","02_cefotaxime","03_cefotaxime","06_ceftriaxone","08_ceftriaxone"
                                   ,"09_ceftriaxone","11_cefotaxime","12_cefotaxime","13_cefotaxime","14_cefotaxime",
                                   "16_ceftriaxone","17_cefotaxime","18_ceftriaxone","20_cefotaxime","22_ceftriaxone",
                                   "23_ceftriaxone","24_ceftriaxone","25_cefotaxime","27_cefotaxime","28_ceftriaxone","34_cefotaxime")) %>%
    mutate(Donor = str_remove(Donor_Treatment, pattern = "_.*"),Treatment = str_remove(Donor_Treatment, pattern = ".*_")
           ,Donor = ifelse(Donor == "11", "No_Ref_11", Donor)) %>%
    select(-Donor_Treatment)

# read signal overview
Read_signal_grouped = read_csv("Table_S2_overview_readsignal.csv")  

  # take only the viral
  Viral_signal = Read_signal_grouped %>% filter(Sequence_set == "viral") %>% 
    filter(mapped_read_count > 1e6, `fraction in %` > 10) %>% 
    select(donor, day, mapped_read_count, fraction=`fraction in %`) %>%
    mutate(Sampling_day = ifelse( grepl(day,pattern="^-"), sprintf("%+04d", day), sprintf("%03d", day))
      ,Donor = ifelse(nchar(donor)==1, sprintf("%02d", donor), donor)
      ,Donor = ifelse(Donor == "11", "No_Ref_11", Donor)) %>% select(-donor, -day)

  
## richness per sample
richness_per_sample = abundance_tbl %>% 
    
    inner_join(Viral_signal, by = c("Donor","Sampling_day")) %>%
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
    select(vOTU_ID,Donor, Sampling_day, rel_abund)  %>% 
    filter(rel_abund > 0) %>%
    group_by(Donor, Sampling_day) %>% summarise(n_vOTUs_sample = n())

  #plot (Figure S5A)
  richness_per_sample %>% 
    ggplot(aes(x=n_vOTUs_sample)) + geom_histogram(col="black", fill = "royalblue", bins = 9)+
    theme_bw()+theme(text = element_text(size =14))+
    xlab("vOTU richness")+ylab("Count")

# get the number of samples per donor, remove the one with n > 3
  n_samples_per_donor = richness_per_sample %>% 
    group_by(Sampling_day) %>% summarise(n=n()) #%>% filter(n>3) %>% ungroup()
                                                  
# richness per day  
  richness_vOTU_per_day =  abundance_tbl %>% 
    inner_join(Viral_signal, by = c("Donor","Sampling_day")) %>%
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
      #keep only vOTUs per samples, that are detected, remove sample of day -7 and day 90
    filter(rel_abund > 0, Sampling_day %in% n_samples_per_donor$Sampling_day) %>% select(-Donor) %>%  
  mutate( Host_class = str_extract(host_taxonomy_lineage, "(?<=c__)[^;]+")
         ,Host_class=ifelse(grepl(Host_class, pattern = "Clostridia"),"Clostridia",Host_class)
         ,Sampling_day = case_when(Sampling_day %in% c("015","030") ~"15+30",Sampling_day %in% c("004","007") ~ "4+7",T ~ Sampling_day)) %>%  
  group_by(Sampling_day, vOTU_ID,Host_class) %>%
  # count only one species per time point (independent in which donor)
  summarise() %>% 
    group_by(Sampling_day,Host_class) %>% summarise(vOTUs_per_day=n()) %>% ungroup()
  
  # take a unified color schmeme
  host_class_colors = data.frame(Class = c("Actinomycetia","Bacilli","Bacteroidia","Clostridia","Coriobacteriia","Gammaproteobacteria"),
                                 Colors = c("red","lightcoral","royalblue","orange","pink","lightgreen"))
  
  # make a plot df
  richness_vOTU_per_day_plot = richness_vOTU_per_day %>% left_join(host_class_colors, by = c("Host_class"="Class")) %>%
    mutate(Colors = case_when( is.na(Host_class) ~ "grey40", is.na(Colors) ~ "grey90", T ~ Colors)
           )
  
  # make a color vector
  category_colors <- setNames(richness_vOTU_per_day_plot$Colors, richness_vOTU_per_day_plot$Host_class)
  
  # plot (Figure 2A)
  richness_vOTU_per_day_plot %>% ggplot(aes(x=Sampling_day, y=vOTUs_per_day, fill = Host_class ))+
  geom_bar(stat = "identity",width = 0.8, col = NA)+xlab("Day")+ylab("Richness vOTUs")+
  scale_fill_manual(values = category_colors)+
  scale_x_discrete(limits = c("-014","-001","4+7","010","15+30","180")
                  ,labels = c("-14","-1","4+7","10","15+30","180"))+
  theme_bw() + 
    theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") +
  guides(color = guide_legend(ncol = 2)) 

  
###### vOTUs shared across donors
  
  ### long table with means and sds
     vOTUs_across_donors = abundance_tbl %>% 
        inner_join(Viral_signal, by = c("Donor","Sampling_day")) %>%
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
       filter(rel_abund > 0) %>% group_by(vOTU_ID, Donor) %>%
       summarise(mean_ab = round(mean(rel_abund),4), sd_ab = round(sd(rel_abund),4)) %>%
                   replace_na( replace = list("sd_ab"=0)) %>% 
       mutate(mean_sd_ab = str_c(mean_ab," +/-",sd_ab)) %>% select(-mean_ab, -sd_ab) %>% 
       group_by(vOTU_ID) %>% mutate(n_shared_donor = n())
     
     # make table: rows = vOTUs, columns = donors
     vOTUs_across_donors_tbl = vOTUs_across_donors %>% 
       mutate(Donor = str_c("Donor", Donor)) %>%
       select(-n_shared_donor) %>% 
       pivot_wider(names_from = Donor, values_from = mean_sd_ab, values_fill = "nD") %>%
       left_join(vOTUs_across_donors %>% distinct(vOTU_ID,n_shared_donor), by = "vOTU_ID") %>% 
       ungroup()
     
     # plot (Figure SB)
     vOTUs_across_donors_tbl %>% 
       mutate(n_shared_donor = sprintf("%02d",n_shared_donor)) %>%
       dplyr::count(n_shared_donor) %>% 
       #mutate(n_counts = sprintf("%04d",n)) %>%
       ggplot(aes(x=n_shared_donor, y = n)) + 
       geom_col(fill = alpha("limegreen", alpha = 0.6))+
       geom_text(aes(label = n), nudge_y = 120)+
       xlab("Found in 'X' donors")+ylab("counts (vOTUs)")+
       theme_bw()+theme(text = element_text(size = 14))
      
     ### on the genus/VC level
     
     # attach to abundance tbl
     VCs_accross_donors = vOTUs_across_donors %>% 
       left_join(vOTU_tbl, by = "vOTU_ID") %>% ungroup() %>%
       select(VClustering, Donor) %>% group_by(VClustering, Donor) %>% 
       mutate(n_vOTUs = n()) %>% distinct() %>%
       mutate(Donor = str_c("Donor_",Donor)) %>%
       group_by(VClustering) %>% mutate(n_shared_donors = n()) 
     
     # plot (Figure S5B)
     VCs_accross_donors %>% group_by(n_shared_donors) %>% filter(VClustering != "OUT/SINGLE/OVER") %>%
       mutate(n_shared_donors = sprintf("%02d",n_shared_donors)) %>% group_by(VClustering) %>%
       slice(1) %>% ungroup() %>% dplyr::count(n_shared_donors) %>%
       ggplot(aes(x=n_shared_donors, y = n)) + 
       geom_col(fill = alpha("tomato", alpha = 0.6))+
       geom_text(aes(label = n), nudge_y = 5)+
       xlab("Found in 'X' donors")+ylab("counts (VCs)")+
       theme_bw()+theme(text = element_text(size = 14))
  

###### Shannon diversity per sample

    # for vOTUs
  shannon_per_sample = abundance_tbl %>% 
    inner_join(Viral_signal, by = c("Donor","Sampling_day")) %>%
    group_by(Donor)  %>%
    mutate(Sample_ID = str_c(Donor, Sampling_day, sep = "_")
           ,min_read_set = min(mapped_read_count)
           # down size to the smallest count per donor
           ,rary_abundance = ab_abundance*(min_read_set/mapped_read_count)
           # normalise to the size
           ,norm_rare_abundance = rary_abundance/size) %>% 
  group_by(Sample_ID) %>% mutate(
           # compute relative abundance (sum = 100) per sample 
           rel_abund  = round(100*norm_rare_abundance/sum(norm_rare_abundance),9)) %>% ungroup() %>% 
      select(vOTU_ID, Sample_ID, rel_abund) %>%
   pivot_wider(names_from = Sample_ID, values_from = rel_abund) %>% 
   column_to_rownames("vOTU_ID") %>% summarise_all(.funs = "diversity") %>% 
   t() %>% as.data.frame() %>% 
   rownames_to_column("Sample_ID") %>% rename(Shannon=V1) %>%
   mutate(Donor = str_remove(Sample_ID, pattern="_.*"),Sampling_day = str_remove(Sample_ID, pattern = ".*_")) %>%
   # add information on AB treatment
   left_join(donor_treatment) %>%
   mutate(AB_treatment = case_when(Sampling_day %in% c("-001") ~ "Before"
                                ,Donor == "No_Ref_11" & Sampling_day == "-014" ~ "Before"
                                ,Sampling_day %in% c("004","007") ~ "After"
                                ,Donor == "13" & Sampling_day == "010" ~ "After", 
                                T ~ "Other"))

    # group donors in increasing and decreasing groups
    shannon_groups = shannon_per_sample %>% filter(AB_treatment != "Other") %>% 
      select(-Sampling_day, -Sample_ID) %>%
      pivot_wider( names_from = AB_treatment, values_from = Shannon) %>%
      mutate(Shannon_change = ifelse(After>Before, "Increase","Decrease"))
    
    table(shannon_groups$Shannon_change)
   
    # normalize shannon to one day
    shannon_of_min_day = shannon_per_sample %>% group_by(Donor) %>% filter(Sampling_day %in% c("-001")) %>% 
       select(Donor,min_shannon=Shannon) %>% ungroup()  
       
        #norm all shannon values
        Norm_shannon = shannon_per_sample %>% 
          inner_join(shannon_of_min_day) %>% 
          mutate(Norm_shannon = Shannon/min_shannon) %>%
          left_join(shannon_groups %>% select(Donor,Shannon_change), by = "Donor")
    
      ## plot norm shannon (Figure 2B, Figure S5C)
        Norm_shannon %>% ggplot(aes(x=Sampling_day, y=Norm_shannon))+
        geom_line(linewidth = 1.5, aes(group = Donor, color = Donor))  + 
        scale_color_manual(values=  colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set3"))(20))+
        geom_point(size = 2, color="grey70") +
        xlab("Day")+ylab("Shannon (norm.)\ndiversity")+
        scale_x_discrete(limits = c("-014","-007","-001","004","007","010","015","030","090","180")
                         ,labels = c("-14","-7","-1","4","7","10","15","30","90","180"))+
        theme_bw() + theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") +
        guides(color = guide_legend(ncol = 2)) +
        facet_wrap(~Shannon_change, ncol = 2) 
          # for figure S5C
          # facet_wrap(~Treatment, ncol = 2) 

    #### test changes of shannon diversity using wilcoxon tests and visualize using boxplots    
        wilc_tbl = Norm_shannon %>% mutate( 
                TimeFrame = case_when(Sampling_day %in% c("-014","-007") ~ "Pre-Before", 
                        #Sampling_day %in% c("030") ~ "Short-After",
                        Sampling_day %in% c("180") ~ "Late-After",
                        T ~ AB_treatment )) %>%
            filter(TimeFrame !="Other") %>%
            select(Donor, TimeFrame, Shannon) %>%
            pivot_wider(names_from = TimeFrame, values_from = Shannon)  %>% drop_na()


# Wilcoxon tests

  # pre-before vs before
  wilcox.test(wilc_tbl$`Pre-Before`, wilc_tbl$Before, paired = TRUE)

  # before vs after
  wilcox.test(wilc_tbl$Before, wilc_tbl$After, paired = TRUE)
  
  # pre-before vs after
  wilcox.test(wilc_tbl$`Pre-Before`, wilc_tbl$After, paired = T)
  
  # pre-before vs late-after
  wilcox.test(wilc_tbl$`Pre-Before`, wilc_tbl$`Late-After`, paired = T)
  
  # Before vs Late-After
  wilcox.test(wilc_tbl$Before, wilc_tbl$`Late-After`, paired = T)
  
  # After vs Late-After
  wilcox.test(wilc_tbl$After, wilc_tbl$`Late-After`, paired = T)


### plot (Figure 2C)
  wilc_tbl %>% pivot_longer(!contains("Donor")) %>% 
    ggplot(aes(x=name, y = value)) + 
    geom_jitter(width = 0.2, size = 2, shape = 19, col = "grey70" )+
    geom_boxplot(fill = alpha(colour = c("grey","deepskyblue","red","grey") ,alpha = 0.5), outliers = F)+
    scale_x_discrete(limits = c("Pre-Before","Before","After","Late-After"), name = "")+
    ylab("Shannon diversity")+
    theme_bw()+theme(text=element_text(size = 18)) + coord_flip()
        
        
### compute the beta diversity between samples using bray-curtis distances, 
 
  ## normalize all samples to one sample (and not per donor as for Shannon)

  Abundance_tbl_bc = abundance_tbl %>%
  mutate(Donor = ifelse(Donor == "11", "No_Ref_11", Donor)) %>% 
  inner_join(Viral_signal, by = c("Donor","Sampling_day")) %>%  
  filter(!(Donor %in% c("29","No_Ref_11"))) %>%
  mutate(min_read_set = min(mapped_read_count)) %>%
  # rarefy to all samples
  mutate(rary_abundance = ab_abundance*(min_read_set/mapped_read_count),rel_rare_abundance = rary_abundance/size
         ,norm_rel_abund  = round(100*rel_rare_abundance/sum(rel_rare_abundance),9) ) %>% 
  ungroup() %>%
  select(vOTU_ID, Donor, Sampling_day, rel_rare_abundance) %>%
  mutate(Sample_ID = str_c(Donor, Sampling_day, sep = "_") ) %>% select(-Donor, -Sampling_day) %>%
  pivot_wider(names_from = Sample_ID, values_from = rel_rare_abundance)

  ## compute BC distances  
  bray_curt_distances = Abundance_tbl_bc %>% column_to_rownames("vOTU_ID") %>% 
  t(.) %>% vegdist(method = "bray", diag = T)
  
  ## make a heatmap 
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
    # remove Donor_29 and Donor 11 (no sample before AB treatment)
    filter(!(Donor %in% c("29","No_Ref_11"))) %>%
    column_to_rownames("Sample_ID")
  anno_colors = list(Donor= setNames( c(colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set3"))(20),"grey50"),levels(anno_col$Donor)))
  
  # plot using pheatmap package (Figure S5F)
  pheatmap( mat = mat_bray_curt
                     ,cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F
                     ,color = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 5,name = "Greys"))(100))   
                     , border_color = NA 
                     ,annotation_colors = anno_colors
                    ,annotation_col = anno_col, annotation_legend = T )  
  
  #BC tbl 
  bc_distances_df = as.data.frame(as.matrix(bray_curt_distances)) %>% rownames_to_column("Sample_ID_From") %>%
    pivot_longer(cols = !contains("Sample_ID_From"), names_to = "Sample_ID_To", values_to = "BC_distance") %>%
    mutate(Donor_From = str_remove(Sample_ID_From, pattern="_.*"),Donor_To = str_remove(Sample_ID_To, pattern="_.*"),
           Donor_From= ifelse(Donor_From == "No", "No_Ref_11", Donor_From), Donor_To= ifelse(Donor_To == "No", "No_Ref_11", Donor_To),
           Day_From = str_remove(Sample_ID_From, pattern=".*_"),Day_To = str_remove(Sample_ID_To, pattern=".*_")) %>% select(-Sample_ID_From, -Sample_ID_To) %>%
    left_join(donor_treatment %>% rename(AB_From=Treatment), by = c("Donor_From"="Donor")) %>% left_join(donor_treatment%>% rename(AB_To=Treatment), by = c("Donor_To"="Donor"))
  
  
  ## group by the categories and get the mean distances
  
  bc_groups_AB = bc_distances_df %>% 
    # AB values
    filter( AB_From == AB_To) %>%
    group_by(AB_From) %>%
    summarise(mean_BC = mean(BC_distance), Category = "Antibiotic") 
  
  bc_groups_Day = bc_distances_df %>% 
    # remove day -7 and day 90 (lack of samples)
    filter( Day_From == Day_To, !Day_From %in% c("-007","007","090"), !Day_To %in% c("-007","007","090")) %>% 
    # Day values
    group_by(Day_From) %>% summarise(mean_BC = mean(BC_distance), Category = "Day")
    
  bc_groups_Donor = bc_distances_df %>% 
    # remove self comparisons
    filter( Donor_From == Donor_To ) %>%
    # Donor values
    group_by(Donor_From) %>% summarise(mean_BC = mean(BC_distance), Category = "Donor")
  
  bc_groups_plot = bind_rows(bc_groups_AB,bc_groups_Day,bc_groups_Donor) %>%
    select(mean_BC, Category)

  # make plot (Figure 2D) 
  bc_groups_plot %>% ggplot(aes(x=Category, y=mean_BC)) + 
    geom_boxplot(fill = "grey80", outliers = FALSE)+geom_jitter(size = 2, width = 0.2)+
    theme_bw()+ theme(text = element_text(size = 14))
    
  
# same for microbial/ bacterial species (msp)
  
  # first get species name or last entry of the taxonomy 
    extract_last_part = function(string) 
      { parts = strsplit(string, ";")[[1]]; last_part = tail(parts, 1); return(last_part) }
  
    # load overview table
      msp_overview = read_tsv("Table_S6_microbial_species_overview.tsv")  %>%
      mutate(Last_tax_entry = sapply(gtdb_classification, extract_last_part))
    
    # load abundance table, count richness per day
      msp_richness_per_day = read_tsv("Table_S7_microbial_species_abundances.tsv") %>%
        left_join(msp_overview, by = c("id_mgs") ) %>%
         filter(relative_abundance > 0, Sampling_day %in% n_samples_per_donor$Sampling_day) %>%
        mutate(Class = str_extract(gtdb_classification, "(?<=c__)[^;]+"),   
              Sampling_day = case_when(Sampling_day %in% c("015","030") ~"15+30",Sampling_day %in% c("004","007") ~ "4+7",T ~ Sampling_day))%>%
          group_by(Sampling_day,id_mgs,Class) %>%
        summarise() %>% 
        #remove sample of day -7
        group_by(Sampling_day,Class) %>% summarise(msp_per_day=n()) %>% ungroup()
      
      #  make table for plot
      msp_richness_per_day_plot = msp_richness_per_day %>% left_join(host_class_colors, by = c("Class")) %>%
      mutate(Colors = case_when( is.na(Class) ~ "grey40", is.na(Colors) ~ "grey90", T ~ Colors))
      
      # plot (Figure 2A)
       msp_richness_per_day_plot %>% ggplot(aes(x=Sampling_day, y=msp_per_day, fill = Class ))+
        geom_bar(stat = "identity",width = 0.8, col = NA)+xlab("Day")+ylab("Richness msp")+
        scale_fill_manual(values = category_colors)+
        scale_x_discrete(limits = c("-014","-001","4+7","010","15+30","180")
                        ,labels = c("-14","-1","4+7","10","15+30","180"))+
        theme_bw() + 
          theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") +
        guides(color = guide_legend(ncol = 2)) 

#### same for bacteria
  
   msp_overview = msp_overview %>% mutate(Last_tax_entry = str_remove(pattern="s__",extract_last_part(gtdb_classification)))
  
   msp_abundance = read_tsv("Table_S7_microbial_species_abundances.tsv") %>%
        left_join(msp_overview, by = c("id_mgs") ) %>%
         filter(relative_abundance > 0) %>%
        mutate(abdundance_metric = ifelse(relative_abundance> 25, Last_tax_entry, "Low_abundance")) 
  
  
  ###### Shannon diversity per sample
  shannon_per_sample_msp = msp_abundance %>% 
    mutate(Sample_ID = str_c(Donor, Sampling_day, sep = "_")) %>% 
    select(id_mgs, Sample_ID, relative_abundance) %>%
    pivot_wider(names_from = Sample_ID, values_from = relative_abundance, values_fill = 0) %>% 
   column_to_rownames("id_mgs") %>% summarise_all(.funs = "diversity") %>% 
   t() %>% as.data.frame() %>% 
   rownames_to_column("Sample_ID") %>% rename(Shannon=V1) %>%
   
   mutate(Donor = str_remove(Sample_ID, pattern="_.*"), Sampling_day = str_remove(Sample_ID, pattern = ".*_")) %>%
   left_join(donor_treatment) %>%
  mutate(AB_treatment = case_when(Sampling_day %in% c("-001") ~ "Before"
                                ,Donor == "11" & Sampling_day == "-014" ~ "Before"
                                ,Sampling_day %in% c("004","007") ~ "After"
                                ,Donor == "13" & Sampling_day == "010" ~ "After", 
                                T ~ "Other"))

    # group donors in increasing and decreasing groups
    shannon_groups_msp = shannon_per_sample_msp %>% 
      filter(AB_treatment != "Other") %>% 
      select(-Sampling_day, -Sample_ID) %>%
      pivot_wider( names_from = AB_treatment, values_from = Shannon) %>%
      mutate(Shannon_change = ifelse(After>Before, "Increase","Decrease"))
    
    table(shannon_groups_msp$Shannon_change)
       
    # normalize shannon to one day
    shannon_of_min_day_msp = shannon_per_sample_msp %>% group_by(Donor) %>% filter(Sampling_day %in% c("-001")) %>% 
       select(Donor,min_shannon=Shannon) %>% ungroup()  
       
        #norm all shannon values
        Norm_shannon_msp = shannon_per_sample_msp %>% 
          inner_join(shannon_of_min_day_msp) %>% 
          mutate(Norm_shannon = Shannon/min_shannon) %>%
          left_join(shannon_groups_msp %>% select(Donor,Shannon_change), by = "Donor")
    
    ## plot norm shannon (Figure S5E)
      Norm_shannon_msp %>% ggplot(aes(x=Sampling_day, y=Norm_shannon))+
      geom_line(linewidth = 1, aes(group = Donor, color = Donor))  + 
      scale_color_manual(values=  colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set3"))(20))+
      geom_point(size = 2, color="grey70") +
      xlab("Day")+ylab("Shannon (norm.)\ndiversity")+
      scale_x_discrete(limits = c("-014","-007","-001","004","007","010","015","030","090","180")
                       ,labels = c("-14","-7","-1","4","7","10","15","30","90","180"))+
      theme_bw() + theme(text = element_text(size = 14), strip.background = element_blank(), legend.position = "right") +
      guides(color = guide_legend(ncol = 2)) +
      facet_wrap(~Shannon_change, ncol = 2)  
      
  ## test for differences
  wilc_tbl_msp = Norm_shannon_msp %>% mutate( 
  TimeFrame = case_when(Sampling_day %in% c("-014","-007") ~ "Pre-Before", 
                        #Sampling_day %in% c("030") ~ "Short-After",
                        Sampling_day %in% c("180") ~ "Late-After",
                        T ~ AB_treatment )) %>%
  filter(TimeFrame !="Other") %>%
  select(Donor, TimeFrame, Shannon) %>%
  pivot_wider(names_from = TimeFrame, values_from = Shannon)  %>% drop_na()
    
####
  
  # wilcoxon tests

  # pre-before vs before
  wilcox.test(wilc_tbl_msp$`Pre-Before`, wilc_tbl_msp$Before, paired = TRUE)

  # before vs after
  wilcox.test(wilc_tbl_msp$Before, wilc_tbl_msp$After, paired = TRUE)
  
  # pre-before vs after
  wilcox.test(wilc_tbl_msp$`Pre-Before`, wilc_tbl_msp$After, paired = T)
  
  # pre-before vs late-after
  wilcox.test(wilc_tbl_msp$`Pre-Before`, wilc_tbl_msp$`Late-After`, paired = T)
  
  # Before vs Late-After
  wilcox.test(wilc_tbl_msp$Before, wilc_tbl_msp$`Late-After`, paired = T)
  
  # After vs Late-After
  wilcox.test(wilc_tbl_msp$After, wilc_tbl_msp$`Late-After`, paired = T)
  
  ### plot (Figure S5D)
  wilc_tbl_msp %>% pivot_longer(!contains("Donor")) %>% 
    ggplot(aes(x=name, y = value)) + 
    geom_jitter(width = 0.2, size = 2, shape = 19, col = "grey70" )+
    geom_boxplot(fill = alpha(colour = c("grey","deepskyblue","red","grey") ,alpha = 0.5), outliers = F)+
    scale_x_discrete(limits = c("Pre-Before","Before","After","Late-After"), name = "")+
    ylab("Shannon diversity")+
    theme_bw()+theme(text=element_text(size = 18))
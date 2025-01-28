####
library(tidyverse)
library(readxl)

### abundances of vOTUs and microbial species

# read signal overview
Read_signal_grouped = read_csv("Table_S2_overview_readsignal.csv")  

  # take only the viral
  Viral_signal = Read_signal_grouped %>% filter(Sequence_set == "viral") %>% filter(mapped_read_count > 1e6, `fraction in %` > 10) %>% 
    select(donor, day, mapped_read_count, fraction=`fraction in %`) %>%
    mutate(Sampling_day = ifelse( grepl(day,pattern="^-"), sprintf("%+04d", day), sprintf("%03d", day))
      ,Donor = ifelse(nchar(donor)==1, sprintf("%02d", donor), donor)) %>% select(-donor, -day)


# All vOTUs over time

    ### path to vOTUs table
    vOTU_tbl = read_excel("Table_S1_vOTUs_final.xlsx") 
    
    #### path to abundance table
    abundance_vOTU_tbl = read_tsv("Table_S4_abundance_tbl_final_pub.tsv") %>%
       mutate( Donor = ifelse(Donor == "11", "No_Ref_11", Donor),vOTU_Donor = str_c(vOTU_ID, Donor, sep="_")) %>%
       inner_join(Viral_signal, by = c("Donor","Sampling_day")) %>%
      group_by(Donor)  %>%
      mutate(Sample_ID = str_c(Donor, Sampling_day, sep = "_")
           # down size to the smallest count per donor
          ,min_read_set = min(mapped_read_count)
          ,rary_abundance = ab_abundunance*(min_read_set/mapped_read_count)
           # normalise to the size
           ,norm_rare_abundance = rary_abundance/size) %>% 
            group_by(Sample_ID) %>% mutate(
           # compute relative abundance (sum = 100) per sample 
           rel_abund  = round(100*norm_rare_abundance/sum(norm_rare_abundance),9)) %>%
      select(vOTU_ID, Sample_ID,Sampling_day,Donor, rel_abund) %>%
      group_by(vOTU_ID,Donor) %>% 
      mutate(max_abdundance = max(rel_abund) ,abdundance_metric= ifelse(max_abdundance> 25, vOTU_ID, "Low_abundance"))
  
      # plot
     abundance_vOTU_tbl %>% filter(rel_abund > 0) %>% 
      group_by(Donor) %>%
      
      ggplot(aes(x=Sampling_day, y = rel_abund, group = vOTU_ID, color = abdundance_metric)) +
      geom_point(size = 3)+
      geom_line(linewidth = 1.1)+
      theme_bw()+ylab("Relative Abundance %")+xlab("Sample/Day")+
      scale_color_manual(values = c(alpha("grey", alpha = 0.2), rev(colorRampPalette(RColorBrewer::brewer.pal(n=8, name = "Set1"))(31))))+
      scale_x_discrete(limits = c("-014","-007","-001","004","007","010","015","030","090","180")
                       ,labels = c("-14","","-1","4","","10","","30","","180"))+
      scale_y_continuous(breaks = c( seq(0,110,20)), limits = c(0,100))+
      theme(legend.position = "none", text = element_text(size = 14), axis.text = element_text(size = 10)
            , strip.background = element_blank(),panel.grid = element_blank()) +
      facet_wrap(~Donor, scales = "fixed", ncol = 7 )
  
     
  # define tbl to collect p-values   
  pvalue_treshold_df = data.frame(NULL)
  
  # get number of samples per day
  samples_per_day = abundance_vOTU_tbl %>% group_by(Donor,Sampling_day) %>% summarise() %>% 
    group_by(Sampling_day) %>% summarise(n_samples_per_day=n()) %>% ungroup()
  
  #for (i in seq(1,50,1)) {
  i = 25
    abundant_vOTUs = abundance_vOTU_tbl %>% filter(norm_rel_abund > i)
  
    # table for test
    abundant_votus_test = abundant_vOTUs %>% 
      left_join(samples_per_day , by = "Sampling_day") %>% 
      group_by(Sampling_day,n_samples_per_day) %>%  summarise(n_ab_vOTUs = n()) %>% 
      mutate(n_votus_norm = round(n_ab_vOTUs/n_samples_per_day,3)) %>% ungroup() %>% 
      # consider only days with at least 3 samples
      filter(n_samples_per_day>3)
    
      # test using poisson distribution
      poisson_model = glm(n_ab_vOTUs ~ Sampling_day, family = poisson(link = "log"), offset = log(n_samples_per_day), data = abundant_votus_test)
      
      # Display summary of the model
      sum_poission = summary(poisson_model)
      
      # Extract coefficients and p-values
      coef_summary = as.data.frame(summary(poisson_model)$coefficients) %>% rownames_to_column("Category") %>% mutate(cutoff=i)
      # add to p-value dataframe
      pvalue_treshold_df = pvalue_treshold_df %>% bind_rows(coef_summary)
  #  }

 # plot pvalue against cut offs of abundance
  pvalue_treshold_df %>% ggplot(aes(x=cutoff,group=Category, y = `Pr(>|z|)`)) +
    geom_point(col="mediumaquamarine")+geom_line(col="mediumaquamarine")+
    geom_line(y=0.05, col = "red", linetype = "dashed")+
    theme_bw()+ylab("pvalue")+xlab("Cutoff, abundance (%)")+
    scale_x_continuous(breaks = seq(0,50,5), limits = c(0,40))+
    scale_y_continuous(breaks = c(0.05, 0.15,0.5))+
    theme(legend.position = "top", text = element_text(size = 14), strip.background = element_blank(),panel.grid = element_blank()
        ,legend.title = element_blank()) + facet_wrap(~Category)

     
 ## microbial species
 
     # define function to extract last part
     extract_last_part = function(string) { parts = strsplit(string, ";")[[1]]; last_part = tail(parts, 1); return(last_part) }
    
      msp_overview = read_tsv("Table_S5_microbial_species_overview.tsv") %>%
      mutate(Last_tax_entry = sapply(gtdb_classification, extract_last_part))
    
      msp_abundance = read_tsv("Table_S6_microbial_species_abundances.tsv") %>%
        left_join(msp_overview, by = c("id_mgs") ) %>%
        group_by(id_mgs,Donor) %>%
        mutate(max_abdundance = max(relative_abundance) ,abdundance_metric= ifelse(max_abdundance> 25, Last_tax_entry, "Low_abundance")) %>%
         filter(relative_abundance > 0) 
    
      # plot
      msp_abundance %>% ggplot(aes(x=Sampling_day, y = relative_abundance, group = Last_tax_entry, color = abdundance_metric)) +
      geom_vline(xintercept = "-001", color = alpha("lightseagreen", alpha =0.5), size =2)+
      geom_point(size = 3)+
      geom_line(linewidth = 1.1)+
      theme_bw()+ylab("Relative Abundance %")+xlab("Sample/Day")+
      scale_color_manual(values = c(alpha("grey", alpha = 0.2), colorRampPalette(RColorBrewer::brewer.pal(n=8, name = "Paired"))(13)))+
      scale_x_discrete(limits = c("-014","-007","-001","004","007","010","015","030","090","180")
                       ,labels = c("-14","","-1","4","","10","","30","","180"))+
      scale_y_continuous(breaks = c( seq(0,110,20)), limits = c(0,100))+
      theme(legend.position = "bottom", text = element_text(size = 14), axis.text = element_text(size = 10)
            , strip.background = element_blank(),panel.grid = element_blank()) +
      facet_wrap(~Donor, scales = "fixed", ncol = 7 )+
        guides(color = guide_legend(nrow =  5), shape = guide_legend(label = F))
      
  
  ### Dynamics of Parabacteroides distasonis_A (msp_0012) and its phages

  votus_infecting_msp12 = vOTU_tbl %>% 
    filter(grepl(host_taxonomy_lineage, pattern="Parabacteroides distasonis_A")) 
  
  msp12_abundances = abundance_vOTU_tbl %>% 
    semi_join(votus_infecting_msp12, by = c("vOTU_ID")) %>%
    select(Species=vOTU_ID, Donor, Sampling_day, rel_abund) %>% mutate(Category = "vOTU") %>%
    bind_rows(msp_abundance %>% filter(id_mgs =="msp_0012") %>%
    select(Species=id_mgs, Donor,Sampling_day, rel_abund=relative_abundance) %>% 
      mutate(Category = "Bacterial_species")) %>%
    group_by(Donor) %>% mutate(max_abund_donor = max(rel_abund)) %>%
    group_by(Species,Donor) %>% mutate(max_abund_species = max(rel_abund)) %>%
    # take species with least with low abundance 
    filter(max_abund_donor>1,max_abund_species>0.1, Donor %in% c("22","24","01","17"))
  
  # plot    
  msp12_abundances %>% ggplot(aes(x=Sampling_day, y=rel_abund, color = Species, group=Species, shape = Category))+
    geom_point(size = 3)+geom_line(size =1.5)+
    facet_wrap(~Donor, nrow = 2)+
    theme_bw()+
    scale_shape_manual(values = c(15,19))+
    scale_color_manual(values = c("grey60",RColorBrewer::brewer.pal(n = 5, name = "Set1")))+
    xlab("Sampling Day")+ylab("Relative abundance %")+ylim(0,80)+
     scale_x_discrete(limits = c("-014","-007","-001","004","007","010","015","030","090","180")
                    ,labels = c("-14","-7","-1","4","7","10","15","30","90","180"))+
     theme(text=element_text(size = 14), strip.background = element_blank()
          ,legend.position = "top",panel.grid =  element_blank())+
          guides(color = guide_legend(ncol =  2), shape = guide_legend(label = F))
  
  
  

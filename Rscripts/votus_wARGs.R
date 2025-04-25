####
library(tidyverse)
library(readxl)

### AMR votus

# AMR vOTUs over time

### path to vOTUs table
vOTU_tbl = read_tsv("Table_S1_vOTUs_final_v2.tsv") 

# read signal overview
Read_signal_grouped = read_csv("Table_S2_overview_readsignal.csv")  

  # take only the viral
  Viral_signal = Read_signal_grouped %>% filter(Sequence_set == "viral") %>% filter(mapped_read_count > 1e6, `fraction in %` > 10) %>% 
    select(donor, day, mapped_read_count, fraction=`fraction in %`) %>%
    mutate(Sampling_day = ifelse( grepl(day,pattern="^-"), sprintf("%+04d", day), sprintf("%03d", day))
      ,Donor = ifelse(nchar(donor)==1, sprintf("%02d", donor), donor)
      ,Donor = ifelse(Donor == "11", "No_Ref_11", Donor)) %>% select(-donor, -day)

#### path to abundance table
abundance_amr_tbl = read_tsv("Table_S5_abundance_tbl_final_pub_v2.tsv") %>%
  # filter for viral signal and compute relative abundance
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
           rel_abund  = round(100*norm_rare_abundance/sum(norm_rare_abundance),9),
          Donor = ifelse(Donor == "11", "No_Ref_11", Donor),vOTU_Donor = str_c(vOTU_ID, Donor, sep="_")) %>%
  # look only on vOTUs with ARGs, and only the ones, that are present in samples
  inner_join(vOTU_tbl %>% filter(grepl(With_ARG,pattern="Yes"))) %>% filter(rel_abund > 0)


# plot vOTUs (Figure 3)
abundance_amr_tbl %>%
  ggplot(aes(x=Sampling_day, y = rel_abund, group = vOTU_Donor, color = Donor)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.5) +
  theme_bw()+ylab("Fraction %")+xlab("Sample/Day")+
  scale_color_manual(values = c(colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set3"))(20),"grey50"))+
  #scale_y_continuous( limits = c(0,0.8), breaks = c(0, 0.2, 0.4, 0.6, 0.8))+
  scale_x_discrete(limits = c("-014","-007","-001","004","007","010","015","030","090","180")
                   ,labels = c("-14","-7","-1","4","7","10","15","30","90","180"))+
  theme(legend.position = "bottom", text = element_text(size = 14), strip.background = element_blank(),panel.grid = element_blank()) +
  guides(color = guide_legend(nrow = 2))+
  facet_wrap(~vOTU_ID, scales = "fixed", nrow = 3)


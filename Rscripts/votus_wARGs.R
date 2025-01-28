####
library(tidyverse)
library(readxl)

### AMR votus

# AMR vOTUs over time

### path to vOTUs table
vOTU_tbl = read_excel("NutriPhage/Projects/predires_phageome_22_volunteers/data/Table_S1_vOTUs_final.xlsx") 

#### path to abundance table
abundance_amr_tbl = read_tsv("NutriPhage/Projects/predires_phageome_22_volunteers/data/Table_S4_abundance_tbl_final.tsv") %>%
  inner_join(vOTU_tbl %>% filter(grepl(With_ARG,pattern="Yes"))) %>%
  mutate( Donor = ifelse(Donor == "11", "No_Ref_11", Donor),vOTU_Donor = str_c(vOTU_ID, Donor, sep="_"))

# plot vOTUs
abundance_amr_tbl %>%
  ggplot(aes(x=Sampling_day, y = norm_rel_abund, group = vOTU_Donor, color = Donor)) +
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


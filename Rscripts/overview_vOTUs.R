### libraries

library(tidyverse)
library(readxl)
library(UpSetR)

# path to vOTUs table
vOTU_tbl = read_excel("Table_S1_vOTUs_final.xlsx") 

### Source plot (Figure 1A)
   vOTU_tbl %>% dplyr::count(source, name = "counts") %>% mutate(Name="source") %>%
      ggplot(aes(x=Name, y=counts, fill=source)) + geom_col(col="black")+
      scale_fill_manual(values = alpha(c("grey40","white","royalblue"), alpha = 0.8))+
      theme_bw()+theme(text = element_text(size = 14))+xlab("")+ylab("#counts")+coord_flip()

### checkV_quality plot (Figure 1B)

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
      main.bar.color = "skyblue",
      sets.bar.color = "orange",
      text.scale = 1.8,
      point.size = 6)

# Sizes and qualities (Figure 1D)

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
    
# Temperate phages (Figure 1E)

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

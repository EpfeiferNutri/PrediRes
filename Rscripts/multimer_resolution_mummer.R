### load required libaries ###

library("data.table")
library("foreach")
library("RColorBrewer")
library("seqinr")
library("tidyverse")

## detect potential repeats by using the blastn output 
col_blastn = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen")

# compare DNA sequences of contigs by blastn (all vs all). Consider only the ones, covering more than 110 %

# read in blastn output
blastn_out = fread("path_to_blast.out"
  , col.names =  col_blastn)%>% 
  #consider only self-hits
  filter(qseqid == sseqid) %>%
  group_by(qseqid,sseqid,qlen,slen) %>% 
  summarise(length_align = sum(length), n_seg = n()) %>%
  mutate(q_cov = round(100*length_align/qlen,1), s_cov = round(100*length_align/slen,1)) %>%
  # filter out sequences that have a self-match larger 10% of own size
  filter(q_cov > 110)

### load in contigs with sequence information (ID, length, sequence)
  contigs_all = read.fasta("path_to_fasta"
                                  , as.string = T, forceDNAtolower = F)

  contigs_all_df = data.frame( contig = unlist(lapply(contigs_all,getName)), contig_size = unlist(lapply(contigs_all,getLength))
                           ,sequence = unlist(lapply(contigs_all,getSequence, as.string = T)) 
                           ); rm(contigs_all)
  
### safe potential multimers as single fasta files for MUMmer comparison

# make a dataframe with all potential contigs, and safe them as single fasta files over a loop
mutlimer_contigs = blastn_out %>% select(contig=qseqid) %>% left_join(contigs_all_df)

 foreach(i = seq_along(mutlimer_contigs$contig)) %do% {

   tmp_seq = mutlimer_contigs %>% filter(contig == mutlimer_contigs$contig[i])
   
   write.fasta( sequences = as.list(tmp_seq$sequence), names = as.list(tmp_seq$contig)
                ,file.out = str_c("path_to_safe_fasta",tmp_seq$contig,".fasta"))
 }

# use exact-tandems functions of MUMmer to detect sequences with tandem repeats larger > 1000bo
 # exact-tandems multimers/file.fasta 1000 > tandems/out.tsv

list_files = dir("path_to_exact_tandems/")

tandems = data.frame(NULL)
foreach(i = seq_along(list_files)) %do% {
#i=1
tmp_file = fread(str_c("path_to_exact_tandems/",list_files[i])) %>%
  mutate(contig = str_remove(list_files[i], pattern = ".tsv$"))
tandems = tandems %>% bind_rows(tmp_file)
}

tandems_to_be_checked = tandems %>% left_join(contigs_all_df) %>% mutate(size_ratio =  contig_size/UnitLen )

### check by using MAFFT homepage (https://mafft.cbrc.jp/alignment/software/) for multimeric sequences

  # remove the the ones without tandem
  tandems_verified = tandems_to_be_checked %>% 
  filter( !contig %in% c("UHGV-1325805_15880")) %>% 
  group_by(contig) %>% 
  # give tandems repeat a number
  mutate(tanden_nr = 1:n())

# Extract one copy of the sequence and remove the tandem one. Then safe as fasta sequence.

foreach(i=seq_along(tandems_verified$contig)) %do% {

tmp_contigs_all = read.fasta(str_c("path_to_fasta/",tandems_verified$contig[i],".fasta")
                                  , as.string = F, forceDNAtolower = F)
tmp_contigs_all_df = data.frame( contig = unlist(lapply(tmp_contigs_all,getName)), contig_size = unlist(lapply(tmp_contigs_all,getLength))
                           ,sequence = unlist(lapply(tmp_contigs_all,getSequence, as.string = F)) ) %>%  
  rowid_to_column("Position") %>% 
  filter(Position > tandems_verified$Start[i], Position < (tandems_verified$Start[i]+tandems_verified$UnitLen[i])) %>% 
  summarise( contig = str_c(tandems_verified$contig[i],"__notan_",tandems_verified$tanden_nr[i],"_",tandems_verified$Start[i],"_to_",(tandems_verified$Start[i]+tandems_verified$UnitLen[i]))
             ,sequence = str_c(sequence, collapse = ""), contig_size = n())
  
write.fasta( sequences = as.list(tmp_contigs_all_df$sequence), names = as.list(tmp_contigs_all_df$contig)
                ,file.out = str_c("path_to_safe",tmp_contigs_all_df$contig,".fasta"))
}

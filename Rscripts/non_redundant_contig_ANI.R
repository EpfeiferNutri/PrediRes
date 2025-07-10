##### libraries #####
library("data.table")
library("igraph")
library("foreach")
library("seqinr")
library("tidyverse")

# de-cluster contigs afer BLASTN + ANI using checkV scripts

# load in all contigs
contigs_all = read.fasta("/path-to/contigs.fna"
                                  , as.string = T, forceDNAtolower = F)

contigs_all_df = data.frame( contig = unlist(lapply(contigs_all,getName)), contig_size = unlist(lapply(contigs_all,getLength))
                           ,sequence = unlist(lapply(contigs_all,getSequence, as.string = T)) 
                           ); rm(contigs_all)

# # read in ANI tsv, computed with ani-calc.py of checkV
 ANI_df = fread("path-to/contig-ANI.tsv")
 
 ANI_processed = ANI_df %>% filter(pid > 95) %>% filter(qname != tname) %>%
   # add sizes
   left_join(contigs_all_df %>% select(qname=contig, qsize = contig_size), by = "qname") %>%
   left_join(contigs_all_df %>% select(tname=contig, tsize = contig_size), by = "tname") %>% 
   # remove hits that cover less than 85% from the shorter contig
   filter( (qcov > 85 & qsize <= tsize) | (tcov > 85 & tsize <= qsize) )
 
 write_tsv(ANI_processed, "path-to-ANI.processed.tsv")

ANI_processed = fread("path-to-ANI.processed.tsv")

### make a graph based on ANI
Edges_df = ANI_processed %>% select(qname, tname, pid, qcov, tcov)

Vert_df = Edges_df %>% select(contig_id=qname) %>% distinct() %>% bind_rows(
  Edges_df %>% select(contig_id=tname) %>% distinct()) %>% distinct() %>% 
  left_join(contigs_all_df, by = c("contig_id"="contig"))

#make graph
graph_network = graph_from_data_frame(vertices = Vert_df, d = Edges_df, directed = F)

# # safe graph
#   write.graph(graph_network, format = "graphml", file = "/path to safe/")

# single linkage clustering
single_linkage_clusters = components(graph_network)

clusters = data.frame(one_link_cluster = single_linkage_clusters$membership) %>% rownames_to_column("contig") %>% 
  group_by(one_link_cluster) %>% mutate(cluster_size = n()) %>% ungroup() %>% 
  left_join(contigs_all_df %>% select(contig, contig_size))

# make a loop to iterively remove already clustered sequence and check if remaining once are clustered
  increment_size = 3e3
  max_contig_size = max(clusters$contig_size)
  #max_contig_size = 10e4
  breaks = seq(3000, max_contig_size , by = increment_size)
  labels = paste(breaks[-1], breaks[-length(breaks)], sep = " - ")
  
  pool_cluster_size_typed = clusters %>% 
    mutate(size_type = cut(contig_size, breaks = breaks, labels = labels, right = FALSE) 
           ,size_type = ifelse( is.na(size_type), str_c(max_contig_size), str_c(size_type)))
  pool_edges = Edges_df
  
# step-wise sub cluster 

  redundant_edges = data.frame(NULL)
  kept_nodes = data.frame(NULL)
  # get the order of the size types
  sizes_types_sorted = unique((pool_cluster_size_typed %>% arrange(desc(contig_size)))$size_type)
  
  #loop
  for(i in seq_along(sizes_types_sorted) ) {
    # take the nodes with the specific size type
    tmp_nodes_clusters = pool_cluster_size_typed %>% filter(size_type==max(sizes_types_sorted[i]))
    
    # take their edges
    tmp_edges = pool_edges %>% semi_join(tmp_nodes_clusters, by = c("qname"="contig")) %>%
      semi_join(tmp_nodes_clusters, by = c("tname"="contig"))
    # make network
    tmp_net = graph_from_data_frame(vertices = tmp_nodes_clusters, d = tmp_edges, directed = F)
    # test if they are clustered
    tmp_clusters = data.frame(tmp_one_link_cluster = components(tmp_net)$membership) %>% rownames_to_column("contig") %>%
      group_by(tmp_one_link_cluster) %>% mutate(cluster_size = n()) %>% ungroup() %>%
      left_join(tmp_nodes_clusters%>% select(contig, one_link_cluster), by = "contig")
    # take largest per cluster
    tmp_one_per_cluster = tmp_clusters %>% left_join(contigs_all_df, by = "contig")%>%
      group_by(tmp_one_link_cluster) %>% arrange(desc(contig_size)) %>% slice(1) %>% ungroup() %>%
      mutate(tmp_one_link_cluster_stp2 = str_c(one_link_cluster,tmp_one_link_cluster,cluster_size,sep="_" ))  
    # keep the largest per cluster  
    kept_nodes = kept_nodes %>% bind_rows(tmp_one_per_cluster)  
    # remove any linked nodes and edges from the network that are smaller since they are redundant
      # get edges to redundant nodes
      tmp_redundant_edges = pool_edges %>% 
                semi_join(tmp_nodes_clusters,by = c("qname"="contig")) %>%
                semi_join( pool_cluster_size_typed %>% filter(size_type!=max(sizes_types_sorted[i])), by = c("tname"="contig"))
      # remove from network
      pool_cluster_size_typed = pool_cluster_size_typed %>% anti_join(tmp_nodes_clusters, by = "contig") %>%
        anti_join(tmp_redundant_edges, by = c("contig"="tname") )
      pool_edges = pool_edges %>% anti_join(tmp_redundant_edges)
      # keep the redundant nodes
      redundant_edges = redundant_edges %>% bind_rows(tmp_redundant_edges)
      # repeat
        }
    
# get non redundant contigs
non_redundant_contigs = contigs_all_df %>% anti_join(clusters,by = c("contig")) %>%
  bind_rows(kept_nodes) %>%  select(contig, contig_size, sequence) 

#safe as fastas
write.fasta(sequences = as.list(non_redundant_contigs$sequence), names = as.list(non_redundant_contigs$contig), 
               file.out = "non_redundant.fasta")  



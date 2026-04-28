##### libraries #####
library("igraph")
library("seqinr")
library("tidyverse")

# de cluster contigs afer BLASTN + ANI using checkV scripts

contig_list =  dir(pattern = "*.fna|*.fasta")
# load in all contigs

contigs_all = read.fasta(file = contig_list, as.string = T, forceDNAtolower = F)
contigs_all_df = data.frame( contig = unlist(lapply(contigs_all,getName)), contig_size = unlist(lapply(contigs_all,getLength))
                           #,sequence = unlist(lapply(contigs_all,getSequence, as.string = T)) 
                           ) %>% filter(contig_size > 3e3)
## read in ANI tsv
  ANI_tbl = dir(pattern = "*.ANI.tsv|*anicalc.tsv") 
  ANI_df = read_tsv(ANI_tbl)

ANI_processed = ANI_df %>% filter(pid > 95) %>% filter(qname != tname) %>%
  # add sizes
  left_join(contigs_all_df %>% select(qname=contig, qsize = contig_size), by = "qname") %>%
  left_join(contigs_all_df %>% select(tname=contig, tsize = contig_size), by = "tname") %>%
  # remove hits that cover less than 85% from the shorter contig
  filter( (qcov > 85 & qsize <= tsize) | (tcov > 85 & tsize <= qsize) )

### make a graph based on ANI
Edges_df = ANI_processed %>% select(qname, tname, pid, qcov, tcov)

Vert_df = Edges_df %>% select(contig_id=qname) %>% distinct() %>% bind_rows(
  Edges_df %>% select(contig_id=tname) %>% distinct()) %>% distinct() %>% 
  left_join(contigs_all_df, by = c("contig_id"="contig")) %>% filter(contig_size > 3e3)

#make graph
graph_network = graph_from_data_frame(vertices = Vert_df, d = Edges_df, directed = F)

# plot(graph_network, vertex.label = NA)
# write.graph(graph_network, format = "graphml" ,file="graph.graphml")

# single linkage clustering
single_linkage_clusters = components(graph_network)

#get clusters
clusters = data.frame(one_link_cluster = single_linkage_clusters$membership) %>% rownames_to_column("contig") %>% 
  group_by(one_link_cluster) %>% mutate(cluster_size = n()) %>% ungroup() %>% 
  left_join(contigs_all_df %>% select(contig, contig_size))

# make a loop to iterively remove already clustered sequence and check if remaining once are clustered
  increment_size = 3e3
  max_contig_size = max(clusters$contig_size)
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
  #i=39
    # take the nodes with the specific size type
    tmp_nodes_clusters = pool_cluster_size_typed %>% filter(size_type==max(sizes_types_sorted[i]))
    
    # take their edges
    tmp_edges = pool_edges %>% semi_join(tmp_nodes_clusters, by = c("qname"="contig")) %>%
     bind_rows(pool_edges %>% semi_join(tmp_nodes_clusters, by = c("tname"="contig")))
    # make a vert df
    tmp_vert = tmp_nodes_clusters %>% select(contig) %>% 
      bind_rows(tmp_edges %>% select(contig=tname), tmp_edges %>% select(contig=qname)) %>% distinct()
    # make network
    tmp_net = graph_from_data_frame(vertices = tmp_vert, d = tmp_edges, directed = F)
    # test if they are clustered
    tmp_clusters = data.frame(tmp_one_link_cluster = components(tmp_net)$membership) %>% rownames_to_column("contig") %>%
      group_by(tmp_one_link_cluster) %>% mutate(cluster_size = n()) %>% ungroup() 
    # take largest per cluster
    tmp_one_per_cluster = tmp_clusters %>% left_join(contigs_all_df, by = "contig")%>%
      group_by(tmp_one_link_cluster) %>% arrange(desc(contig_size)) %>% slice(1) %>% ungroup()  
    # keep the largest per cluster  
    kept_nodes = kept_nodes %>% bind_rows(tmp_one_per_cluster)  
    # remove any linked nodes and edges from the network that are smaller since they are redundant
      
  # get redundant nodes
      tmp_redundant_nodes = bind_rows(pool_edges %>% semi_join(tmp_one_per_cluster,by = c("qname"="contig")) %>% select(contig=tname), 
                                      pool_edges %>% semi_join(tmp_one_per_cluster,by = c("tname"="contig")) %>% select(contig=qname)) %>% distinct()
      
  # remove from network
      # from the nodes
      pool_cluster_size_typed = pool_cluster_size_typed %>% anti_join(bind_rows(tmp_one_per_cluster,tmp_redundant_nodes), by = "contig") 
      # and the edges
      pool_edges = pool_edges %>% 
        # bind_rows(
        anti_join( bind_rows(tmp_redundant_nodes,tmp_one_per_cluster), by = c("qname"="contig")) %>%
        anti_join( bind_rows(tmp_redundant_nodes,tmp_one_per_cluster), by = c("tname"="contig"))
      # keep the nodes that are now not connected anymore
      tmp_unredundant_nodes = tmp_nodes_clusters %>% anti_join(bind_rows(tmp_redundant_nodes,tmp_one_per_cluster), by = "contig")
      pool_cluster_size_typed = pool_cluster_size_typed %>% bind_rows(tmp_unredundant_nodes)
      # repeat
      }

kept_nodes = kept_nodes %>% bind_rows(pool_cluster_size_typed %>% distinct(contig))   
  
# make the contigs tbl again, but add sequence information
contigs_all_df = data.frame( contig = unlist(lapply(contigs_all,getName)), contig_size = unlist(lapply(contigs_all,getLength))
                           ,sequence = unlist(lapply(contigs_all,getSequence, as.string = T))) %>%
  filter(contig_size > 3e3)

# get non redundant contigs
non_redundant_contigs = contigs_all_df %>% anti_join(clusters,by = c("contig")) %>%
  bind_rows(kept_nodes) %>%  select(contig) %>% left_join(contigs_all_df) 

# get the cluster
cluster_sequences = ANI_processed %>% select(contig=qname, tname) %>%
  semi_join(kept_nodes, by = c("contig")) %>%
  bind_rows(ANI_processed %>% select(contig=tname,tname=qname) %>% semi_join(kept_nodes,by = c("contig"))) %>%
  distinct() %>%  rename(kept_node=contig, member=tname) %>%
  group_by(kept_node) %>% mutate(n_related_sequences=n()) %>% 
  group_by(member) %>% mutate(n_related_kept_nodes=n()) %>%
  ungroup()

#safe as fastas
write.fasta(sequences = as.list(non_redundant_contigs$sequence), names = as.list(non_redundant_contigs$contig), 
               file.out = "non_redundant_contigs.fasta")

## safe sequences that cluster together
cluster_sequences %>% write_tsv("clusters.tsv")

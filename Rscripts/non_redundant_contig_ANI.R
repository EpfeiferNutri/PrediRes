##### libraries #####
library("data.table")
library("igraph")
library("foreach")
library("seqinr")
library("tidyverse")

# de-cluster contigs afer BLASTN + ANI using checkV scripts

# load in all contigs
contigs_all = read.fasta("/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/non_viral_v2/removing_redundancy/Non_viral_pre_cluster.fna"
                                  , as.string = T, forceDNAtolower = F)

contigs_all_df = data.frame( contig = unlist(lapply(contigs_all,getName)), contig_size = unlist(lapply(contigs_all,getLength))
                           ,sequence = unlist(lapply(contigs_all,getSequence, as.string = T)) 
                           ); rm(contigs_all)

# # read in ANI tsv
# ANI_df = fread("/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/non_viral_v2/removing_redundancy/non_viral.ANI.tsv")
# 
# ANI_processed = ANI_df %>% filter(pid > 95) %>% filter(qname != tname) %>%
#   # add sizes
#   left_join(contigs_all_df %>% select(qname=contig, qsize = contig_size), by = "qname") %>%
#   left_join(contigs_all_df %>% select(tname=contig, tsize = contig_size), by = "tname") %>% 
#   # remove hits that cover less than 85% from the shorter contig
#   filter( (qcov > 85 & qsize <= tsize) | (tcov > 85 & tsize <= qsize) )
# 
# write_tsv(ANI_processed, "/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/non_viral_v2/removing_redundancy/non_viral.ANI.processed.tsv")

ANI_processed = fread("/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/non_viral_v2/removing_redundancy/non_viral.ANI.processed.tsv")

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
               file.out = "/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/non_viral_v2/non_redudant_non_viral_final_v2.fasta")  



### check again checkV clustering

# load in all contigs
contigs_all = read.fasta("/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/cleaning_redundancy_vOTUs/viral_contigs_ALL_pooled_no_checkv.fna"
                                  , as.string = T, forceDNAtolower = F)

contigs_all_df = data.frame( contig = unlist(lapply(contigs_all,getName)), contig_size = unlist(lapply(contigs_all,getLength))
                           ,sequence = unlist(lapply(contigs_all,getSequence, as.string = T)) 
                           ); rm(contigs_all)

# read ANI tbl
ANI = fread("/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/cleaning_redundancy_vOTUs/ALL_pooled_viral_nocheckV.ANI.tsv") %>%
  # add sizes
  filter(pid > 95) %>% filter(qname != tname) %>%
   left_join(contigs_all_df %>% select(qname=contig, qsize = contig_size), by = "qname") %>%
   left_join(contigs_all_df %>% select(tname=contig, tsize = contig_size), by = "tname") %>% 
#   # remove hits that cover less than 85% from the shorter contig
   filter( (qcov > 85 & qsize <= tsize) | (tcov > 85 & tsize <= qsize) )


# checkV clusters
clusters = fread("/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/cleaning_redundancy_vOTUs/ALL_pooled_viral_nocheckV.clusters.tsv" 
                 , col.names = c("cluster","member"), header = F) %>% rowid_to_column("cluster_ID")
clusters_members = clusters %>% separate_rows(member, sep = ",") %>% group_by(cluster) %>% mutate(n_member = n()) %>% ungroup()

# vOTUs (own clustering)
agg_clusters= fread("/Volumes/gem-calc/Users/Eugen2/Metaviromes/Phageome_CdH/vOTUs/vOTUs_ID_contig.tsv")

####
differences = clusters %>% anti_join(agg_clusters, by = c("cluster"="contig")) %>% 
  separate_rows(member, sep = ",") %>% group_by(cluster) %>% mutate(n_member = n()) %>% ungroup()

singletons = differences %>% filter(n_member == 1)

# are singletons, singletons?
singletons_ANI = ANI %>% inner_join(singletons, by = c("qname"="member") ) %>% rename(q_cluster_ID = cluster_ID) %>%
  left_join(clusters_members %>% select(t_cluster_ID=cluster_ID, member), by = c("tname"="member"))

# examples
#ex 1
singletons_ANI_ex  = c("2644","857")

#all 
all_examples = singletons_ANI %>% select(cluster_ID = t_cluster_ID) %>% bind_rows(singletons_ANI %>% select(cluster_ID = q_cluster_ID)) %>% distinct()

# all singletons
singletons_all = clusters_members %>% filter(n_member == 1)

singletons_ANI_all = ANI %>% inner_join(singletons_all, by = c("qname"="member") ) %>% rename(q_cluster_ID = cluster_ID) %>%
  left_join(clusters_members %>% select(t_cluster_ID=cluster_ID, member), by = c("tname"="member"))

all_examples = singletons_ANI_all %>% select(cluster_ID = q_cluster_ID) %>% 
  bind_rows(singletons_ANI %>% select(cluster_ID = t_cluster_ID)) %>% distinct()

### example graph
nodes_exp = all_examples %>% left_join(clusters_members, by = c("cluster_ID")) %>% distinct() %>%
  select(member,cluster_ID, n_member) %>% left_join(contigs_all_df %>% select(member=contig,contig_size), by = c("member"))

edges_exp = ANI %>% semi_join(nodes_exp, by = c("qname"="member")) %>% 
  semi_join(nodes_exp, by = c("tname"="member"))

#make graph
      net_example = graph_from_data_frame(vertices = nodes_exp, d = edges_exp, directed = F)

      # safe graph
      write.graph(net_example, format = "graphml", file = "~/Desktop/cluster_ANI_95_cov85_very_all_example_checkV_single_new.graphml")

# ### old script
# 
# 
#   # step-wise sub cluster
#   cluster_types = clusters %>%  
#     #mutate(max_size = max(contig_size)) %>% ungroup() %>%
#     mutate(Cluster_Type = case_when(contig_size > 2e4 ~ ">20kb", contig_size > 1e4 ~ "10-20kb", contig_size > 5e3 ~"5-10kb", T ~ "<5kb"))%>%
#     group_by(one_link_cluster,Cluster_Type) %>% mutate(n_types = n()) %>% ungroup()
#     
#     # >20 kb clusters
#     nodes_clusters_20k = cluster_types %>% filter(Cluster_Type==">20kb")
#   
#     edges_cl20 = Edges_df %>% semi_join(nodes_clusters_20k, by = c("qname"="contig")) %>%
#       semi_join(nodes_clusters_20k, by = c("tname"="contig"))
#     
#     cl20_net = graph_from_data_frame(vertices = nodes_clusters_20k, d = edges_cl20, directed = F)
#     
#     clusters_20 = data.frame(one_link_cluster_stp2_20 = components(cl20_net)$membership) %>% rownames_to_column("contig") %>%
#       group_by(one_link_cluster_stp2_20) %>% mutate(cluster_size = n()) %>% ungroup() %>%
#       left_join(nodes_clusters_20k%>% select(contig, one_link_cluster), by = "contig")
#     
#     one_per_cluster20 = clusters_20 %>% left_join(viral_contigs_all_df, by = "contig")%>%
#       group_by(one_link_cluster_stp2_20) %>% arrange(desc(contig_size)) %>% slice(1) %>% ungroup() %>%
#       mutate(one_link_cluster_stp2 = str_c(one_link_cluster_stp2_20,"_20"))
# 
#     # cases that cluster in between >20 w/ 10-20kb
#     l20_to_lower_edges = Edges_df %>% inner_join(nodes_clusters_20k %>% select(contig,qCluster_Type=Cluster_Type), by = c("qname"="contig")) %>%
#       inner_join(cluster_types %>% filter(Cluster_Type != ">20kb") %>% select(contig,one_link_cluster, tCluster_Type=Cluster_Type, n_types), by = c("tname"="contig"))  
#     
#     # 10-20 kb clusters
#     nodes_clusters_10k = cluster_types %>% filter( Cluster_Type == "10-20kb") %>%
#       # and all contigs have hit large contigs
#       anti_join(l20_to_lower_edges, by = c("contig"="tname"))
#     
#     edges_cl10 = Edges_df %>% semi_join(nodes_clusters_10k, by = c("qname"="contig")) %>%
#       semi_join(nodes_clusters_10k, by = c("tname"="contig"))
#     
#     cl10_net = graph_from_data_frame(vertices = nodes_clusters_10k, d = edges_cl10, directed = F)
#     
#     clusters_10 = data.frame(one_link_cluster_stp2_10 = components(cl10_net)$membership) %>% rownames_to_column("contig") %>%
#       group_by(one_link_cluster_stp2_10) %>% mutate(cluster_size = n()) %>% ungroup() %>%
#       left_join(nodes_clusters_10k%>% select(contig, one_link_cluster), by = "contig")
#     
#     one_per_cluster10 = clusters_10 %>% left_join(viral_contigs_all_df, by = "contig")%>%
#       group_by(one_link_cluster_stp2_10) %>% arrange(desc(contig_size)) %>% slice(1) %>% ungroup()%>%
#       mutate(one_link_cluster_stp2 = str_c(one_link_cluster_stp2_10,"_10"))
#     
#     # cases in between 10 to 5-10 kb
#     l10_to_lower_edges = Edges_df %>% inner_join(nodes_clusters_10k %>% select(contig,qCluster_Type=Cluster_Type), by = c("qname"="contig")) %>%
#       inner_join(cluster_types %>% filter( Cluster_Type != "10-20kb") %>% 
#                    select(contig,one_link_cluster, tCluster_Type=Cluster_Type, n_types), by = c("tname"="contig"))  
#     
#     # 5-10 kb clusters
#     nodes_clusters_5k = cluster_types %>% filter( Cluster_Type == "5-10kb") %>%
#       anti_join(l20_to_lower_edges, by = c("contig"="tname")) %>%
#       anti_join(l10_to_lower_edges, by = c("contig"="tname")) 
#     
#     edges_cl5 = Edges_df %>% semi_join(nodes_clusters_5k, by = c("qname"="contig")) %>%
#       semi_join(nodes_clusters_5k, by = c("tname"="contig"))
#     
#     cl5_net = graph_from_data_frame(vertices = nodes_clusters_5k, d = edges_cl5, directed = F)
#     
#     clusters_5 = data.frame(one_link_cluster_stp2_5 = components(cl5_net)$membership) %>% rownames_to_column("contig") %>%
#       group_by(one_link_cluster_stp2_5) %>% mutate(cluster_size = n()) %>% ungroup() %>%
#       left_join(nodes_clusters_5k%>% select(contig, one_link_cluster,cluster_size_old = cluster_size), by = "contig")
#     
#     one_per_cluster5 = clusters_5 %>% left_join(viral_contigs_all_df, by = "contig")%>%
#       group_by(one_link_cluster_stp2_5) %>% arrange(desc(contig_size)) %>% slice(1) %>% ungroup()%>%
#       mutate(one_link_cluster_stp2 = str_c(one_link_cluster_stp2_5,"_05"))
#     
#     l5_to_lower_edges = Edges_df %>% inner_join(nodes_clusters_5k %>% select(contig,qCluster_Type=Cluster_Type), by = c("qname"="contig")) %>%
#       inner_join(cluster_types %>% filter( Cluster_Type != "5-10kb") %>% select(contig,one_link_cluster, tCluster_Type=Cluster_Type, n_types), by = c("tname"="contig"))  
#     
#     #< 5kb
#     nodes_clusters_3k = cluster_types %>% filter(Cluster_Type == "<5kb") %>% 
#        anti_join(l20_to_lower_edges, by = c("contig"="tname")) %>%
#        anti_join(l10_to_lower_edges, by = c("contig"="tname")) %>%
#        anti_join(l5_to_lower_edges, by = c("contig"="tname"))
#     
#     edges_cl3 = Edges_df %>% semi_join(nodes_clusters_3k, by = c("qname"="contig")) %>%
#       semi_join(nodes_clusters_3k, by = c("tname"="contig"))
#     
#     cl3_net = graph_from_data_frame(vertices = nodes_clusters_3k, d = edges_cl3, directed = F)
#     
#     clusters_3 = data.frame(one_link_cluster_stp2_3 = components(cl3_net)$membership) %>% rownames_to_column("contig") %>%
#       group_by(one_link_cluster_stp2_3) %>% mutate(cluster_size = n()) %>% ungroup() %>%
#       left_join(nodes_clusters_3k %>% select(contig, one_link_cluster,cluster_size_old = cluster_size), by = "contig")
#     
#     one_per_cluster3 = clusters_3 %>% left_join(viral_contigs_all_df, by = "contig")%>%
#       group_by(one_link_cluster_stp2_3) %>% arrange(desc(contig_size)) %>% slice(1) %>% ungroup()%>%
#       mutate(one_link_cluster_stp2 = str_c(one_link_cluster_stp2_3,"_03"))
#     
# ###
# one_virus_per_cluster = bind_rows(one_per_cluster20,one_per_cluster10, one_per_cluster5, one_per_cluster3) %>% 
#   select(contig, contig_size, one_link_cluster_stp2)  
# 
# # overview_clusters
# clustering_stp2 = bind_rows(  
#   clusters_20 %>% mutate(one_link_cluster_stp2_type = str_c("20_",one_link_cluster_stp2_20))
#   ,clusters_10 %>% mutate(one_link_cluster_stp2_type = str_c("10_",one_link_cluster_stp2_10))
#   ,clusters_5 %>% mutate(one_link_cluster_stp2_type = str_c("05_",one_link_cluster_stp2_5))
#   ,clusters_3 %>% mutate(one_link_cluster_stp2_type = str_c("03_",one_link_cluster_stp2_3)) ) %>%
#   select(contig, one_link_cluster, one_link_cluster_stp2_type)
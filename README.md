Files and scripts to reproduce figures and analysis from 
# Antibiotic perturbation of the human gut phageome preserves its individuality and promotes blooms of virulent phages
https://doi.org/10.1101/2025.02.07.636470

![Présentation1](https://github.com/user-attachments/assets/076b8d6f-e53b-487e-8569-b7fa3988c21b)


 <br><br>

**Supplementary Tables folder contains all required tables.**
 <br>
 
**`Table S1`: Overview of viral species (=vOTUs)**
* vOTU_ID - Specifier/ID for the vOTUs
* Novel species - Yes/No, based on ANI 95% and AF 85% clustering with viral species of UHGV
* VC - Viral cluster, as given by vConTACT v2. vOTUs that were assigned to seveal clusters (overlap) have several VC indicated (separated by ;). Outlier and singletons were put together to OUT/SINGLE 
* VClustering - Same as VC, but overlaps were added together to OUT/SINGLE/OVER  
* size - contig and genomes sizes in bp
* Novel_genus - Yes/No. A novel genus is defined to contain only novel species.
* dominant - Yes/No. Yes is attributed, if found in at least one donor with an abundance level >25%.
* source - One of three approaches used to obtain the contig: SpadesMetaV (assembly using Spades with --metaviral), SpadesFractionR )random subsampling of lowly covered readsets and assembly with Spades (default)), UHGV (contigs of UHGV that were covered by the reads at least to 50% and a sequencing depth of 1)
* checkV quality - Quality indexes of checkV (Complete, High-, Medium-, Low-quality, and not-determined)
* viral-lifecycle - Prediction of the lifecycle by BACPHLIP for genomes with at least High-quality. A temperate state was assigned, when prediction (ranges between 0 and 1) was >0.5 towards temperate. Otherwise it was set as virulent.
* Related_to_PP - Yes or No,  relation to a phage-plasmids was determined by wGRR or vConTACT v2
* With_ARG - No: no antibiotic resistance gene (ARG) was detected. Yes: ARGs with vOTUs (predicted with AMRfinderPlus), and the encoding gene
* viral_taxonomy_lineage - Taxonomy was predicted by geNomad
* host_taxonomy - Host and their taxonomy were predicted by iPHoP. For vOTUs recovered from UHGV, the already predicted host was taken.
* geNomad, checkV (if classed medium, high or complete), VIBRANT, viralVerify: 1: predicted or analyzed by one of these tools, 0: not detected or determined
   
**`Table S2`: Overview of read and sample sets** 
* Sequence_set - Category of sequences (viral, non-viral, proviral, UHGV, UHGG, unampped). All reads sets were mapped on theses types and the fraction was computed
* mapped_read_count - Absolute number of mapped reads
* mill_reads_per_sample - Number of total reads in the sample (in Millions)
* fraction in % = Fraction of reads mapped to the type of sequences
* donor - Donor ID of the read set. Donor 29 was excluded from the analysis.
* day - Day or sample ID (-14,-7,-1,4,7,10,15,30,90,180) of the read set.

**`Table S3`: Overview of related P-Ps that were retrieved from PMID:38378896**
* Name of P-P as stated in NCBI - Name entry of NCBI RefSeq
* Accession (RefSeq) - Accession ID of RefSeq
* Number of related vOTUs - Sum of the number of vOTUs detected by the two methods (wGRR and vConTACT v2)
* Comment to P-P state - Some phage-plasmids experimentally proven to follow this life cycle, others are suggested and for some any evidence is lacking. Here an additional comment with reference to litature (if available) is left.
* Method_detected (wGRR-only, vContact v2 only, or by the two)	] - Specification of the method that was used to detect the relation to the phage-plasmid.
* vOTUs - List of vOTUs related to the phage-plasmids. Separated by semi-colon.

**`Table S4`: vOTUs across donors** 
* vOTU_ID - as in Table S1
* Donor01 - Donor34: mean and sd of relative abundances in %, nD: not determinable

**`Table S5`: Abundance table for all vOTUs**
*  vOTU_ID and Size - As in Table S1
*  Donor  - ID of the volunteer who gave the sample
*  Sampling_Day - Time point of the sample
*  ab_abundance - Absolute counts of reads mapped on the vOTU 
  
**`Table S6`: Overview of microbial species detected in PMID: 38468305** 
* id_mgs - species ID
* in_n_donor - Number in how many donors this species was detected
* donorID_detected - ID of donors in which the species was detected
* gtdb_classification - taxonomy of the species according to the Genome Taxonomy Database
  
**`Table S7`: Relative abundance of microbial species**
* id_mgs - as in Table S5
* relative_abundance - computed with METEOR. Relative abundance was computed on the mean of the most frequent signature genes
* Donor - Donor_ID
* Sampling_day - Time point (day) of the sample

 <br><br>

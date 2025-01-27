Files and scripts to reproduce figures and analysis from 
# Third-generation cephalosporin antibiotics induce phage bursts in the human gut microbiome
 ([Manuscript Link])

![image](https://github.com/user-attachments/assets/aaa421e5-425c-4da1-9556-1fa6fb3548fe)

**Data and Required Files:**

**In the folder Supplementary Tables all required tables were placed.**
 
**Table S1: Overview of the viral species (=vOTUs)**
* vOTU_ID - Specifier/ID for the vOTUs,
* size - contig and genomes sizes in bp
* source - One of three approaches used to obtain the contig: SpadesMetaV (assembly using Spades with --metaviral), SpadesFractionR )random subsampling of lowly covered readsets and assembly with Spades (default)), UHGV (contigs of UHGV that were covered by the reads at least to 50% and a sequencing depth of 1)
* checkV quality - Quality indexes of checkV (Complete, High-, Medium-, Low-quality, and not-determined)
* viral-lifecycle - Prediction of the lifecycle by BACPHLIP for genomes with at least High-quality. A temperate state was assigned, when prediction (ranges between 0 and 1) was >0.5 towards temperate. Otherwise it was set as virulent.
* found_in_n_donors - Number of donors, in whom the vOTUs was detected
* donorID_sample - ID of the donor in who it was detect and respective sample point (in backets)
* Related_to_PP - Yes or No,  relation to a phage-plasmids was determined by wGRR or vConTACT v2
* With_ARG - No: no antibiotic resistance gene (ARG) was detected. Yes: ARGs with vOTUs (predicted with AMRfinderPlus), and the encoding gene
* viral_taxonomy_lineage - Taxonomy was predicted by geNomad
* host_taxonomy - Host and their taxonomy were predicted by iPHoP. For vOTUs recovered from UHGV, the already predicted host was taken.
* geNomad, checkV (at least medium quality), VIBRANT, viralVerify: 1: predicted or analyzed by one of these tools, 0: not detected or determined
   
**Table S2: Overview of read and sample sets** 
* Sequence_set - Category of sequences (viral, non-viral, proviral, UHGV, UHGG, unampped). All reads sets were mapped on theses types and the fraction was computed
* mapped_read_count - Absolute number of mapped reads
* mill_reads_per_sample - Number of total reads in the sample (in Millions)
* fraction in % = Fraction of reads mapped to the type of sequences
* donor - Donor ID of the read set. Donor 29 was excluded from the analysis.
* day - Day or sample ID (-14,-7,-1,4,7,10,15,30,90,180) of the read set.

**Table S3: Overview of related P-Ps that were retrieved from PMID:38378896**
* Name of P-P as stated in NCBI - Name entry of NCBI RefSeq
* Accession (RefSeq) - Accession ID of RefSeq
* Number of related vOTUs - Sum of the number of vOTUs detected by the two methods (wGRR and vConTACT v2)
* Comment to P-P state - Some phage-plasmids experimentally proven to follow this life cycle, others are suggested and for some any evidence is lacking. Here an additional comment with reference to litature (if available) is left.
* Method_detected (wGRR-only, vContact v2 only, or by the two)	] - Specification of the method that was used to detect the relation to the phage-plasmid.
* vOTUs - List of vOTUs related to the phage-plasmids. Separated by semi-colon.

* **Table S4: Abundance table for all vOTUs**
*  vOTU_ID and Size - As in Table S1
*  Donor  - ID of the volunteer who gave the sample
*  Sampling_Day - Time point of the sample
*  rary_abundance - Counts of reads were down-sized to the smallest read set per donor, and used only to compute the Shannon index
*  
*  
* **Table S5:** [Brief description of the supporting file]
* 
* 
* **Table S6:** [Brief description of the supporting file]
* 
* 

**RScripts to reproduce the figures and analysis were placed in the folder Rscripts**

* **[overview_vOTUs]:** Generates Figure [Figure Number] 
* **[diversity_analysis_vOTUs]:** Generates Figure [Figure Number] 
* **[votus_wARGs]:** Generates Figure [Figure Number] 
* **[abundance_over_time]:** Generates Figure [Figure Number]
* **[non_redundant_contig_ANI]:** Used to dereplicate redundant sequences (to define vOTUs) using a network approach, and a species definition of 95% ANI over 85% AF. Details are listed in the Methods part under * Clustering viral contigs into viral Operational Taxonomical Units (vOTUs). The multifasta sequence file, blastn and anicalc.py (from CheckV) are required.
* **[multimer_resolution_mummer]:** This script is used to detect and solve concatemeric sequences. It requires input from blastn, and a MUMmer 4.0 comparison to detect tandem repeats.

# All R packages that are required for the analysis are listed in the respective script files.

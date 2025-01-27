Files and scripts to reproduce figures and analysis from 
# Third-generation cephalosporin antibiotics induce phage bursts in the human gut microbiome
 ([Manuscript Link])

![image](https://github.com/user-attachments/assets/aaa421e5-425c-4da1-9556-1fa6fb3548fe)

**Data and Required Files:**

* **Table S1:** Overview of the viral species (=vOTUs).
* Table contains:
  vOTU_ID - Specifier/ID for the vOTUs,
  size - contig and genomes sizes in bp
  source - One of three approaches used to obtain the contig: SpadesMetaV (assembly using Spades with --metaviral), SpadesFractionR )random subsampling of lowly covered readsets and assembly with Spades (default)), UHGV (contigs of UHGV that were covered by the reads at least to 50% and a sequencing depth of 1)
  checkV quality - Quality indexes of checkV (Complete, High-, Medium-, Low-quality, and not-determined)
  viral-lifecycle - Prediction of the lifecycle by BACPHLIP for genomes with at least High-quality. A temperate state was assigned, when prediction (ranges between 0 and 1) was >0.5 towards temperate. Otherwise it was set as virulent.
 found_in_n_donors - Number of donors, in whom the vOTUs was detected
donorID_sample - ID of the donor in who it was detect and respective sample point (in backets)
Related_to_PP - Yes or No,  relation to a phage-plasmids was determined by wGRR or vConTACT v2
With_ARG - No: no antibiotic resistance gene (ARG) was detected. Yes: ARGs with vOTUs (predicted with AMRfinderPlus), and the encoding gene
viral_taxonomy_lineage - Taxonomy was predicted by geNomad
host_taxonomy - Host and their taxonomy were predicted by iPHoP. For vOTUs recovered from UHGV, the already predicted host was taken.
geNomad, checkV (at least medium quality), VIBRANT, viralVerify: 1: predicted or analyzed by one of these tools, 0: not detected or determined

* **Table S2:** [Brief description of the supporting file]
* **Table S3:** [Brief description of the supporting file]
* **Table S4:** [Brief description of the data file and its contents]
* **Table S5:** [Brief description of the supporting file]
* **Table S6:** [Brief description of the supporting file]

**R Scripts:**

* **[overview_vOTUs]:** Generates Figure [Figure Number] 
* **[diversity_analysis_vOTUs]:** Generates Figure [Figure Number] 
* **[votus_wARGs]:** Generates Figure [Figure Number] 
* **[abundance_over_time]:** Generates Figure [Figure Number]
* **[non_redundant_contig_ANI]:** Used to dereplicate redundant sequences (to define vOTUs) using a network approach, and a species definition of 95% ANI over 85% AF. Details are listed in the Methods part under * Clustering viral contigs into viral Operational Taxonomical Units (vOTUs). The multifasta sequence file, blastn and anicalc.py (from CheckV) are required.
* **[multimer_resolution_mummer]:** This script is used to detect and solve concatemeric sequences. It requires input from blastn, and a MUMmer 4.0 comparison to detect tandem repeats.

# All R packages that are required for the analysis are listed in the respective script files.

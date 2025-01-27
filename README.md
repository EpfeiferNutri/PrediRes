Files and scripts to reproduce figures and analysis from 
# Third-generation cephalosporin antibiotics induce phage bursts in the human gut microbiome
 ([Manuscript Link])

![image](https://github.com/user-attachments/assets/aaa421e5-425c-4da1-9556-1fa6fb3548fe)


This repository contains all R scripts necessary to reproduce the figures presented in the manuscript "[Manuscript Title]" ([Manuscript Link]). 

**Data and Required Files:**

* **Table S1:** [Link to data file or instructions to obtain it] 
    * Description: [Brief description of the data file and its contents]
* **Table S2::** [Link to file or instructions to obtain it]
    * Description: [Brief description of the supporting file]
* **Table S3:** [Link to file or instructions to obtain it]
    * Description: [Brief description of the supporting file]
* **Table S4:** [Link to data file or instructions to obtain it] 
    * Description: [Brief description of the data file and its contents]
* **Table S5::** [Link to file or instructions to obtain it]
    * Description: [Brief description of the supporting file]
* **Table S6:** [Link to file or instructions to obtain it]
    * Description: [Brief description of the supporting file]

**R Scripts:**

* **[overview_vOTUs]:** Generates Figure [Figure Number] 
* **[diversity_analysis_vOTUs]:** Generates Figure [Figure Number] 
* **[votus_wARGs]:** Generates Figure [Figure Number] 
* **[abundance_over_time]:** Generates Figure [Figure Number]
* **[non_redundant_contig_ANI]:** Used to dereplicate redundant sequences (to define vOTUs) using a network approach, and a species definition of 95% ANI over 85% AF. Details are listed in the Methods part under * Clustering viral contigs into viral Operational Taxonomical Units (vOTUs). The multifasta sequence file, blastn and anicalc.py (from CheckV) are required.
* **[multimer_resolution_mummer]:** This script is used to detect and solve concatemeric sequences. It requires input from blastn, and a MUMmer 4.0 comparison to detect tandem repeats.

# All R packages that are required for the analysis are listed in the respective script files.

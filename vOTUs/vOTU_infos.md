**Nucleotide and protein fasta sequences of the vOTUs are placed in a zenodo repository: `10.5281/zenodo.14758483`**

* Link to the repository:
  
https://zenodo.org/records/14758483?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjYxMmUyNDNmLWYyYmMtNGFjNC05ODZiLWFiMmQ2ZTc1YmM3MSIsImRhdGEiOnt9LCJyYW5kb20iOiJhNGIzYzAxMjQ2MDg4MDkxMDM2N2E3YjI2MzJkZjRmYSJ9.5QuWvikEPoIcYzcVqZKbwihgxErV_OVtfviCRtzP0PIi4Gf-q6ybLoT4cCTP_UsOkUUUNj8T19XNSuTiKtJfpg

<br>

* `amrfinder_out`: AMRfinderPlus output on geNomad-predicted protein sequences.
* `gggenomes_votus_amrfinder.R`: R script for plotting vOTU with ARGs  (Figure S6).
<br>
 
**In folder `homology_to_P-Ps`:**
  * `g2g_pp_votus.tsv` - output table of vConTACT v2 (P-Ps with vOTUs)
  * `wGRR_vOTUs_pp.zip` - compressed wGRR network table (P-Ps to vOTUs)
  * `supplemental_file_1_pps_votus_overview` - summary of the the wGRR and vConTACT v2 analysis

<br>

Folder `comparisont_to_UHGV`  
  * meta information (`metadata_UHGV.zip`, zipped) on the UHGV database (last accessed Jan'24). This information was used to generate Figures S2AB with variations of ggplot2 functions listed in `Rscripts/overview_vOTUs.R`.
  * `vOTUs_from_UHGV.tsv` table indicating which vOTUs originated from UHGV. UHGV sequences had to be (on average) covered at least 50%, and with a depth of 1 (see Methods).  

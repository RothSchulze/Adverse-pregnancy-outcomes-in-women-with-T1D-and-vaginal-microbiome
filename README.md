# Adverse-pregnancy-outcomes-in-women-with-T1D-and-vaginal-microbiome

This repository contains a file describing the specific commands used to process the bacterial and fungal raw data from manuscript "Adverse pregnancy outcomes in women with type 1 diabetes are associated with multiple alterations in the vaginal microbiome" using Qiime2 version 2024.10 and the mapping files used for this. The raw sequences can be found in NCBI SRA https://www.ncbi.nlm.nih.gov/sra (project number PRJNA1037097). It also contains files with the R code used for analysing the samples along with R objects and metadata used as input. 

## QIIME2 files

-> 2025 QIIME2 version 2024.10 Vaginal project: Document containing the description and commands of the Qiime2 pipeline used for processing the raw sequences until obtaining taxonomic vs samples table. Also with this representative sequences and phylogenetic trees were obtain. In this document it is also explained how the qiime2 files were transformed into biom files for input to phyloseq in R.

-> The following metadata files were used in that pipeline:

  -> Mapping_file_Seq_8-Run.csv
  -> Mapping_file_Run_17.csv
  -> Mapping_file_Run_8_17.csv
  -> Mapping_file_Vaginal_ITS1_sequencing.csv
  -> Mapping_file_Vaginal_ITS_sequencing_S.csv

## R code files and objects

The R codes are contained in 4 different sequntial files. 

**1) Phyloseq_Mother_data_formating_Vaginal_Microbiome_runs_8_17_2025.R**:

Contains the code used to format the results from Qimme2 for the vaginal microbiome and filtering of low abundance sequences in which the input files are:

-> Vaginal_merged-tables8_17G_wTaxMet.biom: This contains the bacterial features with taxonomic classification

-> Mapping_file_Run_8_17.csv: This contains the initial sequence metadata

-> rooted-tree_QVaginal8_17.nwk: This contains the bacterial phylogenetic tree generated in qimme2.

-> rep-seqs-dada2_QVaginal8_17.fasta: Fasta file with representative sequences

-> Metadata_8_17_Merged2PET1D.csv: Contains the metadata used to merged semples that were sequenced in duplicate.

Resulting output file: Bacterial_Vaginal_Features_Phyloseq_Obj_Filtered1.RData
Also: rep-seqs-Bacterial_VCMCFHS_AFilt1.fasta

**1) Phyloseq_Mother_data_formating_Vaginal_Mycobiome_2025.R**:

Contains the code used to format the results from Qimme2 for the vaginal mycobiome and filtering of low abundance sequences in which the input files are:

-> table_ITS1_R1_wTaxMet.biom: his contains the bacterial features with taxonomic classification and initial metadata

-> rooted_tree_ITS1_R1.nwk: This contains the fungal phylogenetic tree generated in qimme2.

-> rep-seqs_ITS1_R1.fasta: Fasta file with representative sequences

Resulting output file: Myco_ITS1_R1_CU_Filter1.RData
Also: rep-seqs-Myco_ITS1_R1_CU_Filter1.fasta

**2) Vaginal_micro-mycobiome_data_formatting.Rmd**:

Contains the code used to format Bacterial_VCMCFHS_AFilt1 in the R object Bacterial_Vaginal_Features_Phyloseq_Obj_Filtered1.RData in order to improve the phylogenetic classification of bacterial taxa using BLASTn and also to format mycobiome data. Files used for this:

-> Bacterial_Vaginal_Features_Phyloseq_Obj_Filtered1.RData: Generated above (filtered bacteria taxa)

-> BLASTn_classification_vaginal_microbiome.csv: Contains improved classification obtained with BLASTn.

-> Myco_ITS1_R1_CU_Filter1.RData: Generated above (filtered fungal taxa)

-> Metadata_8_17_Merged2PET1D.csv: Metadata file used in this step to obtain sample names.

Resulting output files:

-> MicrobiomeObjects.RData: Filtered bacterial taxa aglomerated to Species, Genus, Family, Order and Phylum taxonomic levels.

-> MycobiomeObjects.RData: Filtered fungal taxa aglomerated to Species Hypothesis, Genus, Family, Order and Phylum taxonomic levels.

**3) T1D analysis vaginal microbiome.Rmd**

Contains the code used to analyse the alpha and beta diversity and differential abundance of the vaginal bacterial and fungal microbiomes in the context of T1D status. Also the effect of HLA status was evaluated here. Input data:

-> MicrobiomeObjects.RData: see above

-> MycobiomeObjects.RData: see above

-> Metadata_8_17_Merged2PET1D.csv: This is the metadata per sample used for the analysis.

**3) Adverse pregnancy outcomes analysis vaginal microbiome.Rmd**

Contains the code used to analyse the alpha and beta diversity and differential abundance of the vaginal bacterial and fungal microbiomes in the context of the pregnancy complications preeclampsia and pre-term delivery.

-> MicrobiomeObjects.RData: see above

-> MycobiomeObjects.RData: see above

-> Metadata_8_17_Merged2PET1D.csv: This is the metadata per sample used for the analysis.


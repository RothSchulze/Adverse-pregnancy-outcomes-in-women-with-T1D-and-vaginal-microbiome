# Adverse-pregnancy-outcomes-in-women-with-T1D-and-vaginal-microbiome

This repository contains files with the code used for analysing the samples from manuscript "Adverse pregnancy outcomes in women with type 1 diabetes are associated with multiple alterations in the vaginal microbiome" in R. In also includes R objects and metadata used as input. 

The R codes are contained in 4 different sequntial files. 

**1) Phyloseq_Mother_data_formating_Vaginal_Microbiome_runs_8_17.R**:

Contains the code used to format the results from Qimme3 for the vaginal microbiome and filtering of low abundance sequences in which the input files are:

-> table200-160_Run8_17_wTaxMet.biom: This contains the bacterial features with taxonomic classification

-> Mapping_file_Run_8_17.csv: This contains the initial sequence metadata

-> rooted-tree_Vaginal.nwk: This contains the bacterial phylogenetic tree generated in qimme3.

-> Metadata_8_17_Merged.csv: Contains the metadata used to merged semples that were sequenced in duplicated.

Resulting output file: Mothers_Vaginal_OTU_Phyloseq_Obj_Filtered1.RData

**2) Vaginal_micro-mycobiome_data_formatting.Rmd**:

Contains the code used to format Mothers_Vaginal_OTU_Phyloseq_Obj_Filtered1.RData in order to improve the phylogenetic classification of bacterial taxa using BLASTn and also to format mycobiome data. Files used for this:

-> Mothers_Vaginal_OTU_Phyloseq_Obj_Filtered1.RData: Generated above (filtered bacteria taxa)

-> BLASTn_classification_vaginal_microbiome2.csv: Contains improved classification obtained with BLASTn.

-> Metadata_8_17_Merged.csv: Contains the metadata used to merged semples that were sequenced in duplicated.

-> rooted-tree_Vaginal.nwk: This contains the phylogenetic tree generated in qimme3.

-> table_ITS1_R1_wTaxMet.biom: This contains the fungal with taxonomic classification

-> rooted-tree_ITS1_R1.nwk: This contains the fungal phylogenetic tree generated in qimme3.

Resulting output files:

-> MicrobiomeObjectsWHLAsamples.RData: Filtered bacterial taxa aglomerated to OTU, Genus, Family, Order and Phylum taxonomic levels.

-> MycobiomeObjects.RData: Filtered fungal taxa aglomerated to OTU, Genus, Family, Order and Phylum taxonomic levels.

**3) T1D analysis vaginal microbiome.Rmd**

Contains the code used to analyse the alpha and beta diversity and differential abundance of the vaginal bacterial and fungal microbiomes in the context of T1D status. Also the effect of HLA status was evaluated here. Input data:

-> MicrobiomeObjectsWHLAsamples.RData: see above

-> MycobiomeObjects.RData: see above

-> Metadata_8_17_Merged2PET1D.csv: This is the metadata per sample used for the analysis.

**3)T1D-analysis-vaginal-microbiome.html**

Contains the knitr file from 3) T1D analysis vaginal microbiome.Rmd

**3) Adverse pregnancy outcomes analysis vaginal microbiome.Rmd**

Contains the code used to analyse the alpha and beta diversity and differential abundance of the vaginal bacterial and fungal microbiomes in the context of the pregnancy complications preeclampsia and pre-term delivery.

-> MicrobiomeObjectsWHLAsamples.RData: see above

-> MycobiomeObjects.RData: see above

-> Metadata_8_17_Merged2PET1D.csv: This is the metadata per sample used for the analysis.

**3)Adverse-pregnancy-outcomes-analysis-vaginal-microbiome.html**

Contains the knitr file from 3) Adverse pregnancy outcomes analysis vaginal microbiome.Rmd

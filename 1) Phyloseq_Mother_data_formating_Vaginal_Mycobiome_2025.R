library("phyloseq")
library("ggplot2")
library("devtools")
library("metagMisc")
library("metagenomeSeq")
library("Biostrings")

setwd("~/Documents/ENDIA Project/Data and analyses/ENDIA Projects analysis/2025 Sequencing processing with qiime2 2024.10/Vaginal project/Diabetologia submission 2025/Code and R objects/") # HERE YOU SET YOUR DIRECTORY WHERE DATA IS CONTAINED

####################################################################
##                                                               ###
##                  IMPORT AND SUBSET DATA                       ###
##               WITH TAXA AND METADATA INFO                     ###
####################################################################

Myco_ITS1_R1 <- import_biom("./table_ITS1_R1_wTaxMet.biom")
colnames(tax_table(Myco_ITS1_R1)) <- c("Kindom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Sp.Hypotheses ")
### Import tree ITS1_R1
Phylogenetic_TreeITS1_R1 <- read_tree("./rooted_tree_ITS1_R1.nwk")
#Import fasta file
RefFasta <- readDNAStringSet("./rep-seqs_ITS1_R1.fasta")
### Making a phyloseq object 
Myco_ITS1_R1_C <- phyloseq(otu_table(Myco_ITS1_R1),  sample_data(Myco_ITS1_R1), tax_table(Myco_ITS1_R1), phy_tree(Phylogenetic_TreeITS1_R1), refseq(RefFasta))

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1427 taxa and 153 samples ]
#sample_data() Sample Data:       [ 153 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 1427 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1427 tips and 1416 internal nodes ]
#refseq()      DNAStringSet:      [ 1427 reference sequences ]

##Only keep fungal sequences
Myco_ITS1_R1_CU <- subset_taxa(Myco_ITS1_R1_C, Kindom == "k__Fungi")

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1427 taxa and 153 samples ]
#sample_data() Sample Data:       [ 153 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 1427 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1427 tips and 1416 internal nodes ]
#refseq()      DNAStringSet:      [ 1427 reference sequences ]

####################################################################
##                                                               ###
##                  ABUNDANCE FILTER Features                    ###
##                                                               ###
####################################################################

Myco_ITS1_R1_CU_dataframe <- t(as.data.frame(otu_table(Myco_ITS1_R1_CU)))

# function to perform pre-filtering
low.count.removal = function(data, percent=0.01)
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

result.filter = low.count.removal(Myco_ITS1_R1_CU_dataframe, percent=0.01)
length(result.filter$keep.otu) 
Keep_Feat <- names(result.filter$keep.otu)

Myco_ITS1_R1_CU_Filter1 <- prune_taxa(Keep_Feat, Myco_ITS1_R1_CU)
Myco_ITS1_R1_CU_Filter1 <-prune_taxa(taxa_sums(Myco_ITS1_R1_CU_Filter1)>=1, Myco_ITS1_R1_CU_Filter1)

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 593 taxa and 153 samples ]
#sample_data() Sample Data:       [ 153 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 593 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 593 tips and 591 internal nodes ]
#refseq()      DNAStringSet:      [ 593 reference sequences ]

save(Myco_ITS1_R1_CU_Filter1, file="./Myco_ITS1_R1_CU_Filter1.RData")
sequences <- refseq(Myco_ITS1_R1_CU_Filter1)
writeXStringSet(sequences, filepath = "./rep-seqs-Myco_ITS1_R1_CU_Filter1.fasta")

  library("phyloseq")
  library("ggplot2")
  library("devtools")
  library("metagMisc")
  library("metagenomeSeq")
  library("Biostrings")
  
  ####################################################################
  ##                                                               ###
  ##                  IMPORT AND SUBSET DATA                       ###
  ##               WITH TAXA AND METADATA INFO                     ###
  ####################################################################
  
  ## Vaginal samples were sequenced in 2 different runs. All the sequences in the different runs were processed separatelly using the same parameters with QIIME2.
  
  setwd("~/Documents/ENDIA Project/Data and analyses/ENDIA Projects analysis/2025 Sequencing processing with qiime2 2024.10/Vaginal project/Diabetologia submission 2025/Code and R objects/") # HERE YOU SET YOUR DIRECTORY WHERE DATA IS CONTAINED
  
  ### Import data Non Filtered 
  Bacterial_V_wT <- import_biom("./Vaginal_merged-tables8_17G_wTaxMet.biom", parseFunction = parse_taxonomy_greengenes) #
  tax_table(Bacterial_V_wT) <- tax_table(Bacterial_V_wT)[, c(1:6,8)]
  colnames(tax_table(Bacterial_V_wT))[1] <- "Kindom"
  
  NewData <- read.table(file="./Mapping_file_Run_8_17.csv", header=TRUE, row.names=1, sep="\t")
  #Import fasta file
  RefFasta <- readDNAStringSet("./rep-seqs-dada2_QVaginal8_17.fasta")
  
  ### Import tree
  Phylogenetic_Tree <- read_tree("./rooted-tree_QVaginal8_17.nwk")
  
  ### Making a phyloseq object 
  Bacterial_VC <- phyloseq(otu_table(Bacterial_V_wT),  sample_data(NewData), tax_table(Bacterial_V_wT), phy_tree(Phylogenetic_Tree), refseq(RefFasta))
  
  #phyloseq-class experiment-level object
  #otu_table()   OTU Table:         [ 4112 taxa and 800 samples ]
  #sample_data() Sample Data:       [ 800 samples by 5 sample variables ]
  #tax_table()   Taxonomy Table:    [ 4112 taxa by 7 taxonomic ranks ]
  #phy_tree()    Phylogenetic Tree: [ 4112 tips and 3992 internal nodes ]
  #refseq()      DNAStringSet:      [ 4112 reference sequences ]
  
  ## Merge the wells based on the SubjectID and trimester (well replicates) ##
  Bacterial_VCM <-merge_samples(Bacterial_VC, "Merging", fun = sum)
  ### Recover sample data after merging
  New_sample_data <- read.table(file="./Metadata_8_17_Merged2PET1D.csv", header=TRUE, sep="\t", row.names=1)
  ### Making a phyloseq object 
  Bacterial_VCMC <- phyloseq(otu_table(Bacterial_VCM),  sample_data(New_sample_data), tax_table(Bacterial_VCM), phy_tree(Bacterial_VCM), refseq(RefFasta))
  
  ##FILTER SEQUENCES CLASSIFIED AS UNASSIGNED
  Bacterial_VCMCFHS <- subset_taxa(Bacterial_VCMC, Kingdom != "Unassigned")
  #Remove QLD-BRI-MATER-084
  Bacterial_VCMCFHS <- subset_samples(Bacterial_VCMCFHS, rownames(sample_data(Bacterial_VCMCFHS)) != "QLD-BRI-MATER-084") ## Because it has MODY
  Bacterial_VCMCFHS <- subset_samples(Bacterial_VCMCFHS, rownames(sample_data(Bacterial_VCMCFHS)) != "NSW-SYD-CHW-011") ## It has contamination
  
  #phyloseq-class experiment-level object
  #otu_table()   OTU Table:         [ 2442 taxa and 333 samples ]
  #sample_data() Sample Data:       [ 333 samples by 24 sample variables ]
  #tax_table()   Taxonomy Table:    [ 2442 taxa by 7 taxonomic ranks ]
  #phy_tree()    Phylogenetic Tree: [ 2442 tips and 2324 internal nodes ]
  #refseq()      DNAStringSet:      [ 2442 reference sequences ]
  
  ####################################################################
  ##                                                               ###
  ##                  ABUNDANCE FILTER Features                    ###
  ##                                                               ###
  ####################################################################
  
  Bacterial_VCMCFHS_dataframe <- as.data.frame(otu_table(Bacterial_VCMCFHS))
  
  # function to perform pre-filtering
  low.count.removal = function(data, percent=0.01)
  {
    keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
    data.filter = data[,keep.otu]
    return(list(data.filter = data.filter, keep.otu = keep.otu))
  }
  
  result.filter = low.count.removal(Bacterial_VCMCFHS_dataframe, percent=0.01)
  length(result.filter$keep.otu) 
  Keep_OTUS <- names(result.filter$keep.otu)
  
  Bacterial_VCMCFHS_AFilt1 <- prune_taxa(Keep_OTUS, Bacterial_VCMCFHS)
  Bacterial_VCMCFHS_AFilt1 <-prune_taxa(taxa_sums(Bacterial_VCMCFHS_AFilt1)>=1, Bacterial_VCMCFHS_AFilt1)
  #phyloseq-class experiment-level object
  #otu_table()   OTU Table:         [ 157 taxa and 333 samples ]
  #sample_data() Sample Data:       [ 333 samples by 24 sample variables ]
  #tax_table()   Taxonomy Table:    [ 157 taxa by 7 taxonomic ranks ]
  #phy_tree()    Phylogenetic Tree: [ 157 tips and 156 internal nodes ]
  #refseq()      DNAStringSet:      [ 157 reference sequences ]
  
  save(Bacterial_VCMCFHS_AFilt1, file="./Bacterial_Vaginal_Features_Phyloseq_Obj_Filtered1.RData")
  sequences <- refseq(Bacterial_VCMCFHS_AFilt1)
  writeXStringSet(sequences, filepath = "./rep-seqs-Bacterial_VCMCFHS_AFilt1.fasta")
  
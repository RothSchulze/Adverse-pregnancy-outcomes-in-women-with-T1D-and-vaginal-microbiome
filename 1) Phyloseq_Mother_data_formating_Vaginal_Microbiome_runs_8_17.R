  library("phyloseq")
  library("ggplot2")
  library("devtools")
  library("metagMisc")
  library("metagenomeSeq")
  
  ####################################################################
  ##                                                               ###
  ##            IMPORT AND SUBSET DATA                             ###
  ##             DATA WITH TAXA INFO                               ###
  ####################################################################
  
  ## Vaginal samples were sequenced in 2 different runs. All the sequences in the different runs were processed separatelly using the same parameters with QIIME2.
  
  setwd("~/Documents/ENDIA Project/Data and analyses/ENDIA Projects analysis/ENDIA_Run_17/") # HERE YOU SET YOUR DIRECTORY WHERE DATA IS CONTAINED
  
  ### Import data Non Filtered but already with wells merged
  Mother_V_wT <- import_biom("./table200-160_Run8_17_wTaxMet.biom", parseFunction = parse_taxonomy_greengenes) #
  tax_table(Mother_V_wT) <- tax_table(Mother_V_wT)[, 1:7]
 
  NewData <- read.table(file="./Mapping_file_Run_8_17.csv", header=TRUE, row.names=1, sep="\t")
  
   ### Import tree
  Phylogenetic_Tree <- read_tree("./rooted-tree_Vaginal.nwk")
  
  ### Making a phyloseq object 
  Mother_VC <- phyloseq(otu_table(Mother_V_wT),  sample_data(NewData), tax_table(Mother_V_wT), phy_tree(Phylogenetic_Tree))
  
  #phyloseq-class experiment-level object
  #otu_table()   OTU Table:         [ 3877 taxa and 798 samples ]
  #sample_data() Sample Data:       [ 798 samples by 8 sample variables ]
  #tax_table()   Taxonomy Table:    [ 3877 taxa by 7 taxonomic ranks ]
  #phy_tree()    Phylogenetic Tree: [ 3877 tips and 3691 internal nodes ]
  
  ## Merge the wells based on the SubjectID and trimester (well replicates) ##
  Mother_VCM <-merge_samples(Mother_VC, "Merging", fun = sum)
  ### Recover sample data after merging
  New_sample_data <- read.table(file="./Metadata_8_17_Merged.csv", header=TRUE, sep="\t", row.names=1)
  ### Making a phyloseq object 
  Mother_VCMC <- phyloseq(otu_table(Mother_VCM),  sample_data(New_sample_data), tax_table(Mother_VCM), phy_tree(Mother_VCM))
  Mother_VCMC <-prune_taxa(taxa_sums(Mother_VCMC)>=1, Mother_VCMC)
  
  ##FILTER SEQUENCES CLASSIFIED AS UNASSIGNED
  
  #>0faa013f28eefe0075f29d086e1a3d78 (Homo sapiens - mitochondrion)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATCTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCTGCTGGCACGCTGCGATGTGC
  #>1e32ae2ce24bb1630642f001cdf8230b (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATTTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCGGCTGACACGCTGCGATGTGC
  #>24644d3b88ea950eb3808da0677e6585 (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATCTAAAACACTCTTTACGCCGGTTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCTGCTGGCACGCTGCGATGTGC
  #>3f10c7baa45d1b1ad223c1335553714d (Homo sapiens)
  #GCACGGATATATTGTGTTAAGTAAAATGGGCTCAAGTTTCACCATGGGAGTATATAAAGGATCAGGATCAAAGTTCAGTTCCTCCAACTGCCCTAACATAAAAACAAAGTGTGGGGCCCGGGCACGGTGGCTCACCCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCAGGCGGATCACCTGAGGTCGGGAGTTCGAGACCAGCCTGGCCAGCATGGTGAAACCCTGTCCCTACTAAAAATGCAAAAA
  #>49400b8d1e87279d5c83721cb2a3bd06 (Homo sapiens)
  #CCCAGTTTGGGTCTTAGTTATTCTGTGTTCAGATGCATTAAAGCCACTTTCGTAGTTTATTTTGTATCAACTGGAGTTTTTTACAACTCACGTGAATTTTAGCTTTATTGAGGGGAATTGATCTAAAACACTCTTTATGCCGGTTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCGGCTGGCACGCTGCGATGTGCA
  #>5dc75f90143aa434749a58b5acc3468a (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATCTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCGGCTGACACGCTGCGATGTGC
  #>79e42ea2548ffe33aac88734e05800aa (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATCTAAAACACTCTTTACGCCGGTTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCGGCTGACACGCTGCGATGTGC
  #>810041e3b8768fc658970f2d1f32aafc (Homo sapiens)
  #CCCAGTTTGGGTCTTAGTTATTCTGTGTTCAGATGCATTAAAGCCACTTTCGTAGTTTATTTTGTATCAACTGGAGTTTTTTACAACTCACGTGAATTTTAGCTTTATTGAGGGGAATTGATCTAAAACACTCTTTATGCCGGTTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCTGCTGACACGCTGCGATGTGCA
  #>8962925f25819e4d87d865143e9a517c (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATTTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCGGCTGGCACGCTGCGATGTGC
  #>9e99c7af2554f7f36548b0990c9fcabc (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGCGGGGGTGATCTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCGGCTGGCACGCTGCGATGTGC
  #>9f3553d7d2c1df60041a1e646a55424f (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATCTAAAACACTCTTTACGCCGGTTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCTGCTGACACGCTGCGATGTGC
  #>af8b688ab204d05db355a9ad84864b7a (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATCTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCTGCTGACACGCTGCGATGTGC
  #>b9c88e582fc19c7079506f1f23c797a4 (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGCGGGGGTGATCTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCTGCTGACACGCTGCGATGTGC
  #>cb6996812a8b778c2e41cc2f28ac43a6 (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATCTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCGGCTGGCACGCTGCGATGTGC
  #>d0ed3e3343d151f9e33cda23a363a07c (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATCTAAAACACTCTTTACGCCGGTTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCGGCTGGCACGCTGCGATGTGC
  #>ecd6799d260cba9a5dc32371f7bf4c9b (Homo sapiens)
  #CCCAGTTTGGGTCTTAGCTATTGTGTGTTCAGATATGTTAAAGCCACTTTCGTAGTCTATTTTGTGTCAACTGGAGTTTTTTACAACTCAGGTGAGTTTTAGCTTTATTGGGGAGGGGGTGATTTAAAACACTCTTTACGCCGGCTTCTATTGACTTGGGTTAATCGTGTTACCGCGGCTGCTGACACGCTGCGATGTGC
  
  Mother_VCMCFHS <- subset_taxa(Mother_VCMC, Kingdom != "Unassigned")
  
  #phyloseq-class experiment-level object
  #otu_table()   OTU Table:         [ 2263 taxa and 334 samples ]
  #sample_data() Sample Data:       [ 334 samples by 4 sample variables ]
  #tax_table()   Taxonomy Table:    [ 2263 taxa by 7 taxonomic ranks ]
  #phy_tree()    Phylogenetic Tree: [ 2263 tips and 2078 internal nodes ]
  
  ####################################################################
  ##                                                               ###
  ##                  ABUNDANCE FILTER Features                    ###
  ##                                                               ###
  ####################################################################
  
  Mother_VCMCFHS_dataframe <- as.data.frame(otu_table(Mother_VCMCFHS))
  
  # function to perform pre-filtering
  low.count.removal = function(data, percent=0.01)
  {
    keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
    data.filter = data[,keep.otu]
    return(list(data.filter = data.filter, keep.otu = keep.otu))
  }
  
  result.filter = low.count.removal(Mother_VCMCFHS_dataframe, percent=0.01)
  length(result.filter$keep.otu) 
  Keep_OTUS <- names(result.filter$keep.otu)
  
  Mother_VCMCFHS_AFilt1 <- prune_taxa(Keep_OTUS, Mother_VCMCFHS)
  
  #phyloseq-class experiment-level object
  #otu_table()   OTU Table:         [ 159 taxa and 334 samples ]
  #sample_data() Sample Data:       [ 334 samples by 4 sample variables ]
  #tax_table()   Taxonomy Table:    [ 159 taxa by 7 taxonomic ranks ]
  #phy_tree()    Phylogenetic Tree: [ 159 tips and 158 internal nodes ]
  
  
  #############################################################
  ##                                                        ###
  ##                  OTUs FROM FEATURES                    ###
  ##                                                        ###
  #############################################################
  
  # tip_glom: Can be used to create a non-trivial OTU Table, if a phylogenetic tree is available. It is use as a means to merge taxa in our dataset that are closely related. 
  # In this case, we specify a threshold patristic distance. Taxa more closely related than this threshold are merged. This is especially useful when a dataset has many taxa 
  # that lack a taxonomic assignment at the level you want to investigate
  # tip_glom agglomerates OTUs at a given height in the tree. All tips of the tree separated by a cophenetic distance smaller than h will be agglomerated into one taxa using merge_taxa.
  
  Mother_VCM_OTU <- tip_glom(Mother_VCMCFHS, h = 0.03)
  
  ## Staying only with OTUs that have at least one sequence
  Mother_VCM_OTU <-prune_taxa(taxa_sums(Mother_VCM_OTU)>=1, Mother_VCM_OTU)
  
  #phyloseq-class experiment-level object
  #otu_table()   OTU Table:         [ 616 taxa and 334 samples ]
  #sample_data() Sample Data:       [ 334 samples by 4 sample variables ]
  #tax_table()   Taxonomy Table:    [ 616 taxa by 7 taxonomic ranks ]
  #phy_tree()    Phylogenetic Tree: [ 616 tips and 615 internal nodes ]
  

  ####################################################################
  ##                                                               ###
  ##                  ABUNDANCE FILTER                             ###
  ##                                                               ###
  ####################################################################
  
  Mother_VCM_OTU_dataframe <- as.data.frame(otu_table(Mother_VCM_OTU))
  
  result.filter = low.count.removal(Mother_VCM_OTU_dataframe, percent=0.01)
  length(result.filter$keep.otu) 
  Keep_OTUS <- names(result.filter$keep.otu)
  
  Mother_VCM_OTU_AFilt1 <- prune_taxa(Keep_OTUS, Mother_VCM_OTU)
  
  Mother_VCM_OTU_AFilt1
  #otu_table()   OTU Table:         [ 89 taxa and 334 samples ]
  #sample_data() Sample Data:       [ 334 samples by 4 sample variables ]
  #tax_table()   Taxonomy Table:    [ 89 taxa by 7 taxonomic ranks ]
  #phy_tree()    Phylogenetic Tree: [ 89 tips and 88 internal nodes ]
  
  save(Mother_VCM_OTU_AFilt1, file="./Mothers_Vaginal_OTU_Phyloseq_Obj_Filtered1.RData")
  
  
  
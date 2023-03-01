### data analysis with microeco package - agricultural waste paper - Guilherme Martins (16/07/2022)

# load the package
library(microeco)
library(magrittr)
library(ggplot2)

# Get and set work directory
getwd()
setwd("C:/Users/gui_l/OneDrive/Mestrado/Artigo 1 - De novo")
path <- "C:/Users/gui_l/OneDrive/Mestrado/Artigo 1 - De novo"
list.files(path)

# load the example data; 16S rRNA gene amplicon sequencing dataset
# metadata table
sample_info_16S <- read.table("sample_info_16S.txt", header = TRUE, row.names = 1, sep = "\t", dec = ",")

# feature table
otu_table_16S <- read.table("otu_table_16S.txt", header = TRUE, row.names = 1, sep = "\t", dec = ",")

# taxonomic assignment table
taxonomy_table_16S <- read.table("taxonomy_table_16S.txt", header = TRUE, row.names = 1, sep = "\t", dec = ",")

# load the environmental data table if it is not in sample table
env_data_16S <- read.table("env_data_16S.txt", header = TRUE, row.names = 1, sep = "\t", dec = ",")

# set.seed is used to fix the random number generation to make the results repeatable
set.seed(123)

# make the plotting background same with the tutorial
theme_set(theme_bw())

# Make sure that the data types of sample_table, otu_table and tax_table are all data.frame format as the following part shows
class(otu_table_16S)
otu_table_16S[1:5, 1:5]

class(taxonomy_table_16S)
taxonomy_table_16S[1:5, 1:3]


#######################################################################
            ### Organize data for further analysis ###


# Remove NAs to proceed with diversity analysis
# make the taxonomic information unified, very important
taxonomy_table_16S %<>% tidy_taxonomy

# make sure that the rownames of sample information table are sample names
class(sample_info_16S)
sample_info_16S[1:5, ]

# make sure that the environmental data are stored as env data
class(env_data_16S)
env_data_16S[1:5, 1:5]

# In R6 class, '$new' is the original method used to create a new object of class
# Let's create a microtable object with more information
dataset <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S)
dataset

# reorder groups
dataset$sample_table$Group %<>% factor(., levels = c("FC", "PL", "CM"))
dataset$sample_table$Type %<>% factor(., levels = c("Fresh", "Composted"))
str(dataset$sample_table)

# remove OTUs which are not assigned in the Kingdom "k__Archaea" or "k__Bacteria"
dataset$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
dataset

# make the OTU and sample information consistent across all files in the dataset object
dataset$tidy_dataset()
print(dataset)

# let's use sample_sums() to check the sequence numbers in each sample
dataset$sample_sums() %>% range

# Rarefy data to reduce bias in diversity index
dataset$rarefy_samples(sample.size = 23475)

# Check the sequence numbers again
dataset$sample_sums() %>% range

# let's calculate the taxa abundance at each taxonomic rank
# use default parameters
dataset$cal_abund()

# return dataset$taxa_abund
class(dataset$taxa_abund)

# show part of the relative abundance at Phylum level
dataset$taxa_abund$Phylum[1:5, 1:5]

# The function save_abund() can be used to save the taxa abundance file to a local place easily.
dataset$save_abund(dirpath = "taxa_abund")




#######################################################################
### Calculate alpha and beta diversity ###

### Alpha diversity ###

# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = FALSE)

# return dataset$alpha_diversity
class(dataset$alpha_diversity)

# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")

# calculate beta diversity

# unifrac = FALSE means do not calculate unifrac metric
# require GUniFrac package installed
dataset$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")




################################################################################


### Composition-based class ###

# taxa_abund list in the object of microtable class must be first calculated
# create trans_abund object


# use 10 Phyla with the highest abundance in the dataset.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 12)


# Relative abundance with two facets at Phylum level
# require package ggh4x, first run install.packages("ggh4x") if not installed
p1 <- t1$plot_bar(others_color = "grey70", facet = "Group", facet2 = "Type", xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)
p1
# save plot with 600 dpi resolution
dev.print(tiff, "relative_abundance_facet.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# show the heatmap with the high abundant genera
# show 30 taxa at Genus level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 30)

# Plot
t1$plot_heatmap(facet = "Group", 
                facet2 = "Type",  
                xtext_keep = FALSE, 
                withmargin = TRUE)


#  save plot with 600 dpi resolution
dev.print(tiff, "heatmap_abundance2.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# The trans_venn class is used for venn analysis, i.e. shared and unique taxa
# To analyze the unique and shared OTUs of groups, we first merge samples according to the "Group" column of sample_table

# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Group")
# dataset1 is a new microtable object
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn() # The integer is OTU number # The percentage data is the sequence number/total sequence number
#  save plot with 600 dpi resolution
dev.print(tiff, "venn_diagram.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# transform the results of venn plot to the traditional feature-sample table, that is, another object of microtable class
dataset1 <- dataset$merge_samples(use_group = "Group")
t1 <- trans_venn$new(dataset1)

# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)


# Venn diagram by group
dataset1 <- dataset$merge_samples(use_group = "Group")
# dataset1 is a new microtable object
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn() # The integer is OTU number # The percentage data is the sequence number/total sequence number
#  save plot with 600 dpi resolution
dev.print(tiff, "venn_diagram2.tiff", compression = "lzw", res=600, height=15, width=15, units="cm")




################################################################################

### Diversity-based class ###

# The data_alpha is used for the following differential test and plotting
t1 <- trans_alpha$new(dataset = dataset, group = "Unified")
# return t1$data_stat
t1$data_stat[1:5, ]

# test the differences among groups using:
# Kruskal-Wallis Rank Sum Test (overall test when groups > 2)
t1$cal_diff(method = "KW")
# return t1$res_diff
t1$res_diff[1:5, ]

# Dunn's Kruskal-Wallis Multiple Comparisons (for paired groups when groups > 2)
t1$cal_diff(method = "KW_dunn")
# return t1$res_diff
t1$res_diff[1:5, ]

# ANOVA with multiple comparisons
t1$cal_diff(method = "anova")
# return t1$res_diff
t1$res_diff[1:5, ]

# plot the mean and se of alpha diversity for each group, and add the anova result
t1$cal_diff(method = "KW")
t1$plot_alpha(measure = "Chao1")
t1$plot_alpha(measure = "ACE")
t1$plot_alpha(measure = "Simpson")
t1$plot_alpha(measure = "Shannon")
t1$plot_alpha(measure = "Observed")


# Clustering plot is also a frequently used method

# plot and compare the group distances
# calculate and plot sample distances within groups
t1$cal_group_distance()
# return t1$res_group_distance
t1$plot_group_distance(distance_pair_stat = TRUE)

# calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$plot_group_distance(distance_pair_stat = TRUE)

# use replace_name to set the label name, group parameter used to set the color
t1$plot_clustering(group = "Unified", replace_name = c("Group", "Type"))
#  save plot with 600 dpi resolution
dev.print(tiff, "clustering.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")



# PERMANOVA is often used in the differential test of distances among groups
# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova

# manova for each paired groups
t1$cal_manova(manova_all = FALSE)
t1$res_manova

# PERMDISP is implemented to check multivariate homogeneity of groups dispersions (variances)
# for the whole comparison and for each paired groups
t1$cal_betadisper()
t1$res_betadisper




##########################################################################

###  Model-based class ###

# the differential test result is stored in the object$res_diff

# ALDex2 analysis at Genus level
t2 <- trans_diff$new(dataset = dataset,
                     method = "ALDEx2_t",
                     group = "Unified",
                     taxa_level = "Genus",
                     ALDEx2_sig = c("we.ep")
)


# metastat analysis at Genus level
#t1 <- trans_diff$new(dataset = dataset, method = "metastat", group = "Unified", taxa_level = "Genus")

# t1$res_diff is the differential test result
t2$res_diff$Taxa

# t1$res_abund is the group abundance

# the metastat can run the comparisons for each paired group
# the user should use select_group to select the required pair

# select_group should be one of groups in t1$res_diff$Comparison
# FC_f vs FC_c differential abundance
fc_dif <- t1$plot_diff_abund(use_number = 1:15,
                   col = c("tomato2", "tomato4"),
                   select_group = "FC_f - FC_c",
                   add_sig = TRUE,
                   use_se = TRUE,
                   add_sig_label = "Significance",
                   coord_flip = FALSE,
                   keep_prefix = FALSE)

fc_dif


# PL_f vs PL_c differential abundance
pl_dif <- t1$plot_diff_abund(use_number = 1:15,
                             col = c("seagreen2", "seagreen4"),
                             select_group = "PL_f - PL_c",
                             add_sig = TRUE,
                             use_se = TRUE,
                             add_sig_label = "Significance",
                             coord_flip = FALSE,
                             keep_prefix = FALSE)

pl_dif

# CM_f vs CM_c differential abundance
cm_dif <- t1$plot_diff_abund(use_number = 1:15,
                             col = c("steelblue2", "steelblue4"),
                             select_group = "CM_f - CM_c",
                             add_sig = TRUE,
                             use_se = TRUE,
                             add_sig_label = "Significance",
                             coord_flip = FALSE,
                             keep_prefix = FALSE)

cm_dif

# Arrange plots on one page
ggarrange(fc_dif, pl_dif, cm_dif,
          labels = c("A", "B", "C"),
          common.legend = FALSE,
          legend = "right",
          ncol = 1, nrow = 3)

#  save plot with 800 dpi resolution
dev.print(tiff, "diff_abundance.tiff", compression = "lzw", res=600, height=30, width=15, units="cm")


# LEfSe combines the non-parametric test and linear discriminant analysis
t1 <- trans_diff$new(dataset = dataset,
                     method = "lefse",
                     group = "Type",
                     alpha = 0.01,
                     lefse_subgroup = NULL
                     )

# see t1$res_diff for the result
lefse_table <- t1$res_diff

lefse_abund <- t1$res_abund

# Write table
write.table(lefse_table, "lefse_table.txt", sep = "\t", row.names = FALSE)



# From v0.8.0, threshold is used for the LDA score selection.
p2 <- t1$plot_diff_bar(threshold = 4)
p2
#  save plot with 600 dpi resolution
dev.print(tiff, "LEfSE_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# ANOVA method and transpose
t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Unified", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:15, col = c("tomato2", "tomato4", "seagreen2", "seagreen4", "steelblue2", "steelblue4"), add_sig = T, coord_flip = F)
#  save plot with 600 dpi resolution
dev.print(tiff, "ANOVA_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=35, units="cm")


# Try to run those examples
# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = dataset, method = "KW", group = "Unified", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:15, col = c("tomato2", "tomato4", "seagreen2", "seagreen4", "steelblue2", "steelblue4"), add_sig = T, coord_flip = F)
#  save plot with 600 dpi resolution
dev.print(tiff, "KW_diff_abundance.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")





#################################################################################

### Explainable class ###

# If there may be some NA in the user's env data,add complete_na = TRUE when creating the trans_env object

# add_data is used to add the environmental data
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)

# show the autocorrelations among variables
# use group parameter to show the distributions of variables and the autocorrelations across groups
t1$cal_autocor(group = "Unified", col = c("tomato2", "tomato4", "seagreen2", "seagreen4", "steelblue2", "steelblue4"))
#  save plot with 600 dpi resolution
dev.print(tiff, "autocorr.tiff", compression = "lzw", res=600, height=40, width=40, units="cm")


# RDA at Genus level
t1$cal_ordination(method = "RDA", taxa_level = "Family")
t1$res_ordination
t1$res_ordination_R2
t1$res_ordination_terms
t1$res_ordination_axis



# Results from RDA model:
# Model: rda(formula = use_data ~ pH + C.org + N.total + N.inorg + P + K + Ca + Mg + S + Cu + Zn, data = env_data)

# Df Variance       F Pr(>F)    
# pH        1 13405157 98.5914  0.001 ***
# C.org     1  4623820 34.0070  0.001 ***
# N.total   1  4385379 32.2533  0.001 ***
# N.inorg   1   219733  1.6161  0.207    
# P         1  2806280 20.6394  0.001 ***
# K         1  1442067 10.6060  0.001 ***
# Ca        1   599719  4.4108  0.009 ** 
# Mg        1    69795  0.5133  0.735    
# S         1   201020  1.4784  0.225    
# Cu        1   536438  3.9454  0.020 *  
# Zn        1   113112  0.8319  0.487    
# Residual  6   815801                   

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# As the main results of RDA are related with the projection and angles between different arrows,
# we adjust the length of the arrow to show them clearly using several parameters.
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)

# t1$res_rda_trans is the transformed result for plotting
rda_1 <- t1$plot_ordination(plot_color = "Group",
                   plot_shape = "Type",
                   col = c("tomato3", "seagreen3", "steelblue3"),
                   plot_type = c("point")
                   )


#  save plot with 600 dpi resolution
dev.print(tiff, "RDA_genus.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


# Mantel test can be used to check significant correlations between environmental variables and distance matrix
t1$cal_mantel(use_measure = "bray")
t1$res_mantel


# perform a correlation heatmap between environmental and taxa data
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)
t1$cal_cor(use_data = "other",
           other_taxa = t2$res_diff$Taxa[1:20],
           cor_method = "spearman",
           p_adjust_method = "fdr",
           by_group = "Group",
           taxa_name_full = FALSE
           )


# plot the correlation results using plot_cor function
# filter phylum that do not have at least one ***
 corrpheat <- t1$plot_cor(color_vector = c("#A50026", "white", "#053061"),
             filter_feature = c("*", "**", "***"),
             pheatmap = TRUE,
             ylab_type_italic = TRUE,
             cluster_ggplot = "col"
             )


#  save plot with 600 dpi resolution
dev.print(tiff, "env_corr.tiff", compression = "lzw", res=600, height=15, width=25, units="cm")


# pheatmap method
# clustering heatmap; require pheatmap package
# Let's take another color pallete: install.packages('pheatmap')
t1$plot_cor(pheatmap = TRUE, filter_feature = c("", "*", "**"), color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))
#  save plot with 600 dpi resolution
dev.print(tiff, "pheatmap_phylum_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


#correlations between environmental variables and taxa for different groups
# calculate correlations for different groups using parameter by_group

# first create trans_diff object as a demonstration
t2 <- trans_diff$new(dataset = dataset, method = "rf", group = "Unified", rf_taxa_level = "Genus")
# then create trans_env object
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)
# use other_taxa to select taxa you need
t1$cal_cor(by_group = "Unified", use_data = "other", cor_method = "spearman", 
           p_adjust_method = "fdr", other_taxa = t2$res_diff$Taxa[1:25])
t1$plot_cor()
#  save plot with 600 dpi resolution
dev.print(tiff, "env_corr.tiff", compression = "lzw", res=600, height=15, width=40, units="cm")


# relationship between environmental factors and alpha diversity
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)
# use add_abund_table parameter to add the extra data table
t1$cal_cor(add_abund_table = dataset$alpha_diversity, cor_method = "spearman")
# Let's try to use ggplot2 with clustering plot
t1$plot_cor(cluster_ggplot = "row")
#  save plot with 600 dpi resolution
dev.print(tiff, "env_div_corr.tiff", compression = "lzw", res=600, height=15, width=20, units="cm")


#######################################################
# Merge plots with ggpubr
library(ggpubr)


### Arrange plots on one page ###
ggarrange(rda_1, corrpheat, 
          labels = c("A", "B"),
          common.legend = FALSE,
          legend = "right",
          ncol = 1, nrow = 2)

#  save plot with 600 dpi resolution
dev.print(tiff, "FIG2.tiff", compression = "lzw", res=600, height=25, width=20, units="cm")

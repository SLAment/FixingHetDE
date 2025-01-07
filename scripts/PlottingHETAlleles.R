#!/usr/bin/env Rscript

# hnwdAlleles: Visualizing hnwd gene alleles
#############################################################################
# =======================================
# S. Lorena Ament-Velasquez
# 2024-08-06
#############################################################################
# ============================
# Load the necessary libraries
# ============================
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

# ============================
# Input files
# ============================
# # setwd(system("pwd", intern = T))
# alleles <- read.table("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/IFBMontpellier/02_RepeatExplorer/reports/WD40_classification_10-11-12-14-30-32-39.txt", header = TRUE) %>% mutate(a_strain = paste0(allele, "_", strain))
# 
# # Output
# assem_pfile <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/IFBMontpellier/02_RepeatExplorer/results/WD40_assemblies_hets.pdf"
# hetE_pfile <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/IFBMontpellier/02_RepeatExplorer/results/WD40_assemblies_hetE.pdf"
# hetD_pfile <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/IFBMontpellier/02_RepeatExplorer/results/WD40_assemblies_hetD.pdf"
# oldVsNew_pfile <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/IFBMontpellier/02_RepeatExplorer/results/WD40_assemblies_oldsVsNew.pdf"

## Snakemake
# Input
alleles <- read.table(snakemake@input$alleles, header = TRUE) %>% mutate(a_strain = paste0(allele, "_", strain))

# Output
assem_pfile <- snakemake@output$HETs
hetD_pfile <- snakemake@output$hetD
hetE_pfile <- snakemake@output$hetE
oldVsNew_pfile <- snakemake@output$oldsVsNew

# ============================
# Functions
# ============================
getREPcolors <- function(repnumberlist, palette = "Set3", mincol = 11){
  # Get the unique values and generate a color palette
  unique_labels <- unique(repnumberlist)
  num_colors <- length(unique_labels)
  color_palette <- brewer.pal(min(mincol, num_colors), palette)
  if (num_colors > mincol) {
    color_palette <- colorRampPalette(color_palette)(num_colors)
  }
  
  # Create a named vector for the palette
  named_palette <- setNames(color_palette, unique_labels)
  named_palette["0"] <- "black"
  return(named_palette)
}

# ============================
# Compare assemblies
# ============================
## Compare long and short-read assemblies
assemblystrains <- c("PaYp", "PaYp-default", "PaYp-default2", "PaYp-allkmers", "PaYp-allkmers2", 
                     "PaZp", "PaZp-default", "PaZp-default2", "PaZp-allkmers", "PaZp-allkmers2", 
                     "PaWa63p", "PaWa63p-default", "PaWa63p-default2", "PaWa63p-allkmers", "PaWa63p-allkmers2",
                     "PaWa137m", "PaWa137m-allkmers")

#### het-d
# Extract the het-d alleles
hetD <- alleles %>% filter(geneid == "het-d")

hetD_named_palette <- getREPcolors(hetD$repnumber) # Get colors for the repeats

hetD_assem <- filter(hetD, strain %in% assemblystrains)
hetD_assem$strain <- factor(hetD_assem$strain, levels = assemblystrains) # Force the levels to be in that order

# Create the plot
hetD_assem_p <- ggplot(hetD_assem, aes(x = position, y = strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", linewidth = 2) +
  scale_fill_manual(values = hetD_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(levels(hetD_assem$strain))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  ggtitle("het-d") +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face="italic"),
        legend.position = "none",
        panel.grid = element_blank())

#### het-e
# Extract the het-e alleles
hetE <- alleles %>% filter(geneid == "het-e")

hetE_named_palette <- getREPcolors(sort(hetE$repnumber), palette = "Spectral") # Get colors for the repeats

hetE_assem <- filter(hetE, strain %in% assemblystrains)
hetE_assem$strain <- factor(hetE_assem$strain, levels = assemblystrains) # Force the levels to be in that order

# Create the plot
hetE_assem_p <- ggplot(hetE_assem, aes(x = position, y = strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", linewidth = 2, alpha = 0.7) +
  scale_fill_manual(values = hetE_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(levels(hetE_assem$strain))) +
  scale_x_continuous(limits = c(0, 13)) +
  theme_minimal() +
  ggtitle("het-e") +
  theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face="italic"),
        legend.position = "none",
        panel.grid = element_blank())

#### het-r
# Extract the het-r alleles
# Unlike het-d and het-e, I won't show the alleles for het-r in the Wa pop because it's beside the point
hetR <- alleles %>% filter(geneid == "het-r", strain %in% assemblystrains)

hetR_named_palette <- getREPcolors(sort(hetR$repnumber), palette = "Set2", mincol = 8) # Get colors for the repeats

hetR_assem <- filter(hetR, strain %in% assemblystrains)
hetR_assem$strain <- factor(hetR_assem$strain, levels = assemblystrains) # Force the levels to be in that order

# Create the plot
hetR_assem_p <- ggplot(hetR_assem, aes(x = position, y = strain, fill = factor(repnumber), label = repnumber)) +
    geom_tile(color = "white", linewidth = 2) +
    scale_fill_manual(values = hetR_named_palette) +
    geom_text(color = "black") +
    scale_y_discrete(limits = rev(levels(hetR_assem$strain))) +
    theme_minimal() +
    scale_x_continuous(limits = c(0, 13)) +
    ggtitle("het-r") +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face="italic"),
          legend.position = "none",
          # axis.text = element_text(size = 14, face = "bold"),
          # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0),
          panel.grid = element_blank())

assem_p <- plot_grid(hetD_assem_p, hetE_assem_p, hetR_assem_p, nrow = 1, align = "hv")

ggsave(plot = assem_p,
       filename = assem_pfile,
       width = 18.1, height = 6.5)
# The plot is meant to be modified later manually in Inkscape
# A couple of Illumina sequences are not complete because the Ns produced an unclassified Cterm
## -------------------
# ============================
# Functional alleles het-e
# ============================
strainshetE <- c("CmEmm", "ChEhDap", "PaYp", "FJ897789", "CoEcp", "CoEfp", "PaZp")
# FJ897789 is missing the first 5 amino acids of the first repeat but those happen 
# to be conserved in all het-e sequences, so it's likely repeat 11 too. The last 
# black repeat happens because the last amino acids in the sequence correspond to 
# the cryptic repeat at the Cterm that doesn't get recognized as such in full 
# sequences.

functional_hetE <- hetE %>% filter(strain %in% strainshetE) 

# Fix levels
functional_hetE$strain <- factor(functional_hetE$strain, levels = strainshetE) # Force the levels to be in that order

functional_hetE_p <- ggplot(functional_hetE, aes(x = position, y = a_strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", size = 2, alpha = 0.7) +
  scale_fill_manual(values = hetE_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(unique(functional_hetE$a_strain))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) 

hetE_e4 <- hetE %>% filter(strain %in% c("Podan2", "CaDam"))
e4_alleles <- ggplot( hetE_e4, aes(x = position, y = a_strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", size = 2, alpha = 0.7) +
  scale_fill_manual(values = hetE_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(levels(as.factor(hetE_e4$a_strain)))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) 

hete_strains <- c("PaWa46p", "PaWa58m", "PaWa87p", "PaWa28m", "PaWa63p", "PaWa137m", "PaTgp", "PaWa100p", "PaWa21m", "PaWa53m")
hetE_Wa <- hetE %>% filter(strain %in% hete_strains)
hetE_Wa$a_strain <- factor(hetE_Wa$a_strain, levels = paste0("e_", hete_strains)) # Force the levels to be in that order

WaE_alleles <- ggplot( hetE_Wa, aes(x = position, y = a_strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", size = 2, alpha = 0.7) +
  scale_fill_manual(values = hetE_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(levels(as.factor(hetE_Wa$a_strain)))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) 

hetEplots <- plot_grid(functional_hetE_p, e4_alleles, WaE_alleles,
                       nrow = 3, align = "hv", rel_heights = c(1.7, 0.65, 2.1))

ggsave(plot = hetEplots,
       filename = hetE_pfile,
       width = 7, height = 9)

# ============================
# Functional alleles het-d
# ============================
strainshetD <- c("ChEhDap", "CsDfp", "PaYp")

functional_hetD <- hetD %>% filter(strain %in% strainshetD) 

# Fix levels
functional_hetD$strain <- factor(functional_hetD$strain, levels = strainshetD) # Force the levels to be in that order

functional_hetD_p <- ggplot(functional_hetD, aes(x = position, y = a_strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", size = 2, alpha = 0.7) +
  scale_fill_manual(values = hetD_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(unique(functional_hetD$a_strain))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) 

hetD_d3 <- hetD %>% filter(strain %in% c("PaZp", "CoEfp", "CmEmm", "Podan2"))
d3_alleles <- ggplot( hetD_d3, aes(x = position, y = a_strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", size = 2, alpha = 0.7) +
  scale_fill_manual(values = hetD_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = (levels(as.factor(hetD_d3$a_strain)))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) 

hetd_strains <- c("PaWa53m", "PaWa137m", "PaTgp", "PaWa87p", "PaWa100p", "PaWa58m", "PaWa21m", "PaWa46p", "PaWa28m", "PaWa63p")
hetD_Wa <- hetD %>% filter(strain %in% hetd_strains)
hetD_Wa$a_strain <- factor(hetD_Wa$a_strain, levels = paste0("d_", hetd_strains)) # Force the levels to be in that order

WaD_alleles <- ggplot( hetD_Wa, aes(x = position, y = a_strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", size = 2, alpha = 0.7) +
  scale_fill_manual(values = hetD_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(levels(as.factor(hetD_Wa$a_strain)))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) 

hetDplots <- plot_grid(functional_hetD_p, d3_alleles, WaD_alleles,
          nrow = 3, align = "hv", rel_heights = c(0.7, 0.9, 2))

ggsave(plot = hetDplots,
       filename = hetD_pfile,
       width = 6.5, height = 7.5)

# ============================
# Comparing alleles with previous publications
# ============================
hetd_YD2 <- c("PaYp", "likelyWrong")
hetD_YD2 <- hetD %>% filter(strain %in% hetd_YD2)

YD2_alleles <- ggplot( hetD_YD2, aes(x = position, y = a_strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", size = 2, alpha = 0.7) +
  scale_fill_manual(values = hetD_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(levels(as.factor(hetD_YD2$a_strain)))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) 

hetE_Ec <- hetE %>% filter(strain %in% c("CoEcp", "AF323582", "AF323583"))
Ec_alleles <- ggplot( hetE_Ec, aes(x = position, y = a_strain, fill = factor(repnumber), label = repnumber)) +
  geom_tile(color = "white", size = 2, alpha = 0.7) +
  scale_fill_manual(values = hetE_named_palette) +
  geom_text(color = "black") +
  scale_y_discrete(limits = rev(levels(as.factor(hetE_Ec$a_strain)))) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 13)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

oldseqsVsNew <- plot_grid(YD2_alleles, Ec_alleles,
          nrow = 2, align = "hv", rel_heights = c(1,1.3))

ggsave(plot = oldseqsVsNew,
       filename = oldVsNew_pfile,
       width = 6.5, height = 2.7)

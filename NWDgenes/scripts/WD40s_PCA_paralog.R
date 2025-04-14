#!/usr/bin/env Rscript

### PCA of the WD40 repeats in the HNWD family
#############################################################################
# =======================================
# Sandra Lorena Ament Velasquez
# 2023-03-08 - 2025-04-14
# =======================================
# https://arftrhmn.net/quick-pca-analysis-from-sequence-alignment-data-in-r/
# https://stats.stackexchange.com/questions/222/what-are-principal-component-scores
# https://thomasadventure.blog/posts/turning-your-ggplot2-code-into-a-function/
# https://builtin.com/data-science/step-step-explanation-principal-component-analysis

# ============================
# Load the necessary libraries
# ============================
library(adegenet) # It loads the package ade4 too
library(ggplot2)
library(dplyr)
library(ape)
library(patchwork)
library(cowplot)

# ============================
# Read input data
# ============================
# meta file that contains sequence features that we are interested to visualize 
# in the PCA. The meta file must have a column containing isolate ids that 
# matched fasta id as well.

## Input
meta <- read.table(snakemake@input$meta, sep='\t', header = T)
dnaseqs <- ape::read.dna(snakemake@input$fasta, format = "fasta")

## Output
pcaall <- snakemake@output$pcaall
pcacanon <- snakemake@output$pcacanon
pcasansr <- snakemake@output$pcasansr
pcahete <- snakemake@output$pcahete
pcahetd <- snakemake@output$pcahetd
paperfigfile <- snakemake@output$paper

# ============================
# Define some functions
# ============================
## General function to make a PCA dataframe for plotting
makePCAfromAlignment <- function(dnaseqs, nf = 10){
  # Make genind object from the DNAbin object
  multisites <- adegenet::DNAbin2genind(dnaseqs, polyThres=0)
  
  # Deal with missing data before making PCA
  # sum(is.na(multisites$tab)) 
  X <- adegenet::tab(multisites, freq = TRUE, NA.method = "mean") # replace NA by the mean allele frequencies
  
  # Calculate the actual Principal Component Analysis
  pca <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = nf)

  # Make a dataframe from the results
  pca.dataset <- data.frame(Sequence = row.names(pca$li),
                            as.data.frame(pca$li)) # the row coordinates i.e. the principal components

  # Re-assign the name of the sequence based on the Sequence variable
  pca.dataset <- base::merge(meta, pca.dataset, by = "Sequence")
  
  return(list(pca.dataset, pca))
}

## Get variance explained by a PC
getvarPCA <- function(pca, axis = 1){
  Eigenvalues <- pca$eig # Note to self: weights = loadings = eigenvalues
  Variance <- Eigenvalues / sum(Eigenvalues) 
  VarAxis <- (100 * signif(Variance[axis], 4)) %>% round(digits=2)
  return(VarAxis)
}

## Make a theme for the PCA
clean_theme <- function() {
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
}

colorgenes = c("het-r" = "#ffd42aff", 
               "het-d" = "#f8766cff", 
               "het-e" = "#00bd5cff" , 
               "hnwd1" = "#00a6ffff", 
               "hnwd3" = "#de950fff", 
               "nwd1" = "#00c1a7ff", 
               "nwd2" = "#64b100ff", 
               "nwdp-2" = "#01b8dcff", 
               "nwd5" = "#b88dffff",
               "nwd6" = "#f06eedff",
               "nwd3" = "#f867cdff")

shapesgenes = c("het-r" = 15, 
               "het-d" = 16, 
               "het-e" = 17 , 
               "hnwd1" = 18, 
               "hnwd3" = 21, 
               "nwd1" = 20, 
               "nwd2" = 24, 
               "nwdp-2" = 19, 
               "nwd5" = 23,
               "nwd6" = 22,
               "nwd3" = 25)

## Plot different combinations of PCs
# The Switch allows you to turn on and off the shape of the points:
# https://stackoverflow.com/questions/22915337/if-else-condition-in-ggplot-to-add-an-extra-layer
plot4pcas <- function(dnaseqs, classi = Gene, Switch = TRUE){
  pcadf <- makePCAfromAlignment(dnaseqs)[[1]]
  pcaraw <- makePCAfromAlignment(dnaseqs)[[2]]
  
  # # Make a list of shapes (more general but stochastic depending on the dataset)
  # nuclassi <- nrow(unique(pcadf[ deparse(substitute(classi)) ]))
  # # shapevals <- seq(1,nuclassi, 1) # i don't like the first shapes
  # shapevals <- c(seq(15, 15+nuclassi, 1), seq(1,nuclassi, 1))
  
  # Aesthetics
  alphaval <- 0.5
  
  p1 <- ggplot(pcadf, aes(Axis1, Axis2, colour= {{classi}})) +  # The !! is necessary to interpret the variable (unquote)
    geom_point({if(Switch)aes(shape={{classi}})}, size=3, alpha=alphaval) +
    scale_shape_manual(values=shapesgenes) +
    scale_color_manual(values = colorgenes) +
    xlab(paste0("PC1 (", getvarPCA(pcaraw, 1), " %)")) +
    ylab(paste0("PC2 (", getvarPCA(pcaraw, 2), " %)")) +
    theme_bw() +
    clean_theme() +
    theme(legend.position="none")
  
  p2 <- ggplot(pcadf, aes(Axis1, Axis3, colour= {{classi}})) + 
    geom_point({if(Switch)aes(shape={{classi}})}, size=3, alpha=alphaval) +
    scale_shape_manual(values=shapesgenes) +
    scale_color_manual(values = colorgenes) +
    xlab(paste0("PC1 (", getvarPCA(pcaraw, 1), " %)")) +
    ylab(paste0("PC3 (", getvarPCA(pcaraw, 3), " %)")) +
    theme_bw() +
    clean_theme() +
    theme(legend.position="none")
  
  p3 <- ggplot(pcadf, aes(Axis2, Axis3, colour= {{classi}})) + 
    geom_point({if(Switch)aes(shape={{classi}})}, size=3, alpha=alphaval) +
    scale_shape_manual(values=shapesgenes) +
    scale_color_manual(values = colorgenes) +
    xlab(paste0("PC2 (", getvarPCA(pcaraw, 2), " %)")) +
    ylab(paste0("PC3 (", getvarPCA(pcaraw, 3), " %)")) +
    theme_bw() +
    clean_theme() 
    # theme(legend.position="none")
  
  p4 <- ggplot(pcadf, aes(Axis3, Axis4, colour= {{classi}})) + 
    geom_point({if(Switch)aes(shape={{classi}})}, size=3, alpha=alphaval) +
    scale_shape_manual(values=shapesgenes) +
    scale_color_manual(values = colorgenes) +
    xlab(paste0("PC3 (", getvarPCA(pcaraw, 3), " %)")) +
    ylab(paste0("PC4 (", getvarPCA(pcaraw, 4), " %)")) +
    theme_bw() +
    clean_theme() +
    theme(legend.position="none")
  
  return((p1 + p2 + p3))
  # return((p1 + p2) / ( p3 + p4))
}

# ---- Make different PCA combinations ---
## All Genes
allfamily <- meta %>% 
  filter(Gene %in% c("het-r", "het-d", "het-e", "hnwd1", "hnwd3", "nwd1", "nwd2", "nwd3", "nwd5", "nwd6", "nwdp-2")) %>% 
  filter(Species == "anserina")
dnaseqs_all <- dnaseqs[which(labels(dnaseqs) %in% allfamily$Sequence),]

ggsave(plot = plot4pcas(dnaseqs_all), 
       filename = pcaall, 
       width = 10, height = 3)

## The het-e relatives
familyhete <- meta %>% filter(Gene %in% c("het-e", "nwd6", "nwdp-2", "hnwd3")) %>% 
  filter(Species %in% c("anserina"))
dnaseqs_familyhete <- dnaseqs[which(labels(dnaseqs) %in% familyhete$Sequence),]

ggsave(plot = plot4pcas(dnaseqs_familyhete), 
       filename = pcahete, 
       width = 10, height = 3)

## The het-d relatives
familyhetd <- meta %>% filter(Gene %in% c("het-d", "hnwd1")) %>% 
  filter(Species %in% c("anserina"))
dnaseqs_familyhetd <- dnaseqs[which(labels(dnaseqs) %in% familyhetd$Sequence),]

plot4pcas(dnaseqs_familyhetd) + geom_point(aes(label=Sequence))
dnaseqs <- dnaseqs_familyhetd

ggsave(plot = plot4pcas(dnaseqs_familyhetd), 
       filename = pcahetd, 
       width = 10, height = 3)

## Put them together for the paper
paperfig <- plot_grid(plot4pcas(dnaseqs_all), plot4pcas(dnaseqs_familyhete), nrow=2, labels=c('b', 'c'), align = "hv")
ggsave(file = paperfigfile, plot = paperfig, width = 10, height = 6)


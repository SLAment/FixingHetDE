#!/usr/bin/env Rscript

### LOGOS of the WD40 repeats in the NWD family
#############################################################################
# =======================================
# Sandra Lorena Ament Velasquez
# 2024-06-18 - 2025-04-14
# =======================================

# ============================
# Load the necessary libraries
# ============================
library(ggplot2)
library(dplyr)
library(ggseqlogo)

# ============================
# Inputs and outputs
# ============================
# Snakemake Input
hete <- snakemake@input$hete
hetd <- snakemake@input$hetd
hetr <- snakemake@input$hetr
hnwd1 <- snakemake@input$hnwd1
hnwd3 <- snakemake@input$hnwd3
nwd1 <- snakemake@input$nwd1
nwd2 <- snakemake@input$nwd2
nwd3 <- snakemake@input$nwd3
nwd5 <- snakemake@input$nwd5
nwd6 <- snakemake@input$nwd6
nwdp2 <- snakemake@input$nwdp2

# Snakemake Output
eachgene <- snakemake@output$eachgene

# ============================
# Functions 
# ============================
## Read the fasta into a vector character for ggseqlogo
fasta2vector <- function(fasta_file){
  fasta_contents <- readLines(fasta_file)
  # Initialize variables
  sequences <- c()
  current_sequence <- ""
  
  # Parse the FASTA file contents
  for (line in fasta_contents) {
    if (startsWith(line, ">")) {
      # If a sequence was being read, add it to the list
      if (current_sequence != "") {
        sequences <- c(sequences, current_sequence)
        current_sequence <- ""
      }
    } else {
      # Concatenate lines to form the sequence
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  # Add the last sequence read to the list
  if (current_sequence != "") {
    sequences <- c(sequences, current_sequence)
  }
  return(sequences)
}


# # Print the character vector to verify
# hetd <- fasta2vector(fasta_file)
# 
# ggseqlogo( hetd, method = 'bits' ) # default
# ggseqlogo( hetd, method = 'prob' )
# ggseqlogo( hetd, col_scheme='chemistry' ) # For amino acids you can pick chemistry, hydrophobicity, clustalx, taylor

# ============================
# Make Logos 
# ============================
# Make a list with all the genes
allseqs <- list("het-e" = fasta2vector(hete), 
                nwd6 = fasta2vector(nwd6),
                "nwdp-2" = fasta2vector(nwdp2),
                hnwd3 = fasta2vector(hnwd3), 
                "het-r" = fasta2vector(hetr), 
                "het-d" = fasta2vector(hetd),
                hnwd1 = fasta2vector(hnwd1),
                nwd1 = fasta2vector(nwd1), 
                nwd2 = fasta2vector(nwd2), 
                nwd3 = fasta2vector(nwd3), 
                nwd5 = fasta2vector(nwd5)
                )

anserina <- as.vector(unlist(allseqs))
eachgene_plot <- ggseqlogo(c(allseqs, list(All = anserina)), ncol = 1) + 
  # theme(strip.text = element_text(size = 25, face = "italic"))
  theme(strip.text = element_text(size = 25), 
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

ggsave(plot = eachgene_plot, 
       filename = eachgene, 
       width = 11, height = 20)

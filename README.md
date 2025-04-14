# FixingHetDE

Here you'll find the code associated with the paper:

Ament-Velásquez et al. (2025) Reconstructing NOD-like receptor alleles with high internal conservation in *Podospora anserina* using long-read sequencing (under review)

----

Here I have three pipelines in different folders:

- HNWDphylogeny: Small pipeline (`HNWDphylogeny.smk`) to produce maximum likelihood phylogenies of different parts of an NWD protein alignment. (The actual figure in the paper was put together in Inkscape).
- JustHetGenes: contains the pipeline `HETwd40explorer.smk`, more or less as in the [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.13.632504v1) version. Here we compare assemblies, alleles, and WD40 repeats of *het-d*, *het-e*, and *het-r*. Most figures in the paper come from here.
- NWDgenes: contains the pipeline `NWD40explorer.smk`, which extends the pipeline above to all NWD genes with HIC WD40 repeats. Here we produced the LOGO and PCA analysis.


They are all [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines that more or less depend on [conda](https://docs.conda.io/en/latest/) environments for plotting.

----

Disclaimer: These scripts and files are provided "as is" and without any express or implied warranties, including, without limitation, the implied warranties of merchantability and fitness for a particular purpose.

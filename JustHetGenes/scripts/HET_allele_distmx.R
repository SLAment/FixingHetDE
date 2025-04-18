#!/usr/bin/env Rscript

### HET_allele_distmx: assigning colors to the WD40 HIC repeats based on physicochemical characteristics
#############################################################################
# ===========================================================================
# Brendan Furneaux
# with minor additions by Sandra Lorena Ament Velasquez
# 2025/01/09-10
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(Biostrings)
library(ggplot2)
library(ggforce)
library(reshape2)
library(dplyr)

# ============================
# Outputs
# ============================
# Dissimilarity heatmaps
hetd_disheatmap <- "results/het-d_disheatmap.png"
hete_disheatmap <- "results/het-e_disheatmap.png"
hetr_disheatmap <- "results/het-r_disheatmap.png"

# To make an explanatory figure of the method
LABhete3D <- "results/LABhete3D.pdf"
uth_heatmap <- "results/UTH.pdf"

# ============================
#### functions ####
# ============================
# add up the pairwise scores (or distances) between two amino acid sequences
# matrix should have row names and column names corresponding to the amino acids
aa_score <- function(aa1, aa2, matrix){
  sum(
    matrix[cbind(as.matrix(aa1), as.matrix(aa2))]
  )
}

# generate a distance/score matrix for a set of amino acid sequences
aa_distmatrix <- function(aa_seq, matrix){
  n <- length(aa_seq)
  distmatrix <- matrix(0, n, n, dimnames = list(names(aa_seq), names(aa_seq)))
  for (i in 1:n){
    for (j in i:n){
      s <- aa_score(aa_seq[[i]], aa_seq[[j]], matrix)
      distmatrix[i, j] <- s
      distmatrix[j, i] <- s
    }
  }
  return(distmatrix)
}

# Fit a LAB color palette to a distance matrix using simulated annealing
#
# inspired by gecos
# https://github.com/biotite-dev/gecos
lab_nmds <- function(
    dist_matrix,
    l_range = c(0, 100), # allowed range for L (luminance)
    a_range = c(-100, 100), # allowed range for A (green-red)
    b_range = c(-100, 100), # allowed range for B (blue-yellow)
    contrast = 100, # higher value gives higher weight to high-contrast palettes
    n_iter = 1000, # number of iterations
    beta = c(1, 500), # inverse temperature schedule
    stepsize = c(10, 0.2), # step size schedule
    seed = NULL, # random seed
    verbose = TRUE
) {
  # find identical sequences
  dupes <- which(dist_matrix == 0, arr.ind = TRUE)
  # for each pair, only keep the one where the smaller index is in col 1
  dupes <- dupes[dupes[,1] <= dupes[,2],]
  # for each sequence, only keep the first time it appears in col 2
  # this retains exactly one row for each input, which is in col 2
  # then the lowest-numbered identical sequence is in col 1
  # (if the sequence is its own least-numbered identical sequence then
  # col 1 == col 2)
  dupes <- dupes[!duplicated(dupes[,2]),]
  
  # figure out how much the rows will move after we remove the duplicates
  # (i.e. those where col 1 != col 2)
  shift <- cumsum(dupes[,1] != dupes[,2])
  
  # remove the duplicates
  matches <-dupes[dupes[,1] != dupes[,2],]
  dedup_matrix <- if (nrow(matches) > 0) {
    dist_matrix[-matches[,2], -matches[,2]]
  } else {
    dist_matrix
  }
  
  # renumber col 1 to match entries in dedup_matrix instead of dist_matri
  dupes[,1] <- dupes[,1] - shift[dupes[,1]]
  
  n <- nrow(dedup_matrix)
  
  # the target is for the distance matrix in LAB space to match the
  # given matrix, up to scaling.
  # The dist object contains only the lower triangle
  target_dist <- as.dist(dedup_matrix)
  target_dist <- target_dist / mean(target_dist)
  
  # initialize the LAB palette
  if (!is.null(seed)) {
    set.seed(seed)
  }
  palette <- colorspace::LAB(
    L = runif(n, l_range[1], l_range[2]),
    A = runif(n, a_range[1], a_range[2]),
    B = runif(n, b_range[1], b_range[2])
  )
  # re-sample the coordinates to be representable in the RGB gamut
  hexcodes <- colorspace::hex(palette)
  while (any(is.na(hexcodes))) {
    n_fix <- sum(is.na(hexcodes))
    palette@coords[is.na(hexcodes),] <-
      matrix(
        c(runif(n_fix, l_range[1], l_range[2]),
          runif(n_fix, a_range[1], a_range[2]),
          runif(n_fix, b_range[1], b_range[2])),
        nrow = n_fix
      )
    hexcodes <- colorspace::hex(palette)
  }
  
  # loss function from Gecos
  # Similar but not identical to the standard NMDS loss function
  score <- function(palette) {
    palette_dist <- dist(palette@coords)
    palette_contrast <- mean(palette_dist)
    palette_dist <- palette_dist / palette_contrast
    return(sum((palette_dist - target_dist)^2) + contrast / palette_contrast)
  }
  
  current_score <- score(palette)
  best_score <- current_score
  best_palette <- palette
  
  # exponential step and temperature schedules
  step_schedule <- stepsize[1] * (stepsize[2] / stepsize[1]) ^ ((seq_len(n_iter) - 1) / (n_iter - 1))
  beta_schedule <- beta[1] * (beta[2] / beta[1]) ^ ((seq_len(n_iter) - 1) / (n_iter - 1))
  
  # optimize the palette
  for (i in seq_len(n_iter)) {
    proposal <- palette
    proposal@coords <- proposal@coords + rnorm(n * 3, 0, step_schedule[i])
    
    # force all proposals into the requested ranges
    proposal@coords[,1] <- pmax(pmin(proposal@coords[,1], l_range[2]), l_range[1])
    proposal@coords[,2] <- pmax(pmin(proposal@coords[,2], a_range[2]), a_range[1])
    proposal@coords[,3] <- pmax(pmin(proposal@coords[,3], b_range[2]), b_range[1])
    
    # re-sample points which are outside the RGB gamut
    proposal_hexcodes <- colorspace::hex(proposal)
    while (any(is.na(proposal_hexcodes))) {
      n_fix <- sum(is.na(proposal_hexcodes))
      which_fix <- which(is.na(proposal_hexcodes))
      proposal@coords[which_fix,] <-
        palette@coords[which_fix,] + rnorm(n_fix * 3, 0, step_schedule[i])
      proposal@coords[which_fix,1] <- pmax(pmin(proposal@coords[which_fix,1], l_range[2]), l_range[1])
      proposal@coords[which_fix,2] <- pmax(pmin(proposal@coords[which_fix,2], a_range[2]), a_range[1])
      proposal@coords[which_fix,3] <- pmax(pmin(proposal@coords[which_fix,3], b_range[2]), b_range[1])
      proposal_hexcodes <- colorspace::hex(proposal)
    }
    proposal_score <- score(proposal)
    if (proposal_score < current_score) {
      # accept if the score improves
      palette <- proposal
      current_score <- proposal_score
      if (verbose) cat("Iteration", i, "Score", score(palette), "\n")
    } else {
      # if the score doesn't improve, accept with probability dependent on the
      # amount of increase and the current temperature
      if (runif(1) < exp((current_score - proposal_score) * beta_schedule[i])) {
        palette <- proposal
        current_score <- proposal_score
      }
    }
    
    # keep track of the best scoring palette, it's possible that it isn't the
    # last one sampled
    if (current_score < best_score) {
      best_palette <- palette
      best_score <- current_score
    }
  }
  
  # add the identical sequences back in
  hexcodes <- colorspace::hex(best_palette)[dupes[,1]]
  names(hexcodes) <- rownames(dist_matrix)
  attr(hexcodes, "score") <- best_score
  attr(hexcodes, "lab") <- best_palette
  
  return(hexcodes)
}

# run lab_nmds several times with different seeds, choosing the result with the
# best score
meta_lab_nmds <- function(dist_matrix, seeds = sample(1:1e6, 10), ...) {
  best_palette <- NULL
  best_score <- Inf
  for (seed in seeds) {
    palette <- lab_nmds(dist_matrix, seed = seed, ...)
    if (attr(palette, "score") < best_score) {
      best_palette <- palette
      best_score <- attr(palette, "score")
    }
  }
  return(best_palette)
}

# --- Lore: heatmap of all repeats ---
repeatheatmap <- function(distmatrix, ang = 0, hj = 0.5){
  # Convert matrix to a data frame in long format
  dism_long <- melt(as.matrix(distmatrix))
  
  # Check if labels have numbers
  contains_numbers <- any(grepl("[0-9]", dism_long$Var1)) || any(grepl("[0-9]", dism_long$Var1))

  if (contains_numbers) {
    # Case 1: Labels like 'e1', 'e2' -> strip prefix and order numerically
    names(dism_long) <- c("R1o", "R2o", "Dissimilarity")
    
    dism_long <- dism_long %>%
      mutate(
        R1 = gsub("[^0-9]", "", R1o),
        R2 = gsub("[^0-9]", "", R2o))
  
    # Make sure the factor order follows the numeric order (R1 has the same levels as R2)
    r1_levels <- unique(dism_long$R1) 
    
    dism_long <- dism_long %>%
      mutate(
        R1 = factor(R1, levels = r1_levels),
        R2 = factor(R2, levels = r1_levels)) %>%
      select(-R1o, -R2o)
  } else {
    # Case 2: Labels like 'M', 'A', etc. 
    names(dism_long) <- c("R1", "R2", "Dissimilarity")
  }
  
  # Plot heatmap
  p <- ggplot(dism_long, aes(x = R1, y = R2, fill = Dissimilarity)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +
    labs(x = "", y = "", fill = "Dissimilarity") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = ang, hjust = hj),
      axis.title = element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
}

# ============================
#### -------- Work it! --------
# ============================
# Amino acid functional distance matrix from Urbina Tang and Higgs
# https://doi.org/10.1007/s00239-005-0051-1
uth <- read.delim("data/UTH_matrix.txt", header=TRUE, row.names=1, sep = " ") |>
  as.matrix()

# nucleotide sequences for all the HIC repeats
repeat_seq <- Biostrings::readDNAStringSet("alignments/WDrepeats_10-11-12-14-30-32-39.fa")
# translate to amino acids
repeat_aa <- Biostrings::translate(repeat_seq)
# load metadata
repeat_metadata <- read.delim(
  "alignments/WDrepeats_10-11-12-14-30-32-39_metadata.txt",
  header=TRUE,
  row.names=1,
  sep = "\t"
)

#### Palette for all HIC repeats ####
# except those with stop codon

stops <- which(Biostrings::letterFrequency(repeat_aa, "*") > 0)
nostop_aa <- repeat_aa[-stops]
nostop_metadata <- repeat_metadata[names(nostop_aa),]

nostop_distmatrix <- aa_distmatrix(nostop_aa, uth)
set.seed(1)
nostop_palette <- meta_lab_nmds(
  dist_matrix = nostop_distmatrix,
  n_iter = 1e4,
  contrast = 100,
  l_range = c(20, 100)
)
nostop_metadata$fill <- nostop_palette

##### show het-d #####

hetd_strains <- c(
  ## reactive
  "ChEhDap", "CsDfp", "PaYp",
  ## non.reactive
  "PaZp", "CoEfp", "CmEmm", "Podan2",
  ## unknown
  "PaWa53m", "PaWa137m", "PaTgp", "PaWa87p", "PaWa100p", "PaWa58m", "PaWa21m", "PaWa46p", "PaWa28m", "PaWa63p"
)

nostop_metadata |>
  dplyr::filter(Strain %in% hetd_strains, Gene == "het-d") |>
  dplyr::mutate(Strain = ordered(Strain, levels = rev(hetd_strains))) |>
  ggplot(aes(x = Position, y = Strain, label = substr(RepeatType, 2, 10), fill = fill)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(color = "black")  +
  scale_fill_identity() +
  coord_fixed() +
  theme_minimal() +
  scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

##### show het-e #####

hete_strains <- c(
  # reactive
  "CmEmm", "ChEhDap", "PaYp", "PaA", "CoEcp", "PaZp", "CoEfp",
  # non.reactive
  "PaSp", "CaDam",
  #unknown
  "PaWa46p", "PaWa58m", "PaWa28m", "PaWa87+", "PaTgp", "PaWa63p", "PaWa137m",
  "PaWa100p", "PaWa21m", "PaWa53m"
)

nostop_metadata |>
  dplyr::filter(Strain %in% hete_strains, Gene == "het-e") |>
  dplyr::mutate(Strain = ordered(Strain, levels = rev(hete_strains))) |>
  ggplot(aes(x = Position, y = Strain, label = substr(RepeatType, 2, 10), fill = fill)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(color = "black")  +
  scale_fill_identity() +
  coord_fixed() +
  theme_minimal() +
  scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

##### show het-r #####

nostop_metadata |>
  dplyr::filter(Gene == "het-r") |>
  ggplot(aes(x = Position, y = Strain, label = substr(RepeatType, 2, 10), fill = fill)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(color = "black")  +
  scale_fill_identity() +
  coord_fixed() +
  theme_minimal() +
  scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

#### using only "variable" positions ####

variable_pos <- c(10, 11, 12, 14, 30, 32, 39)

variable_aa <- Biostrings::AAMultipleAlignment(nostop_aa)

colmask(variable_aa, invert = TRUE) <-
  IRanges::IRanges(start = variable_pos, width = 1L) |>
  IRanges::reduce()

variable_aa <- methods::as(variable_aa, "AAStringSet") |>
  unique()

names(variable_aa) <- sub("_.+", "", names(variable_aa))

variable_distmatrix <- aa_distmatrix(variable_aa, uth)


#### het-d ####

hetd_variants <- which(startsWith(names(variable_aa), "d"))
hetd_variable_aa <- variable_aa[hetd_variants]
hetd_variable_distmatrix <- variable_distmatrix[hetd_variants, hetd_variants]

set.seed(2)
hetd_palette <- meta_lab_nmds(
  dist_matrix = hetd_variable_distmatrix,
  n_iter = 1e4,
  contrast = 100,
  l_range = c(40, 80),
  a_range = c(20, 100),
  b_range = c(-50, 50)
)

nostop_metadata |>
  dplyr::filter(Strain %in% hetd_strains, Gene == "het-d") |>
  dplyr::mutate(Strain = ordered(Strain, levels = rev(hetd_strains))) |>
  ggplot(aes(x = Position, y = Strain, label = substr(RepeatType, 2, 10), fill = RepeatType)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(color = "black")  +
  scale_fill_manual(values = hetd_palette) +
  coord_fixed() +
  theme_minimal() +
  scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

# heatmap of all het-d repeats
ggsave(plot = repeatheatmap(hetd_variable_distmatrix),
       filename = hetd_disheatmap,
       width = 4.2, height = 5)

#### het-e ####

hete_variants <- which(startsWith(names(variable_aa), "e"))
hete_variable_aa <- variable_aa[hete_variants]
hete_variable_distmatrix <- variable_distmatrix[hete_variants, hete_variants]

set.seed(3)
hete_palette <- meta_lab_nmds(
  dist_matrix = hete_variable_distmatrix,
  n_iter = 1e4,
  contrast = 100,
  l_range = c(40, 90)
)

nostop_metadata |>
  dplyr::filter(Strain %in% hete_strains, Gene == "het-e") |>
  dplyr::mutate(
    Position = ifelse(Strain == "PaA", Position + 1, Position),
    Strain = ordered(Strain, levels = rev(hete_strains))
  ) |>
  ggplot(aes(x = Position, y = Strain, label = substr(RepeatType, 2, 10), fill = RepeatType)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(color = "black")  +
  scale_fill_manual(values = hete_palette) +
  coord_fixed() +
  theme_minimal() +
  scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

# heatmap of all het-e repeats
ggsave(plot = repeatheatmap(hete_variable_distmatrix),
       filename = hete_disheatmap,
       width = 6, height = 6)

##### color-space plot #####
# manual 3D plot

# functions to project a 3-dimensional point to projected x- and y-coordinates
# in a 2D plot using basis vectors
rotate_x <- function(x, y, z) {
  x * ix + y *iy + z * iz
}

rotate_y <- function(x, y, z) {
  x * jx + y *jy + z * jz
}

ix = 1
iy = -1/sqrt(2)
iz = 0

jx = 0
jy = -1/sqrt(2)
jz = 1

# Plot the points
(heteLAB <- tibble::enframe(hete_palette) |>
  dplyr::bind_cols(attr(hete_palette, "lab")@coords) |>
  dplyr::arrange(B) |>
  ggplot(aes(
    x = rotate_x(A, B, L-40), ## offset L by 40 to keep it compact
    y = rotate_y(A, B, L-40),
    color = name,
    fill = name,
    label = substr(name, 2, 5)
  )) +
    # draw the axes
  annotate(
    "segment",
    x = rep(rotate_x(-40, -30, 0), 3),
    y = rep(rotate_y(-40, -30, 0), 3),
    xend = c(
      rotate_x(40, -30, 0),
      rotate_x(-40, 45, 0),
      rotate_x(-40, -30, 60)
    ),
    yend = c(
      rotate_y(40, -30, 0),
      rotate_y(-40, 45, 0),
      rotate_y(-40, -30, 40)
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.8
  ) +
    # draw axis labels
  annotate(
    "text",
    x = c(
      rotate_x(40, -30, 3),
      rotate_x(-40, 45, 3),
      rotate_x(-40, -30, 42)
    ),
    y = c(
      rotate_y(40, -30, 3),
      rotate_y(-40, 45, 3),
      rotate_y(-40, -30, 42)
    ),
    label = c("A", "B", "L"),
    size = 7,
    color = "black"
  ) +
    # draw foreshortened "shadow" points in AB-plane
  geom_ellipse(
    aes(x0 = rotate_x(A, B, 0), y0 = rotate_y(A, B, 0), a = 1, b = 1/2,
        angle = 0)) +
    # draw linking lines
  geom_segment(
    aes(xend = rotate_x(A, B, 0), yend = rotate_y(A, B, 0)),
    linetype = "dashed"
  ) +
    # draw selected HET repeats
  geom_point(
    size = 10,
    data = ~dplyr::filter(., name %in% c("e1", "e17", "e24")),
    color = "black",
    shape = 22,
    stroke = 1
  ) +
    # draw the rest of the HET repeats
  geom_point(
    size = 7,
    shape = 15,
    data = ~dplyr::filter(., !name %in% c("e1", "e17", "e24"))
  ) +
    # draw repeat labels
  geom_text(
    aes(size = ifelse(name %in% c("e1", "e17", "e24"), 5, 3.5)),
    color = "black") +
    # theme
  scale_size_identity(guide = "none") +
  scale_color_manual(values = hete_palette, guide = NULL,
                     aesthetics = c("color", "fill")) +
  theme_void() +
  coord_fixed() )

# Save the 3D plot for the schematic figure
ggsave(plot = heteLAB,
       filename = LABhete3D,
       width = 5, height = 5)

#### het-r ####

hetr_variants <- which(startsWith(names(variable_aa), "r"))
hetr_variable_aa <- variable_aa[hetr_variants]
hetr_variable_distmatrix <- variable_distmatrix[hetr_variants, hetr_variants]

set.seed(1)
hetr_palette <- meta_lab_nmds(
  dist_matrix = hetr_variable_distmatrix,
  n_iter = 1e4,
  contrast = 100,
  l_range = c(60, 90),
  a_range = c(-50, 50),
  b_range = c(0, 100)
)

nostop_metadata |>
  dplyr::filter(Gene == "het-r") |>
  ggplot(aes(x = Position, y = Strain, label = substr(RepeatType, 2, 10), fill = RepeatType)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(color = "black")  +
  scale_fill_manual(values = hetr_palette) +
  coord_fixed() +
  theme_minimal() +
  scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

# heatmap of all het-r repeats
ggsave(plot = repeatheatmap(hetr_variable_distmatrix),
       filename = hetr_disheatmap,
       width = 4.2, height = 5)


#### translate palettes ####

sed_palette_swap <- function(old_palette, new_palette, file) {
  old_palette <- tibble::enframe(old_palette, value = "old") |>
    dplyr::mutate(
      name = readr::parse_number(name),
      old = tolower(old)
    )
  new_palette <- tibble::enframe(new_palette, value = "new") |>
    dplyr::mutate(
      name = readr::parse_number(name),
      new = tolower(new)
    )
  
  dplyr::inner_join(old_palette, new_palette, by = "name") |>
    glue::glue_data("s/fill:{old}/fill:{new}/") |>
    writeLines(file)
}

hetR_mod_palette = c(
  "1" = "#66c2a5",
  "2" = "#cf9c76",
  "3" = "#cf948c",
  "4" = "#969dca",
  "5" = "#d58ec4",
  "6" = "#c6b18b",
  "7" = "#b7d84c",
  "8" = "#f6d832",
  "9" = "#efcc6b",
  "10" = "#d5be9d",
  "11" = "#b3b3b3",
  "16" = "#cccccc"
)

sed_swap_file <- file("figures/colorswap.sed", "w")

# Named palettes defined in PlottingHETAlleles.R

sed_palette_swap(hetD_named_palette, hetd_palette, sed_swap_file)
sed_palette_swap(hetE_named_palette, hete_palette, sed_swap_file)
sed_palette_swap(hetR_named_palette, hetr_palette, sed_swap_file)
sed_palette_swap(hetR_mod_palette, hetr_palette, sed_swap_file)

close(sed_swap_file)

# To modify existing SVG files:
# sed -f figures/colorswap.sed {old_file.svg} > {new_file.svg}

# To remove opacity:
# sed -ri 's/fill-opacity:[0-9.]+;/fill-opacity:1;/g' {new_file.svg}

# --- Lore: heatmap of the amino acid matrix for illustration ---
# heatmap of all het-e repeats
ggsave(plot = repeatheatmap(uth), 
       filename = uth_heatmap, 
       width = 3.8, height = 4)


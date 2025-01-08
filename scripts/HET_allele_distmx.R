library(Biostrings)
library(ggplot2)

# add up the pairwise scores (or distances) between two amino acid sequences
# matrix should have row names and column names corresponding to the amino acids
aa_score <- function(aa1, aa2, matrix){
  score <- 0
  for (i in 1:nchar(aa1)){
    score <- score + matrix[letter(aa1, i), letter(aa2, i)]
  }
  return(score)
}

# generate a distance/score matrix for a set of amino acid sequences
aa_distmatrix <- function(aa_seq, matrix){
  n <- length(aa_seq)
  distmatrix <- matrix(0, n, n, dimnames = list(names(aa_seq), names(aa_seq)))
  for (i in 1:n){
    for (j in 1:n){
      distmatrix[i, j] <- aa_score(aa_seq[[i]], aa_seq[[j]], matrix)
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
  decode <- which(dist_matrix == 0, arr.ind = TRUE)
  decode <- decode[decode[,1] <= decode[,2],]
  decode <- decode[!duplicated(decode[,2]),]

  # distance matrix not including any duplicates
  shift <- cumsum(decode[,1] != decode[,2])
  matches <-decode[decode[,1] != decode[,2],]
  dedup_matrix <- dist_matrix[-matches[,2], -matches[,2]]
  decode[,1] <- decode[,1] - shift[decode[,1]]

  n <- nrow(dedup_matrix)
  target_dist <- as.dist(dedup_matrix)
  target_dist <- target_dist / mean(target_dist)

  #initialize the LAB palette
  if (!is.null(seed)) {
    set.seed(seed)
  }
  palette <- colorspace::LAB(
    L = runif(n, l_range[1], l_range[2]),
    A = runif(n, a_range[1], a_range[2]),
    B = runif(n, b_range[1], b_range[2])
  )
  # force the coordinates to be representable in the RGB gamut
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

  score <- function(palette) {
    palette_dist <- dist(palette@coords)
    palette_contrast <- mean(palette_dist)
    palette_dist <- palette_dist / palette_contrast
    return(sum((palette_dist - target_dist)^2) + contrast / palette_contrast)
  }

  current_score <- score(palette)
  best_score <- current_score
  best_palette <- palette

  step_schedule <- stepsize[1] * (stepsize[2] / stepsize[1]) ^ ((seq_len(n_iter) - 1) / (n_iter - 1))
  beta_schedule <- beta[1] * (beta[2] / beta[1]) ^ ((seq_len(n_iter) - 1) / (n_iter - 1))

  # optimize the palette
  for (i in seq_len(n_iter)) {
    proposal <- palette
    proposal@coords <- proposal@coords + rnorm(n * 3, 0, step_schedule[i])
    proposal@coords[,1] <- pmax(pmin(proposal@coords[,1], l_range[2]), l_range[1])
    proposal@coords[,2] <- pmax(pmin(proposal@coords[,2], a_range[2]), a_range[1])
    proposal@coords[,3] <- pmax(pmin(proposal@coords[,3], b_range[2]), b_range[1])
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
      palette <- proposal
      current_score <- proposal_score
      if (verbose) cat("Iteration", i, "Score", score(palette), "\n")
    } else {
      if (runif(1) < exp((current_score - proposal_score) * beta_schedule[i])) {
        palette <- proposal
        current_score <- proposal_score
      }
    }
    if (current_score < best_score) {
      best_palette <- palette
      best_score <- current_score
    }
  }

  # add the identical sequences back in
  hexcodes <- colorspace::hex(best_palette)[decode[,1]]
  names(hexcodes) <- rownames(dist_matrix)

  return(hexcodes)
}

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

# Choose only het-e

hete_strains <- c(
  # reactive
  "CmEmm", "ChEhDap", "PaYp", "PaA", "CoEcp", "PaZp", "CoEfp",
  # non.reactive
  "PaSp", "CaDam",
  #unknown
  "PaWa46p", "PaWa58m", "PaWa28m", "PaWa87+", "PaTgp", "PaWa63p", "PaWa137m",
  "PaWa100p", "PaWa21m", "PaWa53m"
  )

hete_metadata <- repeat_metadata[repeat_metadata$Strain %in% hete_strains & repeat_metadata$Gene == "het-e",]

# don't include the illumina ones
hete_metadata <- hete_metadata[!grepl("default|allkmer|L28125", row.names(hete_metadata)),]

hete_aa <- repeat_aa[row.names(hete_metadata),]

# calculate distance matrix
hete_distmatrix <- aa_distmatrix(hete_aa, uth)

# plot the distance matrix
# hete_distmatrix |>
#   tibble::as_tibble(rownames = "rep1") |>
#   tidyr::pivot_longer(-rep1, names_to = "rep2", values_to = "dist") |>
#   ggplot(aes(x = rep1, y = rep2, fill = dist)) +
#   geom_tile()

# embed the distance matrix in 3-dimensional space
rand.seed <- 2

# using standard NMDS
# here the problem is that there is no guarantee that all the points will be
# possible to convert to RGB, so everything needs to be "squished" into a
# rectangular region that works.  This limits the possible contrast

# hete_nmds <- vegan::metaMDS(hete_distmatrix, distance = "euclidean", k = 3)
# # scale to fit into CIE-LAB color space
# hete_nmds$points[,1] <- (hete_nmds$points[,1] - min(hete_nmds$points[,1])) / (max(hete_nmds$points[,1]) - min(hete_nmds$points[,1])) * 35 + 30
# hete_nmds$points[,2] <- (hete_nmds$points[,2] - min(hete_nmds$points[,2])) / (max(hete_nmds$points[,2]) - min(hete_nmds$points[,2])) * 70 - 20
# hete_nmds$points[,3] <- (hete_nmds$points[,3] - min(hete_nmds$points[,3])) / (max(hete_nmds$points[,3]) - min(hete_nmds$points[,3])) * 70 - 20
# hete_metadata$fill <- colorspace::hex(colorspace::LAB(hete_nmds$points))

# using hand-coded LAB NMDS
# this version is a bit slower, but it enforces that the colors fit nicely in LAB
# during every iteration, so it can more completely utilize the color space.
hete_metadata$fill <- lab_nmds(
  dist_matrix = hete_distmatrix,
  seed = rand.seed,
  n_iter = 1e4,
  contrast = 100,
  l_range = c(40, 80)
)

# plot
hete_metadata |>
  dplyr::mutate(
    Position = Position + ifelse(Strain == "PaA", 1, 0),
    Strain = ordered(Strain, levels = rev(hete_strains))
    ) |>
  ggplot(aes(x = Position, y = Strain, label = RepeatType, fill = fill)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text() +
  scale_fill_identity() +
  coord_fixed() +
  theme_minimal() +
  scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

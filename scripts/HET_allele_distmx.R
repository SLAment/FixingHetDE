library(Biostrings)
library(ggplot2)

#### functions ####

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

system.time(nostop_distmatrix <- aa_distmatrix(nostop_aa, uth))
set.seed(1)
nostop_palette <- meta_lab_nmds(
  dist_matrix = nostop_distmatrix,
  n_iter = 1e4,
  contrast = 100,
  l_range = c(40, 80)
)
nostop_metadata$fill <- nostop_palette

#### show het-d ####

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

#### show het-e ####

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

#### show het-r ####

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

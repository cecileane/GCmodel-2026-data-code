################################################################################
## Install packages if necessary

# remotes::install_github("lamho86/phylolm")
# requires v2.7.0 or higher

# remotes::install_github("larabreithaupt/seastaR")
# remotes::install_github("pbastide/seastaR", ref = "correct_sigma2_inference", force = TRUE)
# Use corrected version, see PR #7: https://github.com/larabreithaupt/seastaR/pull/7

################################################################################
## Load packages
library(ape)
library(phylolm)
library(seastaR)
library(tidyverse)
library(ggthemes)

library(here)
result_dir <- here("tomato_analysis", "results")
dir.create(result_dir, showWarnings=FALSE)

columnwidth <- 3 # in inches # 1 point / 72 = 1 inch
twocolumnwidth <- 6

################################################################################
## Time tree from
## https://github.com/mhibbins/genetreepruningalg/blob/e2627a711e088844f19c680262e05c433ac81a45/analyses/tomato_analysis/make_tomato_test_input_files.py#L13
tomato_tree <- read.tree(file = here("tomato_analysis", "tomato_tree.phy"))
plot(tomato_tree)

################################################################################
## Floral traits from Hibbins et al 2023
# tomato_traits <- read.csv("https://raw.githubusercontent.com/larabreithaupt/seastaR/refs/heads/main/analyses/flower_morphometrics2.csv")
tomato_traits <- read.csv(here("tomato_analysis", "flower_morphometrics2.csv"))

################################################################################
## Format data as in
## https://raw.githubusercontent.com/larabreithaupt/seastaR/refs/heads/main/analyses/tomato_seastar_analysis.R

# convert from years to coales units
tomato_tree[["edge.length"]] <- tomato_tree[["edge.length"]] / 400000

# subset phylo object into high and low ILS knot groups:
# High ILS knot: S. galapagense 0436, S. cheesmaniae 3124, S. pimpinellifolium 1269
# Low ILS knot: S. pennellii 3778, S. pennellii 0716, S. pimpinellifolium 1589

high_ils_species <- c("S.galapagense", "S.cheesmaniae", "S.pim1269")
low_ils_species <- c("S.pen3778", "S.pen0716", "S.pim1589")
pruned_tree_low<-drop.tip(tomato_tree,tomato_tree$tip.label[-match(low_ils_species, tomato_tree$tip.label)])
pruned_tree_high<-drop.tip(tomato_tree,tomato_tree$tip.label[-match(high_ils_species, tomato_tree$tip.label)])

# get species tree covar matrix
low_ils_C <- ape::vcv(pruned_tree_low)
high_ils_C <- ape::vcv(pruned_tree_high)

#format traits and get full var covar matrices
low_ils_Cstar <- get_full_matrix(pruned_tree_low)
high_ils_Cstar <- get_full_matrix(pruned_tree_high)

low_ils_accession <- c("LA1589", "LA3778", "LA0716")
high_ils_accession <- c("LA0436", "LA3124", "LA1269")

low_ils_filtered_traits <- dplyr::filter(tomato_traits, AccessionID %in% low_ils_accession)
high_ils_filtered_traits <- dplyr::filter(tomato_traits, AccessionID %in% high_ils_accession)

# calculate mean within each accession
low_ils_mean_acc <- low_ils_filtered_traits %>%
  dplyr::group_by(AccessionID) %>%
  summarize_at(vars(SE_PlantMean:SE_bin), list(name = mean))
high_ils_mean_acc <- high_ils_filtered_traits %>%
  dplyr::group_by(AccessionID) %>%
  summarize_at(vars(SE_PlantMean:SE_bin), list(name = mean))

# get specific traits
low_ils_traits <- low_ils_mean_acc[ ,c(1, 3:5)]
high_ils_traits <- high_ils_mean_acc[ ,c(1, 3:5)]

# reorder rows of trait df to match covar matrix
low_ils_traits <- low_ils_traits %>% dplyr::slice(match(low_ils_accession, AccessionID))
high_ils_traits <- high_ils_traits %>% dplyr::slice(match(high_ils_accession, AccessionID))

# estimate sigma^2 value

#LOW ILS
low_ils_CD_C_s2 <- sigma2_inference(low_ils_C, unlist(low_ils_traits[, 2]))
low_ils_AL_C_s2 <- sigma2_inference(low_ils_C, unlist(low_ils_traits[, 3]))
low_ils_SL_C_s2 <- sigma2_inference(low_ils_C, unlist(low_ils_traits[, 4]))

low_ils_CD_Cstar_s2 <- sigma2_inference(low_ils_Cstar, unlist(low_ils_traits[, 2]))
low_ils_AL_Cstar_s2 <- sigma2_inference(low_ils_Cstar, unlist(low_ils_traits[, 3]))
low_ils_SL_Cstar_s2 <- sigma2_inference(low_ils_Cstar, unlist(low_ils_traits[, 4]))

#HIGH ILS
high_ils_CD_C_s2 <- sigma2_inference(high_ils_C, unlist(high_ils_traits[, 2]))
high_ils_AL_C_s2 <- sigma2_inference(high_ils_C, unlist(high_ils_traits[, 3]))
high_ils_SL_C_s2 <- sigma2_inference(high_ils_C, unlist(high_ils_traits[, 4]))

high_ils_CD_Cstar_s2 <- sigma2_inference(high_ils_Cstar, unlist(high_ils_traits[, 2]))
high_ils_AL_Cstar_s2 <- sigma2_inference(high_ils_Cstar, unlist(high_ils_traits[, 3]))
high_ils_SL_Cstar_s2 <- sigma2_inference(high_ils_Cstar, unlist(high_ils_traits[, 4]))

################################################################################
## Format data to use phylolm

# reorder rows of trait df to match tree
low_ils_accession_species <- data.frame(accession = c("LA3778", "LA0716", "LA1589"),
                                        species = c("S.pen3778", "S.pen0716", "S.pim1589"))
high_ils_accession_species <- data.frame(accession = c("LA0436", "LA3124", "LA1269"),
                                         species = c("S.galapagense", "S.cheesmaniae", "S.pim1269"))

low_ils_traits <- low_ils_traits %>% dplyr::slice(match(low_ils_accession_species[match(pruned_tree_low$tip.label, low_ils_accession_species$species), ]$accession, AccessionID))
low_ils_traits <- data.frame(low_ils_traits)
rownames(low_ils_traits) <- low_ils_accession_species[match(low_ils_traits$AccessionID, low_ils_accession_species$accession), ]$species

high_ils_traits <- high_ils_traits %>% dplyr::slice(match(high_ils_accession_species[match(pruned_tree_high$tip.label, high_ils_accession_species$species), ]$accession, AccessionID))
high_ils_traits <- data.frame(high_ils_traits)
rownames(high_ils_traits) <- high_ils_accession_species[match(high_ils_traits$AccessionID, high_ils_accession_species$accession), ]$species

# check that phylolm and cstar give the same results
stopifnot(isTRUE(all.equal(high_ils_AL_C_s2,
                           phylolm(AL_PlantMean_name ~ 1, phy = pruned_tree_high, data = high_ils_traits, REML = TRUE)$sigma2,
                           check.attributes = FALSE, check.class = FALSE)))
stopifnot(isTRUE(all.equal(sigma2_inference(vcv(pruned_tree_high), high_ils_traits[["AL_PlantMean_name"]]),
                           phylolm(high_ils_traits[["AL_PlantMean_name"]] ~ 1, phy = pruned_tree_high, data = high_ils_traits, REML = TRUE)$sigma2,
                           check.attributes = FALSE, check.class = FALSE)))
stopifnot(isTRUE(all.equal(sigma2_inference(low_ils_C, low_ils_traits$CD_PlantMean_name),
                           phylolm(low_ils_traits[["CD_PlantMean_name"]] ~ 1, phy = pruned_tree_low, data = low_ils_traits, REML = TRUE)$sigma2,
                           check.attributes = FALSE, check.class = FALSE)))
stopifnot(isTRUE(all.equal(sigma2_inference(get_full_matrix(pruned_tree_high), high_ils_traits[["AL_PlantMean_name"]]),
                           high_ils_AL_Cstar_s2,
                           check.attributes = FALSE, check.class = FALSE)))

################################################################################
## Construct a table with results
estim_sigma2 <- function(trait, ILS, method) {
  if (ILS == "low ILS") {
    tree <- pruned_tree_low
    trait_data <- low_ils_traits
  } else {
    tree <- pruned_tree_high
    trait_data <- high_ils_traits
  }
  trait_name <- paste0(trait, "_PlantMean_name")
  if (method == "BM") {
    return(phylolm(trait_data[[trait_name]] ~ 1, phy = tree, data = trait_data, REML = TRUE, model = "BM")$sigma2)
  }
  if (method == "BMws") {
    return(phylolm(trait_data[[trait_name]] ~ 1, phy = tree, data = trait_data, REML = TRUE, model = "BM", measurement_error = TRUE)$sigma2)
  }
  if (method == "ILS") {
    return(phylolm(trait_data[[trait_name]] ~ 1, phy = tree, data = trait_data, REML = TRUE, model = "GC")$sigma2)
  }
  if (method == "ILS 0") {
    return(phylolm(trait_data[[trait_name]] ~ 1, phy = tree, data = trait_data, REML = TRUE, model = "GC",
                   starting.value = list(lambda_GC = 0), upper.bound = list(lambda_GC = 0), lower.bound = list(lambda_GC = 0))$sigma2)
  }
  if (method == "ILS 1") {
    return(phylolm(trait_data[[trait_name]] ~ 1, phy = tree, data = trait_data, REML = TRUE, model = "GC",
                   starting.value = list(lambda_GC = 1), upper.bound = list(lambda_GC = 1), lower.bound = list(lambda_GC = 1))$sigma2)
  }
  if (method == "seastaR") {
    return(sigma2_inference(get_full_matrix(tree), trait_data[[trait_name]])[1, 1])
  }
}

fit_sigma2 <- NULL
for (trait in c("CD", "AL", "SL")) {
  for (ILS in c("low ILS", "high ILS")) {
    for (method in c("BM", "seastaR", "ILS 1", "ILS 0")) {
      tmp <- data.frame(trait = trait,
                        ILS = ILS,
                        method = method,
                        sigma2 = estim_sigma2(trait, ILS, method))
      fit_sigma2 <- rbind(fit_sigma2, tmp)
    }
  }
}

fit_sigma2$method <- factor(fit_sigma2$method, levels = c("seastaR", "ILS 1", "ILS 0", "BM"))
levels(fit_sigma2$method) <- c("seastaR", "GC \u03bb=1", "GC \u03bb=0", "BM")
fit_sigma2$trait <- factor(fit_sigma2$trait, levels = c("CD", "AL", "SL"))
levels(fit_sigma2$trait) <- c("corolla diameter", "anther length", "stigma length")
p <- ggplot(fit_sigma2, aes(x=ILS, y=sigma2, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "triplet", y = "estimated σ²", fill = "method") +
  # facet_grid(rows = vars(trait), scales = "free_y") +
  scale_y_sqrt() +
  facet_wrap(vars(trait), scales = "free") +
  theme_bw() +
  scale_fill_viridis_d(option = "B", direction = -1) +
  theme(text=element_text(size=unit(9, 'pt')))
  # scale_fill_manual(
  #   values = c(seastaR = "#F0E442", 'ILS 1' = "#E69F00", 'ILS 0' = "#D55E00", BM = "#000000")
  # )
p

ggsave(filename = file.path(result_dir, "tomato_high_low.pdf"),
       plot = p,
       width = 1.6*columnwidth,
       height = 0.8*columnwidth,
       unit = "in",
       device = grDevices::cairo_pdf)

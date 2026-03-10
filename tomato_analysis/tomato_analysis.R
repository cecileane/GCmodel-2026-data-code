################################################################################
## Install packages if necessary

# remotes::install_github("lamho86/phylolm")
# requires v2.7.0 or higher

################################################################################
## Load packages
library(ape)
library(phylolm)
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
# convert from years to coales units
tomato_tree[["edge.length"]] <- tomato_tree[["edge.length"]] / 400000

################################################################################
## Floral traits from Hibbins et al 2023
# tomato_traits <- read.csv("https://raw.githubusercontent.com/larabreithaupt/seastaR/refs/heads/main/analyses/flower_morphometrics2.csv")
tomato_traits <- read.csv(here("tomato_analysis", "flower_morphometrics2.csv"))

################################################################################
## Format data - alt - keep only one accession number per species

tomato_traits$sp_accession <- paste0(tomato_traits$Sp_ID, "_", tomato_traits$AccessionID)

## accession for each species
aasp <- sapply(unique(tomato_traits$Sp_ID), function(x) summary(factor(subset(tomato_traits, Sp_ID == x)$AccessionID)))
aasp
# when there are several accession per species, we chose the one with the largest number of replicates.
# this happens for the following species
aasp[sapply(aasp, function(x) length(x) > 1)]

# tip names vs traits
poptospecies = data.frame(
  population = c( # not included : "S.huaylasense", "S.tuberosum"
    "S.galapagense",  "S.cheesmaniae",  "S.pim1269",         "S.pim1589",
    "S.neorickii",    "S.arcanum",      "S.peruvianum",      "S.corneliomulleri",
    "S.chilense",     "S.habrochaites", "S.pen3778",         "S.pen0716"),
  species_id  = c( # not included : "SCER", "SCHM", "SLYC"
    "SGAL_LA0436",    "SCHE_LA3124",   "SPIM_LA1269",         "SPIM_LA1589",
    "SNEO_LA1321",    "SARC_LA0385",   "SPER_LA0153",         "SCOR_LA0107",
    "SCHL_LA1971",    "SHAB_LA1777",   "SPEN_LA3778",         "SPEN_LA0716")
)

# trim tree to data
tomato_tree_trim <- keep.tip(tomato_tree, poptospecies$population)
tomato_tree_trim$tip.label[match(poptospecies$population, tomato_tree_trim$tip.label)] <- poptospecies$species_id
plot(tomato_tree)
plot(tomato_tree_trim)

# trim data to tree
tomato_traits_trim <- tomato_traits[tomato_traits$sp_accession %in% poptospecies$species_id, ]
tomato_traits_trim <- tomato_traits_trim[do.call(c, sapply(tomato_tree_trim$tip.label, function(x) grep(x, tomato_traits_trim$sp_accession))), ]
tomato_traits_trim$sample_id <- make.unique(tomato_traits_trim$sp_accession)
rownames(tomato_traits_trim) <- tomato_traits_trim$sample_id
# add replicates
tomato_tree_trim_rep <- add.individuals(tomato_tree_trim, tomato_traits_trim, species_id = "sp_accession", sample_id = "sample_id")
plot(tomato_tree_trim_rep)
summary(factor(tomato_traits_trim$sp_accession))

################################################################################
## Plot data

low_ils_accession <- c("LA1589", "LA3778", "LA0716")
high_ils_accession <- c("LA0436", "LA3124", "LA1269")

okabeito <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(PhylogeneticEM)
tree_plot <- tomato_tree_trim_rep
ntaxa <- length(tree_plot$tip.label)
tip_names <- tomato_tree_trim_rep$tip.label
tip_species <- cumsum(!duplicated(sub("\\.[1-9]$", "", tip_names)))
# tip colors
tip_cols <- rep(c("grey30", "grey60"), 6)
tip_cols <- rep(tip_cols, times = table(tip_species))
tip_cols[do.call(c, lapply(low_ils_accession, function(x) grep(x, tree_plot$tip.label)))] <- okabeito[2]
tip_cols[do.call(c, lapply(high_ils_accession, function(x) grep(x, tree_plot$tip.label)))] <- okabeito[3]
# trait names
plot_matrix <- as.matrix(tomato_traits_trim[, c("CD_PlantMean", "AL_PlantMean", "SL_PlantMean")])
colnames(plot_matrix) <- c("CD", "AL", "SL")
# tip labels
tree_plot$tip.label[-c(2, 5, 9, 13, 16, 20, 24, 27, 30, 33, 36, 39)] <- ""
tree_plot$tip.label <- sub("\\.[1-9]$", "", tree_plot$tip.label)
# remove trait intercept
trait_intercept <- floor(apply(plot_matrix, 2, min))
trait_max <- floor(apply(plot_matrix, 2, max))
plot_matrix <- plot_matrix - rep(1, nrow(plot_matrix)) %*% t(trait_intercept)

pdf(file = file.path(result_dir, "tomato_tree_plot.pdf"),
    width = 1.25*twocolumnwidth,
    height = 1.25*1.5*columnwidth)
par(mar = c(0, 0.1, 0, 0))
# plot
plot(params_BM(p = 3), tree_plot, t(plot_matrix),
     edge.width = 3, show.tip.label = TRUE,
     automatic_colors = FALSE, color_characters = tip_cols, label_cex = 0.85,
     show_axis_traits = FALSE, axis_cex = 0.85)
# tip labels
tiplabels(pch = 19, col = tip_cols)
# trait axis
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
pos_last_tip <- max(lastPP$xx)
size_labels <- max(strwidth(tree_plot$tip.label, cex = 0.85))
h_p <- max(ape::node.depth.edgelength(tree_plot))
x.lim.max <- h_p + 3 * h_p / 3 + h_p / 4
available_x <- x.lim.max - pos_last_tip - size_labels
ell <- available_x / (3 +  (3 + 1) / 3)# length for the plot of one character
offset <- ell / 3
for (t in 1:3){
  eccart_g <- offset
  mult <- ell / (apply(plot_matrix, 2, max)[t] - apply(plot_matrix, 2, min)[1])
  axis(1, at = pos_last_tip + eccart_g + c(0, mult*(trait_max - trait_intercept)[t]),
       labels = c(trait_intercept[t], trait_max[t]),
       pos = -ntaxa/10 + ntaxa/15,
       cex.axis = 0.85,
       padj = -0.75,
       las = 0)
  pos_last_tip <- pos_last_tip + ell + offset
}
text(max(lastPP$xx) + available_x / 2, -4, "mm", pos = 1, cex = 0.85)
# axis time
axisPhylo(pos = -ntaxa/10 + ntaxa/15, cex.axis = 0.85, padj = -0.75)
text(3, -4, "coalescent units", pos = 1, cex = 0.85)
# legend low/high ils
legend(0, 5,
       c("low ILS triplet", "high ILS triplet"),
       col = okabeito[2:3], pch = 19,
       box.lwd = 0, border = "white", ncol = 2, cex = 0.85, pt.cex = 1)
dev.off()


################################################################################
## Function to get within-population variation expected from the GC model
sigma2_within_gc <- function(fit) {
  if (fit$model != "GC") return(0.0)
  return(fit$sigma2 * (1 + (fit$optpar - 1) * exp(-fit$mean.tip.height)))
}


################################################################################
## Construct a table with results
traitnames = c("CD_PlantMean","AL_PlantMean","SL_PlantMean")

res <- data.frame(
  trait = rep(traitnames, each = 4 * 2),
  method = rep(c("GC 1", "GC 0", "GC", "BM", "GC 1 wsp", "GC 0 wsp", "GC wsp", "BM wsp"), 3),
  methodMod = rep(c("GC", "GC", "GC", "BM"), 3 * 2),
  lambda = rep(c(1, 0, NA, NA), 3 * 2),
  measurement_error = rep(rep(c(FALSE, TRUE), each = 4), 3),
  sigma2 = rep(NA, 3*4*2),
  s2_wsp = rep(NA, 3*4*2),
  aic = rep(NA, 3*4*2)
)

# remove BM without measurement error (no sense : all samples are the same)
res <- res[!(res$method == "BM" & res$measurement_error == FALSE), ]

for (r in seq_len(nrow(res))) {
  if (!is.na(res[r, "lambda"])) {
    fit <- phylolm(tomato_traits_trim[[res[r, "trait"]]] ~ 1,
                   phy = tomato_tree_trim_rep, data = tomato_traits_trim, REML = TRUE,
                   model = res[r, "methodMod"], measurement_error = res[r, "measurement_error"],
                   starting.value = list(lambda_GC = res[r, "lambda"]),
                   upper.bound = list(lambda_GC = res[r, "lambda"]),
                   lower.bound = list(lambda_GC = res[r, "lambda"]))
    corAIC <- 2 # phylolm does not know that lambda is fixed
  } else {
    fit <- phylolm(tomato_traits_trim[[res[r, "trait"]]] ~ 1,
                   phy = tomato_tree_trim_rep, data = tomato_traits_trim, REML = TRUE,
                   model = res[r, "methodMod"], measurement_error = res[r, "measurement_error"])
    corAIC <- 0
  }
  res[r, "sigma2"] <- fit$sigma2
  res[r, "s2_wsp"] <- fit$sigma2_error + sigma2_within_gc(fit)
  res[r, "aic"] <- fit$aic - corAIC
}

res$wAIC <- exp(-(res$aic - tapply(res$aic, res$trait, min)[res$trait])/2)
res$wAIC <- res$wAIC / tapply(res$wAIC, res$trait, sum)[res$trait]
res$isbestaic <- res$aic == tapply(res$aic, res$trait, min)[res$trait]
write.csv(res, file.path(result_dir, "phylolm_wholetree_withrep.csv"))


res$method <- factor(res$method, levels = c("GC 1", "GC 0", "GC", "GC 1 wsp", "GC 0 wsp", "GC wsp", "BM wsp"))
levels(res$method) <- c("GC \u03bb=1", "GC \u03bb=0", "GC", "BM wsp", "GC \u03bb=1 wsp", "GC \u03bb=0 wsp", "GC wsp")
res$trait <- factor(res$trait, levels = traitnames)
levels(res$trait) <- c("corolla diameter", "anther length", "stigma length")
res$measurement_error <- factor(res$measurement_error, levels = c(FALSE, TRUE))
levels(res$measurement_error) <- c("no wsp var", "with wsp var")

## wAIC
p <- ggplot(res, aes(x=method, y=wAIC, fill=isbestaic)) +
  geom_col() +
  scale_fill_manual(values = c("grey50", "black")) + guides(fill = "none") +
  labs(x = "", y = "wAIC") +
  facet_grid(cols = vars(trait), scales = "free") +
  scale_x_discrete(labels = c(expression(GC~(λ==1)),
                              expression(GC~(λ==0)),
                              expression(GC),
                              expression(BM+σ[w]^2),
                              expression(GC~(λ==1)+σ[w]^2),
                              expression(GC~(λ==0)+σ[w]^2),
                              expression(GC+σ[w]^2))) +
  theme_bw() +
  theme(text=element_text(size=unit(9, 'pt')),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=unit(6, 'pt')))
p

ggsave(filename = file.path(result_dir, "tomato_wAIC.pdf"),
       plot = p,
       width = 1.6*columnwidth,
       height = 0.8*columnwidth,
       unit = "in",
       device = grDevices::cairo_pdf)

## sigma2 and sigma2_error
resbis <- tidyr::pivot_longer(res, cols = c("s2_wsp", "sigma2"), names_to = "var", values_to = "value")
resbis$var <- factor(resbis$var, levels = c("sigma2", "s2_wsp"))
levels(resbis$var) <- c("estimated σ²", "estimated total \n within-pop. variance")

p <- ggplot(resbis, aes(x=method, y=value, color=measurement_error)) +
  geom_point() +
  scale_colour_manual(values = c("black", "grey50")) + guides(colour = "none") +
  labs(x = "", y = "") +
  ggh4x::facet_grid2(cols = vars(trait), rows = vars(var), scales = "free_y", switch = "y", independent = "y") +
  scale_x_discrete(labels = c(expression(GC~(λ==1)),
                              expression(GC~(λ==0)),
                              expression(GC),
                              expression(BM+σ[w]^2),
                              expression(GC~(λ==1)+σ[w]^2),
                              expression(GC~(λ==0)+σ[w]^2),
                              expression(GC+σ[w]^2))) +
  theme_bw() +
  scale_y_sqrt() +
  theme(text=element_text(size=unit(9, 'pt')),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=unit(6, 'pt')),
        strip.placement = "outside", strip.background.y = element_blank())
p

ggsave(filename = file.path(result_dir, "tomato_variances.pdf"),
       plot = p,
       width = 1.8*columnwidth,
       height = 1.1*columnwidth,
       unit = "in",
       device = grDevices::cairo_pdf)

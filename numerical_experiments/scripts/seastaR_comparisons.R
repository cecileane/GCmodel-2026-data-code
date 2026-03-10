library(ape) # v5.8-1
# remotes::install_github("larabreithaupt/seastaR") # v0.1.0
library(seastaR)
library(phytools) # required by seastaR. v2.3-0

# remotes::install_github("lamho86/phylolm") # v2.7.0
library(phylolm)

library(ggplot2)
library(cowplot)

library(here)
result_dir <- here("numerical_experiments", "results")
run_date_today <- paste0(format(Sys.time(), "%Y-%m-%d"))
run_date <- run_date_today
# to load data calculated earlier, and fine-tune figures
run_date <- "2026-03-03"
result_dir_date <- file.path(result_dir, run_date)
dir.create(result_dir, showWarnings=FALSE)
dir.create(result_dir_date, showWarnings=FALSE)

columnwidth <- 3 # in inches # 1 point / 72 = 1 inch
twocolumnwidth <- 6

################################################################################
## vcv_gc function
################################################################################

#' @title Variance matrix for the GC model
#'
#' @description
#' Computes the un-scaled variance matrix for the Gaussian Coalescent model,
#' using function `phylolm::transf.branch.lengths`.
#' Requires `phylolm` version 2.7.0 or above.
#'
#' @param tree a phylogenetic tree of class ape::phylo.
#' Branch lengths must be in coalescent units.
#' @param lambda the ratio between the ancestral population variance and the
#' evolutionary rate variance.
#'
#' @return The un-scaled variance matrix of the tree.
#'
vcv_gc <- function(tree, lambda) {
  tree_trans <- phylolm::transf.branch.lengths(tree, model = "GC", parameters = list(lambda_GC = lambda))$tree
  return(vcv(tree_trans))
}

################################################################################
## Nested caterpillar tree
# To fine-tune the figure without re-running calculations:
# 1. adjust `run_date` above, on/near line 17
# 2. skip the code until the "Plot the (co)variances" section
################################################################################

#' @title Construct nested caterpillar tree
#'
#' @description
#' Construct a caterpillar tree with total height `Ttot`,
#' one outgroup `C`,
#' and a regular caterpillar clade with MRCA at time `t` from the root.
#' Clade `(A,B)` is the most nested in the caterpillar.
#'
#' @param n the number of species in the tree
#' @param t time between the root and the MRCA of the nested caterpillar tree
#' @param Ttot total height of the tree
#'
#' @return a phylogenetic tree with n tips
#'
get_nested_ladder_tree <- function(n, t, Ttot) {
  tree_ladder <- stree(n, type = "left")
  Tclade <- Ttot - t
  tree_ladder$edge.length <- rep(Tclade/(n-1), 2*n-2)
  tree_ladder$edge.length[2*(1:(n-2)) + 1] <- Tclade * (1 - 0:(n-3)/(n-1))
  tree_ladder$edge.length[2*n-2] <- Tclade * (1 - (n-3)/(n-1))
  tree_ladder$edge.length[2] <- t
  tree_ladder$edge.length[1] <- Ttot
  if (n > 3) {
    tree_ladder$tip.label <- c("C", paste0("B'", seq_len(n-3)), "B", "A")
  } else {
    tree_ladder$tip.label <- c("C", "B", "A")
  }
  return(tree_ladder)
}

## plot the tree
plot(get_nested_ladder_tree(n = 20, t = 5, Ttot = 10))
edgelabels(text = "T", edge = 1)
edgelabels(text = "t", edge = 2)

#' @title Get variance matrices for the triplet tree
#'
#' @description
#' Construct a nested caterpillar tree with `get_nested_ladder_tree`,
#' and get its associated variance matrix with `ape::vcv`, `seastaR` and `vcv_gc`.
#'
#' @param n the number of species in the tree
#' @param t time between the root and the MRCA of the nested caterpillar tree
#' @param Ttot total height of the tree
#' @param lambda variance ratio parameter for `vcv_gc`
#'
#' @return a data frame with the variances and covariance of traits at tips A, and C.
#'
get_var_covar <- function(n, t, Ttot, lambda) {
  tree_ladder <- get_nested_ladder_tree(n, t, Ttot)
  Ctree <- vcv(tree_ladder)
  Cstar <- seastaR::get_full_matrix(tree_ladder)
  Cils <- vcv_gc(tree_ladder, lambda)
  res <- data.frame(n = n, t = t, Ttot = Ttot, lambda = lambda,
                    var = c(Ctree[1, 1], Ctree[n, n], Ctree[1, n],
                            Cstar[1, 1], Cstar[n, n], Cstar[1, n],
                            Cils[1, 1], Cils[n, n], Cils[1, n]),
                    tip = rep(c("C", "A", "AC"), 3),
                    method = rep(c("Ctree", "Cstar", "Cils"), each = 3))
  return(res)
}


## computations
Ttot_all <- c(3, 30)
n_all <- c(3, 10, 20)
t_all <- seq(0, 1, 0.1)
lambda_all <- seq(0, 2, 0.1)

all_var <- NULL
for (Ttot in Ttot_all) {
  for (n in n_all) {
    for (lambda in lambda_all) {
      all_var_tmp <- lapply(t_all, function(tt) get_var_covar(n = n, t = tt * Ttot, Ttot = Ttot, lambda = lambda))
      all_var_tmp <- do.call(rbind, all_var_tmp)
      all_var <- rbind(all_var, all_var_tmp)
    }
  }
}

all_var_all_lambda <- all_var
write.table(all_var_all_lambda, file = file.path(result_dir_date, "all_var_A_C_AC.csv"))

## Plot the (co)variances - lambda fixed to 1

if(run_date != run_date_today){ # then load prior computations
  csv_filename = file.path(result_dir_date, "all_var_A_C_AC.csv")
  if(file.exists(csv_filename)){
    all_var_all_lambda <- read.table(file = csv_filename)
    cat("loaded prior data in 'all_var_all_lambda':",
        nrow(all_var_all_lambda), "rows,", ncol(all_var_all_lambda), "columns\n")
  } else {
    cat("no loading of prior data: no file\n", csv_filename, "\n", sep="")
  }
}

if ("lambda" %in% colnames(all_var_all_lambda)){
  all_var <- subset(all_var_all_lambda, lambda == 1)
}

p1 <-
  ggplot(
    subset(all_var, tip %in% c("A", "C") & method %in% c("Cstar")),
    aes(x = t / Ttot, y = var, color = as.factor(n))) +
  facet_grid(
    rows = vars(factor(Ttot, c(30, 3))), cols = vars(tip),
    scales = "free",
    labeller = labeller(.rows=function(x){paste("T =",x)})) +
  geom_line(
    data = subset(all_var, tip %in% c("A", "C") & method %in% c("Cils")),
    aes(group = interaction(method, n)),
    linetype = "solid",
    color = "black", linewidth=0.4) +
  geom_line(
    data = subset(all_var, tip %in% c("A", "C") & method %in% c("Ctree")),
    aes(group = interaction(method, n)),
    linetype = "solid",
    color = gray(0.5), linewidth=0.4) +
  geom_line(aes(group = interaction(method, n)),
    linetype = 6, linewidth=0.8) +
  scale_color_viridis_d(
    begin = 0.2, end = 0.8, option = "B",
    guide = guide_legend(override.aes=list(size=2), title="n taxa", order=1)) +
  ylab("variance") +
  xlab("clade height (t) / root height (T)") +
  scale_x_continuous(breaks = c(0,.5,1), labels=c("0","0.5","1"),
                     minor_breaks=seq(0,1,0.1)) +
  theme_bw() +
  theme(strip.text.y = element_blank(), # remove facet label: same as in p2
        panel.grid.major = element_line(linewidth=0.1 , color=grey(0.9)),
        panel.grid.minor = element_line(linewidth=0.05, color=grey(0.9)))

p2 <- ggplot(
    subset(all_var, tip %in% c("AC") & method %in% c("Cstar")),
    aes(x = t / Ttot, y = var, color = as.factor(n))) +
  facet_grid(
    rows = vars(factor(Ttot, c(30, 3))), cols = vars(tip),
    scales = "free",
    labeller = labeller(.rows=function(x){paste("T =",x)})) +
  geom_line(
    data = subset(all_var, tip %in% c("AC") & method %in% c("Cils")),
    aes(group = interaction(method, n)),
    linetype = "solid",
    color = "black", linewidth=0.4) +
  geom_line(
    data = subset(all_var, tip %in% c("AC") & method %in% c("Ctree")),
    aes(group = interaction(method, n)),
    linetype = "solid",
    color = gray(0.5), linewidth=0.4) +
  geom_line(aes(group = interaction(method, n)),
    linetype = 6, linewidth=0.8) +
  scale_color_viridis_d(
    begin = 0.2, end = 0.8, option = "B",
    guide = guide_legend(override.aes=list(size=2), title="n taxa", order=1)) +
  ylab("covariance") +
  xlab("t/T") +
  scale_x_continuous(breaks = c(0,.5,1), labels=c("0","0.5","1"),
                     minor_breaks=seq(0,1,0.1)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth=0.1 , color=grey(0.9)),
    panel.grid.minor = element_line(linewidth=0.05, color=grey(0.9)))

legend <- get_legend( # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p <- cowplot::plot_grid(p1 + theme(legend.position="none"),
                        p2 + theme(legend.position="none"),
                        legend,
                        nrow = 1,
                        rel_widths = c(1.7, 1, 0.5))
p

ggsave(filename = file.path(result_dir_date, "triplet_caterpillar_comparisons.pdf"),
       plot = p,
       width = twocolumnwidth,
       height = columnwidth,
       unit = "in")

## Plot - varying lambda

cils_vars <- subset(all_var_all_lambda,
    n == 3 & t == 1.5 & Ttot == 3 & method != "Cstar")

p1 <-
  ggplot(subset(cils_vars, tip %in% c("A", "C") & method == "Cils"),
         aes(x = lambda, y = var)) +
  facet_grid(
    cols = vars(tip),
    scales = "free") +
  geom_line(
    data = subset(cils_vars, tip %in% c("A", "C") & method %in% c("Ctree")),
    aes(group = interaction(method)),
    linetype = "solid",
    color = gray(0.5), linewidth=0.2) +
  geom_line(aes(group = interaction(method))) +
  ylab("variance") +
  xlab(expression(lambda)) +
  theme_bw() +
  theme(strip.text.y = element_blank(),
        panel.grid.major = element_line(linewidth=0.1 , color=grey(0.9)),
        panel.grid.minor = element_line(linewidth=0.05, color=grey(0.9)))

p2 <-
  ggplot(subset(cils_vars, tip %in% c("AC") & method == "Cils"),
         aes(x = lambda, y = var)) +
  facet_grid(
    cols = vars(tip),
    scales = "free") +
  geom_line(
    data = subset(cils_vars, tip %in% c("AC") & method %in% c("Ctree")),
    aes(group = interaction(method)),
    linetype = "solid",
    color = gray(0.5), linewidth=0.2) +
  geom_line(aes(group = interaction(method))) +
  ylab("variance") +
  xlab(expression(lambda)) +
  theme_bw() +
  theme(strip.text.y = element_blank(),
        panel.grid.major = element_line(linewidth=0.1 , color=grey(0.9)),
        panel.grid.minor = element_line(linewidth=0.05, color=grey(0.9)))

p <- cowplot::plot_grid(p1 + theme(legend.position="none"),
                        p2 + theme(legend.position="none"),
                        nrow = 1,
                        rel_widths = c(1.7, 1, 0.5))
p

ggsave(filename = file.path(result_dir_date, "triplet_caterpillar_lambda.pdf"),
       plot = p,
       width = twocolumnwidth,
       height = columnwidth,
       unit = "in")

################################################################################
## mixed setting with ladder - four tips
# To fine-tune the figure without re-running calculations:
# 1. adjust `run_date` above, on/near line 14
# 2. skip the code until the "Plot the covariances" section
################################################################################

#' @title Construct nested caterpillar tree
#'
#' @description
#' Construct a caterpillar tree with total height `Ttot`,
#' one outgroup `(C,D)`,
#' and a regular caterpillar clade with MRCA at time `t` from the root.
#' Clade `(A,B)` is the most nested in the caterpillar.
#' Both clades `(A,B)` and `(C,D)` are at distance `9*Ttot/10` from the root,
#' and have the same covariance under the classical BM.
#'
#' @param n the number of species in the tree
#' @param t time between the root and the MRCA of the nested caterpillar tree
#' @param Ttot total height of the tree
#'
#' @return a phylogenetic tree with n tips
#'
get_nested_ladder_tree <- function(n, t, Ttot) {
  Tclade <- Ttot - t
  nclade <- n - 1
  tree_ladder <- stree(nclade, type = "left")
  if (nclade > 3) {
    tree_ladder$edge.length <- rep((Tclade - Ttot / 10)/(nclade-3), 2*nclade-2)
    tree_ladder$edge.length[2*(1:(nclade-2)) + 1] <- Tclade - (Tclade - Ttot / 10) * 0:(nclade-3)/(nclade-3)
    tree_ladder$edge.length[2*nclade-2] <- Ttot / 10
    tree_ladder$edge.length[2] <- t
    tree_ladder$edge.length[1] <- Ttot
    tree_ladder$tip.label <- c("C", paste0("B'", seq_len(nclade-3)), "B", "A")
  } else {
    tree_ladder$edge.length[4] <- Ttot / 10
    tree_ladder$edge.length[3] <- Ttot / 10
    tree_ladder$edge.length[2] <- 9 * Ttot / 10
    tree_ladder$edge.length[1] <- Ttot
    tree_ladder$tip.label <- c("C", "B", "A")
  }
  tree_ladder <- phytools::bind.tip(tree_ladder, tip.label = "D",
    edge.length = Ttot / 10, where = 1, position = Ttot / 10)
  return(tree_ladder)
}

## plot the tree
plot(get_nested_ladder_tree(20, 5, 10))
edgelabels(text = "9T/10", edge = 1)
edgelabels(text = "T/10", edge = 38)
edgelabels(text = "T/10", edge = 2)
edgelabels(text = "t", edge = 4)

#' @title Get variance matrices for the triplet tree
#'
#' @description
#' Construct a nested caterpillar tree with `get_nested_ladder_tree`,
#' and get its associated variance matrix with `ape::vcv`, `seastaR` and `vcv_gc`.
#'
#' @param n the number of species in the tree
#' @param t time between the root and the MRCA of the nested caterpillar tree
#' @param Ttot total height of the tree
#' @param lambda variance ratio parameter for `vcv_gc`
#'
#' @return a data frame with the variances and covariance of traits at tips A, and C.
#'
get_var_covar <- function(n, t, Ttot, lambda) {
  tree_ladder <- get_nested_ladder_tree(n, t, Ttot)
  Ctree <- vcv(tree_ladder)
  Cstar <- seastaR::get_full_matrix(tree_ladder)
  Cils <- vcv_gc(tree_ladder, lambda)
  res <- data.frame(
    n = n, t = t, Ttot = Ttot, lambda = lambda,
    var = c(Ctree[1, 1], Ctree[n, n], Ctree[1, n], Ctree[n-1, n], Ctree[1, 2],
            Cstar[1, 1], Cstar[n, n], Cstar[1, n], Cstar[n-1, n], Cstar[1, 2],
            Cils[1, 1], Cils[n, n], Cils[1, n], Cils[n-1, n], Cils[1, 2]),
    tip = rep(c("C", "A", "AC", "AB", "CD"), 3),
    method = rep(c("Ctree", "Cstar", "Cils"), each = 5))
  return(res)
}


## computations
Ttot_all <- c(3, 30)
n_all <- c(4, 10, 20)
t_all <- c(seq(0, 0.8, 0.1), 0.89999) # replace t=.9 by 0.89999 because
# C*(A,B) = NaN if t/T = 0.9: the tree has 0-length branches in that case
lambda_all <- seq(0, 2, 0.1)

all_var <- NULL
for (Ttot in Ttot_all) {
  for (n in n_all) {
    for (lambda in lambda_all) {
      all_var_tmp <- lapply(t_all, function(tt) get_var_covar(n = n, t = tt * Ttot, Ttot = Ttot, lambda = lambda))
      all_var_tmp <- do.call(rbind, all_var_tmp)
      all_var <- rbind(all_var, all_var_tmp)
    }
  }
}

all_var_all_lambda <- all_var
write.table(all_var_all_lambda, file = file.path(result_dir_date, "all_var_A_C_AC_AB_CD.csv"))

## Plot the covariances - lambda = 1

if (run_date != run_date_today){ # then load prior computations
  csv_filename = file.path(result_dir_date, "all_var_A_C_AC_AB_CD.csv")
  if (file.exists(csv_filename)){
    all_var_all_lambda <- read.table(file = csv_filename)
    cat("loaded prior data in 'all_var':",
      nrow(all_var_all_lambda), "rows,", ncol(all_var_all_lambda), "columns\n")
  } else {
    cat("no loading of prior data: no file\n", csv_filename, "\n", sep="")
  }
}

if ("lambda" %in% colnames(all_var_all_lambda)){
  all_var <- subset(all_var_all_lambda, lambda == 1)
} else {
  all_var <- all_var_all_lambda
}

# customize tickmarks to make them similar for T=3 versus 30
cov_grid.major.breaks = c(1,2,2.7, 10,20,27) # for T=3 then T=30
# expand limits less than default: to avoid 2.7 on the plot for T=30.
cov_grid.minor.breaks = c(.5,1.5,2.5, 5,15,25) # otherwise halfway between major

p <- ggplot(
    subset(all_var, tip %in% c("AB", "CD") & method %in% c("Cstar")),
    aes(x = t / Ttot, y = var, color = as.factor(n))) +
  facet_grid(
    rows = vars(factor(Ttot, c(30, 3))), cols = vars(tip),
    scales = "free",
    labeller = labeller(.rows=function(x){paste("T =",x)})) +
  geom_line(
    data = subset(all_var, tip %in% c("AB", "CD") & method %in% c("Ctree")),
    aes(group = interaction(method, n)),
    linetype = "solid",
    color = gray(0.5), linewidth=0.4) +
  geom_line(
    data = subset(all_var, tip %in% c("AB", "CD") & method %in% c("Cils")),
    aes(group = interaction(method, n)),
    linetype = "solid",
    color = "black", linewidth=0.4) +
  geom_line(
    data = subset(all_var, tip %in% c("AB", "CD")),
    aes(group = interaction(method, n)),
    linetype = 6, linewidth=0.8) +
  scale_color_viridis_d(
    begin = 0.2, end = 0.8, option = "B",
    guide = guide_legend(override.aes=list(size=2), title="n taxa", order=1)) +
  ylab("covariance") +
  xlab("clade height (t) / root height (T)") +
  scale_x_continuous(breaks = c(0,.5,.9), labels=c("0","0.5","0.9"),
                     minor_breaks=seq(0,0.9,0.1), limits=c(0,0.9)) +
  scale_y_continuous(expand=expansion(mult = c(0.025,.05)),
                     breaks = cov_grid.major.breaks,
                     minor_breaks=cov_grid.minor.breaks) +
  theme_bw() + theme(
    panel.grid.major = element_line(linewidth=0.1 , color=grey(0.9)),
    panel.grid.minor = element_line(linewidth=0.05, color=grey(0.9)))
p

ggsave(filename = file.path(result_dir_date, "quadruplet_caterpillar_comparisons.pdf"),
       plot = p,
       width = twocolumnwidth,
       height = columnwidth,
       unit = "in")

## Plot the covariances - lambda varying

cils_vars <- subset(all_var_all_lambda, n == 4 & t == 1.5 & Ttot == 3 & method != "Cstar")

# customize tickmarks to make them similar for T=3 versus 30
cov_grid.major.breaks = c(1,2,2.7, 10,20,27) # for T=3 then T=30
# expand limits less than default: to avoid 2.7 on the plot for T=30.
cov_grid.minor.breaks = c(.5,1.5,2.5, 5,15,25) # otherwise halfway between major

p <- ggplot(subset(cils_vars, tip %in% c("AB", "CD") & method == "Cils"),
            aes(x = lambda, y = var)) +
  facet_grid(
    cols = vars(tip),
    scales = "free",
    labeller = labeller(.rows=function(x){paste("T =",x)})) +
  geom_line(
    data = subset(cils_vars, tip %in% c("AB", "CD") & method %in% c("Ctree")),
    linetype = "solid", # "dotted"
    color = gray(0.5), linewidth=0.2) +
  geom_line() +
  scale_color_viridis_d(
    end = 0.8, option = "B",
    guide = guide_legend(override.aes=list(size=2), title="n taxa", order=1)) +
  ylab("covariance") +
  xlab(expression(lambda)) +
  theme_bw() + theme(
    panel.grid.major = element_line(linewidth=0.1 , color=grey(0.9)),
    panel.grid.minor = element_line(linewidth=0.05, color=grey(0.9)))
p

ggsave(filename = file.path(result_dir_date, "quadruplet_caterpillar_lambda.pdf"),
       plot = p,
       width = twocolumnwidth,
       height = columnwidth,
       unit = "in")

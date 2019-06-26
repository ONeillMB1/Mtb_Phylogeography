library(ggplot2)
library(reshape2)
library(cowplot)
library(viridis)
library(akima)
library(dplyr)
library(scales)

#Functions to generate plots
sfs <- function(dir_path) {
  observed_SFS <- read.table(paste(dir_path, "observedSFS.txt", sep=""))
  neutral_SFS <- read.table(paste(dir_path, "neutralModelSFS.txt", sep=""))
  expansion_SFS <- read.table(paste(dir_path, "expansionModelSFS.txt", sep=""))
  growth_SFS <- read.table(paste(dir_path, "growthModelSFS.txt", sep=""))
  
  SFS <- rbind(data.frame(Model="Observed",
                          Count=as.numeric(as.character(observed_SFS$V1)),
                          Frequency=as.numeric(row.names(observed_SFS))),
               data.frame(Model="Neutral", 
                          Count=as.numeric(as.character(neutral_SFS$V1)),
                          Frequency=as.numeric(row.names(neutral_SFS))),
               data.frame(Model = "Expansion", 
                          Count=as.numeric(as.character(expansion_SFS$V1)), 
                          Frequency=as.numeric(row.names(expansion_SFS))),
               data.frame(Model = "Growth", 
                          Count=as.numeric(as.character(growth_SFS$V1)), 
                          Frequency=as.numeric(row.names(growth_SFS))))
  
  sfs.p <- ggplot(SFS, aes(x=Frequency, y=Count, fill=Model)) +
    geom_bar(stat="identity", position="dodge", show.legend = TRUE) + 
    #scale_fill_brewer(palette="Paired") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(.8,.6),
          panel.background = element_blank())
  
  return(sfs.p)
}

expansion_params <- function(dir_path) {
  expansion_params <- read.table(paste(dir_path, "dadi_expansionModel.txt", sep=""), header = TRUE)
  expansion_params <- setNames(expansion_params, c("v","t","LL"))
  params_df <- melt(expansion_params[,1:2], variable.name = "Parameter", value.name = "Value")
  expansion.nu.plot <- ggplot(data = params_df[1:100,], aes(x=Parameter, y=Value)) + geom_boxplot() +
    xlab("") +
    theme(text = element_text(size=16),
          axis.text.x = element_text(angle=90, vjust=1, size = 16))
  
  expansion.tau.plot <- ggplot(data = params_df[101:200,], aes(x=Parameter, y=Value)) + geom_boxplot() + 
    xlab("") +
    theme(text = element_text(size=16),
          axis.text.x = element_text(angle=90, vjust=1, size = 16))
  
  ex_params.p <- plot_grid(expansion.nu.plot, expansion.tau.plot, align = "h")
  
  return(ex_params.p)
}

expansion_params_corr <- function(dir_path) {
  expansion_params <- read.table(paste(dir_path, "dadi_expansionModel.txt", sep=""), header = TRUE)
  expansion_params <- setNames(expansion_params, c("v","t","LL"))
  dat <- cor(as.matrix(expansion_params[,c("v","t", "LL")]), use = "pairwise.complete.obs", method = "pearson")
  dat.m = melt(dat)
  
  exp.params.corr.p <- ggplot(dat.m, aes(x=Var2,y=Var1)) +
    geom_tile(aes(fill= value)) +
    geom_text(aes(fill = value, label = round(value, 2)), size=6) +
    xlab(NULL) + ylab(NULL) +
    scale_fill_gradient2(high="#b8b8b8", mid="white", low="#b8b8b8", midpoint=0, limits =c(-1,1), guide=FALSE) +
    theme_bw() +
    theme(text = element_text(size=16),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major= element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          panel.grid.minor= element_blank())
  
  return(exp.params.corr.p)
}

expansion_AIC <- function(dir_path, ll) {
  expansion_params <- read.table(paste(dir_path, "dadi_expansionModel.txt", sep=""), header = TRUE)
  expansion_params <- setNames(expansion_params, c("v","t","LL"))
  Expansion <- expansion_params$LL
  ll <- cbind(ll, Expansion)
  avg_expansion_ll <- median(expansion_params$LL)
  cat("Average LL: ", avg_expansion_ll, "\n")
  avg_expansion_AIC <- 2 * 2 - 2 * avg_expansion_ll
  cat("Average AIC: ", avg_expansion_AIC, "\n")
  max_expansion_ll <- max(expansion_params$LL)
  cat("Max LL: ", max_expansion_ll, "\n")
  max_expansion_AIC <- 2 * 2 - 2 * max_expansion_ll
  cat("Max AIC: ", max_expansion_AIC, "\n")
  max_expansion_ll_idx <- which(expansion_params$LL == max_expansion_ll)
  max_expansion_nu <- expansion_params$v[max_expansion_ll_idx]
  cat("v: ", max_expansion_nu, "\n")
  max_expansion_tau <- expansion_params$t[max_expansion_ll_idx]
  cat("t: ", max_expansion_tau, "\n")
  
  return(ll)
}

expansion_ls <- function(dir_path, no_bins) {
  param_grid_dadi <- read.table(paste(dir_path, "likelihood_grid_expansion.txt", sep=""), header=TRUE)
  param_grid_dadi_sub <- param_grid_dadi[
    param_grid_dadi$nu >= 
      min(param_grid_dadi[param_grid_dadi$LL > quantile(param_grid_dadi$LL, probs = c(0.25)), 'nu']) &
    param_grid_dadi$nu <= 
      max(param_grid_dadi[param_grid_dadi$LL > quantile(param_grid_dadi$LL, probs = c(0.25)), 'nu']) &
    param_grid_dadi$T >= 
      min(param_grid_dadi[param_grid_dadi$LL > quantile(param_grid_dadi$LL, probs = c(0.25)), 'T']) &
    param_grid_dadi$T <= 
      max(param_grid_dadi[param_grid_dadi$LL > quantile(param_grid_dadi$LL, probs = c(0.25)), 'T']),]
  ls.p <- ggplot(param_grid_dadi_sub, aes(x=nu, y=T, z=LL)) +
    geom_tile(aes(fill = LL)) + 
    stat_contour(colour = "black", bins=no_bins) +
    scale_fill_viridis(option = "inferno") +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
    ylab(expression(tau)) +
    xlab(expression(nu)) +
    theme_bw() +
    theme(text = element_text(size=14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="right",
          legend.title=element_blank(),
          legend.text = element_text(size=12))
}

growth_params <- function(dir_path) {
  growth_params <- read.table(paste(dir_path, "dadi_growthModel.txt", sep=""), header = TRUE)
  growth_params <- setNames(growth_params, c("v","t","LL"))
  params_df <- melt(growth_params[,1:2], variable.name = "Parameter", value.name = "Value")
  growth.nu.plot <- ggplot(data = params_df[1:100,], aes(x=Parameter, y=Value)) + geom_boxplot() +
    xlab("") +
    theme(text = element_text(size=16),
          axis.text.x = element_text(angle=90, vjust=1, size = 16))
  
  growth.tau.plot <- ggplot(data = params_df[101:200,], aes(x=Parameter, y=Value)) + geom_boxplot() + 
    xlab("") +
    theme(text = element_text(size=16),
          axis.text.x = element_text(angle=90, vjust=1, size = 16))
  
  ex_params.p <- plot_grid(growth.nu.plot, growth.tau.plot, align = "h")
  
  return(ex_params.p)
}

growth_params_corr <- function(dir_path) {
  growth_params <- read.table(paste(dir_path, "dadi_growthModel.txt", sep=""), header = TRUE)
  growth_params <- setNames(growth_params, c("v","t","LL"))
  dat <- cor(as.matrix(growth_params[,c("v","t", "LL")]), use = "pairwise.complete.obs", method = "pearson")
  dat.m = melt(dat)
  
  exp.params.corr.p <- ggplot(dat.m, aes(x=Var2,y=Var1)) +
    geom_tile(aes(fill= value)) +
    geom_text(aes(fill = value, label = round(value, 2)), size=6) +
    xlab(NULL) + ylab(NULL) +
    scale_fill_gradient2(high="#b8b8b8", mid="white", low="#b8b8b8", midpoint=0, limits =c(-1,1), guide=FALSE) +
    theme_bw() +
    theme(text = element_text(size=16),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major= element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          panel.grid.minor= element_blank())
  
  return(exp.params.corr.p)
}

growth_AIC <- function(dir_path, ll) {
  growth_params <- read.table(paste(dir_path, "dadi_growthModel.txt", sep=""), header = TRUE)
  growth_params <- setNames(growth_params, c("v","t","LL"))
  growth <- growth_params$LL
  ll <- cbind(ll, growth)
  avg_growth_ll <- median(growth_params$LL)
  cat("Average LL: ", avg_growth_ll, "\n")
  avg_growth_AIC <- 2 * 2 - 2 * avg_growth_ll
  cat("Average AIC: ", avg_growth_AIC, "\n")
  max_growth_ll <- max(growth_params$LL)
  cat("Max LL: ", max_growth_ll, "\n")
  max_growth_AIC <- 2 * 2 - 2 * max_growth_ll
  cat("Max AIC: ", max_growth_AIC, "\n")
  max_growth_ll_idx <- which(growth_params$LL == max_growth_ll)
  max_growth_nu <- growth_params$v[max_growth_ll_idx]
  cat("v: ", max_growth_nu, "\n")
  max_growth_tau <- growth_params$t[max_growth_ll_idx]
  cat("t: ", max_growth_tau, "\n")
  
  return(ll)
}

growth_ls <- function(dir_path, no_bins) {
  param_grid_dadi <- read.table(paste(dir_path, "likelihood_grid_growth.txt", sep=""), header=TRUE)
  param_grid_dadi_sub <- param_grid_dadi[
    param_grid_dadi$nu >= 
      min(param_grid_dadi[param_grid_dadi$LL > quantile(param_grid_dadi$LL, probs = c(0.25)), 'nu']) &
      param_grid_dadi$nu <= 
      max(param_grid_dadi[param_grid_dadi$LL > quantile(param_grid_dadi$LL, probs = c(0.25)), 'nu']) &
      param_grid_dadi$T >= 
      min(param_grid_dadi[param_grid_dadi$LL > quantile(param_grid_dadi$LL, probs = c(0.25)), 'T']) &
      param_grid_dadi$T <= 
      max(param_grid_dadi[param_grid_dadi$LL > quantile(param_grid_dadi$LL, probs = c(0.25)), 'T']),]
  ls.p <- ggplot(param_grid_dadi_sub, aes(x=nu, y=T, z=LL)) +
    geom_tile(aes(fill = LL)) + 
    stat_contour(colour = "black", bins=no_bins) +
    scale_fill_viridis(option = "inferno") +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
    ylab(expression(tau)) +
    xlab(expression(nu)) +
    theme_bw() +
    theme(text = element_text(size=14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="right",
          legend.title=element_blank(),
          legend.text = element_text(size=12))
}

ll_comp <- function(ll) {
  ll$a <- NULL
  comp.p <- ggplot(data = melt(ll, variable.name = "Model", value.name = "LL"), aes(x = Model, y = LL)) + 
    geom_boxplot(aes(color = Model)) + 
    theme(text = element_text(size=16))
  
  return(comp.p)
}

comb_data <- function(dir_path) {
  growth_params <- read.table(paste(dir_path, "dadi_growthModel.txt", sep=""), header = TRUE)
  growth_params <- setNames(growth_params, c("v","t","LL"))
  expansion_params <- read.table(paste(dir_path, "dadi_expansionModel.txt", sep=""), header = TRUE)
  expansion_params <- setNames(expansion_params, c("v","t","LL"))
  df <- rbind(data.frame(model="expansion", 
                        v = expansion_params$v, 
                        t = expansion_params$t, 
                        LL = expansion_params$LL)#,
              #data.frame(model="growth", 
                        #v = growth_params$v, 
                        #t = growth_params$t, 
                        #LL = growth_params$LL)
              )
  df.m <- melt(df, id.vars=c("model"))
  p <- ggplot(df.m) + 
    facet_wrap(~variable, scales='free_y') +
    scale_y_continuous(breaks=pretty_breaks()) +
    geom_boxplot(aes(model, value)) +
    xlab("") +
    ylab("")
  
  return(p)
}


#Data

#All/OldWorld
all_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/OldWorld/"
all.ll <- data.frame(a=rep(NA, times = 100))
all.sfs <- sfs(all_dir)

all.exp.params <- expansion_params(all_dir)
all.exp.params.corr <- expansion_params_corr(all_dir)
all.ll <- expansion_AIC(all_dir, all.ll)
all.exp.ls <- expansion_ls(all_dir, 100)

all.growth.params <- growth_params(all_dir)
all.growth.params.corr <- growth_params_corr(all_dir)
all.ll <- growth_AIC(all_dir, all.ll)
all.growth.ls <- growth_ls(all_dir, 100)

all.comp <- ll_comp(all.ll)
all.comb.p <- comb_data(all_dir)

#Lineage 1
lin1_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin1/"
lin1.ll <- data.frame(a=rep(NA, times = 100))
lin1.sfs <- sfs(lin1_dir)

lin1.exp.params <- expansion_params(lin1_dir)
lin1.exp.params.corr <- expansion_params_corr(lin1_dir)
lin1.ll <- expansion_AIC(lin1_dir, lin1.ll)
lin1.exp.ls <- expansion_ls(lin1_dir, 100)

lin1.growth.params <- growth_params(lin1_dir)
lin1.growth.params.corr <- growth_params_corr(lin1_dir)
lin1.ll <- growth_AIC(lin1_dir, lin1.ll)
lin1.growth.ls <- growth_ls(lin1_dir, 100)

lin1.comp <- ll_comp(lin1.ll)
lin1.comb.p <- comb_data(lin1_dir)

l1 <- plot_grid(
  lin1.sfs,
  lin1.comb.p,
  lin1.exp.ls,
  ncol=3
)


#Lineage 2
lin2_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin2/"
lin2.ll <- data.frame(a=rep(NA, times = 100))
lin2.sfs <- sfs(lin2_dir)

lin2.exp.params <- expansion_params(lin2_dir)
lin2.exp.params.corr <- expansion_params_corr(lin2_dir)
lin2.ll <- expansion_AIC(lin2_dir, lin2.ll)
lin2.exp.ls <- expansion_ls(lin2_dir, 100)

lin2.growth.params <- growth_params(lin2_dir)
lin2.growth.params.corr <- growth_params_corr(lin2_dir)
lin2.ll <- growth_AIC(lin2_dir, lin2.ll)
lin2.growth.ls <- growth_ls(lin2_dir, 100)

lin2.comp <- ll_comp(lin2.ll)

#Lineage 3
lin3_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin3/"
lin3.ll <- data.frame(a=rep(NA, times = 100))
lin3.sfs <- sfs(lin3_dir)

lin3.exp.params <- expansion_params(lin3_dir)
lin3.exp.params.corr <- expansion_params_corr(lin3_dir)
lin3.ll <- expansion_AIC(lin3_dir, lin3.ll)
lin3.exp.ls <- expansion_ls(lin3_dir, 100)

lin3.growth.params <- growth_params(lin3_dir)
lin3.growth.params.corr <- growth_params_corr(lin3_dir)
lin3.ll <- growth_AIC(lin3_dir, lin3.ll)
lin3.growth.ls <- growth_ls(lin3_dir, 100)

lin3.comp <- ll_comp(lin3.ll)

#Lineage 4
lin4_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin4/"
lin4.ll <- data.frame(a=rep(NA, times = 100))
lin4.sfs <- sfs(lin4_dir)

lin4.exp.params <- expansion_params(lin4_dir)
lin4.exp.params.corr <- expansion_params_corr(lin4_dir)
lin4.ll <- expansion_AIC(lin4_dir, lin4.ll)
lin4.exp.ls <- expansion_ls(lin4_dir, 100)

lin4.growth.params <- growth_params(lin4_dir)
lin4.growth.params.corr <- growth_params_corr(lin4_dir)
lin4.ll <- growth_AIC(lin4_dir, lin4.ll)
lin4.growth.ls <- growth_ls(lin4_dir, 100)

lin4.comp <- ll_comp(lin4.ll)

#Lineage 5
lin5_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin5/"
lin5.ll <- data.frame(a=rep(NA, times = 100))
lin5.sfs <- sfs(lin5_dir)

lin5.exp.params <- expansion_params(lin5_dir)
lin5.exp.params.corr <- expansion_params_corr(lin5_dir)
lin5.ll <- expansion_AIC(lin5_dir, lin5.ll)
lin5.exp.ls <- expansion_ls(lin5_dir, 100)

lin5.growth.params <- growth_params(lin5_dir)
lin5.growth.params.corr <- growth_params_corr(lin5_dir)
lin5.ll <- growth_AIC(lin5_dir, lin5.ll)
lin5.growth.ls <- growth_ls(lin5_dir, 100)

lin5.comp <- ll_comp(lin5.ll)

#Lineage 6
lin6_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin6/"
lin6.ll <- data.frame(a=rep(NA, times = 100))
lin6.sfs <- sfs(lin6_dir)

lin6.exp.params <- expansion_params(lin6_dir)
lin6.exp.params.corr <- expansion_params_corr(lin6_dir)
lin6.ll <- expansion_AIC(lin6_dir, lin6.ll)
lin6.exp.ls <- expansion_ls(lin6_dir, 100)

lin6.growth.params <- growth_params(lin6_dir)
lin6.growth.params.corr <- growth_params_corr(lin6_dir)
lin6.ll <- growth_AIC(lin6_dir, lin6.ll)
lin6.growth.ls <- growth_ls(lin6_dir, 100)

lin6.comp <- ll_comp(lin6.ll)

#Lineage 7
lin7_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin7/"
lin7.ll <- data.frame(a=rep(NA, times = 100))
lin7.sfs <- sfs(lin7_dir)

lin7.exp.params <- expansion_params(lin7_dir)
lin7.exp.params.corr <- expansion_params_corr(lin7_dir)
lin7.ll <- expansion_AIC(lin7_dir, lin7.ll)
lin7.exp.ls <- expansion_ls(lin7_dir, 100)

lin7.growth.params <- growth_params(lin7_dir)
lin7.growth.params.corr <- growth_params_corr(lin7_dir)
lin7.ll <- growth_AIC(lin7_dir, lin7.ll)
lin7.growth.ls <- growth_ls(lin7_dir, 100)

lin7.comp <- ll_comp(lin7.ll)


#Full
ow_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/OldWorld//"
ow.ll <- data.frame(a=rep(NA, times = 100))
ow.sfs <- sfs(ow_dir)

ow.exp.params <- expansion_params(ow_dir)
ow.exp.params.corr <- expansion_params_corr(ow_dir)
ow.ll <- expansion_AIC(ow_dir, ow.ll)
ow.exp.ls <- expansion_ls(ow_dir, 100)

ow.growth.params <- growth_params(ow_dir)
ow.growth.params.corr <- growth_params_corr(ow_dir)
ow.ll <- growth_AIC(ow_dir, ow.ll)
ow.growth.ls <- growth_ls(ow_dir, 100)

ow.comp <- ll_comp(ow.ll)


#Plot

comp <- plot_grid(
  lin1.comp,
  lin2.comp,
  lin3.comp,
  lin4.comp,
  lin5.comp,
  lin6.comp,
  lin7.comp,
  ncol=2,
  align="v"
)
comp


mod <- plot_grid(
  lin1.exp.params, lin1.growth.params, lin1.comp,
  lin2.exp.params, lin2.growth.params, lin2.comp,
  lin3.exp.params, lin3.growth.params, lin3.comp,
  lin4.exp.params, lin4.growth.params, lin4.comp,
  lin5.exp.params, lin5.growth.params, lin5.comp,
  lin6.exp.params, lin6.growth.params, lin6.comp,
  lin7.exp.params, lin7.growth.params, lin7.comp,
  ncol=3,
  align="v"
)
mod




#################
lin7_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin7/"
lin7.exp.ls <- expansion_ls(lin7_dir, 100)
lin7.growth.ls <- growth_ls(lin7_dir, 100)
lin7.comb.p <- comb_data(lin7_dir)

l7 <- plot_grid(
  lin7.comb.p,
  lin7.exp.ls,
  lin7.growth.ls,
  ncol=3
)
l7

lin6_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin6/"
lin6.exp.ls <- expansion_ls(lin6_dir, 100)
lin6.growth.ls <- growth_ls(lin6_dir, 100)
lin6.comb.p <- comb_data(lin6_dir)

l6 <- plot_grid(
  lin6.comb.p,
  lin6.exp.ls,
  lin6.growth.ls,
  ncol=3
)
l6

lin5_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin5/"
lin5.exp.ls <- expansion_ls(lin5_dir, 100)
lin5.growth.ls <- growth_ls(lin5_dir, 100)
lin5.comb.p <- comb_data(lin5_dir)

l5 <- plot_grid(
  lin5.comb.p,
  lin5.exp.ls,
  lin5.growth.ls,
  ncol=3
)
l5


lin4_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin4/"
lin4.exp.ls <- expansion_ls(lin4_dir, 100)
lin4.growth.ls <- growth_ls(lin4_dir, 100)
lin4.comb.p <- comb_data(lin4_dir)

l4 <- plot_grid(
  lin4.comb.p,
  lin4.exp.ls,
  lin4.growth.ls,
  ncol=3
)
l4

lin3_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin3/"
lin3.exp.ls <- expansion_ls(lin3_dir, 100)
lin3.growth.ls <- growth_ls(lin3_dir, 100)
lin3.comb.p <- comb_data(lin3_dir)

l3 <- plot_grid(
  lin3.comb.p,
  lin3.exp.ls,
  lin3.growth.ls,
  ncol=3
)
l3


lin2_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin2/"
lin2.exp.ls <- expansion_ls(lin2_dir, 100)
lin2.growth.ls <- growth_ls(lin2_dir, 100)
lin2.comb.p <- comb_data(lin2_dir)

l2 <- plot_grid(
  lin2.comb.p,
  lin2.exp.ls,
  lin2.growth.ls,
  ncol=3
)
l2


lin1_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/lin1/"
lin1.exp.ls <- expansion_ls(lin1_dir, 100)
lin1.growth.ls <- growth_ls(lin1_dir, 100)
lin1.comb.p <- comb_data(lin1_dir)

l1 <- plot_grid(
  lin1.comb.p,
  lin1.exp.ls,
  lin1.growth.ls,
  ncol=3
)
l1

ggsave(paste("figs/dadi/", format(Sys.time(), "%y-%m-%d_"), "SI_lin1_dadi", sep = ""), plot=l1, width=9, height=3)


ow_dir <- "C:/Users/Mary/PepLab/Mtb_Phylogeography_v2/data/170404_dadi/OldWorld/"
ow.exp.ls <- expansion_ls(ow_dir, 100)
ow.growth.ls <- growth_ls(ow_dir, 100)
ow.comb.p <- comb_data(ow_dir)

ow <- plot_grid(
  ow.comb.p,
  ow.exp.ls,
  ow.growth.ls,
  ncol=3
)
ow

#expansion LS
lin1.exp.ls <- expansion_ls(lin1_dir, 100)
lin2.exp.ls <- expansion_ls(lin2_dir, 100)
lin3.exp.ls <- expansion_ls(lin3_dir, 100)
lin4.exp.ls <- expansion_ls(lin4_dir, 100)
lin5.exp.ls <- expansion_ls(lin5_dir, 100)
lin6.exp.ls <- expansion_ls(lin6_dir, 100)
lin7.exp.ls <- expansion_ls(lin7_dir, 100)


exp.ls.p <- plot_grid(ow.exp.ls,
          lin1.exp.ls, 
          lin2.exp.ls,
          lin3.exp.ls,
          lin4.exp.ls,
          lin5.exp.ls,
          lin6.exp.ls,
          lin7.exp.ls,
          ncol=4)

#growth LS
lin1.gro.ls <- growth_ls(lin1_dir, 100)
lin2.gro.ls <- growth_ls(lin2_dir, 100)
lin3.gro.ls <- growth_ls(lin3_dir, 100)
lin4.gro.ls <- growth_ls(lin4_dir, 100)
lin5.gro.ls <- growth_ls(lin5_dir, 100)
lin6.gro.ls <- growth_ls(lin6_dir, 100)
lin7.gro.ls <- growth_ls(lin7_dir, 100)


gro.ls.p <- plot_grid(lin1.gro.ls, 
                      lin2.gro.ls,
                      lin3.gro.ls,
                      lin4.gro.ls,
                      lin5.gro.ls,
                      lin6.gro.ls,
                      lin7.gro.ls,
                      ncol=4)

#params exp
ow.comb.p <- comb_data(ow_dir)
lin1.comb.p <- comb_data(lin1_dir)
lin2.comb.p <- comb_data(lin2_dir)
lin3.comb.p <- comb_data(lin3_dir)
lin4.comb.p <- comb_data(lin4_dir)
lin5.comb.p <- comb_data(lin5_dir)
lin6.comb.p <- comb_data(lin6_dir)
lin7.comb.p <- comb_data(lin7_dir)

params.exp.p <- plot_grid(ow.comb.p,
                      lin1.comb.p, 
                      lin2.comb.p,
                      lin3.comb.p,
                      lin4.comb.p,
                      lin5.comb.p,
                      lin6.comb.p,
                      lin7.comb.p,
                      ncol=4)



ggplot(data = melt(l2.ll, variable.name = "Model", value.name = "LL"), aes(x = Model, y = LL)) + 
  geom_boxplot(aes(color = Model)) + 
  theme(text = element_text(size=16))
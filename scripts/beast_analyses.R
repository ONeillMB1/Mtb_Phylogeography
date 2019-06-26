#source("https://www.bioconductor.org/biocLite.R")

#biocLite("ggtree")
#devtools::install_github(c("hadley/ggplot2", "GuangchuangYu/ggtree"))

require(ggplot2)
require(ggtree)
require(cowplot)
require(scales)
require(splitstackshape)
require(reshape2)
require(dplyr)
require(viridis)

#function to generate image files
exportPlot <- function(gplot, filename, width=2, height=1.5) {
  ggsave(paste(filename,'.pdf',sep=""), gplot, width=width, height=height)
  postscript(file=paste(filename,'.eps',sep=""), width=width, height=height)
  print(gplot)
  dev.off()
  png(file=paste(filename,'.png', sep=""), width=width*100, height=height*100)
  print(gplot)
  dev.off()
}

#code to make axis pretty
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # remove the exponent for zero 
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
  l <- gsub("e\\+","e",l) 
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

#scale factor
scale_factor <- function(d, tmrca) {
  sf <- tmrca/max(d$height_median, na.rm=TRUE)
  return(sf)
}

#modify tree tip labs
mod.tip.labs <- function(tree, d) {
  tre_tips <- subset(d, isTip == TRUE, c(label, Geography))
  tre_tips$samp <- sapply(strsplit(tre_tips$label, "_"), "[[", 1)
  tre_tips$iso2 <- sapply(strsplit(tre_tips$label, "_"), "[[", 2)
  tre_tips$lin <- sapply(strsplit(tre_tips$label, "_"), "[[", 4)
  tre_tips$newName <- paste(tre_tips$samp, tre_tips$lin, tre_tips$iso2, tre_tips$Geography, sep="_")
  tree@phylo$tip.label <- tre_tips$newName
  return(tree)
}

#modify tree dat file
mod.tree.dat <- function(d, sf){
  d$Geography.set.prob = gsub("\\{|\\}", "", d$Geography.set.prob)
  d$Geography.set = gsub("\\{|\\}", "", d$Geography.set)
  #Extract parent node info
  d$parent_Geography <- d[match(d$parent, d$node), 'Geography']
  d$parent_Geography.prob <- d[match(d$parent, d$node), 'Geography.prob']
  d$parent_height_median <- d[match(d$parent, d$node), 'height_median']
  d$parent_posterior <- d[match(d$parent, d$node), 'posterior']
  #modify node heights with scale factor
  d.mod <- d
  d.mod$height_median <- d.mod$height_median * sf
  d.mod$parent_height_median <- d.mod$parent_height_median * sf
  #determine if migration has occurred 
  d.mod$mig <- ifelse(d.mod$Geography != d.mod$parent_Geography, "Yes", "No")
  #remove low confidence nodes
  d.mod <- subset(d.mod, d.mod$posterior >= 0.8 | is.na(d.mod$posterior))
  d.mod <- subset(d.mod, d.mod$parent_posterior >= 0.8)
  return(d.mod)
}

#Convert mod tree data to migration data
lin.mig.dat <- function(d.mod) {
  rr <- range(d.mod[,c('height_median', 'parent_height_median')])
  bb <- seq(floor(rr[1]),ceiling(rr[2]))
  mm <- apply(d.mod[,c('parent_height_median', 'height_median')],1,function(x) bb %in% seq(floor(x['height_median']),ceiling(x["parent_height_median"])))

  cl <- apply(d.mod[,c('parent_height_median', 'height_median')],1,function(x) bb %in% seq(floor(x['height_median']),ceiling(x["parent_height_median"])))
  l <- apply(cl,1,sum)
  clun <- apply(d.mod[d.mod$mig == "Yes",c('parent_height_median', 'height_median')],1,function(x) bb %in% seq(floor(x['height_median']),ceiling(x["parent_height_median"])))
  lun <- apply(clun,1,sum)
  
  gr <- apply(cl,1,function(x) length(unique((as.vector(as.matrix(d.mod[x,c('parent_Geography', 'Geography')]))))))
  
  d.mig <- cbind(data.frame(bb), data.frame(l), data.frame(lun), data.frame(gr))
  
  return(d.mig)
}
lin567.mig.dat <- function(d.mod) {
  rr <- range(d.mod[,c('height_median', 'parent_height_median')])
  bb <- seq(floor(rr[1]),ceiling(rr[2]))
  mm <- apply(d.mod[,c('parent_height_median', 'height_median')],1,function(x) bb %in% seq(floor(x['height_median']),ceiling(x["parent_height_median"])))
  
  cl <- apply(d.mod[,c('parent_height_median', 'height_median')],1,function(x) bb %in% seq(floor(x['height_median']),ceiling(x["parent_height_median"])))
  l <- apply(cl,1,sum)
  #clun <- apply(d.mod[d.mod$mig == "Yes",c('parent_height_median', 'height_median')],1,function(x) bb %in% seq(floor(x['height_median']),ceiling(x["parent_height_median"])))
  #lun <- apply(clun,1,sum)
  
  gr <- apply(cl,1,function(x) length(unique((as.vector(as.matrix(d.mod[x,c('parent_Geography', 'Geography')]))))))
  
  d.mig <- cbind(data.frame(bb), data.frame(l), data.frame(gr))
  
  return(d.mig)
}

#migration plot
gen.mig.p <- function(d.mig) {
  mig.p <- ggplot(d.mig[d.mig$l > 4,]) + 
  geom_line(aes(bb, lun/l), colour = "black") + 
  theme(legend.position = "none") + 
  xlab("") + 
  ylab("Migration Events") + 
  scale_x_reverse(limits=c(2000, 0), 
                    breaks=seq(500,2000, by=500)) + 
  #geom_hline(yintercept = 0.0, linetype="dotted") +
  geom_ribbon(aes(x=bb, ymax=((lun/l) + 1/l), ymin=ifelse(((lun/l) - 1/l)<0, 0, ((lun/l) - 1/l))), alpha=.2) +
  scale_fill_manual(values=c("gray")) + 
  scale_y_continuous(limits=c(-0.25,.75), breaks = c(0, .25, .50)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(linetype="dotted", colour="grey"),
        #panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank()) 
  return(mig.p)
}
gen.full.mig.p <- function(d.mig) {
  mig.p <- ggplot(d.mig[d.mig$l > 4,]) + 
    geom_line(aes(bb, lun/l), colour = "black") + 
    theme(legend.position = "none") + 
    xlab("") + 
    ylab("Migration Events") + 
    scale_x_reverse(limits=c(max(d.mig[d.mig$l>4, 'bb']), 0), 
                    breaks=seq(500,max(d.mig[d.mig$l>4, 'bb']), by=500)) + 
    #geom_hline(yintercept = 0.0, linetype="dotted") +
    geom_ribbon(aes(x=bb, ymax=((lun/l) + 1/l), ymin=((lun/l) - 1/l)), alpha=.2) +
    scale_fill_manual(values=c("gray")) + 
    scale_y_continuous(limits=c(-0.25,.75), breaks = c(0, .25, .50)) +
    theme_bw() +
    theme(panel.grid.major.y = element_line(linetype="dotted", colour="grey"),
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()) #,
  return(mig.p)
} 

#range expansion plot
gen.gr.p <- function(d.mig) {
  gr.p <- ggplot(d.mig[d.mig$l > 4,]) + 
    geom_line(aes(bb, gr), colour = "black") + 
    theme(legend.position = "none") + 
    xlab("") + 
    ylab("Regions") + 
    scale_x_reverse(limits=c(2000, 0), 
                    breaks=seq(500,2000, by=500)) + 
    scale_y_continuous(limits=c(0,16), breaks = c(4,8,12)) +
    theme_bw() +
    theme(panel.grid.major.y = element_line(linetype="dotted", colour="grey"),
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()) 
  return(gr.p)
}
gen.full.gr.p <- function(d.mig) {
  gr.p <- ggplot(d.mig[d.mig$l > 4,]) + 
    geom_line(aes(bb, gr), colour = "black") + 
    theme(legend.position = "none") + 
    xlab("") + 
    ylab("Regions") + 
    scale_x_reverse(limits=c(max(d.mig[d.mig$l>4, 'bb']), 0), 
                    breaks=seq(500,max(d.mig[d.mig$l>4, 'bb']), by=500)) + 
    scale_y_continuous(limits=c(0,16), breaks = c(4,8,12)) +
    theme_bw() +
    theme(panel.grid.major.y = element_line(linetype="dotted", colour="grey"),
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()) 
  return(gr.p)
}

#bsp
gen.bsp.p <- function(d.bsp, sf, d.mig) {
  #gt <- 62/(365/7) #62 weeks -> years
  d.bsp.scaled <- d.bsp * sf
  bsp.p <- ggplot(subset(d.bsp.scaled, d.bsp.scaled$Time <= max(d.mig[d.mig$l>4, 'bb']))) + 
    #geom_line(aes(Time, Median/gt)) + 
    geom_line(aes(Time, Median)) + 
    #geom_ribbon(aes(x=Time, ymax=Upper/gt, ymin=Lower/gt), alpha=.2) + 
    geom_ribbon(aes(x=Time, ymax=Upper, ymin=Lower), alpha=.2) + 
    scale_color_manual(values = c("gray")) + 
    scale_fill_manual(values = c("gray")) + 
    xlab("") + 
    ylab("Ne * Gen Time") + 
    scale_x_reverse(limits=c(2000, 0), 
                    breaks=seq(500,2000, by=500)) + 
    scale_y_continuous(labels=fancy_scientific) +
    theme_bw() +
    theme(panel.grid.major.y = element_line(linetype="dotted", colour="grey"),
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()#,
          #axis.title.y = element_text(angle = 0)
          )
  return(bsp.p)
}  
gen.full.bsp.p <- function(d.bsp, sf, d.mig) {
  #gt <- 62/(365/7) #62 weeks -> years
  d.bsp.scaled <- d.bsp * sf
  bsp.p <- ggplot(subset(d.bsp.scaled, d.bsp.scaled$Time <= max(d.mig[d.mig$l>4, 'bb']))) + 
    #geom_line(aes(Time, Median/gt)) + 
    geom_line(aes(Time, Median)) + 
    #geom_ribbon(aes(x=Time, ymax=Upper/gt, ymin=Lower/gt), alpha=.2) + 
    geom_ribbon(aes(x=Time, ymax=Upper, ymin=Lower), alpha=.2) + 
    scale_color_manual(values = c("gray")) + 
    scale_fill_manual(values = c("gray")) + 
    xlab("") + 
    ylab("Ne * Gen Time") + 
    scale_x_reverse(limits=c(max(d.mig[d.mig$l>4, 'bb']), 0), 
                    breaks=seq(500,max(d.mig[d.mig$l>4, 'bb']), by=500)) + 
    scale_y_continuous(labels=fancy_scientific) +
    theme_bw() +
    theme(panel.grid.major.y = element_line(linetype="dotted", colour="grey"),
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()#,
          #axis.title.y = element_text(angle = 0)
          )
  return(bsp.p)
}  

#BSP plot (log scale)
gen.bsp.log.p <- function(d.bsp, sf, d.mig) {
  #gt <- 62/(365/7) #62 weeks -> years
  d.bsp.scaled <- d.bsp * sf
  bsp.log.p <- ggplot(subset(d.bsp.scaled, d.bsp.scaled$Time <= max(d.mig[d.mig$l>4, 'bb']))) + 
    #geom_line(aes(Time, Median/gt)) + 
    geom_line(aes(Time, Median)) + 
    #geom_ribbon(aes(x=Time, ymax=Upper/gt, ymin=Lower/gt), alpha=.2) + 
    geom_ribbon(aes(x=Time, ymax=Upper, ymin=Lower), alpha=.2) + 
    scale_color_manual(values = c("gray")) + 
    scale_fill_manual(values = c("gray")) + 
    xlab("") + 
    ylab("Ne * Gen Time") + 
    scale_x_reverse(limits=c(2000, 0), 
                    breaks=seq(500,2000, by=500)) + 
    scale_y_continuous(trans = 'log10',
                     limits = c(40, 300000),
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid.major.y = element_line(linetype="dotted", colour="grey"),
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()#,
          #axis.title.y = element_text(angle = 0)
          )
  return(bsp.log.p)
}
gen.full.bsp.log.p <- function(d.bsp, sf, d.mig) {
  #gt <- 62/(365/7) #62 weeks -> years
  d.bsp.scaled <- d.bsp * sf
  bsp.log.p <- ggplot(subset(d.bsp.scaled, d.bsp.scaled$Time <= max(d.mig[d.mig$l>4, 'bb']))) + 
    #geom_line(aes(Time, Median/gt)) + 
    geom_line(aes(Time, Median)) + 
    #geom_ribbon(aes(x=Time, ymax=Upper/gt, ymin=Lower/gt), alpha=.2) + 
    geom_ribbon(aes(x=Time, ymax=Upper, ymin=Lower), alpha=.2) + 
    scale_color_manual(values = c("gray")) + 
    scale_fill_manual(values = c("gray")) + 
    xlab("") + 
    ylab("Ne * Gen Time") +  
    scale_x_reverse(limits=c(max(d.mig[d.mig$l>4, 'bb']), 0), 
                    breaks=seq(500,max(d.mig[d.mig$l>4, 'bb']), by=500)) + 
    scale_y_continuous(trans = 'log10',
                       limits = c(40, 250000),
                       breaks = trans_breaks('log10', function(x) 10^x),
                       labels = trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid.major.y = element_line(linetype="dotted", colour="grey"),
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()#,
          #axis.title.y = element_text(angle = 0)
          )
  return(bsp.log.p)
}

#Pie Chart Tree
gen.pie.tree <- function(tree, d, sf) {
  #subset the tips
  tre_tips <- subset(d, isTip == TRUE, c(label, Geography))
  #split the tips by UN
  tre_iso = split(tre_tips$label, tre_tips$Geography)
  tre.g <- groupOTU(tree, tre_iso)
  tre.g@phylo$edge.length <- tre.g@phylo$edge.length * sf

  t <- ggtree(tre.g, mrsd="2000-01-01", ladderize=FALSE) + geom_tiplab(aes(color=group)) + theme_tree2()
  
  d <- d[d$isTip == FALSE,]
  d <- d[,c('node', 'Geography.set', 'Geography.set.prob', 'posterior')]
  d <- d[d$posterior >= 0.8,] #remove nodes with posterior prob limit < 0.8 

  d$Geography.set.prob = gsub("\\{|\\}", "", d$Geography.set.prob)
  d$Geography.set = gsub("\\{|\\}", "", d$Geography.set)

  dlong = cSplit(d, splitCols = c("Geography.set", "Geography.set.prob"), sep = c(",", ","), direction = "long")
  dlong <- data.frame(dlong)
  dlong <- dlong[,-4]

  d.wide = as.data.frame(dcast(dlong, node ~ Geography.set))
  d.wide[is.na(d.wide)] <- 0
 
  pies <- nodepie(d.wide, cols=2:length(d.wide), alpha=.6)
  tp <- inset(t, pies, width=0.05)
  tp
  return(tp)
}

#Make migration matrix
make.mig.mat <- function(rates, bayes) {
  names(rates) <- c('statistic', 'mean', 'stdErr', 'median', 'hpdLower', 'hpdUpper', 'ESS', '50hpdLower', '50hpdUpper')
  geo.rates <- rates[grep("Geography.rates", rates$statistic),]
  geo.rates$statistic <- as.character(geo.rates$statistic)
  geo.rates$FROM <- sapply(strsplit(geo.rates$statistic, "[.]"), "[[", 3)
  geo.rates$TO <- sapply(strsplit(geo.rates$statistic, "[.]"), "[[", 4)
  
  geo.rates$comb <- paste(geo.rates$FROM, geo.rates$TO, sep="-")
  bayes$comb <- paste(bayes$FROM, bayes$TO, sep="-")
  
  rate.dat <- merge(geo.rates, bayes, by = "comb")
  rate.dat$mean[rate.dat$BAYES_FACTOR <= 5] <- NA
  rate.dat$median[rate.dat$BAYES_FACTOR <= 5] <- NA

  return(rate.dat)
}

#Migration Matrix Plot
gen.mig.mat.p <- function(rates.comb){  
  mig.mat.p <- ggplot(rates.comb, aes(x=TO.x, y=FROM.x, fill=median)) + 
    geom_tile(colour="gray") +
    scale_fill_viridis(na.value = "transparent") +
    #xlab(NULL) + ylab(NULL) +
    theme_bw(base_size=12) +
    xlab("") +
    ylab("") +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_text(angle=45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) 
  return(mig.mat.p)
}


#####################################################

### Full tree
#Read in the MCC tree
tree <- read.beast("data/OldWorld_GTR_G_BSP_DiscreteUN_MCC.tree")

#Group tree by lineages

tree <- groupClade(tree, node=MRCA(tree, tree@phylo$tip.label[grep("L1", tree@phylo$tip.label)]), "L1")
tree <- groupClade(tree, node=MRCA(tree, tree@phylo$tip.label[grep("L2", tree@phylo$tip.label)]), "L2")
tree <- groupClade(tree, node=MRCA(tree, tree@phylo$tip.label[grep("L3", tree@phylo$tip.label)]), "L3")
tree <- groupClade(tree, node=MRCA(tree, tree@phylo$tip.label[grep("L4", tree@phylo$tip.label)]), "L4")
tree <- groupClade(tree, node=MRCA(tree, tree@phylo$tip.label[grep("L5", tree@phylo$tip.label)]), "L5")
tree <- groupClade(tree, node=MRCA(tree, tree@phylo$tip.label[grep("L6", tree@phylo$tip.label)]), "L6")
tree <- groupClade(tree, node=MRCA(tree, tree@phylo$tip.label[grep("L7", tree@phylo$tip.label)]), "L7")

#Read full BSP data
full.bsp <- read.table("data/OldWorld_GTR_G_BSP_DiscreteUN_BSP.csv", skip=1, header=T)

#Convert tree data to dataframe
full <- fortify(tree)

#Modify the tip labels to include UN regions
tree <- mod.tip.labs(tree, full)

#Generate dataframe with updated tip labels
full <- fortify(tree)

#Define the scale factor based on the Mummy Mtb rate
full.sf <-  1/(5e-8*(3838249/60787))

#Make data readable
full.mod <- mod.tree.dat(full, full.sf)

#Drop Unkown samps
full.mod <- subset(full.mod, full.mod$Geography != "Unknown")

#Migration for the full tree
full.mig <- lin.mig.dat(full.mod)

#migration plots
full.mig.p <- gen.full.mig.p(full.mig)
full.gr.p <- gen.full.gr.p(full.mig)

#bsp plots
full.bsp.p <- gen.full.bsp.p(full.bsp,full.sf,full.mig)

#bsp plots log scale
full.bsp.log.p <- gen.full.bsp.log.p(full.bsp,full.sf,full.mig)

#range exp plot
full.gr.p <- gen.full.gr.p(full.mig)

#Figure - Old World
fig.full <- plot_grid(
  full.bsp.log.p,
  full.mig.p,
  full.gr.p,
  ncol=1,
  align="v"
)
fig.full

#save the figure
#ggsave(paste("figs/bsp_mig_range/", format(Sys.time(), "%y-%m-%d_"), "full__bsp_mig_range.pdf", sep = ""), plot=fig.full, width=9, height=6.75)

#pie chart trees
full.pie.p <- gen.pie.tree(tree, full, full.sf)

#Save the pie chart trees
#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "full.pie.pdf", sep = ""), plot=full.pie.p, width=14, height=79.3*1.5, limitsize=FALSE)

#Extract migration info for each lineage from the full tree
full.l1.mod <- full.mod[full.mod$L1 == "1",]
full.l1.mig <- lin.mig.dat(full.l1.mod)
full.l1.mig.p <- gen.mig.p(full.l1.mig)
full.l1.gr.p <- gen.gr.p(full.l1.mig)

full.l2.mod <- full.mod[full.mod$L2 == "1",]
full.l2.mig <- lin.mig.dat(full.l2.mod)
full.l2.mig.p <- gen.mig.p(full.l2.mig)
full.l2.gr.p <- gen.gr.p(full.l2.mig)

full.l3.mod <- full.mod[full.mod$L3 == "1",]
full.l3.mig <- lin.mig.dat(full.l3.mod)
full.l3.mig.p <- gen.mig.p(full.l3.mig)
full.l3.gr.p <- gen.gr.p(full.l3.mig)

full.l4.mod <- full.mod[full.mod$L4 == "1",]
full.l4.mig <- lin.mig.dat(full.l4.mod)
full.l4.mig.p <- gen.mig.p(full.l4.mig)
full.l4.gr.p <- gen.gr.p(full.l4.mig)

full.l5.mod <- full.mod[full.mod$L5 == "1",]
full.l5.mig <- lin567.mig.dat(full.l5.mod)

full.l6.mod <- full.mod[full.mod$L6 == "1",]
full.l6.mig <- lin567.mig.dat(full.l6.mod)

full.l7.mod <- full.mod[full.mod$L7 == "1",]
full.l7.mig <- lin567.mig.dat(full.l7.mod)

#TMRCA's for individual lineages
l1.mrca <- full[full$node == MRCA(tree, tree@phylo$tip.label[grep("L1", tree@phylo$tip.label)]), 'height_median']*full.sf
l2.mrca <- full[full$node == MRCA(tree, tree@phylo$tip.label[grep("L2", tree@phylo$tip.label)]), 'height_median']*full.sf
l3.mrca <- full[full$node == MRCA(tree, tree@phylo$tip.label[grep("L3", tree@phylo$tip.label)]), 'height_median']*full.sf
l4.mrca <- full[full$node == MRCA(tree, tree@phylo$tip.label[grep("L4", tree@phylo$tip.label)]), 'height_median']*full.sf
l5.mrca <- full[full$node == MRCA(tree, tree@phylo$tip.label[grep("L5", tree@phylo$tip.label)]), 'height_median']*full.sf
l6.mrca <- full[full$node == MRCA(tree, tree@phylo$tip.label[grep("L6", tree@phylo$tip.label)]), 'height_median']*full.sf
l7.mrca <- full[full$node == MRCA(tree, tree@phylo$tip.label[grep("L7", tree@phylo$tip.label)]), 'height_median']*full.sf

##BSP
#Read BSP data
l1.bsp.m <- read.table("data/lin1_UN_BSP.data", skip=1, header=T)
l2.bsp.m <- read.table("data/lin2_UN_BSP.data", skip=1, header=T)
l3.bsp.m <- read.table("data/lin3_UN_BSP.data", skip=1, header=T)
l4.bsp.m <- read.table("data/lin4_UN_BSP.data", skip=1, header=T)
l5.bsp.m <- read.table("data/lin5_UN_BSP.data", skip=1, header=T)
l6.bsp.m <- read.table("data/lin6_UN_BSP.data", skip=1, header=T)
l7.bsp.m <- read.table("data/lin7_UN_BSP.data", skip=1, header=T)


#Read ind lin trees
l1.tree.m <- read.beast("data/lin1_MCC.tree")
l2.tree.m <- read.beast("data/lin2_MCC.tree")
l3.tree.m <- read.beast("data/lin3_MCC.tree")
l4.tree.m <- read.beast("data/lin4_MCC.tree")
l5.tree.m <- read.beast("data/lin5_MCC.tree")
l6.tree.m <- read.beast("data/lin6_MCC.tree")
l7.tree.m <- read.beast("data/lin7_MCC.tree")

l1.m <- fortify(l1.tree.m)
l2.m <- fortify(l2.tree.m)
l3.m <- fortify(l3.tree.m)
l4.m <- fortify(l4.tree.m)
l5.m <- fortify(l5.tree.m)
l6.m <- fortify(l6.tree.m)
l7.m <- fortify(l7.tree.m)

l1.sf.m <- scale_factor(l1.m, l1.mrca)
l2.sf.m <- scale_factor(l2.m, l2.mrca)
l3.sf.m <- scale_factor(l3.m, l3.mrca)
l4.sf.m <- scale_factor(l4.m, l4.mrca)
l5.sf.m <- scale_factor(l5.m, l5.mrca)
l6.sf.m <- scale_factor(l6.m, l6.mrca)
l7.sf.m <- scale_factor(l7.m, l7.mrca)

#bsp plots
l1.bsp.p.m <- gen.bsp.p(l1.bsp.m,l1.sf.m,full.l1.mig)
l2.bsp.p.m <- gen.bsp.p(l2.bsp.m,l2.sf.m,full.l2.mig)
l3.bsp.p.m <- gen.bsp.p(l3.bsp.m,l3.sf.m,full.l3.mig)
l4.bsp.p.m <- gen.bsp.p(l4.bsp.m,l4.sf.m,full.l4.mig)
l5.bsp.p.m <- gen.bsp.p(l5.bsp.m,l5.sf.m,full.l5.mig)
l6.bsp.p.m <- gen.bsp.p(l6.bsp.m,l6.sf.m,full.l6.mig)
l7.bsp.p.m <- gen.bsp.p(l7.bsp.m,l7.sf.m,full.l7.mig)

#bsp plots log scale
l1.bsp.log.p.m <- gen.bsp.log.p(l1.bsp.m,l1.sf.m,full.l1.mig)
l2.bsp.log.p.m <- gen.bsp.log.p(l2.bsp.m,l2.sf.m,full.l2.mig)
l3.bsp.log.p.m <- gen.bsp.log.p(l3.bsp.m,l3.sf.m,full.l3.mig)
l4.bsp.log.p.m <- gen.bsp.log.p(l4.bsp.m,l4.sf.m,full.l4.mig)
l5.bsp.log.p.m <- gen.bsp.log.p(l5.bsp.m,l5.sf.m,full.l5.mig)
l6.bsp.log.p.m <- gen.bsp.log.p(l6.bsp.m,l6.sf.m,full.l6.mig)
l7.bsp.log.p.m <- gen.bsp.log.p(l7.bsp.m,l7.sf.m,full.l7.mig)

fig <- plot_grid(
  l1.bsp.log.p.m,
  l2.bsp.log.p.m,
  l3.bsp.log.p.m,
  l4.bsp.log.p.m,
  full.l1.mig.p,
  full.l2.mig.p,
  full.l3.mig.p,
  full.l4.mig.p,
  full.l1.gr.p,
  full.l2.gr.p,
  full.l3.gr.p,
  full.l4.gr.p,
  ncol=4,
  align="v"
)


supfig <- plot_grid(
  full.bsp.log.p,
  l5.bsp.log.p.m,
  full.mig.p,
  l6.bsp.log.p.m,
  full.gr.p,
  l7.bsp.log.p.m,
  ncol=2,
  align="v"
)
supfig

#save figure
#ggsave(paste("figs/bsp_mig_range/", format(Sys.time(), "%y-%m-%d_"), "lin_bsp_mig_range.pdf", sep = ""), plot=fig, width=18, height=6)

#ggsave(paste("figs/bsp_mig_range/", format(Sys.time(), "%y-%m-%d_"), "SI_bsp_mig_range.pdf", sep = ""), plot=supfig, width=9, height=6)


#pie chart trees
names(l1.m) <- gsub("geography", "Geography", names(l1.m))
names(l2.m) <- gsub("geography", "Geography", names(l2.m))
names(l3.m) <- gsub("geography", "Geography", names(l3.m))
names(l4.m) <- gsub("geography", "Geography", names(l4.m))

l1.pie.p <- gen.pie.tree(l1.tree.m, l1.m, l1.sf.m)
l2.pie.p <- gen.pie.tree(l2.tree.m, l2.m, l2.sf.m)
l3.pie.p <- gen.pie.tree(l3.tree.m, l3.m, l3.sf.m)
l4.pie.p <- gen.pie.tree(l4.tree.m, l4.m, l4.sf.m)

#Save the pie chart trees
#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "lin1.mean.pie.pdf", sep = ""), plot=l1.pie.p, width=14, height=14.8*1.5)
#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "lin2.mean.pie.pdf", sep = ""), plot=l2.pie.p, width=14, height=30*1.5)
#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "lin3.mean.pie.pdf", sep = ""), plot=l3.pie.p, width=14, height=10.8*1.5)
#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "lin4.mean.pie.pdf", sep = ""), plot=l4.pie.p, width=14, height=23.7*1.5)

#country-level pie charts
l1.co.tree <- read.beast("data/lin1_iso2_MCC.tree")
l3.co.tree <- read.beast("data/lin3_iso2_MCC.tree")
l2.co.tree <- read.beast("data/lin2_iso2_MCC.tree")
l4.co.tree <- read.beast("data/lin4_iso2_MCC.tree")

l1.co <- fortify(l1.co.tree)
names(l1.co) <- gsub("location", "Geography", names(l1.co))
l3.co <- fortify(l3.co.tree)
names(l3.co) <- gsub("location", "Geography", names(l3.co))
l2.co <- fortify(l2.co.tree)
names(l2.co) <- gsub("location", "Geography", names(l2.co))
l4.co <- fortify(l4.co.tree)
names(l4.co) <- gsub("location", "Geography", names(l4.co))

#Get UN region for tip labels
meta <- read.table("data/TableS1_StrainInfo.txt", header=T, sep='\t', na.strings="NA", stringsAsFactors=FALSE)

#Modify tip labels
mod.co.tip.labs <- function(tree, d) {
  tre_tips <- subset(d, isTip == TRUE, c(label, Geography))
  tre_tips$samp <- sapply(strsplit(tre_tips$label, "_"), "[[", 1)
  tre_tips$iso2 <- meta[match(tre_tips$samp, meta$BioSample), 'Iso2']
  tre_tips$lin <- paste("L", meta[match(tre_tips$samp, meta$BioSample), 'Lineage'], sep="")
  tre_tips$UN <- meta[match(tre_tips$samp, meta$BioSample), 'UN']
  tre_tips$newName <- paste(tre_tips$samp, tre_tips$lin, tre_tips$iso2, tre_tips$UN, sep="_")
  tree@phylo$tip.label <- tre_tips$newName
  return(tree)
}

l1.co.tree <- mod.co.tip.labs(l1.co.tree, l1.co)
l1.co <- fortify(l1.co.tree)
names(l1.co) <- gsub("location", "Geography", names(l1.co))

l2.co.tree <- mod.co.tip.labs(l2.co.tree, l2.co)
l2.co <- fortify(l2.co.tree)
names(l2.co) <- gsub("location", "Geography", names(l2.co))

l3.co.tree <- mod.co.tip.labs(l3.co.tree, l3.co)
l3.co <- fortify(l3.co.tree)
names(l3.co) <- gsub("location", "Geography", names(l3.co))

l4.co.tree <- mod.co.tip.labs(l4.co.tree, l4.co)
l4.co <- fortify(l4.co.tree)
names(l4.co) <- gsub("location", "Geography", names(l4.co))

#scale factor
l1.co.sf <- scale_factor(l1.co, l1.mrca)
l3.co.sf <- scale_factor(l3.co, l3.mrca)
l2.co.sf <- scale_factor(l2.co, l2.mrca)
l4.co.sf <- scale_factor(l4.co, l4.mrca)

#mod data frame
l1.co.mod <- mod.tree.dat(l1.co, l1.co.sf)
l3.co.mod <- mod.tree.dat(l3.co, l3.co.sf)
l2.co.mod <- mod.tree.dat(l2.co, l2.co.sf)
l4.co.mod <- mod.tree.dat(l4.co, l4.co.sf)

l1.co.pie.p <- gen.pie.tree(l1.co.tree, l1.co, l1.co.sf)
l3.co.pie.p <- gen.pie.tree(l3.co.tree, l3.co, l3.co.sf)
l2.co.pie.p <- gen.pie.tree(l2.co.tree, l2.co, l2.co.sf)
l4.co.pie.p <- gen.pie.tree(l4.co.tree, l4.co, l4.co.sf)

#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "lin1.co.pie.pdf", sep = ""), plot=l1.co.pie.p, width=14, height=14.8*1.5)
#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "lin3.co.pie.pdf", sep = ""), plot=l3.co.pie.p, width=14, height=10.8*1.5)
#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "lin2.co.pie.pdf", sep = ""), plot=l2.co.pie.p, width=14, height=30*1.5)
#ggsave(paste("figs/MCC/", format(Sys.time(), "%y-%m-%d_"), "lin4.co.pie.pdf", sep = ""), plot=l4.co.pie.p, width=14, height=23.7*1.5)


#Migration rate matrix
rates <- as.matrix(read.table("data/OldWorld_GTR_G_BSP_DiscreteUN_rates.txt",header=T, sep="\t", row.names=1, na.strings=c("-","NA")))
sigs <- read.table("data/OldWorld_GTR_G_BSP_DiscreteUN_sig.txt", header=T, sep="\t")

rates.mod <- melt(rates)
names(rates.mod) <- c("var1", "var2", "median")

rates.mod.2 <- subset(rates.mod, !is.na(rates.mod$median))

rates.mod.2$comb <- paste(rates.mod.2$var2, rates.mod.2$var1, sep="-")
sigs$comb <- paste(sigs$FROM, sigs$TO, sep="-")

rate.dat <- merge(rates.mod.2, sigs, by = "comb")
rate.dat.c <- rate.dat
rate.dat.c$median[rate.dat.c$BAYES_FACTOR <= 5] <- NA

#####
#generate a network file to be exported for relative rate figure
net <- rate.dat.c[!is.na(rate.dat.c$median), c('var1', 'median', 'var2')]
net$median <- net$median*2

#write.table(net, file = "data/full_net.sif", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", row.names=FALSE, col.names=FALSE, qmethod = "double")
#####

#Lineage-specific Migration rate matrices
#Read the data
l1.bayes <- read.table("data/Lineage1.SpreaD3.output.txt", header=T, sep="\t")
l1.rates <- read.table("data/Lineage1.Geography.rates.txt", skip = 2, header=F, sep="\t")

l2.bayes <- read.table("data/Lineage2.Spread3.output.txt", header=T, sep="\t")
l2.rates <- read.table("data/Lineage2.Geography.rates.txt", skip = 2, header=F, sep="\t")

l3.bayes <- read.table("data/Lineage3.SpreaD3.output.txt", header=T, sep="\t")
l3.rates <- read.table("data/Lineage3.Geography.rates.txt", skip = 2, header=F, sep="\t")

l4.bayes <- read.table("data/Lineage4.SpreaD3.output.txt", header=T, sep="\t")
l4.rates <- read.table("data/Lineage4.Geography.rates.txt", skip = 2, header=F, sep="\t")

#Make the df for plotting
l1.comb <- make.mig.mat(l1.rates, l1.bayes)
l2.comb <- make.mig.mat(l2.rates, l2.bayes)
l3.comb <- make.mig.mat(l3.rates, l3.bayes)
l4.comb <- make.mig.mat(l4.rates, l4.bayes)

net1 <- l1.comb[!is.na(l1.comb$median), c('TO.x', 'median', 'FROM.x')]
net1$median <- net1$median*2
#write.table(net1, file = "data/lin1_net.sif", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", row.names=FALSE, col.names=FALSE, qmethod = "double")

net2 <- l2.comb[!is.na(l2.comb$median), c('TO.x', 'median', 'FROM.x')]
net2$median <- net2$median*2
#write.table(net2, file = "data/lin2_net.sif", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", row.names=FALSE, col.names=FALSE, qmethod = "double")

net3 <- l3.comb[!is.na(l3.comb$median), c('TO.x', 'median', 'FROM.x')]
net3$median <- net3$median*2
#write.table(net3, file = "data/lin3_net.sif", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", row.names=FALSE, col.names=FALSE, qmethod = "double")

net4 <- l4.comb[!is.na(l4.comb$median), c('TO.x', 'median', 'FROM.x')]
net4$median <- net4$median*2
#write.table(net4, file = "data/lin4_net.sif", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", row.names=FALSE, col.names=FALSE, qmethod = "double")

#Plot mig mat
l1.mig.mat.p <- gen.mig.mat.p(l1.comb)
l2.mig.mat.p <- gen.mig.mat.p(l2.comb)
l3.mig.mat.p <- gen.mig.mat.p(l3.comb)
l4.mig.mat.p <- gen.mig.mat.p(l4.comb)
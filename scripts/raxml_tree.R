#Maximum likelihood tree (SI)
require(ggplot2)
require(ggtree)
require(cowplot)
require(scales)
require(splitstackshape)
require(reshape2)
require(dplyr)
require(viridis)
require(ape)
require(phytools)
require(phangorn)


#Read the RAxML tree
rax.tree <- read.tree("data/RAxML.tree")

#Read meta info for samples
meta <- read.table("data/Dataset1.txt", header=T, sep='\t', na.strings="NA", stringsAsFactors=FALSE)

#devise groups
raxUN <- split(rax.tree$tip.label, meta[match(lapply(strsplit(as.character(rax.tree$tip.label), "_"), "[", 1), meta$BioSample), 'UN'])

#midpoint root
raxmpt <- midpoint(rax.tree)

#group tree by UN subregions
raxmpt <- groupOTU(raxmpt, raxUN)

#Get data to color branch by UN region
dat <- fortify(rax.tree)
dat$UN <- meta[match(lapply(strsplit(as.character(dat$label), "_"), "[", 1), meta$BioSample), 'UN']

dat$UN <- as.character(dat$UN)

index <- dat$isTip == FALSE
dat$UN[index] <- "internal"

dat$UN <- as.factor(dat$UN)

p <- ggtree(raxmpt, layou="circular") + geom_treescale(width=3e-05, offset=-10) 

q <- p %<+% dat + aes(color=(UN)) + theme(legend.position="right") + 
  scale_color_manual(name = "UN Geographic Region", values=c('#F8766D','#E88526','#D39200','#B79F00','#666666','#93AA00','#5EB300','#00BF74','#00C19F','#00BFC4','#00B9E3','#00ADFA','#619CFF','#AE87FF'), labels = c('Central Asia','Eastern Africa','Eastern Asia', 'Eastern Europe', 'internal', 'Melanesia', 'Micronesia', 'South Eastern Asia', 'Southern-Africa', 'Southern Asia', 'Southern Europe', 'Western Africa', 'Western Asia', 'Western Europe')) 




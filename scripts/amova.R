# AMOVA

# Load packages
library(pegas)
library(ape)
library(poppr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

# Set the directory
setwd("~/Mtb_Phylogeography")

# Read the SNP alignment
Mtb <- read.FASTA("data/snps_50per_552strains.fa") #read the fasta

# Computer pairwise distance matrix
genetic.dists = dist.dna(Mtb) 
gm <- as.matrix(genetic.dists)
rownames(gm) <- sapply(strsplit(rownames(gm), "_"), "[[", 1)
colnames(gm) <- sapply(strsplit(colnames(gm), "_"), "[[", 1)

# Load meta data file containing the AMOVA bins
meta <- read.table("data/TableS1.txt", header=T, sep="\t", na.strings="NA")

# Remove any isolates for which there is no meta info
samps <- intersect(rownames(gm),meta$BioSample)

gm <- gm[rownames(gm) %in% samps, colnames(gm) %in% samps]
gd <- as.dist(gm)

# Ordering for categories
samp_ord <- labels(gd)

# Generate the categories
l <- as.factor(meta[match(as.factor(samp_ord), meta$BioSample), 'Lineage']) #Lineage
lev1 <- meta[match(as.factor(samp_ord), meta$BioSample), 'Level1'] #Botanical Level 1
un <- meta[match(as.factor(samp_ord), meta$BioSample), 'UN'] #UN region

# Convert the DNA bin to a genind object
Mtb.genind <- DNAbin2genind(Mtb)

# Add categories
sample <- as.factor(samp_ord)
strata(Mtb.genind) <- data.frame(sample)
addStrata(Mtb.genind) <- data.frame(l)
addStrata(Mtb.genind) <- data.frame(lev1)
addStrata(Mtb.genind) <- data.frame(un)

# Convert the genind object to genclone object
Mtb.genclone <- as.genclone(Mtb.genind)

subsamps <- read.table("data/biosamples_subsampling.txt", header=F, sep="\t", stringsAsFactors = FALSE)
subsamps <- data.frame(t(subsamps))
colnames(subsamps) <- paste("sub", seq(1,10), sep="")

# Read in the names of random subsamples of UN subregions (n=20 per region)
subsamps <- read.table("data/biosamples_subsampling.txt", header=F, sep="\t", stringsAsFactors = FALSE)
subsamps <- data.frame(t(subsamps))
colnames(subsamps) <- paste("sub", seq(1,10), sep="")
head(subsamps)

# Create function to generate amovas and sigs
perform_amova <- function(genclone) {
  set.seed(1999)
  
  lev1.am <- poppr.amova(genclone, ~lev1)
  lev1.am.sig <- randtest(lev1.am, nrepet = 1000)
  
  un.am <- poppr.amova(genclone, ~un)
  un.am.sig <- randtest(un.am, nrepet = 1000)
  
  return(list(lev1.am, lev1.am.sig, un.am, un.am.sig))
}

# Run AMOVA for entire collection and each lineage with enough isolates
# lineages 1-4
amovas.list <- c()
for (i in 1:4) {
  subsample <- subsamps[,i]
  subX <- Mtb.genclone %>% setPop(~l) %>% popsub(paste(i)) %>% missingno("loci") %>% missingno("geno")
  subX.amovas <- perform_amova(subX)
  amovas.list[[i]] <- subX.amovas
}
names(amovas.list) <- paste("L", seq(1,4), sep="")

# Old world collection
amovas.list[["OldWorld"]]

# Old World collection and random sub-samples at the UN level
OW.amovas.list <- c()
for (i in 1:10) {
  subsample <- subsamps[,i]
  subX <- Mtb.genclone %>% setPop(~sample) %>% popsub(sublist = subsample) %>% missingno("loci") %>% missingno("geno")
  subX.amovas <- perform_amova(subX)
  OW.amovas.list[[i]] <- subX.amovas
}
names(OW.amovas.list) <- paste("sub", seq(1,10), sep="")

OW.amovas.list[["full"]] <- perform_amova((Mtb.genclone %>% missingno("loci") %>% missingno("geno")))

# Extract information of interest
OW <- rbind(rbind(names(OW.amovas.list),
            sapply(seq(1,11), function(x) OW.amovas.list[[x]][[3]]$componentsofcovariance$`%`),
            sapply(seq(1,11), function(x) OW.amovas.list[[x]][[4]]$pvalue)) %>% t() %>% data.frame() %>% dplyr::mutate(Test="UN"),
            rbind(names(OW.amovas.list),
            sapply(seq(1,11), function(x) OW.amovas.list[[x]][[1]]$componentsofcovariance$`%`),
            sapply(seq(1,11), function(x) OW.amovas.list[[x]][[2]]$pvalue)) %>% t() %>% data.frame() %>% dplyr::mutate(Test="Bot")) 
names(OW) <- c("Sample", "Between", "Within", "Total", "pvalue", "Test")


LL <- rbind(rbind(names(amovas.list),
            sapply(seq(1,4), function(x) amovas.list[[x]][[3]]$componentsofcovariance$`%`),
            sapply(seq(1,4), function(x) amovas.list[[x]][[4]]$pvalue)) %>% t() %>% data.frame() %>% dplyr::mutate(Test="UN"),
            rbind(names(amovas.list),
            sapply(seq(1,4), function(x) amovas.list[[x]][[1]]$componentsofcovariance$`%`),
            sapply(seq(1,4), function(x) amovas.list[[x]][[2]]$pvalue)) %>% t() %>% data.frame() %>% dplyr::mutate(Test="Bot"))
names(LL) <- c("Sample", "Between", "Within", "Total", "pvalue", "Test")


# And Plot
FigS3C <- ggplot() + 
  geom_boxplot(data=OW[grep("sub", OW$Sample),], aes(x=Test, y=as.numeric(as.character(Between)), fill=Test), alpha=0.5) +
  scale_fill_manual(values=c("#FF9999", "#99CCFF")) +
  geom_point(data=OW[grep("sub", OW$Sample),], aes(x=Test, y=as.numeric(as.character(Between)))) +
  geom_point(data=OW[OW$Sample=="full",], aes(x=Test, y=as.numeric(as.character(Between))), color="red") +
  theme_bw() +
  xlab("") +
  ylab("Variance Between Regions") +
  theme(panel.grid = element_blank(),
        legend.position="none")
FigS3C

# Save results to file
write.table(rbind(LL, OW), file="data/AMOVA_results.txt", quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, qmethod="double")

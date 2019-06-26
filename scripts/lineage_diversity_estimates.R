# Effect of subsampling on diversity estimates

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)

# Set the directory
setwd("~/Mtb_Phylogeography")

# Read in the data
ss <- read.table("/Users/moneill/Mtb_Phylogeography/data/lineage_diversity_estimates.txt", header=T, stringsAsFactors = F, na.strings=c("nan", "None", "-inf"))
str(ss)

ss <- tidyr::separate(ss, Alignment, into=c("sample", "lineage"), remove=FALSE)
table(ss$sample)

# Add the estimates from the Full Alignment (calculated previously)
fullSamp_all <- data.frame("fullSamp_all", "fullSamp", "all", 2.13e-03, 2.80e-04, NA)
names(fullSamp_all) <- colnames(ss)
ss <- rbind(ss, fullSamp_all)
table(ss$sample)

# Data wrangling
ss$lineage <- gsub("lin", "Lineage ", ss$lineage)
ss$lineage <- gsub("all", "Old World Collection", ss$lineage)

ss.m <- melt(ss, id.vars=c("Alignment", "sample", "lineage"))
ss.m$col <- as.factor(ifelse(substr(ss.m$sample, 1, 3) == "sub", 0, 1))

ss.m$type <- as.factor(substr(ss.m$sample, 1, 3))

means <- ss.m %>% 
  group_by(type, variable, lineage) %>%
  summarise(mean=mean(value))

means.wide <- spread(means, type, mean)

full <- filter(ss.m, variable != "TajimasD" & type == "ful") %>% select(variable, lineage, sample, value)
comp <- filter(ss.m, variable != "TajimasD" & type == "sub") %>% select(variable, lineage, sample, value)
comp$full <- full[match(paste(comp$variable, comp$lineage, sep="-"), paste(full$variable, full$lineage, sep="-")), 'value']
   
sumComp <- comp %>% 
  group_by(variable, lineage) %>%
  summarise(mean=mean(value), min=min(value), max=max(value))
sumComp$full <- full[match(paste(sumComp$variable, sumComp$lineage, sep="-"), paste(full$variable, full$lineage, sep="-")), 'value']

# Plot!
FigS3A <- ggplot(filter(sumComp, variable=="Theta")) +
  geom_point(aes(x=full, y=mean, color=lineage)) +
  geom_errorbar(aes(x=full, ymin=min, ymax=max, color=lineage)) +
  scale_color_manual(values = c("pink", "blue", "purple", "red", "darkred", "darkgreen", "#FFCC33", "black")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab("Complete Sample") +
  ylab("Subsample") +
  ggtitle(expression("Theta"))

FigS3B <- ggplot(filter(sumComp, variable=="Pi")) +
  geom_point(aes(x=full, y=mean, color=lineage)) +
  geom_errorbar(aes(x=full, ymin=min, ymax=max, color=lineage)) +
  scale_color_manual(values = c("pink", "blue", "purple", "red", "darkred", "darkgreen", "#FFCC33", "black")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab("Complete Sample") +
  ylab("Subsample") +
  ggtitle(expression("Pi"))

plot_grid(FigS3A, FigS3B)


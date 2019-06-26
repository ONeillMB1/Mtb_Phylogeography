#Decline in Diversity

require(ape)
require(PopGenome)
require(ggplot2)
require(reshape2)

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

dec <- read.table("data/decline_in_diversity.txt", header=T, sep="\t")

dec.p <- ggplot(dec) +
  facet_wrap(~Lineage) +
  geom_point(aes(x=MeanDist, y=Pi)) +
  scale_y_continuous(labels=fancy_scientific) +
  ylab(expression(pi)) +
  xlab("Average Distance from Addis Ababa") +
  geom_point(aes(x=MedianDist, y=Pi, color = "red")) +
  theme(legend.position = "none")
dec.p  

#ggsave(paste("figs/decline_in_diversity/", format(Sys.time(), "%y-%m-%d_"), "dec_in_div.pdf", sep = ""), plot=dec.p, width=8, height=8)

lin1.mod.mean = lm(Pi~MeanDist, data=dec[dec$Lineage=="Lin1",])
summary(lin1.mod.mean)

lin1.mod.median = lm(Pi~MedianDist, data=dec[dec$Lineage=="Lin1",])
summary(lin1.mod.mean)

lin2.mod.mean = lm(Pi~MeanDist, data=dec[dec$Lineage=="Lin2",])
summary(lin2.mod.mean)

lin2.mod.median = lm(Pi~MedianDist, data=dec[dec$Lineage=="Lin2",])
summary(lin2.mod.mean)

lin3.mod.mean = lm(Pi~MeanDist, data=dec[dec$Lineage=="Lin3",])
summary(lin3.mod.mean)

lin3.mod.median = lm(Pi~MedianDist, data=dec[dec$Lineage=="Lin3",])
summary(lin3.mod.mean)

lin4.mod.mean = lm(Pi~MeanDist, data=dec[dec$Lineage=="Lin4",])
summary(lin4.mod.mean)

lin4.mod.median = lm(Pi~MedianDist, data=dec[dec$Lineage=="Lin4",])
summary(lin4.mod.mean)

lin5.mod.mean = lm(Pi~MeanDist, data=dec[dec$Lineage=="Lin5",])
summary(lin5.mod.mean)

lin5.mod.median = lm(Pi~MedianDist, data=dec[dec$Lineage=="Lin5",])
summary(lin5.mod.mean)

lin6.mod.mean = lm(Pi~MeanDist, data=dec[dec$Lineage=="Lin6",])
summary(lin6.mod.mean)

lin6.mod.median = lm(Pi~MedianDist, data=dec[dec$Lineage=="Lin6",])
summary(lin6.mod.mean)

lin7.mod.mean = lm(Pi~MeanDist, data=dec[dec$Lineage=="Lin7",])
summary(lin7.mod.mean)

lin7.mod.median = lm(Pi~MedianDist, data=dec[dec$Lineage=="Lin7",])
summary(lin7.mod.mean)

all.mod.mean = lm(Pi~MeanDist, data=dec[dec$Lineage=="All",])
summary(all.mod.mean)

all.mod.median = lm(Pi~MedianDist, data=dec[dec$Lineage=="All",])
summary(all.mod.mean)

# Mantel tests

require(geosphere)
require(ape)
require(vegan)
require(ggplot2)

# Set the directory
setwd("~/Mtb_Phylogeography")

# Load meta data file 
meta <- read.table("data/TableS1.txt", header=T, sep="\t", na.strings="NA", stringsAsFactors = FALSE)
meta <- meta[meta$Country != "Kiribati",] #Not on trade routes
samps <- meta$BioSample

############# Great Circle Distance between Isolates
# Extract coordinates
x <- meta[,c('sample_lon', 'sample_lat')]
rownames(x) <- meta$BioSample
y <- as.matrix(x)

# Calculate great circle distance between isolates
dvse <- apply(y, 1, FUN=function(X) distVincentyEllipsoid(X, y))
rownames(dvse) <- colnames(dvse)

#order the rows and columns
dvse <- dvse[order(rownames(dvse)), order(colnames(dvse))]

############# Distance Via Trade Routes
#Load trade route distances (calculated from NearPost asignments)
trade <- read.table("data/TradeRouteDist_from_NearPosts.txt", header=T, row.names=1, sep='\t')
dm.trade <- as.matrix(trade)

#define rules
calcNearPostDist <- function(lon1, lat1, hub1, nearLon1, nearLat1, lon2, lat2, hub2, nearLon2, nearLat2) {
  if (hub1 == hub2) {
    dd = distVincentyEllipsoid(c(lon1,lat1), c(lon2,lat2))/1000
  } else {
    dd = distVincentyEllipsoid(c(lon1,lat1), c(nearLon1,nearLat1))/1000 + distVincentyEllipsoid(c(lon2,lat2), c(nearLon2,nearLat2))/1000
  }
  return(dd)
}

#calc matrix function
near_dist_mat = function(d) {
  mat = matrix(0, ncol = nrow(d), nrow = nrow(d))
  for(i in 1:nrow(mat)) {
    for(j in 1:nrow(mat)) {
      mat[i,j] = calcNearPostDist(lon1 = d[i,'sample_lon'], lat1 = d[i,'sample_lat'], hub1 = d[i,'NearPost'], nearLon1 = d[i, 'NearPost_lon'], nearLat1 = d[i, 'NearPost_lat'], lon2 = d[j,'sample_lon'], lat2 = d[j,'sample_lat'], hub2 = d[j,'NearPost'], nearLon2 = d[j,'NearPost_lon'], nearLat2 = d[j, 'NearPost_lat'])
    }
  }
  return(mat)
}

nd.m <- near_dist_mat(meta)
rownames(nd.m) <- meta$BioSample
colnames(nd.m) <- meta$BioSample

nd.m <- nd.m[order(rownames(nd.m)), order(colnames(nd.m))]

#order matrices
nd.m <- nd.m[order(rownames(nd.m)), order(colnames(nd.m))]
dm.trade <- dm.trade[order(rownames(dm.trade)), order(colnames(dm.trade))]

identical(rownames(nd.m), rownames(dm.trade))
identical(colnames(nd.m), colnames(dm.trade))
identical(rownames(nd.m), colnames(dm.trade))

#Add the two together for the final distances via trade routes
trade.tot <- (nd.m + dm.trade) * 1000

########### Distances via Waypoints
#Waypoint coordinates
x$Cairo <- apply(y, 1, FUN=function(X) distVincentyEllipsoid(X, c(31,30)))
x$Istanbul <- apply(y, 1, FUN=function(X) distVincentyEllipsoid(X, c(28,41)))
x$PhnomPenh <- apply(y, 1, FUN=function(X) distVincentyEllipsoid(X, c(104,11)))
x$AddisAbada <- apply(y, 1, FUN=function(X) distVincentyEllipsoid(X, c(38.74,9.03)))

x$UN <- meta[match(rownames(x), meta$BioSample), 'UN']
x$continent <- sapply(strsplit(as.character(x$UN), "-"), tail, 1)

dCI <- distVincentyEllipsoid(c(31,30), c(28,41)) #Cairo to Istanbul
dCP <- distVincentyEllipsoid(c(31,30), c(104,11)) #Cairo to Phnom Penh
dIP <- distVincentyEllipsoid(c(28,41), c(104,11)) #Istanbul Phnom Penh

calcWayDistMODIFIED <- function(lon1, lat1, con1, lon2, lat2, con2) {
  if (con1 == con2) {
    dd = distVincentyEllipsoid(c(lon1,lat1), c(lon2,lat2))
  } else if (con1 == "Africa" & con2 == "Europe") {
    afr = distVincentyEllipsoid(c(lon1,lat1), c(31,30)) #dist from isolate to Cairo
    eur = distVincentyEllipsoid(c(lon2,lat2), c(28,41)) #dist from isolate to Istanbul
    dd = afr + dCI + eur
  } else if (con1 == "Europe" & con2 == "Africa") {
    afr = distVincentyEllipsoid(c(lon2,lat2), c(31,30)) #dist from isolate to Cairo
    eur = distVincentyEllipsoid(c(lon1,lat1), c(28,41)) #dist from isolate to Istanbul
    dd = afr + dCI + eur
  } else if (setequal(c("Africa", "Asia"), c(con1, con2))) {
    s1 = distVincentyEllipsoid(c(lon1,lat1), c(31,30)) #dist from isolate to Cairo
    s2 = distVincentyEllipsoid(c(lon2,lat2), c(31,30)) #dist from isolate to Cairo
    dd = s1 + s2
  } else if (con1 == "Africa" & con2 == "Melanesia") {
    afr = distVincentyEllipsoid(c(lon1,lat1), c(31,30)) #dist from isolate to Cairo
    mel = distVincentyEllipsoid(c(lon2,lat2), c(104,11)) #dist from isolate to Phnom Penh
    dd = afr + dCP + mel
  } else if (con1 == "Melanesia" & con2 == "Africa") {
    afr = distVincentyEllipsoid(c(lon2,lat2), c(31,30)) #dist from isolate to Cairo
    mel = distVincentyEllipsoid(c(lon1,lat1), c(104,11)) #dist from isolate to Phnom Penh
    dd = afr + dCP + mel
  } else if (setequal(c("Europe", "Asia"), c(con1, con2))) {
    s1 = distVincentyEllipsoid(c(lon1,lat1), c(28,41)) #dist from isolate to Cairo
    s2 = distVincentyEllipsoid(c(lon2,lat2), c(28,41)) #dist from isolate to Cairo
    dd = s1 + s2
  } else if (con1 == "Europe" & con2 == "Melanesia") {
    eur = distVincentyEllipsoid(c(lon1,lat1), c(28,41)) #dist from isolate to Istanbul
    mel = distVincentyEllipsoid(c(lon2,lat2), c(104,11)) #dist from isolate to Phnom Penh
    dd = eur + dIP + mel
  } else if (con1 == "Melanesia" & con2 == "Europe") {
    eur = distVincentyEllipsoid(c(lon2,lat2), c(28,41)) #dist from isolate to Istanbul
    mel = distVincentyEllipsoid(c(lon1,lat1), c(104,11)) #dist from isolate to Phnom Penh
    dd = eur + dIP + mel
  } else if (setequal(c("Asia", "Melanesia"), c(con1, con2))) {
    s1 = distVincentyEllipsoid(c(lon1,lat1), c(104,11)) #dist from isolate to Phnom Penh
    s2 = distVincentyEllipsoid(c(lon2,lat2), c(104,11)) #dist from isolate to Phnom Penh
    dd = s1 + s2
  }
  return(dd)
}

distance_function = function(d) {
  mat = matrix(0, ncol = nrow(d), nrow = nrow(d))
  for(i in 1:nrow(mat)) {
    for(j in 1:nrow(mat)) {
      mat[i,j] = calcWayDistMODIFIED(lon1 = d[i,'sample_lon'], lat1 = d[i,'sample_lat'], con1 = d[i,'continent'], lon2 = d[j,'sample_lon'], lat2 = d[j,'sample_lat'], con2 = d[j,'continent'])
    }
  }
  return(mat)
}


waypoint.mat <- distance_function(d = x)
rownames(waypoint.mat) <- rownames(x)
colnames(waypoint.mat) <- rownames(x)

#order matrix
waypoint.mat <- waypoint.mat[order(rownames(waypoint.mat)), order(colnames(waypoint.mat))]

######################## Comparisons

distances <- data.frame(dvse=as.vector(dvse[upper.tri(dvse)]), trade=as.vector(trade.tot[upper.tri(trade.tot)]), waypoint=as.vector(waypoint.mat[upper.tri(waypoint.mat)]))

dvse.waypoint.p <- ggplot(distances) + 
  geom_point(aes(dvse, waypoint))

dvse.trade.p <- ggplot(distances) + 
  geom_point(aes(dvse, trade))

trade.waypoint.p <- ggplot(distances) + 
  geom_point(aes(trade, waypoint))

###################### Genetic data

Mtb <- read.FASTA("data/Mantel_OldWorld_snpAln_50per_noKI.fasta") 
genetic.dists = dist.dna(Mtb) #pairwise distance matrix computation
gm <- as.matrix(genetic.dists)
rownames(gm) <- sapply(strsplit(rownames(gm), "_"), "[[", 1)
colnames(gm) <- sapply(strsplit(colnames(gm), "_"), "[[", 1)

gm <- gm[order(rownames(gm)), order(rownames(gm))]

###################### Genetics and physical distance

corr.data <- data.frame(dvse=as.vector(dvse[upper.tri(dvse)]), trade=as.vector(trade.tot[upper.tri(trade.tot)]), waypoint=as.vector(waypoint.mat[upper.tri(waypoint.mat)]), genetic=as.vector(gm[upper.tri(gm)]))

dvse.genetic.p <- ggplot(corr.data) +
  geom_point(aes(dvse, genetic))

dvse.genetic.mod = lm(dvse~genetic, data=corr.data)
summary(dvse.genetic.mod)

trade.genetic.p <- ggplot(corr.data) +
  geom_point(aes(trade, genetic))

trade.genetic.mod = lm(trade~genetic, data=corr.data)
summary(trade.genetic.mod)

waypoint.genetic.p <- ggplot(corr.data) + 
  geom_point(aes(waypoint, genetic))

waypoint.genetic.mod = lm(waypoint~genetic, data=corr.data)
summary(waypoint.genetic.mod)


############ Mantel
man.dvse.st.log <- mantel(gm, decostand(log(dvse+1), "standardize"), permutations=9999)
man.trade.st.log <- mantel(gm, decostand(log(trade.tot+1), "standardize"), permutations=9999)
man.waypoint.st.log <- mantel(gm, decostand(log(waypoint.mat+1), "standardize"), permutations=9999)

man.dvse.st <- mantel(gm, decostand(dvse, "standardize"), permutations=9999)
man.trade.st <- mantel(gm, decostand(trade.tot, "standardize"), permutations=9999)
man.waypoint.st <- mantel(gm, decostand(waypoint.mat, "standardize"), permutations=9999)

#### Subsets
subsamps <- read.table("data/biosamples_subsampling.txt", header=F, sep="\t", stringsAsFactors = FALSE)
subsamps <- data.frame(t(subsamps))
colnames(subsamps) <- paste("sub", seq(1,10), sep="")


ind_lin_mantel <- function(meta, l, gcd, td, wd, gm) {
  #lin.vec = as.vector(subset(meta, Lineage == l, BioSample))
  lin.vec = as.vector(data.frame(meta[,l]))
  names(lin.vec) <- c("BioSample")
  
  gcd <- gcd[rownames(gcd) %in% lin.vec$BioSample, colnames(gcd) %in% lin.vec$BioSample]
  td <- td[rownames(td) %in% lin.vec$BioSample, colnames(td) %in% lin.vec$BioSample]
  wd <- wd[rownames(wd) %in% lin.vec$BioSample, colnames(wd) %in% lin.vec$BioSample]
  gm <- gm[rownames(gm) %in% lin.vec$BioSample, colnames(gm) %in% lin.vec$BioSample]
  
  man.gcd <- mantel(decostand(gcd, "standardize"), gm, permutations=9999)
  man.td <- mantel(decostand(td, "standardize"), gm, permutations=9999)
  man.wd <- mantel(decostand(wd, "standardize"), gm, permutations=9999)
  
  tests <- list(gcd = man.gcd, td = man.td, wd = man.wd)
  return(tests)
}

ind_lin_mantel_nonSt <- function(meta, l, gcd, td, wd, gm) {
  #lin.vec = as.vector(subset(meta, Lineage == l, BioSample))
  lin.vec = as.vector(data.frame(meta[,l]))
  names(lin.vec) <- c("BioSample")
  
  gcd <- gcd[rownames(gcd) %in% lin.vec$BioSample, colnames(gcd) %in% lin.vec$BioSample]
  td <- td[rownames(td) %in% lin.vec$BioSample, colnames(td) %in% lin.vec$BioSample]
  wd <- wd[rownames(wd) %in% lin.vec$BioSample, colnames(wd) %in% lin.vec$BioSample]
  gm <- gm[rownames(gm) %in% lin.vec$BioSample, colnames(gm) %in% lin.vec$BioSample]
  
  man.gcd <- mantel(gcd, gm, permutations=9999)
  man.td <- mantel(td, gm, permutations=9999)
  man.wd <- mantel(wd, gm, permutations=9999)
  
  tests <- list(gcd = man.gcd, td = man.td, wd = man.wd)
  return(tests)
}

subMantelsStLog = list()
for (i in seq(1,10)) {
  subMantelsStLog[[i]] = list()
  sub.man <- ind_lin_mantel(subsamps, i, log(dvse+1), log(trade.tot+1), log(waypoint.mat+1), gm)
  subMantelsStLog[[i]] = sub.man
}

subMantelsSt = list()
for (i in seq(1,10)) {
  subMantelsSt[[i]] = list()
  sub.man <- ind_lin_mantel(subsamps, i, dvse, trade.tot, waypoint.mat, gm)
  subMantelsSt[[i]] = sub.man
}

subMantels = list()
for (i in seq(1,10)) {
  subMantels[[i]] = list()
  sub.man <- ind_lin_mantel_nonSt(subsamps, i, dvse, trade.tot, waypoint.mat, gm)
  subMantels[[i]] = sub.man
}

d1 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantels[[x]][["gcd"]]$statistic)))
d1$test <- "gcd"
d2 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantels[[x]][["td"]]$statistic)))
d2$test <- "td"
d3 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantels[[x]][["wd"]]$statistic)))
d3$test <- "wd"  

d4 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantelsSt[[x]][["gcd"]]$statistic)))
d4$test <- "gcd"
d5 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantelsSt[[x]][["td"]]$statistic)))
d5$test <- "td"
d6 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantelsSt[[x]][["wd"]]$statistic)))
d6$test <- "wd"  

d7 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantelsStLog[[x]][["gcd"]]$statistic)))
d7$test <- "gcd"
d8 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantelsStLog[[x]][["td"]]$statistic)))
d8$test <- "td"
d9 <- data.frame(statistic=do.call(rbind, lapply(seq(1,10), function(x) subMantelsStLog[[x]][["wd"]]$statistic)))
d9$test <- "wd"  

mantels <- data.frame(rbind(gcd=d1$statistic,td=d2$statistic,wd=d3$statistic,
                            gcdST=d4$statistic,tdST=d5$statistic,wdST=d6$statistic,
                            gcdSTLOG=d7$statistic,tdSTLOG=d8$statistic,wdSTLOG=d9$statistic))
names(mantels) <- paste("sub", seq(1,10), sep="")
mantels$test <- c("gcd", "td", "wd", "gcd", "td", "wd","gcd", "td", "wd")
mantels$cor <- c("none", "none", "none", "standardize", "standardize", "standardize","standardize(log(x+1))", "standardize(log(x+1))", "standardize(log(x+1))")
mantels.m <- melt(mantels, id.vars=c("test", "cor"))

mantels.m$test <- as.factor(mantels.m$test)
mantels.m$test <- factor(mantels.m$test,levels(mantels.m$test)[c(2,1,3)])

stP <- ggplot(rbind(mantels.m, fulMant) %>% filter(cor == "standardize")) +
  #facet_wrap(~cor) +
  geom_point(aes(x=test, y= value, group=variable, color=variable)) +
  geom_line(aes(x=test, y= value, group=variable, color=variable)) +
  scale_color_manual(values=c(rep("black", 10), "red")) +
  #geom_point(data=fulMant, aes(x=test, y=statistic), color="red") +
  #geom_line(data=fulMant, aes(x=test, y=statistic, group="1"), color="red") +
  theme_bw() +
  xlab("distance calculation") +
  ylab("Mantel Statistic") +
  theme(panel.grid = element_blank(),
        legend.position="none")



fulMant <- data.frame(rbind(c(statistic=man.dvse$statistic, test="gcd", cor="none"),
                            c(statistic=man.trade$statistic, test="td", cor="none"),
                            c(statistic=man.waypoint$statistic, test="wd", cor="none"),
                            c(statistic=man.dvse.st$statistic, test="gcd", cor="standardize"),
                            c(statistic=man.trade.st$statistic, test="td", cor="standardize"),
                            c(statistic=man.waypoint.st$statistic, test="wd", cor="standardize"),
                            c(statistic=man.dvse.st.log$statistic, test="gcd", cor="standardize(log(x+1))"),
                            c(statistic=man.trade.st.log$statistic, test="td", cor="standardize(log(x+1))"),
                            c(statistic=man.waypoint.st.log$statistic, test="wd", cor="standardize(log(x+1))")
                            ))
fulMant$statistic <- as.numeric(as.character(fulMant$statistic))
fulMant$variable <- "fullData"
names(fulMant)[1] <- "value"







mantels$full <- fulMant[match(mantels$test, fulMant$test), 'statistic']

ggplot(mantels) 
  geom_point(aes(x=full, y=statistic)) +
  #geom_errorbar(aes(x=full, ymin=min, ymax=max, color=lineage)) +
  scale_color_manual(values = c("pink", "blue", "purple", "red", "darkred", "darkgreen", "#FFCC33", "black")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab("Complete Sample") +
  ylab("Subsample") +
  ggtitle(expression("Pi"))
#lin2.man <- ind_lin_mantel(meta, 2, dvse, trade.tot, waypoint.mat, gm)
#lin3.man <- ind_lin_mantel(meta, 3, dvse, trade.tot, waypoint.mat, gm)
#lin4.man <- ind_lin_mantel(meta, 4, dvse, trade.tot, waypoint.mat, gm)
#lin5.man <- ind_lin_mantel(meta, 5, dvse, trade.tot, waypoint.mat, gm)
#lin6.man <- ind_lin_mantel(meta, 6, dvse, trade.tot, waypoint.mat, gm)
#lin7.man <- ind_lin_mantel(meta, 7, dvse, trade.tot, waypoint.mat, gm)


Runs<-readRDS("R_objects/Runs_PRF.rds")

library(NetCoMi)

RunsT <- get_taxadf(Runs, taxlevel = 7, type = "species")

taxa_names(RunsT) <- gsub(taxa_names(RunsT), pattern = "s__", replacement = "")

Y2021T <- subset_samples(RunsT, Year == "2021")
Y2022T <- subset_samples(RunsT, Year == "2022")

net_seas_p <- netConstruct(Y2021T,Y2022T,
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 50),
                           zeroMethod = "pseudo",
                           normMethod = "clr",
                           measure = "pearson",
                           verbose = 0)

netprops1 <- netAnalyze(net_seas_p, 
                        clustMethod = "cluster_fast_greedy")

nclust <- as.numeric(max(names(table(netprops1$clustering$clust1))))


my_cols<- c("#bc80bd",
           "#fb8072",
           "#80b1d3",
           "#fdb462",
           "#b3de69",
           "#fccde5",
           "#ffed6f")

set.seed(102)
plot(netprops1,
     layout="spring",
     nodeFilter = "clustMin",
     nodeFilterPar = 5,
     nodeSize = "eigenvector",
     #groupNames = c("2021", "2022"),
     sameLayout = TRUE,
     colorVec = my_cols,
     nodeTransp = 10,
     nodeSizeSpread = 4, 
     cexNodes = 0.5,
     #cexLabels = 7,
     negDiffCol = TRUE,
     posCol = "dodgerblue4",
     negCol = "red4",
     edgeWidth = 0.3,
     borderCol = "gray40",
     #edgeFilter="threshold",
     #edgeFilterPar=0.01,
     edgeTranspLow = 50, 
     edgeTranspHigh = 30,
     highlightHubs = TRUE
     )


legend(-1, 0.2, cex = 1.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = FALSE)

##
##    ACE.R
##
##    Ancestral Character State Reconstruction of anuran TE's
##    for Rubio et al. 202X
##    Created 06/05/2025     J. Rader
##    Modified 
##
################################################################################

## Housekeeping
# Set the working directory
setwd('C:/Users/Owner/Dropbox/Frogs/Genome_Size_Adam_Stuckert/') # Jonathan laptop
setwd('C:/Users/rader/Dropbox/Frogs/Genome_Size_Adam_Stuckert/') # New laptop (europa)

# Load required packages and functions
require(phytools)
require(plyr)
require(OUwie)
library(caper)


## Read in the data
frog.dat <- read.csv("data/reshaped_with_species_final_Manuscript_version.csv")

frog.dat$Species <- gsub(" ", "_", frog.dat$Species)

## Read in the tree
frogtree <- read.tree('data/altered_TreePL-Rooted_Anura_bestTree.tre')

drop <- setdiff(frogtree$tip.label, frog.dat$Species)
smalltree <- drop.tip(frogtree, drop)

drop2 <- setdiff(frog.dat$Species, smalltree$tip.label) # check that this is empty

# Order the data to match the tree
frog.dat <- frog.dat[order(frog.dat$Species),][rank(smalltree$tip.label),]


## Ancestral character state reconstructions

genome.anc <- anc.ML(smalltree, setNames(frog.dat$Genome.size, frog.dat$Species), maxit=10000, model="BM", CI=TRUE)

tecount.anc <- anc.ML(smalltree, setNames(frog.dat$total.TEs.counts, frog.dat$Species), maxit=10000, model="BM", CI=TRUE)
telength.anc <- anc.ML(smalltree, setNames(frog.dat$total.TEs.length, frog.dat$Species), maxit=10000, model="BM", CI=TRUE)
tepercent.anc <- anc.ML(smalltree, setNames(frog.dat$total.TEs.percentage_final, frog.dat$Species), maxit=10000, model="BM", CI=TRUE)
ty3.anc <- anc.ML(smalltree, setNames(frog.dat$total.TEs.percentage_final, frog.dat$Species), maxit=10000, model="BM", CI=TRUE)


edges <- smalltree$edge
internal_nodes <- unique(edges[edges[,2] > Ntip(smalltree), , drop=FALSE])
internal_edges <- edges[edges[,1] > Ntip(smalltree) & edges[,2] > Ntip(smalltree), , drop=FALSE]

check_overlap <- function(ci1, ci2) {
  return(!(ci1[2] < ci2[1] || ci2[2] < ci1[1]))
}

# Find adjacent parent and child nodes, look for overlapping confidence intervals
# non-overlap is indicative of a significant evolutionary transition
genome.cis <- genome.anc$CI95
for (i in 1:nrow(internal_edges)) {
  parent <- internal_edges[i,1]
  child <- internal_edges[i,2]
  
  ci_parent <- genome.cis[as.character(parent), ]
  ci_child <- genome.cis[as.character(child), ]
  
  overlap <- check_overlap(ci_parent, ci_child)
  cat("Nodes", parent, "and", child, "overlap:", overlap, "\n")
}

percent.cis <- tepercent.anc$CI95
for (i in 1:nrow(internal_edges)) {
  parent <- internal_edges[i,1]
  child <- internal_edges[i,2]
  
  ci_parent <- percent.cis[as.character(parent), ]
  ci_child <- percent.cis[as.character(child), ]
  
  overlap <- check_overlap(ci_parent, ci_child)
  cat("Nodes", parent, "and", child, "overlap:", overlap, "\n")
}

count.cis <- tecount.anc$CI95
for (i in 1:nrow(internal_edges)) {
  parent <- internal_edges[i,1]
  child <- internal_edges[i,2]
  
  ci_parent <- count.cis[as.character(parent), ]
  ci_child <- count.cis[as.character(child), ]
  
  overlap <- check_overlap(ci_parent, ci_child)
  cat("Nodes", parent, "and", child, "overlap:", overlap, "\n")
}

length.cis <- telength.anc$CI95
for (i in 1:nrow(internal_edges)) {
  parent <- internal_edges[i,1]
  child <- internal_edges[i,2]
  
  ci_parent <- length.cis[as.character(parent), ]
  ci_child <- length.cis[as.character(child), ]
  
  overlap <- check_overlap(ci_parent, ci_child)
  cat("Nodes", parent, "and", child, "overlap:", overlap, "\n")
}


# Prepare a data frame to store results
results <- data.frame(
  parent = integer(),
  child = integer(),
  overlap = logical(),
  stringsAsFactors = FALSE
)

# Loop through internal edges and record overlap
for (i in 1:nrow(internal_edges)) {
  parent <- internal_edges[i,1]
  child <- internal_edges[i,2]
  
  ci_parent <- genome.cis[as.character(parent), ]
  ci_child <- genome.cis[as.character(child), ]
  
  overlap <- check_overlap(ci_parent, ci_child)
  
  results[i, ] <- list(parent, child, overlap)
}

cor.test(genome.anc$ace, tepercent.anc$ace, method="pearson")
cor.test(genome.anc$ace, tecount.anc$ace, method="pearson")
cor.test(genome.anc$ace, telength.anc$ace, method="pearson")



## Phylogenetic signal analysis

genome.lambda <- phylosig(smalltree, setNames(frog.dat$Genome.size, frog.dat$Species), method = "lambda", nsim = 1000, test = TRUE)
genome.k <- phylosig(smalltree, setNames(frog.dat$Genome.size, frog.dat$Species), method = "K", nsim = 1000, test = TRUE)

count.lambda <- phylosig(smalltree, setNames(frog.dat$total.TEs.counts, frog.dat$Species), method = "lambda", nsim = 1000, test = TRUE)
count.k <- phylosig(smalltree, setNames(frog.dat$total.TEs.counts, frog.dat$Species), method = "K", nsim = 1000, test = TRUE)

length.lambda <- phylosig(smalltree, setNames(frog.dat$total.TEs.length, frog.dat$Species), method = "lambda", nsim = 1000, test = TRUE)
length.k <- phylosig(smalltree, setNames(frog.dat$total.TEs.length, frog.dat$Species), method = "K", nsim = 1000, test = TRUE)

percent.lambda <- phylosig(smalltree, setNames(frog.dat$total.TEs.percentage_final, frog.dat$Species), method = "lambda", nsim = 1000, test = TRUE)
percent.k <- phylosig(smalltree, setNames(frog.dat$total.TEs.percentage_final, frog.dat$Species), method = "K", nsim = 1000, test = TRUE)

signal.table <- cbind(rbind(genome.lambda$lambda, percent.lambda$lambda, count.lambda$lambda, length.lambda$lambda),
                      rbind(genome.lambda$P, percent.lambda$P, count.lambda$P, length.lambda$P),
                      rbind(genome.k$K, percent.k$K, count.k$K, length.k$K),
                      rbind(genome.k$P, percent.k$P, count.k$P, length.k$P))

signal.table <- round(signal.table, 3)

row.names(signal.table) <- c("Genome Length", "TE Percent", "TE Count", "TE Length")
colnames(signal.table) <- c("Lambda", "P_lambda", "K", "P_K")

write.csv(signal.table, file="data/signal_table.csv", quote = FALSE)



## FIGURES ---------------------------------------------------------------------


## Ancestral character state plots

anctree <- smalltree
anctree$edge.length <- rep(1, length(anctree$edge.length))

pdf(file="figures/tepercent_anc_tree.pdf", width=8, height=10, pointsize=12)

# First, plot egg size
plotTree(anctree, offset=2, ftype="i")

# Get coordinates of the plotted tree
coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Extract X and Y coordinates of tip labels
tip_coords <- data.frame(
  tip_label = anctree$tip.label,
  x = coords$xx[1:length(anctree$tip.label)],  # X-coordinates of tips
  y = coords$yy[1:length(anctree$tip.label)]   # Y-coordinates of tips
)

# set point size based on egg size
# point scaling based on mean mass

mincex <- 0.5
maxcex <- 5

slope <- (maxcex-mincex)/((max(frog.dat$Genome.size))-(min(frog.dat$Genome.size)))
intercept <- mincex - ((min(frog.dat$Genome.size))*(slope))

frog.dat$genome.cex <- (frog.dat$Genome.size*slope)+intercept

points(tip_coords$x, tip_coords$y, pch=21, cex=frog.dat$genome.cex, col="white", bg="springgreen2")

nodelabels(text="", as.numeric(names(genome.anc$ace)), pch=21,
           cex=(genome.anc$ace*slope+intercept), frame="none",
           col="white", bg="springgreen2")

# Build a legend
points(rep(2,4),c(60, 57, 54, 51), cex=(seq(min(frog.dat$Genome.size),max(frog.dat$Genome.size),
                                              (max(frog.dat$Genome.size)-min(frog.dat$Genome.size))/4))*slope+intercept, 
       pch=21, col="white", bg="springgreen2")

text(rep(1.5,4),c(60, 57, 54, 51), 
     formatC(seq(min(frog.dat$Genome.size),max(frog.dat$Genome.size),
                 (max(frog.dat$Genome.size)-min(frog.dat$Genome.size))/3), 
             format = "e", digits = 2), pos=2)

# Now plot adult body size
par(new=TRUE)
plotTree(anctree, offset=2, ftype="i")

# Get coordinates of the plotted tree
coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Extract X and Y coordinates of tip labels
tip_coords <- data.frame(
  tip_label = anctree$tip.label,
  x = coords$xx[1:length(anctree$tip.label)],  # X-coordinates of tips
  y = coords$yy[1:length(anctree$tip.label)]   # Y-coordinates of tips
)

# set point size based on egg size
# point scaling based on mean mass

mincex <- 0.5
maxcex <- 5

slope <- (maxcex-mincex)/((max(frog.dat$total.TEs.percentage_final))-(min(frog.dat$total.TEs.percentage_final)))
intercept <- mincex - ((min(frog.dat$total.TEs.percentage_final))*(slope))

frog.dat$percent.cex <- (frog.dat$total.TEs.percentage_final*slope)+intercept

points(tip_coords$x, tip_coords$y, pch=21, cex=frog.dat$percent.cex, col="white", bg="#ac6aff")

nodelabels(text="", as.numeric(names(tepercent.anc$ace)), pch=21,
           cex=(tepercent.anc$ace*slope+intercept), frame="none",
           col="white", bg="#ac6aff")

# Build a legend
points(rep(2,4),c(60, 57, 54, 51), cex=(c(20, 40, 60, 80)*slope+intercept), pch=21,
       col="white", bg="#ac6aff")

text(rep(2.5,4),c(60, 57, 54, 51), c(20, 40, 60, 80), pos=4)

dev.off()


pdf(file="figures/genomesVtes.pdf", width=8, height=4, pointsize=12)

par(mfrow=c(1,3),
    mar=c(0,0,0,0),
    oma=c(0.1,0.1,0.1,0.1),
    plt=c(0.1, 0.9, 0.25, 0.95))
par(xpd=NA)


plot(tepercent.anc$ace, genome.anc$ace, pch=21, bg="white", col="#00C094",
     xlim=c(min(frog.dat$total.TEs.percentage_final), max(frog.dat$total.TEs.percentage_final)),
     ylim=c(min(frog.dat$Genome.size), max(frog.dat$Genome.size)))
points(frog.dat$total.TEs.percentage_final, frog.dat$Genome.size, pch=21, col="white", bg="#00C094")



plot(tecount.anc$ace, genome.anc$ace, pch=21, bg="white", col="cadetblue4", 
     xlim=c(min(frog.dat$total.TEs.counts), max(frog.dat$total.TEs.counts)),
     ylim=c(min(frog.dat$Genome.size), max(frog.dat$Genome.size)))
points(frog.dat$total.TEs.counts, frog.dat$Genome.size, pch=21, col="white", bg="cadetblue4")



plot(telength.anc$ace, genome.anc$ace, pch=21, bg="white", col="darkolivegreen4", 
     xlim=c(min(frog.dat$total.TEs.length), max(frog.dat$total.TEs.length)),
     ylim=c(min(frog.dat$Genome.size), max(frog.dat$Genome.size)))
points(frog.dat$total.TEs.length, frog.dat$Genome.size, pch=21, col="white", bg="darkolivegreen4")


dev.off()


## Phylosig figure

pdf(file="figures/phylosig.pdf", width=8, height=8, pointsize=12)

par(mfrow=c(4, 2), xpd=NA)

plot.phylosig(genome.lambda)
plot.phylosig(genome.k)
plot.phylosig(percent.lambda)
plot.phylosig(percent.k)
plot.phylosig(count.lambda)
plot.phylosig(count.k)
plot.phylosig(length.lambda)
plot.phylosig(length.k)


dev.off()


#### Currently unusued code

edge_colors <- rep("black", nrow(smalltree$edge))  # default color for all edges

for (i in 1:nrow(results)) {
  parent <- results$parent[i]
  child <- results$child[i]
  overlap <- results$overlap[i]
  
  # Find the row index of the edge in the tree's edge matrix
  edge_index <- which(smalltree$edge[,1] == parent & smalltree$edge[,2] == child)
  
  if (length(edge_index) == 1) {
    edge_colors[edge_index] <- ifelse(overlap, "forestgreen", "firebrick")
  }
}

plot(smalltree, edge.color=edge_colors, cex=1)
title("Edge Colors by CI Overlap")
legend("topleft", legend=c("Overlap", "No Overlap"), col=c("forestgreen", "firebrick"), lwd=2)


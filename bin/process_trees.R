library(ape)
library(ggmap)
library(phyloseq)
library(tidyverse)

contig = commandArgs(trailingOnly=TRUE)[[1]]
tree <- ape::read.tree(paste0(contig,".tree"))

# Optionally set an outgroup.
# tree <- root(tree,outgroup = "outgroup", resolve.root = T)

treeSegs <- phyloseq::tree_layout(
                                phyloseq::phy_tree(tree),
                                ladderize = T
                                )

treeSegs$edgeDT <- treeSegs$edgeDT  %>% 
                   dplyr::mutate(edge.length = 
                                    ifelse(edge.length < 0, 0, edge.length)
                                 , xright = xleft + edge.length
                                 )
edgeMap = aes(x = xleft, xend = xright, y = y, yend = y)
vertMap = aes(x = x, xend = x, y = vmin, yend = vmax)
labelMap <- aes(x = xright+0.0001, y = y, label = OTU)

ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) + 
  geom_segment(vertMap, data = treeSegs$vertDT) +
  geom_text(labelMap, data = dplyr::filter(treeSegs$edgeDT, !is.na(OTU)), na.rm = TRUE, hjust = -0.05) +
  ggmap::theme_nothing() + 
  scale_x_continuous(limits = c(
    min(treeSegs$edgeDT$xleft)-0.15,
    max(treeSegs$edgeDT$xright)+0.15
  ),
  expand = c(0,0)) +
  labs(title = contig)

ggsave(paste0(contig,".svg"), height = 45, width = 10)
ggsave(paste0(contig,".png"), height = 45, width = 10)  
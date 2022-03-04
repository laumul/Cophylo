data <- read.table("results_summary_DTL.txt", header = TRUE)

mean(data[["duplications"]])

boxplot(data$duplications, data$transfers, data$losses, data$speciations, 
        col = c("light blue",  "light green","light pink", "orange"),
        outline = FALSE)

        



###################
### Will's code ###
###################


library(phytools)

spT <- read.tree("h_mcc.tre")

#traits ordered by edge number in tree and ladderized
traits_edge <- read.csv2("final2.csv")

as.numeric(traits_edge$Duplications)

#assign traits to vectors. removed ancestral branch as not displayed
dup_edge <- as.numeric(traits_edge$Duplications)[2:length(as.numeric(traits_edge$Duplications))]
ht_edge <- as.numeric(traits_edge$Transfers)[2:length(as.numeric(traits_edge$Transfers))]
loss_edge <- as.numeric(traits_edge$Losses)[2:length(as.numeric(traits_edge$Losses))]
cosp_edge <- as.numeric(traits_edge$copies)[2:length(as.numeric(traits_edge$copies))]


# set up gradient color scales
n <- 1000

colfunc <- function(x){colorRampPalette(c("blue", x))}
gradient <- colfunc("red")(n)

max_edge <- max(c(dup_edge,ht_edge,loss_edge))

dup_color = gradient[n*max(dup_edge)/max_edge]
ht_color = gradient[n*max(ht_edge)/max_edge]
loss_color = gradient[n*max(loss_edge)/max_edge]


### Figure 4 ###
#only host switching 

pdf("host_switching.pdf")

#windows()
plotBranchbyTrait(ladderize(spT), ht_edge, mode="edges", palette=colfunc(ht_color), xlim=c(0,1.003), no.margin=FALSE)
title(main="Host Switching")
dev.off()




#for multiframe plotting
pdf("evoprobybranch.pdf")
par(mfrow = c(1,3))

plotBranchbyTrait(ladderize(spT), dup_edge, mode="edges", palette=colfunc(dup_color), xlim=c(0,0.25), no.margin=FALSE)
title(main="Parasite Speciation")

plotBranchbyTrait(ladderize(spT), ht_edge, mode="edges", palette=colfunc(ht_color), xlim=c(0,1.25), no.margin=FALSE)
title(main="Host Switching")

plotBranchbyTrait(ladderize(spT), loss_edge, mode="edges", palette=colfunc(loss_color), xlim=c(0,0.75), no.margin=FALSE)
title(main="Parasite Extinction")

dev.off()


pdf("parasite_extinction.pdf")


plotBranchbyTrait(ladderize(spT), loss_edge, mode="edges", palette=colfunc(loss_color), xlim=c(0,0.75), no.margin=FALSE)
title(main="Parasite Extinction")

dev.off()




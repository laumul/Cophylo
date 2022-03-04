library(ape)
setwd("/Users/mulvey/Desktop")
library(phytools)

tr<- read.nexus("badets2.nexus")


plot(tr)
nodelabels()

tr25<-extract.clade(tr,25)
plotTree(tr25)
tr25$tip.label
dtips<-sample(tr25$tip.label,15)
dtips
dt<-drop.tip(tr25,dtips)
plotTree(dt)


write.tree(tr17,"badets.tree")


tip<- c("Wetapolystoma_Almae", "Polystoma_gallieni" ,  "Polystoma_spRV1"  )

tr17<-drop.tip(tr25,tip)



setwd("Desktop/FAU/Thesis_analysis/Paper/ALE_Simulations/")

library(ggplot2)

data <- read.table("probs.csv")

data <- subset(data, select = -(V1))

names(data)[names(data) == "V2"] <- "Speciation"
names(data)[names(data) == "V3"] <- "Host_Switching"
names(data)[names(data) == "V4"] <- "Extinction"
names(data)[names(data) == "V5"] <- "Cospeciation"

clrs<- c("light blue","light green", "light pink", "orange")
irg<- c(0, 6.57, 3.97, 13.4 )
nms<-names(data)
nms[2]

#i=4
plots<-list()

pdf("Simulated_data.pdf")
for (i in 1:length(nms) ){
  it <- nms[i]
 p<- hist(data[,i], col = clrs[i], xlim = c (0,20), main = it)
  abline(v=irg[i], col="black",lty = 1, lwd = 3)
  plots[[i]]<-p
}

multiplot(plotlist=plots, cols =4)
dev.off()



pv<-data$Speciation>irg[1]
pv<-data$Host_Switching> irg[2]
pv<-data$Extinction>irg[3]
table(pv)
###number I got higher than random
pv<-data$Cospeciation<irg[4]
table(pv)
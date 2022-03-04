setwd("Desktop/FAU/Thesis_analysis/Paper/ALE_Simulations/")

library(TreeSim)
library(geiger)
library(ape)
library(phytools)
n<-17
lambda <- 1
mu <- 0.9
frac <-0.6
numbsim<-1000



t<-sim.bd.taxa(n,numbsim,lambda,mu,frac,complete=FALSE,stochsampling=FALSE)

tree<-"tree"
dir.create("output")
temdir<-"output/"

#i=2

names<- c("Ranarugosa_Diplorchis_ranae", "Buforegularis_Eupolystoma_Alluaudi",
          "Schismadermacarens_Eupolystoma_sp1", "Litoriagracilenta_Parapolystoma_sp",
          "Physalaemuscuvieri_Polystoma_cuvieri", "Ptychadenaanchietae_Polystoma_dawiekokis",
          "Rhacophorusmaximus_Polystoma_indicum", "Ranatemporaria_Polystoma_integerrimum",
          "Phrynohyasvenulosa_Polystoma_lopezromani", "Hyperoliusviridiflavus_Polystoma_marmorati",
          "Smiliscabaudinii_Polystoma_naevius", "Hylaversicolor_Polystoma_nearcticum", 
          "Leptodactylusmystaceus_Polystoma_sp2", "Hylacinerea_Polystoma_spHC",
          "Rhacophorusarboreus_Polystoma_spRA", "Rhacophorusomeimontis_Polystoma_spRO",
          "Strongylopusfasciatus_Polystoma_testimagna")


run <- (1:1000)

for (i in 1:length(run)) {
  name <- i
  myfile <- paste0(temdir,tree,i,".tre")
  t[[i]]$tip.label <- names
  write.tree(t[[i]], myfile)
  
}
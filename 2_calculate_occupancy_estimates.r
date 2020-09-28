library(MASS)
require(sp)
require(rgdal)
require(maps)
require(doBy)
require(ggplot2)
require(gridExtra)
require(lubridate)
library(chron)
library(dplyr)
library(chron)
library(geosphere)
library(reshape2)
library(car)
library(data.table)
library(sparta)
library(doParallel)
library(foreach)
nbcl=24
cl<-makeCluster(nbcl) 
registerDoParallel(cl)
nbit=65000
burnin=50000


tab=as.data.frame(fread("BD_Finale_20180914.txt",sep="\t",header=T))
tab$Annee=as.numeric(substr(as.character(tab$Date),1,4))
tab$Mois=as.numeric(substr(as.character(tab$Date),5,6))
tab$Jour=as.numeric(substr(as.character(tab$Date),7,8))
tab$Jour.de.collecte=yday(as.Date(paste(tab$Jour,tab$Mois,tab$Annee,sep="/"),format ="%d/%m/%Y"))
tab$week=week(as.Date(paste(tab$Jour,tab$Mois,tab$Annee,sep="/"),format ="%d/%m/%Y"))
tab=subset(tab,Mois!=0 & Jour!=0 & !is.na(Jour.de.collecte))
tab=subset(tab,Annee>=1950 & Annee<2017)
tab$species[which(tab$species=="Hylaeus bipunctatus")]="Hylaeus signatus"
tab$species[which(tab$species=="Osmia fulviventris")]="Osmia niveata"
tab$species[which(tab$species=="Xylocopa violacea")]="Xylocopa violaceae"
tab$species[which(tab$species=="Hoplosmia spinulosa")]="Osmia spinulosa"
tab$species[which(tab$species=="Chalicodoma ericetorum")]="Megachile ericetorum"
tab$species[which(tab$species=="Panurgus calcaratus calcaratus")]="Panurgus calcaratus"

#FAIRE LA MATRICE
library(plyr)
tab$Longitude2=round_any(tab$Longitude,0.01)
tab$Latitude2=round_any(tab$Latitude,0.01)
tab$dates=as.Date(paste(tab$Jour,tab$Mois,tab$Annee,sep="/"),format ="%d/%m/%Y")
tab=subset(tab,!is.na(dates))
objformat=formatOccData(gsub(" ","_",tab$species,fixed=T), paste(tab$Latitude2,tab$Longitude2,sep="-"),tab$dates,includeJDay = T)
mat=as.data.frame(fread("matrice_pres_abs.txt",sep="\t",header=T))
liste=data.frame(species=unique(tab$species),trend_bay=NA,trend_fre=NA,trend_fre_prob=NA)


#(nrow(liste)-nbcl+1):nrow(liste)
P1=foreach(i=1:nrow(liste),.combine=rbind)%dopar%{
library(sparta)
library(doBy)
model2=occDetFunc(taxa_name =gsub(" ","_",paste(liste$species[i]),fixed=T), occDetdata =objformat$occDetdata, spp_vis = objformat$spp_vis, n_iterations = nbit,
  nyr = 2, burnin = burnin, thinning = 3, n_chains = 2,
  write_results = T, output_dir = getwd(), modeltype = c('ranwalk','halfcauchy','jul_date'),seed = 5)
}

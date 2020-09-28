setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/data")
library(ggplot2)
library(mgcv)
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
library(rgbif)
library(reshape2)
library(car)
library(data.table)
library(lme4)
library(pastecs)
library(fitdistrplus)
library(rgbif)

#CHARGER LISTE ESPECES CIBLES
tab=as.data.frame(fread("BD_Finale_pour phenological_shifts.txt",sep="\t",header=T))

#CHARGER DONNEES GBIF
hym=fread("gbif_hymenos.csv")
hym=hym[,c("gbifid","species","day","month","year","countrycode","locality","decimallatitude","decimallongitude","specieskey","identifiedby",
"institutioncode","order","family","genus")]
names(hym)=c("Reference","Espèce","Jour","Mois","Annee","Pays","locality","Latitude","Longitude","ID","Collecteur","Source","ORDRE","FAMILLE","GENRE")

#CHECKER SYNONIMIE:
library(rgbif)
liste=data.frame(species=setdiff(tab$species,hym$Espèce))
liste$Espèce=NA
liste$Espèce[which(liste$species=="Seladonia confusa")]="Halictus confusus"
liste$Espèce[which(liste$species=="Seladonia tumulorum")]="Halictus tumulorum"
for(i in 1:nrow(liste)){
if(is.na(liste$Espèce[i])){
liste$Espèce[i]=name_backbone(name=liste[i,"species"],rank='species', kingdom='animals')$species
}}

tab=merge(tab,liste,by="species",all.x=T)
tab$Espèce[which(is.na(tab$Espèce))]=tab$species[which(is.na(tab$Espèce))]

liste=data.frame(species=unique(tab$species),Espèce=unique(tab$Espèce))
library(plyr)
#ARRONDIR COORDONNEES POUR FAIRE DES MAILLES
hym$Latitude=round_any(hym$Latitude,0.01)
hym$Longitude=round_any(hym$Longitude,0.01)
#CONSERVER LES DONNEES HYMENOS D'ESPECES NON-CIBLES SOUS UNE SEULE APPELATION POUR SIMPLIFIER:
hym$Espèce[which(!(hym$Espèce %in% liste$Espèce))]="autre"
#FAIRE LA MATRICE
mat=dcast(hym,Latitude+Longitude~Espèce,fun.aggregate=length,value.var="Reference")
#CALCULER LA PRESSION D'ECHANTILLONNAGE PAR MAILLE (NB DE DONNEES)
mat$sampl=apply(mat[,-c(1:2)],1,sum)
#AJOUTER DONNEES DE TEMP POUR CHAQUE MAILLES
if(!requireNamespace("sf")) install.packages("sf")
if(!requireNamespace("rdryad")) devtools::install_github("ropensci/rdryad")
chelsa_raster = "wc2.0_bio_30s_01.tif"
temp_moy <-  raster::raster(chelsa_raster)
chelsa_raster = "wc2.0_bio_30s_04.tif"
conti <-  raster::raster(chelsa_raster)
sf_pts <- sf::st_transform(sf::st_as_sf(mat[,c("Longitude","Latitude")], coords = c("Longitude", "Latitude"), crs = 4326), raster::projection(temp_moy))
mat$temp_moy=raster::extract(temp_moy, as(sf_pts, 'Spatial'))
mat$conti=raster::extract(conti, as(sf_pts, 'Spatial'))/100
ggplot(data=mat,aes(x=Longitude,y=Latitude,color=temp_moy))+geom_point()

mat=as.data.frame(mat)
liste=data.frame(Latitude=mat$Latitude[which(is.na(mat$temp_moy))],Longitude=mat$Longitude[which(is.na(mat$temp_moy))],temp_moy=NA,conti=NA)
euclid_dist=function(x,y){
sqrt((x[,1]-y[,1])^2+(x[,2]-y[,2])^2)}
#AJOUTER DONNEES TEMP AUX DONNEES LEGEREMENT EN MER:
for(i in 1:nrow(liste)){
bidon=subset(mat,Latitude>=liste$Latitude[i]-0.1 & Latitude<=liste$Latitude[i]+0.1 & Longitude<=liste$Longitude[i]+0.1 & Longitude>=liste$Longitude[i]-0.1 & !is.na(temp_moy))
if(nrow(bidon)>0){
mat[which(mat$Latitude==liste$Latitude[i] & mat$Longitude==liste$Longitude[i]),c("temp_moy","conti")]=bidon[which.min(euclid_dist(liste[i,c("Latitude","Longitude")],bidon[,c("Latitude","Longitude")])),c("temp_moy","conti")]
}
}
#CACLUL STI ET SCI
library(dplyr)
liste=data.frame(species=unique(tab$species),Espèce=unique(tab$Espèce))
liste$temp_index=NA
liste$temp_var=NA
liste$nbdatagbif=NA
for(i in 145:nrow(liste)){
#STI
liste$temp_index[i]=sum((mat[,paste(liste$Espèce[i])]/mat$sampl)*mat$temp_moy,na.rm=T)/sum((mat[,paste(liste$Espèce[i])]/mat$sampl),na.rm=T)
#SCI
liste$temp_var[i]=sum((mat[,paste(liste$Espèce[i])]/mat$sampl)*mat$conti,na.rm=T)/sum((mat[,paste(liste$Espèce[i])]/mat$sampl),na.rm=T)
liste$nbdatagbif[i]=sum(mat[,paste(liste$Espèce[i])])
}

ggplot(data=liste,aes(x=conti,y=temp_moy))+geom_point()

write.table(liste,"temperature_indexes_by_species.txt",row.names=F,sep="\t")

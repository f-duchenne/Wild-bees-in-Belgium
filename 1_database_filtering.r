################################################################################################
#CHARGER LES PACKAGES UTILES:
setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/scripts finaux/article/data")
library(lubridate)
library(chron)
library(dplyr)
library(chron)
library(reshape2)
library(data.table)

tab=as.data.frame(fread("BD_Finale_20180914.txt",sep="\t",header=T))
tab$ID=paste0("ID",1:nrow(tab))
traits=as.data.frame(fread("traits_full_new_liste_MG.txt",sep="\t",header=T))
tab$Annee=as.numeric(substr(as.character(tab$Date),1,4))
tab$Mois=as.numeric(substr(as.character(tab$Date),5,6))
tab$Jour=as.numeric(substr(as.character(tab$Date),7,8))
tab$Jour.de.collecte=yday(as.Date(paste(tab$Jour,tab$Mois,tab$Annee,sep="/"),format ="%d/%m/%Y"))
tab$week=week(as.Date(paste(tab$Jour,tab$Mois,tab$Annee,sep="/"),format ="%d/%m/%Y"))
tab=subset(tab,Mois!=0 & Jour!=0 & !is.na(Jour.de.collecte))
tab=subset(tab,Annee>=1950 & Annee<2017)
setwd(dir="C:/Users/Francois/Documents/land use change/elevation")
w=raster("elevation1x1_new.tif")
sf_pts=st_as_sf(tab,coords=c("Longitude","Latitude"))%>% st_set_crs(4326) %>% st_transform(st_crs(w)) #en faire un spatial dataframe
tab$Altitude=raster::extract(w, as(sf_pts, 'Spatial'))


#Correting some mistakes
tab$species[which(tab$species=="Hylaeus bipunctatus")]="Hylaeus signatus"
tab$species[which(tab$species=="Osmia fulviventris")]="Osmia niveata"
tab$species[which(tab$species=="Xylocopa violacea")]="Xylocopa violaceae"
tab$species[which(tab$species=="Hoplosmia spinulosa")]="Osmia spinulosa"
tab$species[which(tab$species=="Chalicodoma ericetorum")]="Megachile ericetorum"
tab$species[which(tab$species=="Lasioglossum sabulosum")]="Lasioglossum monstrificum"

#Homogenize sex:
tab$Sexe[which(tab$Sexe=="f")]="F"
tab$Sexe[which(tab$Sexe=="FE")]="F"
tab$Sexe[which(tab$Sexe=="MA")]="M"
tab$Sexe[which(tab$Sexe=="m")]="M"
tab$Sexe[which(tab$Sexe=="w")]="W"
tab$Sexe[which(tab$Sexe=="IM")]="immature"
tab$Sexe[which(tab$Sexe=="I")]="roi"
tab$Sexe[which(tab$Sexe=="0")]=NA
tab$Sexe[which(tab$Sexe=="4")]="oeufs"
tab$Sexe[which(tab$Sexe=="1")]="M"
tab=subset(tab,Sexe=="F" | Sexe=="M" | Sexe=="W" | is.na(Sexe))

tab$cat=NA
tab$cat[which(tab$Annee<=1980 & tab$Annee>=1950)]="A"
tab$cat[which(tab$Annee>=1980)]="B"
tab$cat[which(tab$Annee>=1990)]="C"
b=summaryBy(Annee~species+cat,data=tab,FUN=length)
b=dcast(b,species~cat,value.var="Annee.length")
b$D=b$A+b$B+b$C
b2=subset(b,D>=30 & !is.na(D))

tab=tab[tab$species %in% b2$species,]
liste=data.frame(species=unique(tab$species))

library(randomForest)
a=0
b=0
filtre_taxo=c("Bombus cryptarum","Bombus terrestris","Bombus lucorum","Bombus magnus","Hylaeus brevicornis","Hylaeus annularis","Andrena strohmella","Andrena pilipes","Andrena minutuloides",
"Lasioglossum monstrificum")
for(i in 1:nrow(liste)){
bidon=subset(tab,tab$species==liste$species[i])
bidon$sexe_inf=NA
if(!(liste$species[i] %in% filtre_taxo)){ 
if(subset(traits,species==liste$species[i])$Sociality=="Solitary"){
bidon$Sexe[which(bidon$Sexe=="W")]=NA}
mi=randomForest(as.factor(Sexe)~Jour.de.collecte+Annee+Altitude+Latitude+Longitude,data=subset(bidon,!is.na(Sexe)))
bidon$sexe_inf=predict(mi,newdata=bidon)
a=a+length(bidon$sexe_inf[which(bidon$Sexe!=bidon$sexe_inf & !is.na(bidon$Sexe))])
b=b+length(bidon$Sexe[which(is.na(bidon$Sexe))])
bidon$sexe_inf[which(!is.na(bidon$Sexe))]=bidon$Sexe[which(!is.na(bidon$Sexe))]
}
if(i==1){res=bidon}else{res=rbind(res,bidon)}
}
#error sex inference:
a/b


res=res[!(res$species %in% filtre_taxo),]

write.table(res,"BD_Finale_pour phenological_shifts.txt",sep="\t",row.names=F)

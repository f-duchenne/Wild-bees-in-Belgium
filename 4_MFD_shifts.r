#CHARGER LES PACKAGES UTILES:
setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/scripts finaux/article/data")
library(mgcv)
library(MASS)
require(doBy)
library(dplyr)
library(reshape2)
library(car)
library(data.table)
library(lme4)
library(MuMIn)


tab=as.data.frame(read.table("BD_Finale_pour phenological_shifts.txt",sep="\t",header=T))
liste=data.frame(species=unique(tab$species))
liste$year_effect=NA
liste$year_err=NA
liste$year_pval=NA
liste$sex_pval=NA
liste$rs=NA
liste$minp=NA
liste$maxp=NA
liste$fquant=NA
liste$lquant=NA
liste$moyp=NA
liste$medp=NA
liste$nbannee=NA
liste$mfd=NA
liste$pheno_length=NA
liste$rechauf=NA
liste$n=NA

for(i in 1:nrow(liste)){
bidon=subset(tab,species==liste$species[i])
model=lmer(Jour.de.collecte~Latitude*Longitude+Altitude+(1|Annee)+(1|sexe_inf)+poly(Annee,2),data=bidon)
if(isSingular(model)){print(i)}
rca=r.squaredGLMM(model)[2]
trend=data.frame(Latitude=mean(bidon$Latitude,na.rm=T),Longitude=mean(bidon$Longitude,na.rm=T),
Altitude=mean(bidon$Altitude,na.rm=T),Annee=1950:2016)

mySumm <- function(.) {
  predict(., newdata=trend, re.form=~(1|Annee),allow.new.levels =T)
}
trend$fit=predict(model, newdata=trend, re.form=~(1|Annee),allow.new.levels =T)
trend$se.fit=summary(bootMer(model, mySumm, nsim=250,re.form=~(1|Annee), type="parametric"))$bootSE
trend$species=liste$species[i]
trend2=trend[trend$Annee %in% unique(bidon$Annee),]

if(i==1){
resquali=trend
resquali2=trend2
}else{
resquali=rbind(resquali,trend)
resquali2=rbind(resquali2,trend2)
}

model2b=lm(fit~Annee,weights=1/se.fit,data=trend2)
summa=summary(model2b)
ano=Anova(model2b)

liste$year_effect[i]=model2b$coef["Annee"]
liste$year_err[i]=summa$coefficient["Annee","Std. Error"]
liste$year_pval[i]=ano["Annee",4]
liste$sex_pval[i]=ano["sexe_inf",4]
liste$rs[i]=summa$adj.r.squared
liste$minp[i]=min(bidon$Annee)
liste$maxp[i]=max(bidon$Annee)
liste$moyp[i]=mean(bidon$Annee)
liste$medp[i]=median(bidon$Annee)
liste$fquant[i]=quantile(bidon$Annee,probs=0.1)
liste$lquant[i]=quantile(bidon$Annee,probs=0.9)
liste$nbannee[i]=length(unique(bidon$Annee))
liste$mfd[i]=mean(bidon$Jour.de.collecte)
liste$pheno_length[i]=sd(bidon$Jour.de.collecte)
liste$rsquared[i]=rca
liste$n[i]=nrow(bidon)
}

trend=as.data.frame(fread("yearly_estimates_of_occupancy.txt",sep="\t",header=T))
trend=trend[trend$species %in% liste$species,]

liste2=merge(resquali,trend,by=c("species","Annee"),all.x=T,all.y=F)
write.table(liste2,"yearly_estimates_of_occupancy_and_mfd.txt",sep="\t",row.names=F)

liste2=merge(resquali2,trend,by=c("species","Annee"),all=T)
write.table(liste2,"yearly_estimates_of_occupancy_and_mfd_only_for_years_withdata.txt",sep="\t",row.names=F)

write.table(liste,"linear_mfd_shifts.txt",sep="\t",row.names=F)
#########################

#figure pour thèse
setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/scripts finaux/article/data")
library(mgcv)
library(MASS)
require(doBy)
library(dplyr)
library(reshape2)
library(car)
library(data.table)
library(lme4)
library(MuMIn)
library(ggeffects)


tab=as.data.frame(read.table("BD_Finale_pour phenological_shifts.txt",sep="\t",header=T))
liste=data.frame(species=unique(tab$species))


i=65

bidon=subset(tab,species==liste$species[i])
model=lmer(Jour.de.collecte~Latitude*Longitude+Altitude+(1|Annee)+(1|sexe_inf)+poly(Annee,2),data=bidon)

mySumm <- function(.) {
  predict(., newdata=trend, re.form=~(1|Annee),allow.new.levels =T)
}


trend=data.frame(Latitude=mean(bidon$Latitude,na.rm=T),Longitude=mean(bidon$Longitude,na.rm=T),
Altitude=mean(bidon$Altitude,na.rm=T),Annee=1950:2016)
trend$fit=predict(model, newdata=trend, re.form=~(1|Annee),allow.new.levels =T)
trend$se.fit=summary(bootMer(model, mySumm, nsim=250,re.form=~(1|Annee), type="parametric"))$bootSE

ggplot()+xlab("Years")+ylab("Day of year")+
geom_jitter(data=bidon,aes(x=Annee,y=Jour.de.collecte),width=0.3,height=0,size=0.7)+theme_bw()+
theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
geom_ribbon(data=trend,aes(x=Annee,ymin=fit-1.96*se.fit, ymax = fit+1.96*se.fit),alpha=0.5,col=NA,
fill="firebrick3")+
geom_line(data=trend,aes(x=Annee,y=fit),col="firebrick3",size=1.2)





setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/scripts finaux/article/data")
library(ggplot2)
library(mgcv)
library(MASS)
require(doBy)
require(gridExtra)
require(lubridate)
library(chron)
library(dplyr)
library(rgbif)
library(reshape2)
library(car)
library(data.table)
library(lme4)
library(RColorBrewer)
library(phia)
library(ggsignif)
library(blme)
library(glmmTMB)


for(j in seq(0.1,0.9,0.1)){
spani=j
liste2=read.table("yearly_estimates_of_occupancy_and_mfd_only_for_years_withdata.txt",sep="\t",header=T)
liste=read.table("linear_mfd_shifts.txt",sep="\t",header=T,na.strings=c("","NA"))
liste2$species=as.character(liste2$species)
liste$species=as.character(liste$species)
liste2$quant_025[liste2$quant_025==0]=1e-16
liste2$quant_975[liste2$quant_975==1]=1-1e-16
liste2[,c("mean2","quant_0252","quant_9752")]=liste2[,c("mean","quant_025","quant_975")]
liste2[which(liste2$rhat>1.1),c("mean2","quant_0252","quant_9752")]=NA

err=function(x){
vec=c()
for(i in 1:length(x)){
vec[i]=sqrt(x[i]^2+x[i-1]^2)
}
return(vec)}

logit_func=function(x){log(x/(1-x))}
for(i in 1:nrow(liste)){
bidon2=subset(liste2,species==liste$species[i])
wci2=logit_func(bidon2$quant_9752)-logit_func(bidon2$quant_0252)
bidon2$occ_derivs=c(NA,diff(logit_func(bidon2$mean2)))
bidon2$occ_derivs_er=err(wci2)
bidon2$pheno_derivs=c(NA,diff(bidon2$fit))
bidon2$pheno_derivs_er=err(bidon2$se.fit)

wci=bidon2$quant_975-bidon2$quant_025
model3=lm(mean~Annee,data=bidon2,weights=1/wci)
liste$trend_effect[i]=model3$coeff[2]
liste$trend_pval[i]=Anova(model3)[1,4]
liste$stat_trend[i]=if(liste$trend_pval[i]>0.05){"stable"}else{if(liste$trend_effect[i]>0){"increase"}else{"decline"}}
liste$stat_year[i]=if(liste$year_pval[i]>0.05){"unaffected"}else{if(liste$year_effect[i]>0){"delay"}else{"advance"}}
if(i==1){res=bidon2}else{res=rbind(res,bidon2)}
}

tabvar=as.data.frame(fread("belgium_variables_SIG.txt",sep="\t",header=T))
newdat=data.frame(Annee=1902:2016)
tabvar=subset(tabvar,Annee>=1900)
model=loess(value~Annee,data=subset(tabvar,variable=="ratio" & !is.na(value)),span=spani)
tabvarb=cbind(newdat,c(NA,diff(predict(model,newdata=newdat))),"ratio")
names(tabvarb)=c("Annee","value","variable")


model=loess(value~Annee,data=subset(tabvar,variable=="temp" & !is.na(value) & Annee>=1902),span=spani)
tabvarc=cbind(newdat,c(NA,diff(predict(model,newdata=newdat))),"temp_trend")
names(tabvarc)=c("Annee","value","variable")


model=loess(value~Annee,data=subset(tabvar,variable=="temp" & !is.na(value) & Annee>=1902),span=spani)
tabvare=cbind(newdat,c(NA,diff(predict(model,newdata=newdat))),"temp")
names(tabvare)=c("Annee","value","variable")
tabvare$value=c(NA,diff(subset(tabvar,variable=="temp" & !is.na(value) & Annee>=1902)$value))


model=loess(value~Annee,data=subset(tabvar,variable=="urban" & !is.na(value)),span=0.2)
tabvarf=cbind(newdat,c(NA,diff(predict(model,newdata=newdat))),"urban")
names(tabvarf)=c("Annee","value","variable")

tabvar=rbind(tabvarb,tabvarc,tabvare,tabvarf)
tabvar=subset(tabvar,Annee>=1950)
tabvar <- tabvar %>% dplyr::group_by(variable) %>% dplyr::mutate(value=scale(value,center=F,scale=T))

tabvar2=dcast(tabvar,Annee~variable,value.var="value")



final=merge(res,tabvar2,by="Annee")
final=merge(final,liste[,c("species","stat_trend","stat_year")],by="species")


final=as.data.frame(final %>% dplyr::group_by(species) %>% 
dplyr::mutate(ndelta=length(which(!is.na(pheno_derivs) & abs(pheno_derivs)<50)),
ndelta2=length(which(!is.na(occ_derivs) & occ_derivs_er<30))))
final$stat_trend=as.factor(final$stat_trend)
final$stat_year=as.factor(final$stat_year)



bidonb=subset(final,!is.na(pheno_derivs) & ndelta>=25 & abs(pheno_derivs)<50)
bidono=subset(final,!is.na(occ_derivs) & ndelta2>=25 & occ_derivs_er<30)



bidono$Annee2=numFactor(bidono$Annee-1950)

model=glmmTMB(occ_derivs~(ratio+temp+temp_trend+urban)*stat_trend+
ou(Annee2+0 | species),data=bidono,weights=(1/occ_derivs_er)^0.2,
control=glmmTMBControl(optCtrl = list(iter.max=10000000, eval.max=10000000)))



sde=c(sqrt(vcov(model)$cond["ratio","ratio"]),
sqrt(vcov(model)$cond["ratio","ratio"]+2*vcov(model)$cond["ratio","ratio:stat_trendstable"]+vcov(model)$cond["ratio:stat_trendstable","ratio:stat_trendstable"]),
sqrt(vcov(model)$cond["ratio","ratio"]+2*vcov(model)$cond["ratio","ratio:stat_trendincrease"]+vcov(model)$cond["ratio:stat_trendincrease","ratio:stat_trendincrease"]),
sqrt(vcov(model)$cond["urban","urban"]),
sqrt(vcov(model)$cond["urban","urban"]+2*vcov(model)$cond["urban","urban:stat_trendstable"]+vcov(model)$cond["urban:stat_trendstable","urban:stat_trendstable"]),
sqrt(vcov(model)$cond["urban","urban"]+2*vcov(model)$cond["urban","urban:stat_trendincrease"]+vcov(model)$cond["urban:stat_trendincrease","urban:stat_trendincrease"]),
sqrt(vcov(model)$cond["temp","temp"]),
sqrt(vcov(model)$cond["temp","temp"]+2*vcov(model)$cond["temp","temp:stat_trendstable"]+vcov(model)$cond["temp:stat_trendstable","temp:stat_trendstable"]),
sqrt(vcov(model)$cond["temp","temp"]+2*vcov(model)$cond["temp","temp:stat_trendincrease"]+vcov(model)$cond["temp:stat_trendincrease","temp:stat_trendincrease"]),
sqrt(vcov(model)$cond["temp_trend","temp_trend"]),
sqrt(vcov(model)$cond["temp_trend","temp_trend"]+2*vcov(model)$cond["temp_trend","temp_trend:stat_trendstable"]+vcov(model)$cond["temp_trend:stat_trendstable","temp_trend:stat_trendstable"]),
sqrt(vcov(model)$cond["temp_trend","temp_trend"]+2*vcov(model)$cond["temp_trend","temp_trend:stat_trendincrease"]+vcov(model)$cond["temp_trend:stat_trendincrease","temp_trend:stat_trendincrease"]))
est=c(summary(model)$coeff$cond["ratio",1],
summary(model)$coeff$cond["ratio",1]+summary(model)$coeff$cond["ratio:stat_trendstable",1],
summary(model)$coeff$cond["ratio",1]+summary(model)$coeff$cond["ratio:stat_trendincrease",1],
summary(model)$coeff$cond["urban",1],
summary(model)$coeff$cond["urban",1]+summary(model)$coeff$cond["urban:stat_trendstable",1],
summary(model)$coeff$cond["urban",1]+summary(model)$coeff$cond["urban:stat_trendincrease",1],
summary(model)$coeff$cond["temp",1],
summary(model)$coeff$cond["temp",1]+summary(model)$coeff$cond["temp:stat_trendstable",1],
summary(model)$coeff$cond["temp",1]+summary(model)$coeff$cond["temp:stat_trendincrease",1],
summary(model)$coeff$cond["temp_trend",1],
summary(model)$coeff$cond["temp_trend",1]+summary(model)$coeff$cond["temp_trend:stat_trendstable",1],
summary(model)$coeff$cond["temp_trend",1]+summary(model)$coeff$cond["temp_trend:stat_trendincrease",1])
dat1=data.frame(est=est,sde=sde,group=rep(c("decline","stable","increase"),4),
varia=rep(c("Agriculture intensification","Urbanization","Inter-annual temp. changes","Temperature trend"),each=3),model="lmer")
dat1$cate="occupancy"
dat1$lwr=dat1$est-1.96*dat1$sde
dat1$upr=dat1$est+1.96*dat1$sde
dat1$signi=">0.05"
dat1$signi[which(dat1$upr<0)]="<0.05"
dat1$signi[which(dat1$lwr>0)]="<0.05"



rm(model)
bidonb$Annee2=numFactor(bidonb$Annee-1950)
model=glmmTMB(pheno_derivs~(ratio+temp+temp_trend+urban)*stat_year+ou(Annee2+0|species),data=bidonb,weights=1/pheno_derivs_er,
control=glmmTMBControl(optCtrl = list(iter.max=100000000, eval.max=100000000)))




sde=c(sqrt(vcov(model)$cond["ratio","ratio"]),
sqrt(vcov(model)$cond["ratio","ratio"]+2*vcov(model)$cond["ratio","ratio:stat_yearunaffected"]+vcov(model)$cond["ratio:stat_yearunaffected","ratio:stat_yearunaffected"]),
sqrt(vcov(model)$cond["ratio","ratio"]+2*vcov(model)$cond["ratio","ratio:stat_yeardelay"]+vcov(model)$cond["ratio:stat_yeardelay","ratio:stat_yeardelay"]),
sqrt(vcov(model)$cond["urban","urban"]),
sqrt(vcov(model)$cond["urban","urban"]+2*vcov(model)$cond["urban","urban:stat_yearunaffected"]+vcov(model)$cond["urban:stat_yearunaffected","urban:stat_yearunaffected"]),
sqrt(vcov(model)$cond["urban","urban"]+2*vcov(model)$cond["urban","urban:stat_yeardelay"]+vcov(model)$cond["urban:stat_yeardelay","urban:stat_yeardelay"]),
sqrt(vcov(model)$cond["temp","temp"]),
sqrt(vcov(model)$cond["temp","temp"]+2*vcov(model)$cond["temp","temp:stat_yearunaffected"]+vcov(model)$cond["temp:stat_yearunaffected","temp:stat_yearunaffected"]),
sqrt(vcov(model)$cond["temp","temp"]+2*vcov(model)$cond["temp","temp:stat_yeardelay"]+vcov(model)$cond["temp:stat_yeardelay","temp:stat_yeardelay"]),
sqrt(vcov(model)$cond["temp_trend","temp_trend"]),
sqrt(vcov(model)$cond["temp_trend","temp_trend"]+2*vcov(model)$cond["temp_trend","temp_trend:stat_yearunaffected"]+vcov(model)$cond["temp_trend:stat_yearunaffected","temp_trend:stat_yearunaffected"]),
sqrt(vcov(model)$cond["temp_trend","temp_trend"]+2*vcov(model)$cond["temp_trend","temp_trend:stat_yeardelay"]+vcov(model)$cond["temp_trend:stat_yeardelay","temp_trend:stat_yeardelay"]))
est=c(summary(model)$coeff$cond["ratio",1],
summary(model)$coeff$cond["ratio",1]+summary(model)$coeff$cond["ratio:stat_yearunaffected",1],
summary(model)$coeff$cond["ratio",1]+summary(model)$coeff$cond["ratio:stat_yeardelay",1],
summary(model)$coeff$cond["urban",1],
summary(model)$coeff$cond["urban",1]+summary(model)$coeff$cond["urban:stat_yearunaffected",1],
summary(model)$coeff$cond["urban",1]+summary(model)$coeff$cond["urban:stat_yeardelay",1],
summary(model)$coeff$cond["temp",1],
summary(model)$coeff$cond["temp",1]+summary(model)$coeff$cond["temp:stat_yearunaffected",1],
summary(model)$coeff$cond["temp",1]+summary(model)$coeff$cond["temp:stat_yeardelay",1],
summary(model)$coeff$cond["temp_trend",1],
summary(model)$coeff$cond["temp_trend",1]+summary(model)$coeff$cond["temp_trend:stat_yearunaffected",1],
summary(model)$coeff$cond["temp_trend",1]+summary(model)$coeff$cond["temp_trend:stat_yeardelay",1])
dat2=data.frame(est=est,sde=sde,group=rep(c("advance","unaffected","delay"),4),
varia=rep(c("Agriculture intensification","Urbanization","Inter-annual temp. changes","Temperature trend"),each=3),model="lmer")
dat2$cate="phenology"
dat2$lwr=dat2$est-1.96*dat2$sde
dat2$upr=dat2$est+1.96*dat2$sde
dat2$signi=">0.05"
dat2$signi[which(dat2$upr<0)]="<0.05"
dat2$signi[which(dat2$lwr>0)]="<0.05"


b=rbind(dat1,dat2)
b$nderivs_pheno=nrow(bidonb)
b$nderivs_occ=nrow(bidono)
b$spani=j


if(j==0.1){bf=b}else{bf=rbind(bf,b)}


}



bf$varia=factor(bf$varia,c("Inter-annual temp. changes","Temperature trend","Urbanization","Agriculture intensification"))
bf$moy=bf$est
ponds1=unique(bidono[,c("species","ndelta2","stat_trend"),]) %>%
dplyr::group_by(stat_trend) %>% dplyr::summarise(n=length(species))
ponds1$cate="occupancy"
names(ponds1)[1]="group"
ponds2=unique(bidonb[,c("species","ndelta","stat_year"),]) %>%
 dplyr::group_by(stat_year) %>% dplyr::summarise(n=length(species))
ponds2$cate="phenology"
names(ponds2)[1]="group"
bf=merge(bf,rbind(ponds1,ponds2),by=c("group","cate"))
bf$stat=factor(bf$group,c("decline","advance","stable","unaffected","increase","delay"))
bf=bf[order(bf$stat),]
bf$stat2=bf$stat
bf$stat=paste0(bf$stat," (n=",bf$n,")")
bf$stat=factor(bf$stat,unique(bf$stat))
bf$signi=factor(bf$signi,c(">0.05","<0.05"))
bf$star="ns"
bf$star[which(bf$pvalinter<0.05)]="*"
bf$star[which(bf$pvalinter<0.01)]="**"
bf$star[which(bf$pvalinter<0.001)]="***"
labo=unique(bf[,c("cate","varia","star")])
labo=labo[order(labo$varia),]

pl1=ggplot(data=subset(bf,cate=="occupancy" & varia!="inter"),aes(x=as.factor(spani),y=moy,col=stat,shape=signi))+
geom_hline(yintercept=0,size=1.2)+
scale_shape_manual(values=c(19,21),guide=F)+scale_shape_manual(values=c(21,19),guide=F,na.value=15,drop=F)+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge(width = 0.50),fill="white")+
theme_bw()+ylab("Standardised effects")+
theme(panel.grid.minor=element_blank(),plot.title=element_text(size=14,face="bold"),
legend.title = element_blank(),axis.title.x=element_blank(),
strip.background = element_blank(),legend.position="bottom")+
scale_color_discrete()+ggtitle("a")+xlab("Maximum time-lag allowed (in years)")+
scale_colour_manual(values=c("darkorchid4","dodgerblue3",
"azure4"))+facet_wrap(~varia,nrow=1)


pl2=ggplot(data=subset(bf,cate=="phenology" & varia!="inter"),aes(x=as.factor(spani),y=moy,col=stat,shape=signi))+
geom_hline(yintercept=0,size=1.2)+
geom_pointrange(aes(ymin=lwr,ymax=upr),position=position_dodge(width = 0.50),fill="white")+
scale_shape_manual(values=c(19,21),guide=F)+scale_shape_manual(values=c(21,19),guide=F,na.value=15,drop=F)+
theme_bw()+ylab("Standardised effects \n")+
theme(panel.grid.minor=element_blank(),plot.title=element_text(size=14,face="bold"),legend.position="bottom",axis.title.x=element_blank(),
strip.background = element_blank(),legend.title = element_blank())+
scale_color_discrete()+ggtitle("b")+xlab("Maximum time-lag allowed (in years)")+
scale_colour_manual(values=c("firebrick4","orange","lemonchiffon3"))+
facet_wrap(~varia,nrow=1)
gridExtra::grid.arrange(pl1,pl2,bottom="Smoothing parameter value",nrow=2)


png(paste0("fig_s5.png"),width=1200,height=1000,res=140)
gridExtra::grid.arrange(pl1,pl2,bottom="Smoothing parameter value",nrow=2)
dev.off();


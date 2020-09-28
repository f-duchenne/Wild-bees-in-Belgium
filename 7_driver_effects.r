setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/scripts finaux/article/data")
library(ggplot2)
library(mgcv)
library(MASS)
require(doBy)
require(ggplot2)
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
library(performance)
library(glmmTMB)
library(optimx)
library(blme)


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

spani=0.5
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


#####################
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



ano1=Anova(model)
ano1$cate="occupancy"

newdat=data.frame(urban=seq(min(bidono$urban),max(bidono$urban),length.out=100),
temp=mean(bidono$temp),ratio=mean(bidono$ratio),stat_trend=rep(c("decline","increase","stable"),each=100),
signi=rep(dat1$signi[dat1$varia=="Urbanization"],each=100),temp_trend=mean(bidono$temp_trend),
delta="Urbanization")
newdat$var=newdat$urban
newdat2=data.frame(urban=mean(bidono$urban),
temp=seq(min(bidono$temp),max(bidono$temp),length.out=100),ratio=mean(bidono$ratio),
stat_trend=rep(c("decline","stable","increase"),each=100),
signi=rep(dat1$signi[dat1$varia=="Inter-annual temp. changes"],each=100),temp_trend=mean(bidono$temp_trend),
delta="Inter-annual temp. changes")
newdat2$var=newdat2$temp
newdat3=data.frame(urban=mean(bidono$urban),
temp=mean(bidono$temp),ratio=seq(min(bidono$ratio),max(bidono$ratio),length.out=10),
stat_trend=rep(c("decline","stable","increase"),each=100),temp_trend=mean(bidono$temp_trend),
signi=rep(dat1$signi[dat1$varia=="Agriculture intensification"],each=100),
delta="Agriculture intensification")
newdat3$var=newdat3$ratio
newdat4=data.frame(urban=mean(bidono$urban),
temp=mean(bidono$temp),ratio=mean(bidono$ratio),stat_trend=rep(c("decline","stable","increase"),each=100),
signi=rep(dat1$signi[dat1$varia=="Temperature trend"],each=100),temp_trend=seq(min(bidono$temp_trend),max(bidono$temp_trend),length.out=100),
delta="Temperature trend")
newdat4$var=newdat4$temp_trend

newdat=rbind(newdat,newdat2,newdat3,newdat4)
newdat$species=NA#as.factor("Bombus lapidarius")
newdat$Annee2=NA#numFactor(30)
newdat$Annee=NA#numFactor(30)
newdat$occ_derivs_er=mean(bidono$occ_derivs_er)
newdat$wei=mean(bidono$wei)
newdat=cbind(newdat,predict(model,newdata=newdat,se.fit=T,allow.new.levels=T))
newdat$star="ns"
newdat$star[which(newdat$pvalinter<0.05)]="*"
newdat$star[which(newdat$pvalinter<0.01)]="**"
newdat$star[which(newdat$pvalinter<0.001)]="***"
newdat$stat_trend=factor(newdat$stat_trend,c("decline","stable","increase"))
newdat$delta=factor(newdat$delta,c("Inter-annual temp. changes","Temperature trend","Urbanization","Agriculture intensification"))
ponds1=unique(bidono[,c("species","ndelta2","stat_trend"),]) %>% dplyr::group_by(stat_trend) %>% dplyr::summarise(n=length(species))
newdat=merge(newdat,ponds1,by=c("stat_trend"))
newdat$stat=paste0(newdat$stat," (n=",newdat$n,")")
newdat=newdat[order(newdat$stat_trend),]
newdat$stat=factor(newdat$stat,unique(newdat$stat))
newdat$signi=factor(newdat$signi,c("<0.05",">0.05"))


pl1=ggplot(data=newdat,aes(x=var,y=fit,col=stat,fill=stat,linetype=signi))+
geom_hline(yintercept=0,size=0.6)+
geom_line(size=1)+geom_ribbon(aes(ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit,x=var),col=NA,alpha=0.2)+
scale_color_manual(values=c("darkorchid4","dodgerblue3",
"azure4"))+
scale_fill_manual(values=c("darkorchid4","dodgerblue3",
"azure4"))+
scale_linetype_manual(values=c("solid","dashed"),guide=F,drop=F)+
#stat_summary(fun.data = n_fun, geom = "text",position=position_dodge(width=0.2),size=5)+
theme_bw()+ylab(paste0("Yearly occupancy changes (logit scale)"))+
theme(panel.grid=element_blank(),strip.background=element_blank(),
plot.title=element_text(size=14,face="bold"),legend.position="right",legend.title=element_blank(),
strip.placement = "outside")+
xlab("Yearly changes")+ggtitle("a")+
facet_wrap(~delta,scales="free_x",nrow=1,strip.position = "bottom")
model2=model



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
dat1=data.frame(est=est,sde=sde,group=rep(c("advance","unaffected","delay"),4),
varia=rep(c("Agriculture intensification","Urbanization","Inter-annual temp. changes","Temperature trend"),each=3),model="lmer")
dat1$cate="occupancy"
dat1$lwr=dat1$est-1.96*dat1$sde
dat1$upr=dat1$est+1.96*dat1$sde
dat1$signi=">0.05"
dat1$signi[which(dat1$upr<0)]="<0.05"
dat1$signi[which(dat1$lwr>0)]="<0.05"


newdat=data.frame(urban=seq(min(bidono$urban),max(bidono$urban),length.out=100),
temp=mean(bidono$temp),ratio=mean(bidono$ratio),stat_year=rep(c("advance","unaffected","delay"),each=100),
signi=rep(dat1$signi[dat1$varia=="Urbanization"],each=100),temp_trend=mean(bidono$temp_trend),
delta="Urbanization")
newdat$var=newdat$urban
newdat2=data.frame(urban=mean(bidono$urban),
temp=seq(min(bidono$temp),max(bidono$temp),length.out=100),ratio=mean(bidono$ratio),
stat_year=rep(c("advance","unaffected","delay"),each=100),
signi=rep(dat1$signi[dat1$varia=="Inter-annual temp. changes"],each=100),temp_trend=mean(bidono$temp_trend),
delta="Inter-annual temp. changes")
newdat2$var=newdat2$temp
newdat3=data.frame(urban=mean(bidono$urban),
temp=mean(bidono$temp),ratio=seq(min(bidono$ratio),max(bidono$ratio),length.out=10),
stat_year=rep(c("advance","unaffected","delay"),each=100),temp_trend=mean(bidono$temp_trend),
signi=rep(dat1$signi[dat1$varia=="Agriculture intensification"],each=100),
delta="Agriculture intensification")
newdat3$var=newdat3$ratio
newdat4=data.frame(urban=mean(bidono$urban),
temp=mean(bidono$temp),ratio=mean(bidono$ratio),stat_year=rep(c("advance","unaffected","delay"),each=100),
signi=rep(dat1$signi[dat1$varia=="Temperature trend"],each=100),temp_trend=seq(min(bidono$temp_trend),max(bidono$temp_trend),length.out=100),
delta="Temperature trend")
newdat4$var=newdat4$temp_trend

newdat=rbind(newdat,newdat2,newdat3,newdat4)
newdat$species=NA#as.factor("Bombus lapidarius")
newdat$Annee2=NA#numFactor(30)
newdat$Annee=NA#numFactor(30)
newdat$pheno_derivs_er=mean(bidono$occ_derivs_er)
newdat$wei=mean(bidonb$wei)
newdat=cbind(newdat,predict(model,newdata=newdat,se.fit=T,allow.new.levels=T))
newdat$star="ns"
newdat$star[which(newdat$pvalinter<0.05)]="*"
newdat$star[which(newdat$pvalinter<0.01)]="**"
newdat$star[which(newdat$pvalinter<0.001)]="***"
newdat$stat_year=factor(newdat$stat_year,c("advance","unaffected","delay"))
newdat$delta=factor(newdat$delta,c("Inter-annual temp. changes","Temperature trend","Urbanization","Agriculture intensification"))
ponds1=unique(bidonb[,c("species","ndelta","stat_year"),]) %>% dplyr::group_by(stat_year) %>% dplyr::summarise(n=length(species))
newdat=merge(newdat,ponds1,by=c("stat_year"))
newdat$stat=paste0(newdat$stat," (n=",newdat$n,")")
newdat=newdat[order(newdat$stat_year),]
newdat$stat=factor(newdat$stat,unique(newdat$stat))
newdat$signi=factor(newdat$signi,c("<0.05",">0.05"))

pl2=ggplot(data=newdat,aes(x=var,y=fit,col=stat,fill=stat,linetype=signi))+
geom_hline(yintercept=0,size=0.6)+
geom_line(size=1)+geom_ribbon(aes(ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit,x=var),col=NA,alpha=0.2)+
scale_color_manual(values=c("firebrick4","orange","lemonchiffon3"))+
scale_fill_manual(values=c("firebrick4","orange","lemonchiffon3"))+
scale_linetype_manual(values=c("solid","dashed"),guide=F,drop=F)+
#stat_summary(fun.data = n_fun, geom = "text",position=position_dodge(width=0.2),size=5)+
theme_bw()+ylab(paste0("Yearly MFD changes (days) \n"))+
theme(panel.grid=element_blank(),strip.background=element_blank(),
plot.title=element_text(size=14,face="bold"),legend.position="right",legend.title=element_blank(),
strip.placement = "outside")+
xlab("Yearly changes")+ggtitle("b")+
facet_wrap(~delta,scales="free_x",nrow=1,strip.position = "bottom")

gridExtra::grid.arrange(pl1,pl2,ncol=1)


png(paste0("fig5.png"),width=1400,height=1000,res=160)
gridExtra::grid.arrange(pl1,pl2,ncol=1)
dev.off();

pdf(paste0("fig5.pdf"),width=10,height=8)
gridExtra::grid.arrange(pl1,pl2,ncol=1)
dev.off();



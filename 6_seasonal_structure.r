setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/scripts finaux/article/data")
require(ggplot2)
require(gridExtra)
library(dplyr)
library(reshape2)
library(car)
library(data.table)
library(RColorBrewer)
library(Hmisc)
library(chron)

liste2=read.table("yearly_estimates_of_occupancy_and_mfd.txt",sep="\t",header=T)
traits=read.table("traits_full_new_liste_MG.txt",sep="\t",header=T)
liste2=merge(liste2,traits,by="species")
bide=data.frame(Annee_group2=cut(1950:2016,seq(1950,2020,10),include.lowest=T,right = FALSE),Annee=seq(1950,2016,1))
bide$Annee_group2=rep(c("1950-1959","1960-1969","1970-1979","1980-1989","1990-1999","2000-2009","2010-2016"),each=10)[1:67]
liste2=merge(liste2,bide,by="Annee")

liste2=liste2 %>%
group_by(Annee_group2,species,Pheno_length) %>%
summarise(abd=mean(mean,na.rm=T),mfd=mean(fit,na.rm=T))
liste2=liste2 %>% group_by(species,Pheno_length) %>%
mutate(abond_ini=abd[Annee_group2=="1950-1959"],
mfd_ini=mfd[Annee_group2=="1950-1959"])


for(i in 1:nrow(liste2)){
abond=dnorm(1:365,liste2$mfd[i],liste2$Pheno_length[i])*liste2$abd[i]
bidon=data.frame(species=liste2$species[i],Annee_group2=liste2$Annee_group2[i],abondance=abond,jour=1:365)
if(i==1){pheno_masse=bidon}else{pheno_masse=rbind(pheno_masse,bidon)}
}

pheno_masse=subset(pheno_masse,!is.na(abondance))
b=pheno_masse %>%
group_by(jour,Annee_group2) %>%
summarise(pres.sum=sum(abondance,na.rm=T),div=vegan::diversity(abondance),
div2=vegan::diversity(abondance,index="simpson"),vari=var(abondance,na.rm=T))
b$Month=paste0(month.abb[month.day.year(b$jour,  c(month = 1, day = 1, year = 1970))$month],".")
moy1=wtd.mean(b$jour[b$Annee_group2=="1950-1959"],b$pres.sum[b$Annee_group2=="1950-1959"])
moy2=wtd.mean(b$jour[b$Annee_group2=="2010-2016"],b$pres.sum[b$Annee_group2=="2010-2016"])
moy1-moy2
pl1=ggplot(b,aes(x=jour,y=pres.sum))+geom_line(aes(color=as.factor(Annee_group2)),size=1.2)+
theme_bw()+xlab("Julian days")+ylab("Total daily occupancy")+labs(color="Years")+
theme(legend.position="right",axis.ticks.x=element_blank(),axis.title.x=element_blank(),
panel.grid.minor = element_blank(),plot.title=element_text(size=14,face="bold"),panel.grid.major = element_blank())+
scale_colour_manual(values=brewer.pal(9,"Greens")[3:9])+ggtitle("a")+
geom_vline(xintercept=c(moy1,moy2),linetype="dashed",color=brewer.pal(9,"Greens")[c(5,9)])+
scale_x_continuous(breaks=seq(15,365,30),labels=unique(b$Month))

for(i in 1:nrow(liste2)){
abond=dnorm(1:365,liste2$mfd[i],liste2$Pheno_length[i])*liste2$abond_ini[i]
bidon=data.frame(species=liste2$species[i],Annee_group2=liste2$Annee_group2[i],abondance=abond,jour=1:365)
if(i==1){pheno_masse=bidon}else{pheno_masse=rbind(pheno_masse,bidon)}
}
library(Hmisc)
pheno_masse=subset(pheno_masse,!is.na(abondance))
b=pheno_masse %>%
group_by(jour,Annee_group2) %>%
summarise(pres.sum=sum(abondance,na.rm=T),div=vegan::diversity(abondance),
div2=vegan::diversity(abondance,index="simpson"),vari=var(abondance,na.rm=T))
b$Month=paste0(month.abb[month.day.year(b$jour,  c(month = 1, day = 1, year = 1970))$month],".")
moy1=wtd.mean(b$jour[b$Annee_group2=="1950-1959"],b$pres.sum[b$Annee_group2=="1950-1959"])
moy2=wtd.mean(b$jour[b$Annee_group2=="2010-2016"],b$pres.sum[b$Annee_group2=="2010-2016"])
moy1-moy2
b$Month=paste0(month.abb[month.day.year(b$jour,  c(month = 1, day = 1, year = 1970))$month],".")
pl1_sans_abond=ggplot(b,aes(x=jour,y=pres.sum))+geom_line(aes(color=as.factor(Annee_group2)),size=1.2)+
theme_bw()+xlab("Julian days")+ylab("Total daily occupancy")+labs(color="Years")+
theme(legend.position="right",axis.ticks.x=element_blank(),axis.title.x=element_blank(),
panel.grid.minor = element_blank(),plot.title=element_text(size=14,face="bold"),panel.grid.major = element_blank())+
scale_colour_manual(values=brewer.pal(9,"Greens")[3:9])+ggtitle("b")+
geom_vline(xintercept=c(moy1,moy2),linetype="dashed",color=brewer.pal(9,"Greens")[c(5,9)])+
scale_x_continuous(breaks=seq(15,365,30),labels=unique(b$Month))



for(i in 1:nrow(liste2)){
abond=dnorm(1:365,liste2$mfd_ini[i],liste2$Pheno_length[i])*liste2$abd[i]
bidon=data.frame(species=liste2$species[i],Annee_group2=liste2$Annee_group2[i],abondance=abond,jour=1:365)
if(i==1){pheno_masse=bidon}else{pheno_masse=rbind(pheno_masse,bidon)}
}


library(Hmisc)
pheno_masse=subset(pheno_masse,!is.na(abondance))
b=pheno_masse %>%
group_by(jour,Annee_group2) %>%
summarise(pres.sum=sum(abondance,na.rm=T),div=vegan::diversity(abondance),
div2=vegan::diversity(abondance,index="simpson"),vari=var(abondance,na.rm=T))
b$Month=paste0(month.abb[month.day.year(b$jour,  c(month = 1, day = 1, year = 1970))$month],".")
moy1=wtd.mean(b$jour[b$Annee_group2=="1950-1959"],b$pres.sum[b$Annee_group2=="1950-1959"])
moy2=wtd.mean(b$jour[b$Annee_group2=="2010-2016"],b$pres.sum[b$Annee_group2=="2010-2016"])
moy1-moy2
b$Month=paste0(month.abb[month.day.year(b$jour,  c(month = 1, day = 1, year = 1970))$month],".")
pl1_sans_mfd=ggplot(b,aes(x=jour,y=pres.sum))+geom_line(aes(color=as.factor(Annee_group2)),size=1.2)+
theme_bw()+xlab("Julian days")+ylab("Total daily occupancy")+labs(color="Years")+
theme(legend.position="right",axis.ticks.x=element_blank(),axis.title.x=element_blank(),
panel.grid.minor = element_blank(),plot.title=element_text(size=14,face="bold"),panel.grid.major = element_blank())+
scale_colour_manual(values=brewer.pal(9,"Greens")[3:9])+ggtitle("c")+
geom_vline(xintercept=c(moy1,moy2),linetype="dashed",color=brewer.pal(9,"Greens")[c(5,9)])+
scale_x_continuous(breaks=seq(15,365,30),labels=unique(b$Month))

grid.arrange(pl1,pl1_sans_abond,pl1_sans_mfd)

pdf(paste0("fig 4.pdf"),width=8,height=10)
grid.arrange(pl1,pl1_sans_abond,pl1_sans_mfd)
dev.off();


####SUpplementary #########
for(i in 1:nrow(liste2)){
abond=dnorm(1:365,liste2$mfd[i],liste2$Pheno_length[i])*liste2$abd[i]
bidon=data.frame(species=liste2$species[i],Annee_group2=liste2$Annee_group2[i],abondance=abond,jour=1:365)
if(i==1){pheno_masse=bidon}else{pheno_masse=rbind(pheno_masse,bidon)}
}

b2=subset(pheno_masse,abondance>0.002) %>%
group_by(jour,Annee_group2) %>%
summarise(pres.sum=sum(abondance,na.rm=T),div=vegan::diversity(abondance),
div2=vegan::diversity(abondance,index="simpson"),vari=var(abondance,na.rm=T),n=length(species))
b2$Month=paste0(month.abb[month.day.year(b2$jour,  c(month = 1, day = 1, year = 1970))$month],".")
pl2=ggplot(subset(b2,jour>=50 & jour<=300),aes(x=jour,y=n))+geom_line(aes(color=as.factor(Annee_group2)),size=1.2)+theme_bw()+xlab("Julian days")+
ylab("Species richness")+labs(color="Years")+theme(legend.position="right",axis.title.x=element_blank(),axis.ticks.x=element_blank(),
panel.grid.minor = element_blank(),plot.title=element_text(size=14,face="bold"),panel.grid.major = element_blank())+
scale_colour_manual(values=brewer.pal(9,"Greens")[3:9])+ggtitle("b")+
scale_x_continuous(breaks=seq(60,300,30),labels=unique(subset(b2,jour>=50 & jour<=300)$Month))


b=pheno_masse %>%
dplyr::group_by(jour,Annee_group2) %>%
dplyr::summarise(pres.sum=sum(abondance,na.rm=T),div=vegan::diversity(abondance),
div2=vegan::diversity(abondance,index="simpson"),vari=var(abondance,na.rm=T))
for(i in seq(0,0.1,0.001)){
bibi=subset(b,pres.sum>=i)
bibi1=bibi %>% group_by("Annee_group2") %>% count(Annee_group2)
bibi1$threshold=i
if(i==0){resi=bibi1}else{resi=rbind(resi,bibi1)}
}
resi=resi[resi$Annee_group2 %in% c("1950-1959","2010-2016"),]
pl3=ggplot(subset(resi,threshold==0.01 | threshold==0.05),aes(x=as.factor(threshold),y=n,color=as.factor(Annee_group2)))+
geom_point(size=1.5,position=position_dodge(width=0.5))+
theme_bw()+xlab("Occupancy threshold")+
ylab("Pollination season duration (in days)")+labs(color="Years")+theme(
panel.grid.minor = element_blank(),plot.title=element_text(size=14,face="bold"),panel.grid.major = element_blank(),
legend.position="none")+
scale_colour_manual(values=brewer.pal(9,"Greens")[c(5,9)])+ggtitle("a")
grid.arrange(pl3,pl2,ncol=2,widths=c(1,2.5))

png(paste0("fig s7.png"),width=1000,height=600,res=130)
grid.arrange(pl3,pl2,ncol=2,widths=c(1,2.5))
dev.off();
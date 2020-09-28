setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/scripts finaux/article/data")
library(ggplot2)
require(doBy)
require(ggplot2)
require(gridExtra)
require(lubridate)
library(chron)
library(dplyr)
library(chron) 
library(reshape2)
library(car)
library(data.table)
library(RColorBrewer)

liste2=read.table("yearly_estimates_of_occupancy_and_mfd.txt",sep="\t",header=T)
liste=read.table("linear_mfd_shifts.txt",sep="\t",header=T,na.strings=c("","NA"))
liste=subset(liste,!is.na(year_effect))
liste2$species=as.character(liste2$species)
liste$species=as.character(liste$species)

for(i in 1:nrow(liste)){
bidon2=subset(liste2,species==liste$species[i])
wci=bidon2$quant_975-bidon2$quant_025
model3=lm(mean~Annee,data=bidon2,weights=1/wci)
liste$trend_effect[i]=model3$coeff[2]
liste$trend_err[i]=summary(model3)$coeff[2,2]
liste$trend_pval[i]=Anova(model3)[1,4]
liste$stat_trend[i]=if(liste$trend_pval[i]>0.05){"stable"}else{if(liste$trend_effect[i]>0){"increase"}else{"decline"}}
liste$stat_year[i]=if(liste$year_pval[i]>0.05){"unaffected"}else{if(liste$year_effect[i]>0){"delay"}else{"advance"}}
}

summary(liste$trend_effect)
library(plotrix)
std.error(liste$trend_effect)
dim(subset(liste,trend_pval<0.05 & trend_effect<0))
summary(liste$trend_effect*66/liste$prob50)
liste2=merge(liste2,liste[,c("species","stat_trend","stat_year")])
liste2=liste2 %>% dplyr::group_by(stat_year) %>% dplyr::mutate(n1=length(unique(species)))
liste2=liste2 %>% dplyr::group_by(stat_trend) %>% dplyr::mutate(n2=length(unique(species)))
liste2$stat_year=paste0(liste2$stat_year," (n=",liste2$n1,")")
liste2$stat_trend=paste0(liste2$stat_trend," (n=",liste2$n2,")")
#PLOT TENDENCES
b1=Rmisc::group.CI(fit~Annee+stat_year,data=liste2)
names(b1)[3:5]=c("upper","mean","lower")
b1$var="phenology"
b1b=Rmisc::group.CI(fit~Annee,data=liste2)
names(b1b)[2:4]=c("upper","mean","lower")
b2=Rmisc::group.CI(mean~Annee+stat_trend,data=liste2)
names(b2)[3:5]=c("upper","mean","lower")
b2$var="occupancy"
b2$stat_trend=factor(b2$stat_trend,c("decline (n=125)","stable (n=30)","increase (n=50)"))
b1$stat_year=factor(b1$stat_year,c("advance (n=83)","unaffected (n=96)","delay (n=26)"))
b2b=Rmisc::group.CI(mean~Annee,data=liste2)
names(b2b)[2:4]=c("upper","mean","lower")

pl1=ggplot(data=b2,aes(x=Annee,y=mean))+
geom_ribbon(aes(ymin=lower,ymax=upper,color=stat_trend,fill=stat_trend),alpha=0.3,colour = NA)+
geom_line(size=1.2,aes(color=stat_trend,fill=stat_trend))+
theme_bw()+xlab("Years")+
ylab("Occupancy probability")+ggtitle("a")+
theme(plot.title=element_text(size=14,face="bold"),panel.grid=element_blank(),legend.position="bottom")+
labs(color="",fill="")+
scale_color_manual(values=c("darkorchid4","dodgerblue3","azure4"))+
scale_fill_manual(values=c("darkorchid4","dodgerblue3","azure4"))+labs(color="",fill="")+
# geom_ribbon(data=b2b,aes(ymin=lower,ymax=upper),alpha=0.3,colour=NA,fill="grey")+
geom_line(data=b2b,aes(ymin=lower,ymax=upper),size=1.2,col="black")

pl2=ggplot(data=b1,aes(x=Annee,y=mean))+
geom_ribbon(aes(ymin=lower,ymax=upper,color=stat_year,fill=stat_year),alpha=0.3,colour = NA)+
geom_line(size=1.2,aes(color=stat_year,fill=stat_year))+
theme_bw()+xlab("Years")+
ylab("MFD (julian days)")+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold"),panel.grid=element_blank(),legend.position="bottom")+
scale_color_manual(values=c("firebrick4","orange","lemonchiffon3"))+
scale_fill_manual(values=c("firebrick4","orange","lemonchiffon3"))+labs(color="",fill="")+
# geom_ribbon(data=b1b,aes(ymin=lower,ymax=upper),alpha=0.3,colour=NA,fill="grey")+
geom_line(data=b1b,aes(ymin=lower,ymax=upper),size=1.2,col="black")

#PLOT VARIABLES
tabvar=as.data.frame(fread("belgium_variables_SIG.txt",sep="\t",header=T))
df_scaled <- tabvar %>% dplyr::group_by(variable) %>% dplyr::mutate(value2=scale(value))

df_scaled$variable=factor(df_scaled$variable,c("temp","urban","ratio","fertilizer"))

biche=subset(df_scaled,Annee>=1930 & variable=="temp")
biche$variable="temp_trend"

pl11=ggplot(subset(df_scaled,Annee>=1930 & variable!="fertilizer"),aes(col=variable,linetype=variable,shape=variable))+
geom_point(data=subset(df_scaled,Annee>=1930 & variable!="fertilizer"),aes(x=Annee,y=value2,color=variable),fill="white")+
stat_smooth(data=subset(df_scaled,Annee>=1930 & variable=="ratio"),size=1,aes(x=Annee,y=value2),se=F,alpha=0.2,span=0.5)+theme_bw()+
stat_smooth(data=subset(df_scaled,Annee>=1930 & variable=="urban"),size=1,aes(x=Annee,y=value2),fill="white",alpha=0.2,span=0.2)+
stat_smooth(data=biche,size=1,aes(x=Annee,y=value2),fill="white",alpha=0.2,span=0.5)+
xlab("Years")+ylab("scaled values")+
scale_colour_manual(values=c("dodgerblue4",brewer.pal(1,"Set1")[1],"firebrick",brewer.pal(1,"Set1")[3]),
labels=c("Agriculture intensity","Temperature","Temperature trend","Total built urban area"))+
theme(legend.position="bottom",plot.title=element_text(size=14,face="bold"),panel.grid=element_blank())+ggtitle("c")+
scale_linetype_manual("",values=c(1,0,1,1),labels=c("Agriculture intensity","Temperature","Temperature trend","Total built urban area"))+
scale_shape_manual("",values=c(16,16,NA,16),labels=c("Agriculture intensity","Temperature","Temperature trend","Total built urban area"))+
labs(color="")+guides(color=guide_legend(nrow=2))


lay=cbind(c(2,3),c(1,1))
gridExtra::grid.arrange(pl11,pl1,pl2,layout_matrix=lay)
png("fig2.png",width=1100,height=800,res=130)
gridExtra::grid.arrange(pl11,pl1,pl2,layout_matrix=lay)
dev.off();

##### FIGURE S6:

pl1=ggplot(liste,aes(x=trend_effect,y=year_effect))+
theme_bw()+xlab("Occupancy linear trends (/year)")+
ylab("MFD linear trends (days/year)")+labs(color="")+
geom_errorbarh(aes(xmin=trend_effect-trend_err,xmax=trend_effect+trend_err),col="grey")+
geom_errorbar(aes(ymin=year_effect-year_err,ymax=year_effect+year_err),col="grey")+
theme(legend.position="right",panel.grid.minor = element_blank(),
plot.title=element_blank(),panel.grid.major = element_blank())+
geom_point()+
geom_hline(yintercept=0,linetype="dashed",color="black")+
geom_vline(xintercept=0,linetype="dashed",color="black")+
guides(color = guide_legend(override.aes = list(linetype = 0)))+ylim(c(-1.4,1.3))+
xlim(c(-0.014,0.012))

library(ggpubr)
library(patchwork)
pl1h=ggplot(liste,aes(x=trend_effect,color=stat_year,fill=stat_year))+
stat_density(alpha=0.4,trim=F,position = "identity",adjust=1)+scale_color_manual(values=c("firebrick4","orange","lemonchiffon3"))+
scale_fill_manual(values=c("firebrick4","orange","lemonchiffon3"))+theme_bw()+
theme(legend.position="top",panel.grid.minor = element_blank(),
plot.title=element_text(size=14,face="bold"),panel.grid.major = element_blank(),axis.text=element_blank(),
axis.title=element_blank(),panel.border=element_blank(),axis.ticks=element_blank())+
labs(fill="",color="")+xlim(c(-0.014,0.012))+ggtitle("")
pl1c=ggplot(liste,aes(x=year_effect,color=stat_trend,fill=stat_trend))+
stat_density(alpha=0.4,trim=F,position = "identity",adjust=1)+scale_color_manual(values=c("darkorchid4","dodgerblue3","azure4"))+
scale_fill_manual(values=c("darkorchid4","dodgerblue3","azure4"))+theme_bw()+
theme(legend.position="right",panel.grid.minor = element_blank(),
plot.title=element_blank(),panel.grid.major = element_blank(),axis.text=element_blank(),
axis.title=element_blank(),panel.border=element_blank(),axis.ticks=element_blank())+
labs(fill="",color="")+xlim(c(-1.4,1.3))+coord_flip()

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(), 
   panel.border = element_blank(),
   panel.background = element_blank(),
   axis.title.x = element_blank(),
   axis.title.y = element_blank(),
   axis.text.x = element_blank(), 
   axis.text.y = element_blank(),
   axis.ticks = element_blank()
     )

pl1h + plot_spacer() + pl1+ pl1c + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1.1, 4))


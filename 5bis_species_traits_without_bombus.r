#PHYLOGENIE
setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/scripts finaux/article/data")
library(ape)
library(picante)
library(doBy)
library(ggtree)
library(phytools)
library(phylotools)
library(caper)
library(phangorn)
library(ggtree)
library(mgcv)
library(car)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

anova.f <- function (object){
	data <- object$data
	tlabels <- attr(terms(object$formula), "term.labels")
	k <- object$k
	n <- object$n
	NR <- length(tlabels) + 1
	rss <- resdf <- rep(NA, NR)
	rss[1] <- object$NSSQ
	resdf[1] <- n - 1
	lm <- object$param["lambda"]
	dl <- object$param["delta"]
	kp <- object$param["kappa"]
	for (i in 1:length(tlabels)) {
		fmla <- as.formula(paste(object$namey, " ~ ", paste(tlabels[1:i], collapse = "+")))
		plm <- pgls(fmla, data, lambda = lm, delta = dl, kappa = kp)
		rss[i + 1] <- plm$RSSQ
		resdf[i + 1] <- (n - 1) - plm$k + 1
	}
	ss <- c(abs(diff(rss)), object$RSSQ)
	df <- c(abs(diff(resdf)), n - k)
	ms <- ss/df
	fval <- ms/ms[NR]
	P <- pf(fval, df, df[NR], lower.tail = FALSE)
	table <- data.frame(df, ss, ms, f = fval, P)
	table[length(P), 4:5] <- NA
	dimnames(table) <- list(c(tlabels, "Residuals"), c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
	structure(table, heading = c("Analysis of Variance Table", sprintf("Sequential SS for pgls: lambda = %0.2f, delta = %0.2f, kappa = %0.2f\n", lm, dl, kp), paste("Response:", deparse(formula(object)[[2L]]))), class = c("anova", "data.frame"))
}


a=read.nexus("tree_final.nexus")
codetab=read.table("access_number.txt",sep="\t",header=T,na.strings="NA")
codetab=unique(codetab[,c("code","species")])
codetab$species=gsub("_"," ",codetab$species)

liste2=read.table("yearly_estimates_of_occupancy.txt",sep="\t",header=T)
liste=read.table("linear_mfd_shifts.txt",sep="\t",header=T)


for(i in 1:nrow(liste)){
bidon2=subset(liste2,species==liste$species[i])
bidon2$Annee=bidon2$Annee-1950
se=bidon2$quant_975-bidon2$quant_025
model=lm(mean~Annee,data=bidon2,weights=1/se)
liste$trend_effect[i]=model$coeff["Annee"]
liste$trend_err[i]=summary(model)$coefficient["Annee","Std. Error"]
liste$trend_pval[i]=Anova(model)["Annee",4]
liste$prob50[i]=bidon2[bidon2$Annee==0,"mean"]
}

resf=merge(liste,codetab,by="species")
traits=read.table("traits_full_new_liste_MG.txt",sep="\t",header=T)
resf=merge(resf,traits,by="species")
resf$stat_trend=">0.05"
resf$stat_trend[which(resf$trend_pval<0.05)]="<0.05"
resf$stat_year=">0.05"
resf$stat_year[which(resf$year_pval<0.05)]="<0.05"
resf$Sociality[resf$Sociality=="Solitary+Primitivelyeusocial"]="Primitively eusocial"
resf$Sociality[resf$Sociality=="Socialparasite"]="Social parasite"
resf$Sociality[resf$Sociality=="primitivelyeusocial"]="Primitively eusocial"
resf$Sociality[resf$Sociality=="primitively eusocial"]="Primitively eusocial"
resf$Sociality[resf$Sociality=="Primitivelyeusocial"]="Primitively eusocial"
resf$Sociality[resf$Sociality=="Cleptoparasite"]="Kleptoparasite"
resf[resf=="Unknown"]=NA
resf[resf==""]=NA
resf=resf[grep("Bombus",resf$species,invert=T),]
resf=subset(resf,!is.na(ITD) & !is.na(Lecty) & !is.na(Overwintering.location))


p3d=comparative.data(a,unique(resf[,c("species","code","trend_effect","year_effect","ITD","mfd","Pheno_length","STI","SCI","Sociality",
"Lecty","Overwintering.location")]),"code",vcv=T,na.omit=F)

model=pgls(year_effect~ITD+mfd+Pheno_length+STI+SCI+Sociality+Lecty+Overwintering.location,data=p3d,lambda='ML')
aicref=model$aic
vec=attr(terms(model$formula),"term.labels")
formu=terms(model$formula)
aicvec=c()
for(i in 1:length(vec)){
new.formula=reformulate(attr(drop.terms(formu,dropx=i),"term.labels"),response="year_effect")
model=pgls(new.formula,data=p3d,lambda='ML')
aicvec[i]=model$aic
}
bidon=data.frame(variable=vec,AIC=aicvec)
bidon=bidon[order(bidon$AIC),]
print(bidon)
ind=min(aicvec)
while(ind<aicref & length(vec)>1){
new.formula=reformulate(attr(drop.terms(formu,dropx=which(attr(formu,'term.labels')==bidon[1,1])),"term.labels"),response="year_effect")
model=pgls(new.formula,data=p3d,lambda='ML')
vec=attr(terms(model$formula),"term.labels")
formu=terms(model$formula)
aicref=model$aic
aicvec=c()
for(i in 1:length(vec)){
new.formula=reformulate(attr(drop.terms(formu,dropx=i),"term.labels"),response="year_effect")
model=pgls(new.formula,data=p3d,lambda='ML')
aicvec[i]=model$aic
}
bidon=data.frame(variable=vec,AIC=aicvec)
bidon=bidon[order(bidon$AIC),]
print(bidon)
ind=min(aicvec)
if(ind>=aicref){model=pgls(reformulate(attr(formu,'term.labels'),response="year_effect"),data=p3d,lambda='ML')}
}


summary(model)

anova.f(model)

model2=model
newdat2=newdat

newdat=data.frame(mfd=mean(p3d$data$mfd),ITD=mean(p3d$data$ITD),mfd=mean(p3d$data$mfd),
Pheno_length=mean(p3d$data$Pheno_length),STI=seq(min(p3d$data$STI),max(p3d$data$STI),0.01),SCI=mean(p3d$data$SCI),
Sociality=rep(unique(p3d$data$Sociality),each=length(seq(min(p3d$data$STI),max(p3d$data$STI),0.01))))
newdat$Overwintering.location="Under ground"
newdat1=newdat
newdat1$Overwintering.location="Above ground"
newdat=rbind(newdat,newdat1)
newdat$fit=predict(model2,newdata=newdat,se=T)
b=summaryBy(fit~STI+Sociality,data=newdat,FUN=mean)
b2=summaryBy(STI~Sociality,data=resf,FUN=c(min,max))
b=merge(b,b2,by="Sociality")
b=subset(b,STI>=STI.min & STI<=STI.max)
resf$Sociality=as.factor(resf$Sociality)
b$Sociality=as.factor(b$Sociality)
resf$Socialty=factor(resf$Sociality,levels =c("Kleptoparasite","Primitively eusocial","Solitary","Social parasite","Social",ordered = T))
b$Socialty=factor(b$Sociality,levels =c("Kleptoparasite","Primitively eusocial","Solitary","Social parasite","Social",ordered = T))
pl1=ggplot()+geom_point(data=resf,aes(x=STI,y=year_effect,col=Sociality))+
theme_bw()+ylab("MFD linear trend (days/year) \n")+ggtitle("c")+
theme(panel.grid.minor=element_blank(),plot.title=element_text(size=14,face="bold"),
legend.position="bottom",legend.key.width=unit(0.2,"cm"),legend.text = element_text(size=10),
strip.background = element_blank(),panel.grid=element_blank())+geom_line(data=b,aes(x=STI,y=fit.mean,col=Sociality))+
xlab("Species Temperature Index (°C)")+scale_color_manual(values=brewer.pal(5,"Set2")[c(1,4,5,2,3)])+
guides(col = guide_legend(nrow = 3, byrow = TRUE,title=""))


model=pgls(trend_effect~ITD+mfd+Pheno_length+STI+SCI+Sociality+Lecty+Overwintering.location,data=p3d,lambda='ML')
aicref=model$aic
vec=attr(terms(model$formula),"term.labels")
formu=terms(model$formula)
aicvec=c()
for(i in 1:length(vec)){
new.formula=reformulate(attr(drop.terms(formu,dropx=i),"term.labels"),response="trend_effect")
model=pgls(new.formula,data=p3d,lambda='ML')
aicvec[i]=model$aic
}
bidon=data.frame(variable=vec,AIC=aicvec)
bidon=bidon[order(bidon$AIC),]
print(bidon)
ind=min(aicvec)
while(ind<aicref & length(vec)>1){
new.formula=reformulate(attr(drop.terms(formu,dropx=which(attr(formu,'term.labels')==bidon[1,1])),"term.labels"),response="trend_effect")
model=pgls(new.formula,data=p3d,lambda='ML')
vec=attr(terms(model$formula),"term.labels")
formu=terms(model$formula)
aicref=model$aic
aicvec=c()
for(i in 1:length(vec)){
new.formula=reformulate(attr(drop.terms(formu,dropx=i),"term.labels"),response="trend_effect")
model=pgls(new.formula,data=p3d,lambda='ML')
aicvec[i]=model$aic
}
bidon=data.frame(variable=vec,AIC=aicvec)
bidon=bidon[order(bidon$AIC),]
print(bidon)
ind=min(aicvec)
if(ind>=aicref){model=pgls(reformulate(attr(formu,'term.labels'),response="trend_effect"),data=p3d,lambda='ML')}
}

summary(model)

newdat=data.frame(mfd=mean(p3d$data$mfd),ITD=seq(min(p3d$data$ITD),max(p3d$data$ITD),0.01),mfd=mean(p3d$data$mfd),
Pheno_length=mean(p3d$data$Pheno_length),STI=mean(p3d$data$STI),SCI=mean(p3d$data$SCI),
Sociality=rep(unique(p3d$data$Sociality),each=length(seq(min(p3d$data$ITD),max(p3d$data$ITD),0.01))))
newdat$fit=predict(model,newdata=newdat,se=T)
b2=summaryBy(ITD~Sociality,data=resf,FUN=c(min,max))
newdat=merge(newdat,b2,by="Sociality")
newdat=subset(newdat,ITD>=ITD.min & ITD<=ITD.max)
resf$Socialty=factor(resf$Sociality,c("Kleptoparasite","Primitively eusocial","Solitary","Social parasite","Social"))
newdat$Socialty=factor(newdat$Sociality,c("Kleptoparasite","Primitively eusocial","Solitary","Social parasite","Social"))
pl2=ggplot(data=resf,aes(x=ITD,y=trend_effect,col=Sociality))+geom_point()+
theme_bw()+ylab("Occupancy linear trend")+ggtitle("b")+
theme(panel.grid.minor=element_blank(),plot.title=element_text(size=14,face="bold"),legend.position="none",
strip.background = element_blank(),panel.grid=element_blank())+
geom_line(data=newdat,aes(x=ITD,y=fit,col=Sociality))+xlab("Intertegular distance (mm)")+
scale_color_manual(values=brewer.pal(5,"Set2")[c(1,4,5,2,3)])


p3d=comparative.data(a,unique(resfb[,c("species","code","trend_effect","year_effect","ITD","stat_year",
"stat_trend","Family","Sociality")]),"code",vcv=T,na.omit=F)

library(ggstance)
library(ggnewscale)

df=data.frame(size=log(p3d$data[,"ITD"]))
rownames(df)=rownames(p3d$data)
nodelabels(cex=0.6)
p=ggtree(p3d$phy)+
geom_cladelabel(node=388,barsize=1,color='black',label="Colletidae",offset=2,
offset.text=5,angle=-90,hjust=0.4,fontsize=2.8)+
geom_cladelabel(node=207,barsize=1,color='black',label="Apidae",offset.text=5,offset=2,
angle=-90,hjust=0.5,fontsize=2.8)+
geom_cladelabel(node=294,barsize=1,color='black',label="Andrenidae",offset.text=5,offset=2,
angle=-90,hjust=0.5,fontsize=2.8)+
geom_cladelabel(node=257,barsize=1,color='black',label="Megachilidae",offset.text=5,offset=2,
angle=-90,hjust=0.5,fontsize=2.8)+
geom_cladelabel(node=400,barsize=1,color='black',label="Melittidae",offset.text=5,offset=2,
angle=-90,hjust=0.5,fontsize=2.8)+
geom_cladelabel(node=342,barsize=1,color='black',label="Halictidae",offset.text=5,offset=2,
angle=-90,hjust=0.5,fontsize=2.8)

df=data.frame(id=rownames(p3d$data),size=log(p3d$data[,"ITD"]),stat_year=p3d$data$stat_year)
###########
p1 <- p %<+% df + geom_tippoint(aes(color=size),size=1)+scale_colour_viridis_c(option="C",na.value="black")+
guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5,breaks=c(-1,0,1)))
df2=data.frame(id=rownames(p3d$data),MFDs=p3d$data[,"year_effect"])
p2 <- facet_plot(p1, panel="MFD linear trends", data=df2, geom=geom_barh, aes(x=MFDs,fill=stat_year),
col=NA,stat='identity') +
scale_fill_manual(values=c("black","darkgrey"))+
theme(legend.position="none",axis.text.x=element_text(size=10),axis.line.x=element_line())+
xlab("days/year")


df3=data.frame(id=rownames(p3d$data),abd=p3d$data[,"trend_effect"],stat_trend=p3d$data$stat_trend)
p2=p2+new_scale_fill()+guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
p3 <- facet_plot(p2, panel="Occupancy linear trends", data=df3, geom=geom_barh, aes(x=abd,fill=stat_trend),col=NA,stat='identity') + 
scale_fill_manual(values=c("black","darkgrey"))+
theme(legend.position="bottom",plot.title=element_text(size=14,face="bold"),axis.line.x=element_line(),
axis.ticks.x=element_line(),axis.title.x=element_text())+
guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5,breaks=c(-1,0,1)),
fill = guide_legend(title.position="top", title.hjust = 0.5))+ggtitle("a")




gridExtra::grid.arrange(p3,pl2,pl1,layout_matrix=rbind(c(1, 2),c(1, 3)),widths=c(1.5,1),heights=c(1,1.3))

pdf("figure_size_phylo.pdf",width=10, height=9)
gridExtra::grid.arrange(p3,pl2,pl1,layout_matrix=rbind(c(1, 2),c(1, 3)),widths=c(1.5,1),heights=c(1,1.3))
dev.off();


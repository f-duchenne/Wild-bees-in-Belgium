################################GRAPHS
library(reshape2)
library(data.table)
library(sparta)

setwd(dir="C:/Users/Francois/Documents/pheno abeilles belges/model_isaac2")
lili=list.files()
lili=lili[grep("rdata",lili)]
lili=lili[grep(".png",lili,invert=T)]
for(i in 1:length(lili)){
load(lili[i])
bibi=summary(out)
bibi$Annee=1950:2016
bibi$species=gsub(".rdata","",gsub("_"," ",lili[i]))
model=gam(mean~s(Annee),data=bibi)
bibi$mean_smooth=predict(model)
bidon=as.data.frame(out$BUGSoutput$summary)
bibi$rhat=bidon[grep("psi.fs",rownames(bidon)),"Rhat"]
bibi$rhat.prop=length(bibi$rhat[bibi$rhat<=1.1])/length(bibi$rhat)
if(i==1){res=bibi}else{res=rbind(res,bibi)}
print(paste(i,length(bibi$rhat[bibi$rhat<=1.1])/length(bibi$rhat),sep=" - "))
}

length(unique(subset(res,rhat.prop>=0.6)$species))

fwrite(res,"trends.txt",sep="\t",row.names=F)
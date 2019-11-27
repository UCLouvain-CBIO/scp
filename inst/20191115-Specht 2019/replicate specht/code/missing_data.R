# Calculate missing data: 

# Type I: missing due to low S/N

frac<-c()
for( X in colnames(ev.matrix.sc.f.n)){
  
  ev.t<-ev.melt[ev.melt$id==X, ]
  
  q.t<-aggregate(quantitation ~ protein, data=ev.t, FUN=median.na)
  
  q.t$quantitation[q.t$quantitation==0]<-NA
  q.t$quantitation[q.t$quantitation==Inf]<-NA
  q.t$quantitation[q.t$quantitation==-Inf]<-NA
  
  frac.t<-length(which(is.na(q.t$quantitation))) / length(q.t$quantitation)
  
  frac<-c(frac, frac.t)
}


# Type II: missing due to disparate observations between sets (ie. identified in some sets, not others)

evg1<-ev.matrix.lowbin0

length(which(evg1==0))/ dim(evg1)[1] / dim(evg1)[2]

length(which(is.na(evg1))) / dim(evg1)[1] / dim(evg1)[2]

evg<-ev.matrix.lowbin0[rownames(ev.matrix.lowbin0)%in%unique(ev.melt$sequence[ev.melt$protein%in%rownames(ev.matrix.sc.f.n)]), colnames(ev.matrix.lowbin0)%in%colnames(ev.matrix.sc.f.n) ]
evg<-collapse_to_protein( filt.mat.rc(evg, na.row, na.col), ev.melt, T)

length(which(evg==0))/ dim(evg)[1] / dim(evg)[2]

length(which(is.na(evg))) / dim(evg)[1] / dim(evg)[2]


frac2<-na.count(t(evg)) / dim(evg)[1]

dvm<-melt(data.frame(frac, frac2)); head(dvm)

ggplot(data = dvm, aes(y=variable, x= value)) + 
  geom_density_ridges(aes(fill=variable)) + 
  scale_fill_manual(values=c("#F3B900", "#313A87")) + 
  theme_pubr() + 
  rremove("legend") + 
  rremove("y.ticks") + 
  rremove("y.text") + 
  ylab("Density") + 
  xlab("Fraction missing data") + 
  font("xylab", size=30) + 
  font("x.text", size=20) 


summary(frac); summary(frac2)

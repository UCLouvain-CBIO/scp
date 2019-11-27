####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Dimensional reduction

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# PC's to display:
PCx<-"PC1"
PCy<-"PC2"

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Data to use: 
mat.sc.imp<-cr_norm_log(matrix.sc.batch)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Calculate missing data per protein, use inverse as weight for PCA: 
prots<-rownames(mat.sc.imp)
weight.v<-c()

for(X in prots){
  
  peps<-ev.melt.pep$sequence[ev.melt.pep$protein==X]
  
  peps_ind<-which(rownames(t3)%in%peps)
  
  mat.t<-t3[peps_ind, ]
  
  weight.t<-NA
  
  if(length(peps_ind)>1) {
  
  weight.t <- 1 - ( sum(na.count(mat.t)) / ( dim(mat.t)[1] * dim(mat.t)[2] ) )
  
  }
  
  if(length(peps_ind)<=1) {
    
    weight.t <- 1 - ( length(which(is.na(mat.t)))/length(mat.t) )
    
  }
  
  
  weight.v<-c(weight.v, weight.t)
  
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Calculate the weighted data matrix: 

X.m <- mat.sc.imp

pca.imp<- t(X.m) %*% diag(weight.v)  %*% X.m / sum(weight.v)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Perform PCA
sc.pca<-prcomp(pca.imp, scale. = T, center = T)

# Percent of variance explained by each principle component
pca_var <- (sc.pca$sdev)^2
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data 
pca.melt <- melt(sc.pca$rotation); colnames(pca.melt)<-c("id","pc","value")
add.cols<-colnames(ev.melt)[4:7]
pca.melt[,add.cols]<-NA

for(X in unique(pca.melt$id)){
  
  pca.melt[pca.melt$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]
  
}


# Re map ... 
pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)

pca.display[,add.cols]<-NA

for(X in unique(pca.display$id)){
  
  pca.display[pca.display$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]
  
}


# Map to TMT channel
pca.display$channel<-c2q[c2q$celltype%in%pca.display$id, "channel"]

# Load in bulk proteomic data to color the PCA plot: 
m.m<-mat.sc.imp[rownames(mat.sc.imp)%in%mono_genes30, ]; dim(m.m)
mac.m<-mat.sc.imp[rownames(mat.sc.imp)%in%mac_genes30, ]; dim(mac.m)

# Color by median value across replicates
int.m<-rowMedians((t(m.m)), na.rm=T)
int.mac<-rowMedians((t(mac.m)), na.rm=T)
pca.display$mono<-int.m
pca.display$mac<-int.mac


# Display celltype
pg1<-ggscatter(pca.display, x =PCx, y = PCy , color="celltype", size = 5, alpha=0.5) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) + 
  font("ylab",size=30) + 
  font("xlab",size=30) + 
  font("xy.text", size=20) + 
  rremove("legend") + 
  scale_color_manual(values = my_colors[2:3]) + 
  annotate("text", x=-0.03, y=0.21,label="Macrophage", color=my_colors[2], size=10)  + 
  annotate("text", x=0.04, y=0.21, label="Monocyte", color=my_colors[3], size=10) + 
  annotate("text", x=0.05-0.01, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  annotate("text", x=0.062-0.014, y=-0.13, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8)

# Adjust color scale: 
varx<-int.m
qs.varx<-quantile(rescale(varx), probs=c(0,0.3,0.5,0.7,1))

# Display celltype
pg2<-ggscatter(pca.display, x =PCx, y = PCy , color="mono", size = 3, alpha=0.7) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) + 
  font("ylab",size=30) + 
  font("xlab",size=30) + 
  font("xy.text", size=20) + 
  scale_color_gradientn(colors=c("purple","yellow2"), values=qs.varx) + 
  #annotate("text", x=-0.04, y=0.21,label="Macrophage", color=my_colors[2], size=10)  + 
 # annotate("text", x=0.04, y=0.21, label="Monocyte", color=my_colors[3], size=10) + 
  #annotate("text", x=0.05-0.01, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  #annotate("text", x=0.062-0.014, y=-0.13, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) + 
  rremove("xylab")+
  rremove("xy.text") + 
  rremove("ticks") + 
  rremove("legend") + 
  ggtitle("Monocyte genes")

# Adjust color scale:
varx<-int.mac
qs.varx<-quantile(rescale(varx), probs=c(0,0.3,0.5,0.7,1))

# Display celltype
pg3<-ggscatter(pca.display, x =PCx, y = PCy , color="mac", size = 3, alpha=0.7) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) + 
  font("ylab",size=30) + 
  font("xlab",size=30) + 
  font("xy.text", size=20) + 
  scale_color_gradientn(colors=c("purple","yellow2"), values=qs.varx) + 
  #annotate("text", x=-0.04, y=0.21,label="Macrophage", color=my_colors[2], size=10)  + 
  #annotate("text", x=0.04, y=0.21, label="Monocyte", color=my_colors[3], size=10) + 
  #annotate("text", x=0.05-0.01, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  #annotate("text", x=0.062-0.014, y=-0.13, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) + 
  rremove("xylab")+
  rremove("xy.text") + 
  rremove("ticks") + 
  ylab("Macrophage-genes") + 
  rremove("legend")+ 
  ggtitle("Macrophage genes")

ggscatter(pca.display, x =PCx, y = PCy , color="mac", size = 3, alpha=0.7) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) + 
  font("ylab",size=30) + 
  font("xlab",size=30) + 
  font("xy.text", size=20) + 
  scale_color_gradientn(colors=c("purple","yellow2"), values=qs.varx) + 
  #annotate("text", x=-0.04, y=0.21,label="Macrophage", color=my_colors[2], size=10)  + 
  #annotate("text", x=0.04, y=0.21, label="Monocyte", color=my_colors[3], size=10) + 
  #annotate("text", x=0.05-0.01, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  #annotate("text", x=0.062-0.014, y=-0.13, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) + 
  rremove("xylab")+
  rremove("xy.text") + 
  rremove("ticks") + 
  ylab("Macrophage-genes") + 
  #rremove("legend")+ 
  ggtitle("Macrophage genes")


pg1 + (pg2 + pg3 + plot_layout(ncol=1, nrow=2, heights=c(1,1))) + plot_layout(ncol=2, nrow=1, widths=c(2,1))



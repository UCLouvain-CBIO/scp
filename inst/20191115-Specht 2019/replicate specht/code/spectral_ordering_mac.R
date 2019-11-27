####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Spectral clustering - macrophage-like cells only

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Get the macrophage-like data only: 
m0.ind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_m0"]
mat.sc.imp.m0<-matrix.sc.batch[, colnames(matrix.sc.batch)%in%m0.ind]
mat.sc.imp.m0<-cr_norm_log(mat.sc.imp.m0)

# Calculate the Laplacian of the correlation matrix (+1 to all values, so that no values are negative)
W <-  1+ cor(mat.sc.imp.m0)
D <- diag( rowSums(W))
L <- D - W 

# Get the eigenvalues and eigenvectors of the Laplacian
eigs<-eigen(L)
vec<-eigs$vectors
val<-eigs$values
plot(val)

# Order the eigenvectors by eigenvalue
vec<-vec[,order(val, decreasing=F)]
vec.m<-vec[,1:2]

# Reformat eigenvector into data frame (taking second eigen vector, first should have eigenvalue 0 and all same values in eigenvector)
vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")

# Plot the eigenvector
pl1<-ggline(vdf, x="num", y="eigen", size=0.001, color="gray60") + theme_pubr() +
  rremove("xy.text") + rremove("xylab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0))


mac_vec<-vec.m[,2]

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Reorder the data matrix
mat.c<-mat.sc.imp.m0[, order(vec.m[,2]) ]

# Record rowMedians of the outer 40 cells on either side of matrix
rm1<-rowMedians(mat.c[,1:40], na.rm=T)
rm2<-rowMedians(mat.c[, (ncol(mat.c)-40):ncol(mat.c) ], na.rm=T)

# Create a score to reorder rows for visualization
mfc<-(rm1-rm2)
names(mfc)<-rownames(mat.sc.imp.m0)
mfc2<-mfc; names(mfc2)<-rownames(matrix.sc.batch); write.csv(mfc2,"fc_2019.csv")
names(mfc)<-1:length(mfc)

mfc.g<-mfc[mfc>0]
mfc.l<-mfc[mfc<0]

mfc.g<-mfc.g[order(mfc.g, decreasing = F)]
mfc.l<-mfc.l[order(mfc.l, decreasing = T)]

mfc.reorder<-c(mfc.g, mfc.l)

names(rm1)<-1:length(rm1)

rm1g<-rm1[names(mfc.g)]
rm1l<-rm1[names(mfc.l)]

rm1g<-rm1g[order(rm1g, decreasing = F)]
rm1l<-rm1l[order(rm1l, decreasing = T)]

rm1.reorder<-c(rm1g,rm1l)

rows.keep<-as.numeric(names(mfc.reorder[abs(mfc.reorder)>quantile(abs(mfc.reorder),probs=0.75)]))

rm.keep<-rm1.reorder[as.numeric(names(rm1.reorder))%in%rows.keep]

# Create a new matrix with re-ordered rows and columns
mat.p<-mat.c[as.numeric(names(rm.keep)),]

# Normalize
mat.p<-cr_norm_log(mat.p)

# mat.p[mat.p>2] <- 2
# mat.p[mat.p< (-2)] <- (-2)

# Reverse order
mat.p<-mat.p[, ncol(mat.p):1]
mac.spectrum.ps<-rownames(mat.p)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Create the color scale

varx<-melt(mat.p)$value
qs.varx<-quantile(rescale(varx), probs=c(0,0.05,0.5,0.95,1))
qs.varx<-quantile(rescale(varx), probs=c(0,0.10,0.5,0.9,1))

my_col2<-c(rgb(0,0,1,0.5),rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),rgb(1,0,0,0.5))
my_col2<-c("blue","blue","white","red","red")
my_col2<-c("blue",rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),"red")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Plot! 

pl2<-ggplot(melt(mat.p), aes(x=Var2, y=Var1, fill=value))+ 
  geom_tile() + 
  scale_fill_gradientn(colors=my_col2, values=qs.varx) + 
  #scale_fill_gradient2(low="blue", mid="white",high="red", midpoint=0) +
  rremove("xy.text") + 
  rremove("ticks") + 
  ylab("Proteins") + 
  xlab("Cells") + 
  font("xylab", size=20)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

dfx<-data.frame(var1<-ev.melt.uniqueID$celltype[match(colnames(mat.c), ev.melt.uniqueID$id)]); colnames(dfx)<-c("var")
dfx$x<-nrow(dfx):1
dfx$y<-1

vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")

vdf$celltype<-dfx$var[order(dfx$x, decreasing=F)]

pl1<-ggbarplot(vdf, x="num", y="eigen", fill="celltype", color="celltype") + theme_pubr() +
  rremove("xy.text") + rremove("xylab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0)) + 
  scale_fill_manual(values=my_colors[c(2,3)])  + 
  scale_color_manual(values=my_colors[c(2,3)]) + 
  rremove("legend") + rremove("axis")

plot_cells<-ggplot(dfx, aes(x=x,y=y, fill=var))+geom_tile() +
  theme_pubr() + 
  rremove("legend") + 
  rremove("xylab") + 
  rremove("axis") + 
  rremove("ticks") + 
  rremove("xy.text") + 
  scale_fill_manual(values=my_colors[c(2,3)]) + 
  theme(plot.margin = margin(0, 1, 0, 5, "cm"))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Extract markers of M1 and M2 polarization and look at their quantitation across spectrum

M2<-read.csv("M2overM1_Geneset.csv")
head(M2)

M1<-read.csv("M1overM2_Geneset.csv")
head(M1)

M2.m<-( mat.c[which(rownames(mat.c)%in%as.character(M2$M2overM12)), ] )
M1.m<-( mat.c[which(rownames(mat.c)%in%as.character(M1$M1overM22)), ] )
m2.v<-rowMedians(t(M2.m), na.rm = T)
m1.v<-rowMedians(t(M1.m), na.rm = T)


m2d<-data.frame(m2.v, 1:length(m2.v))
colnames(m2d)<-c("value", "index")

m1d<-data.frame(m1.v, 1:length(m1.v))
colnames(m1d)<-c("value", "index")

colnames(M2.m)<-1:ncol(M2.m)
melt(M2.m)



rhom2<-cor(m2d$value, m2d$index)
rhom1<-cor(m1d$value, m1d$index)

sp1<-ggscatter(m1d, x="index", y="value", fill="black") + theme_pubr() +
  rremove("xy.text") + rremove("xlab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0)) + 
  rremove("legend") + rremove("axis") + geom_smooth(method = "lm") + geom_text(x=200, y=0.25, label=paste0("r = ",round(rhom1,2)), size=7) + 
  ylab("M1-markers") + font("ylab", size=20)


sp2<-ggscatter(m2d, x="index", y="value", fill="black") + theme_pubr() +
  rremove("xy.text") + rremove("xlab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0)) + 
  rremove("legend") + rremove("axis") + geom_smooth(method = "lm") +  geom_text(x=50, y=0.37, label=paste0("r = ",round(rhom2,2)), size=7) + 
  ylab("M2-markers") + font("ylab", size=20) 

# One point per cell: 
pl1 + pl2 + sp1 + sp2 + plot_layout(ncol=1, nrow=4, heights = c(1,5,2,2))


#####################

colnames(M2.m)<-1:ncol(M2.m)
m2d<-melt(M2.m)

colnames(M1.m)<-1:ncol(M1.m)
m1d<-melt(M1.m)

head(m1d)


rhom2<-cor(m2d$value, m2d$Var2)
rhom1<-cor(m1d$value, m1d$Var2)

sp1<-ggscatter(m1d, x="Var2", y="value", fill="black") + theme_pubr() +
  rremove("xy.text") + rremove("xlab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0)) + 
  rremove("legend") + rremove("axis") + geom_smooth(method = "lm") + geom_text(x=200, y=1.15, label=paste0("r = ",round(rhom1,2)), size=7) + 
  ylab("M1-markers") + font("ylab", size=20) + ylim(c(-1.2,1.2))


sp2<-ggscatter(m2d, x="Var2", y="value", fill="black") + theme_pubr() +
  rremove("xy.text") + rremove("xlab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0)) + 
  rremove("legend") + rremove("axis") + geom_smooth(method = "lm") +  geom_text(x=50, y=0.8, label=paste0("r = ",round(rhom2,2)), size=7) + 
  ylab("M2-markers") + font("ylab", size=20)  + ylim(c(-1,1))

# Multiple data points per cell
pl1 + pl2 + sp1 + sp2 + plot_layout(ncol=1, nrow=4, heights = c(1,5,2,2))



###################

se<-function(x){ sd(x, na.rm=T)/sqrt(length(x))}


qs<-quantile(m1d$Var2, probs=seq(0,1,0.1))[-1]

m1d$quantile<-NA

for(i in round(qs[order(qs, decreasing = T)],0)){
  
  print(i)
  m1d$quantile[m1d$Var2%in%1:i]<-i
  
}

qs<-quantile(m2d$Var2, probs=seq(0,1,0.1))[-1]

m2d$quantile<-NA

for(i in round(qs[order(qs, decreasing = T)],0)){
  
  print(i)
  m2d$quantile[m2d$Var2%in%1:i]<-i
  
}

####
sdf<-data.frame( aggregate(value~quantile, data=m1d, FUN = median), aggregate(value~quantile, data=m1d, FUN = se) )

colnames(sdf)<-c("quantile","median","quantile1","se")

sdf$quantile1<-sdf$quantile1-13


rhom1<-cor(sdf$quantile1, sdf$median)

sp1<-ggscatter(sdf, x="quantile1",y="median") + geom_smooth(method = "lm") + 
  geom_errorbar(data = sdf, aes(x = quantile1, y = median, ymin = median - se, ymax = median + se)) + 
  rremove("xy.text") + 
  rremove("xlab") + 
  rremove("ticks")+
  xlim(c(0, ncol(mat.p)))+
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(breaks=c(-0.1,0,0.15), labels=c("-100%", "0%", "100%") ) + 
  rremove("legend")  +  #geom_text(x=200, y=0.1, label=paste0("r = ",round(rhom1,2)), size=7) + 
  ylab("M1-genes") + font("ylab", size=20) #+ geom_hline(yintercept=0, color="black", size=1)


sdf<-data.frame( aggregate(value~quantile, data=m2d, FUN = median), aggregate(value~quantile, data=m1d, FUN = se) )


colnames(sdf)<-c("quantile","median","quantile1","se")
rhom2<-cor(sdf$quantile1, sdf$median)

sdf$quantile1<-sdf$quantile1-13

sp2<-ggscatter(sdf, x="quantile1",y="median") + geom_smooth(method = "lm") + 
  geom_errorbar(data = sdf, aes(x = quantile1, y = median, ymin = median - se, ymax = median + se)) + 
  rremove("xy.text") + 
  rremove("xlab") + 
  rremove("ticks")+
  xlim(c(0, ncol(mat.p)))+
  scale_x_continuous(expand = c(0, 0)) + 
  #scale_y_continuous(breaks=c(-0.1,0,0.15), labels=c("-100%", "0%", "100%") ) + 
  rremove("legend")  +  #geom_text(x=50, y=0.1, label=paste0("r = ",round(rhom2,2)), size=7) + 
  ylab("M2-genes") + font("ylab", size=20) #+ geom_hline(yintercept=0, color="black", size=1)


#####

# As a summary point and standard error over bins of 26 cells:
pl1 + pl2 + sp1 + sp2 + plot_layout(ncol=1, nrow=4, heights = c(1,5,1.5,1.5))


#M1$M1overM21[M1$M1overM22%in% rownames(M1.m)]
#M2$M2overM11[M2$M2overM12%in% rownames(M2.m)]

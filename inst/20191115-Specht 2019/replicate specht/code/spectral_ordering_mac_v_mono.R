####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Spectral clustering

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Calculate the Laplacian of the correlation matrix (+1 to all values, so that no values are negative)
W <- 1 + cor(matrix.sc.batch)
D <- diag(rowSums(W))
L <- D - W 

# Get the eigenvalues and eigenvectors of the Laplacian
eigs<-eigen(L)
vec<-eigs$vectors
val<-eigs$values

# Sanity check
plot(val)

# Order the eigenvectors by eigenvalue
vec<-vec[,order(val, decreasing=F)]
vec.m<-vec[,1:2]

# Reformat eigenvector into data frame (taking second eigen vector, first should have eigenvalue 0 and all same values in eigenvector)
vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")

# Plot the eigenvector
pl1<-ggline(vdf, x="num", y="eigen", size=0.001, color="gray60") +
  theme_pubr() +
  rremove("xy.text") +
  rremove("xylab") +
  rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0))

# Record re-ordered eigenvector
mac_mono_vec<-vec.m[,2]

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

 
# Reorder the data matrix
mat.c<-mat.sc.imp[, order(vec.m[,2]) ]

# Record rowMedians of the outer 40 cells on either side of matrix
rm1<-rowMedians(mat.c[,1:40], na.rm=T)
rm2<-rowMedians(mat.c[, (ncol(mat.c)-40):ncol(mat.c) ], na.rm=T)

# Create a score to reorder rows for visualization
mfc<-(rm1-rm2)
names(mfc)<-rownames(mat.sc.imp)
mfc2<-mfc; names(mfc2)<-rownames(matrix.sc.batch)
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

rows.keep<-as.numeric(names(mfc.reorder[abs(mfc.reorder)>quantile(abs(mfc.reorder),probs=0.8)]))

rm.keep<-rm1.reorder[as.numeric(names(rm1.reorder))%in%rows.keep]


# Create a new matrix with re-ordered rows and columns
mat.p<-mat.c[as.numeric(names(rm.keep)),]

# Normalize
mat.p<-cr_norm_log(mat.p)

# mat.p[mat.p>2] <- 2
# mat.p[mat.p< (-2)] <- (-2)

# Reverse order to have monocytes on the left
mat.p<-mat.p[, ncol(mat.p):1]

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

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Plot! 

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
  scale_fill_manual(values=my_colors[c(2,3)]) 

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Final figure

pl1 + pl2 + plot_layout(ncol=1, nrow=2, heights = c(2,7))

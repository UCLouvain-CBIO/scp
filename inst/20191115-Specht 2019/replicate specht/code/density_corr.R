# Looking at heterogeniety within clusters of macrophage-like and monocyte cells

# Grab cells
m0.ind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_m0"]
u.ind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_u"]

# Create matrices
mat.u<-matrix.sc.batch[, colnames(matrix.sc.batch)%in%u.ind]
mat.m0<-matrix.sc.batch[, colnames(matrix.sc.batch)%in%m0.ind]

# Calculate all pairwise correlations
m0cor<-c(cor(mat.m0, use = "pairwise.complete.obs"))
ucor<-c(cor(mat.u, use = "pairwise.complete.obs"))

# Combine into data frame: 
ucor<-c(ucor, rep(NA, 57672))
dfxxx<-data.frame(ucor, m0cor)
dxm<-melt(dfxxx)

# Plot! 
ggplot(data = dxm, aes(y=variable, x= value)) + 
  geom_density_ridges(aes(fill=variable)) + 
  xlim(c(-0.2,0.8)) + 
  scale_fill_manual(values=my_colors[c(3,2)]) + 
  theme_pubr() + 
  rremove("legend") + 
  rremove("y.ticks") + 
  rremove("y.text") + 
  ylab("Density") + 
  xlab("Correlation") + 
  font("xylab", size=30) + 
  font("x.text", size=20) 


ggplot(data = dxm, aes(y=variable, x= value)) + 
  geom_density_ridges(aes(fill=variable)) + 
  xlim(c(-0.2,0.8)) + 
  scale_fill_manual(values=my_colors[c(3,2)]) + 
  theme_pubr() + 
  rremove("legend") + 
  rremove("x.ticks") + 
  rremove("x.text") + 
  ylab("Density") + 
  xlab("Correlation") + 
  font("xylab", size=30) + 
  font("y.text", size=20) + 
  coord_flip() 

# Are the means of the distributions the same? 
t.test(m0cor[-which(m0cor==1)], ucor[-which(ucor==1)],  alternative = c("two.sided"))$p.value

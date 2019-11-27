

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Calculate coefficient of variation for proteins with at least 6 peptides
cv1<-prot.cv(evc.matrix.sc, ev.melt, 6)

# Function for calculating the 30th percentile
q30<-function(x){quantile(x, probs=0.30, na.rm = T)}

# Grab indicies corresponding to each cell type and control wells
sc0ind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_0"]
m0ind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_m0"]
uind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_u"]

# Heuristic values for filtering... adjust as needed: 
rRI_thres<- -2.5
cv_thres<- 0.43
contam_thres<- -1.3

# Plot limits for consistency between the dotplots and density distributions, redefine as needed: 
ylim.t = c(0.27, 0.63)
xlim.t = c(-3, -1.2) 

# Calculate statistics on the relative reporter ion (relative to carrier channel) and coefficient of variation of relative quantitation (normalized
# to reference) of each protein for each single cell or control well. 
df<-data.frame(
  cv<-apply(cv1, 2, FUN=median, na.rm=T),
  rRI<-apply(log10(evc.matrix.sc.lb), 2, FUN=q30),
  rRIm<-apply(log10(evc.matrix.sc.lb), 2, FUN=median.na) )

# Name the columns of the data frame
colnames(df)<-c("median_cv","q30_rRI","median_rRI")

# Assign values <-3 to one bin
df$median_rRI[df$median_rRI==-Inf]<- -3
df$median_rRI[df$median_rRI==Inf]<- -3
df$q30_rRI[df$q30_rRI==-Inf]<- -3
df$q30_rRI[df$q30_rRI==Inf]<- -3

# Assign cell types to each well
df$celltype<-NA
df$celltype[colnames(cv1)%in%sc0ind]<-"  control-wells  "
df$celltype[colnames(cv1)%in%m0ind]<-"  macrophages  "
df$celltype[colnames(cv1)%in%uind]<-"  monocytes  "

# Assign batches
df$batch<-1
df$batch<-ev.melt.uniqueID$digest[ match(colnames(cv1), ev.melt.uniqueID$id )]
q
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Get the cells that pass the assigned heutristic filters (from above)

pass_filter<-c()

for(i in 1:nrow(df)){
  # TODO there is a typo here there should be '(df$median_rRI[i] < contam_thres)'
  if( (df$median_cv[i] < cv_thres) & (df$q30_rRI[i] > rRI_thres) & (df$median_rRI[i] < contam_thres)){
    
    pass_filter<-c(pass_filter, i)
    
  }
  
}


wells_pass_filter<-colnames(cv1)[pass_filter]


# Use the reference-normalized matrix (not "evc", which is the column normalized matrix, only used for scoring): 
ev.matrix.sc.f<-ev.matrix.sc[,wells_pass_filter]


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Percent of m0 / u cells that fall in 4th quadrant
pct.passed<-ncol(ev.matrix.sc.f) / length(c(m0ind,uind))
pct.passed0<-length(which(colnames(ev.matrix.sc.f)%in%sc0ind)) / length(sc0ind)
pct.passedm0<-length(which(colnames(ev.matrix.sc.f)%in%m0ind)) / length(m0ind)
pct.passedu<-length(which(colnames(ev.matrix.sc.f)%in%uind)) / length(uind)



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Plotting! 


my_colors<-c("black","#048ABF","#FF5733")


p3<-ggscatter(data=df, x="median_rRI", y="median_cv", color="celltype",alpha=0.5,
              margin.params = list(fill = "celltype", color = "black", size = 0.2, alpha=0.3)) +
  xlim(xlim.t)+
  ylim(ylim.t)+
  ylab("Coefficient of variation\n") +
  xlab(expression('Relative RI, log'[10])) + 
  geom_hline(yintercept=cv_thres, linetype="dashed", color = "gray",size=1.2) + 
  geom_vline(xintercept=rRI_thres, linetype="dashed", color = "gray",size=1.2) + 
  rremove("legend") +
  font("xlab", size=20) + 
  font("ylab", size=20) + 
  font("xy.text", size=15) + 
  scale_color_manual(values=my_colors) +
  annotate("text",x=-1.25, y=0.28, label=paste0(round(pct.passed0,2)*100, "%"), fontface =2,size=7, color=my_colors[1]) + 
  annotate("text",x=-1.25, y=0.336, label=paste0(round(pct.passedm0,2)*100, "%"), fontface =2,size=7, color=my_colors[2]) + 
  annotate("text",x=-1.25, y=0.4, label=paste0(round(pct.passedu,2)*100, "%"), fontface =2,size=7, color=my_colors[3]) +
    scale_x_continuous(breaks=c(-3,-2.5,-2,-1.5), labels = c(paste0("\u2264"," -3"), "-2.5", "-2","-1.5")) +
  ylim(ylim.t) +
  xlim(xlim.t) + 
coord_flip() 


p5<-ggplot(data=df, aes(x = median_rRI)) + 
  geom_density(aes(y=..density..,fill=celltype), alpha=0.7, bw=density(df$median_rRI)$bw/2, position="dodge") +
  theme_pubr() +
  rremove("ylab") +
  rremove("xlab") +
  #rremove("x.axis") +
  rremove("x.ticks") +
  rremove("x.text") +
  rremove("y.axis") +
  rremove("y.ticks") +
  rremove("y.text") +
  rremove("legend") +
  scale_fill_manual(values=my_colors) +
  font("legend.text",size=15) +
  font("title", size=20) + 
  #scale_fill_manual(values=c(rgb(0,0,0,0.9),rgb(0,0,1,0.7),rgb(1,0,0,0.7) )) + 
  scale_fill_manual(values=my_colors)+
  theme(plot.title = element_text(hjust = 0.5, vjust=1)) + 
  xlim(xlim.t) + 
coord_flip() + 
  rremove("x.axis")



p6<-ggplot(data=df, aes(x = median_cv)) + 
  geom_density(aes(y=..density..,fill=celltype), alpha=0.7, bw=density(df$median_cv)$bw, position="dodge") +
  theme_pubr() +
  rremove("ylab") +
  rremove("xlab") +
  rremove("x.axis") +
  rremove("x.ticks") +
  rremove("x.text") +
  #rremove("y.axis") +
  rremove("y.ticks") +
  rremove("y.text") +
  rremove("legend") +
  scale_fill_manual(values=my_colors) +
  font("legend.text",size=15) +
  font("title", size=20) +
  ggtitle("Scoring sample preparation success") +
  scale_fill_manual(values=my_colors)+
  #scale_fill_manual(values=c(rgb(0,0,0,0.9),rgb(0,0,1,0.7),rgb(1,0,0,0.7) )) +
  xlim(ylim.t) + 
  ylim(c(0,25)) + 
  rremove("y.axis")


pL<-get_legend(p5 +  theme(legend.position="right") +rremove("legend.title"))

p6 + pL + p3 + p5 + plot_layout(ncol=2, nrow=2, widths=c(5,1.9), heights=c(1.5,5))





####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Which genes up or down regulated in mono vs. mac in bulk proteomic data? 

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Take the 100 cell TMT11-plex data, put on log2 scale and renormalize: 

p100<-collapse_to_protein(log2(ev.matrix.100.n), ev.melt, T)
p100<-cr_norm_log(p100)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Calculate fold-change between cell types and significance of that fold-change using the distributions of the protein values in 
# each of those cell types: 

# Initialize variables to store protein, p-value, fold-change:
pc<-c()
pval<-c()
fc<-c()

# Iterate over every protein: 
for(X in rownames(p100)){
  
  # Separate data according to cell type
  m0.id<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="m0_100"]
  u.id<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="u_100"]
  
  # As long as there is data for both cell types: 
  if( !all(is.na(p100[rownames(p100)==X, colnames(p100)%in%m0.id])) & !all(is.na(p100[rownames(p100)==X, colnames(p100)%in%u.id])) ){
    
    fc.t<-NA
    pval.t<-NA
    
    pval.t<-t.test( p100[rownames(p100)==X, colnames(p100)%in%m0.id], 
             p100[rownames(p100)==X, colnames(p100)%in%u.id], 
             alternative = "two.sided" )$p.value
    
    fc.t<- mean(p100[rownames(p100)==X, colnames(p100)%in%m0.id], na.rm = T) - 
      mean(p100[rownames(p100)==X, colnames(p100)%in%u.id], na.rm = T)
    
    pc<-c(pc, X); pval<-c(pval, pval.t); fc<-c(fc, fc.t)
  
  }
  
}

# Organize data
df100<-data.frame(pc, pval, fc)

# Take significantly changing proteins
df100<-df100[df100$pval<0.01, ]
df100<-df100[order(df100$fc, decreasing=T),]


# Take top 30 most differentiatial proteins, up and down: 
mono_genes<-df100[df100$fc<0, ]
mono_genes<-mono_genes[order(mono_genes$fc, decreasing=F), ]
mono_genes30<-mono_genes$pc[1:30]
mac_genes<-df100[df100$fc>0, ]
mac_genes30<-mac_genes$pc[1:30]

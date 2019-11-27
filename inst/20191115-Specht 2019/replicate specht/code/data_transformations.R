####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#  Normalization 
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# One step: 
#ev.matrix.sc.f.n<- collapse_to_protein(log2(cr_norm(filt.mat.rc(ev.matrix.sc.f, na.row, na.col))), ev.melt, T )


# Perform normalizations / transformations in multiple steps with visual sanity checks: 
b.t<-"FD"
xlim.t<-c(-2,2)
par(mfrow=c(3,3))

# Original data, normalized to reference channel, filtered for failed wells: 
t0<-ev.matrix.sc.f
hist(c(t0), breaks=b.t, xlim=xlim.t)

# Column then row normalize by median or mean (see source functions): 
t1<-cr_norm(t0)
hist(c(t1), breaks=b.t, xlim=xlim.t)

# Filter for missing data: 
t2<-filt.mat.rc(t1, na.row, na.col)
hist(c(t2), breaks=b.t, xlim=xlim.t)

# Log2 transform: 
t3<-log2(t2)
hist(c(t3), breaks=b.t, xlim=xlim.t)

# Collapse to protein level by median: 
t4<-collapse_to_protein(t3, ev.melt, T)
hist(c(t4), breaks=b.t, xlim=xlim.t)

# Re-column and row normalize: 
t4b<-cr_norm_log(t4)
hist(c(t4b), breaks=b.t, xlim=xlim.t)

# Assign to a final variable name: 
ev.matrix.sc.f.n<-t4b

# Perform similar operations for the 10 and 100-cell data (10-cell data not used in publication)
ev.matrix.10.n<-( cr_norm(ev.matrix.10) )
ev.matrix.100.n<-( cr_norm(ev.matrix.100) )

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#  Imputation 
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## Single cells
imp.input<-ev.matrix.sc.f.n

sc.imp <- hknn(imp.input, k.t)

t5<-sc.imp

####################################################################################################################################

## 10 cells
imp.input<-ev.matrix.10.n

# Filter out rows then columns with High missing data: 
imp.input<-filt.mat.rc(imp.input, na.row, na.col)

# Impute missing data
c10.imp <- hknn(imp.input, k.t)

####################################################################################################################################

## 100 cells
imp.input<-ev.matrix.100.n

# Filter out rows then columns with High missing data: 
imp.input<-filt.mat.rc(imp.input, na.row, na.col)

# Impute missing data
c100.imp <- hknn(imp.input, k.t)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#  Batch correction 
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Single cells
# Define the batches and model:
batch.covs <- ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id %in% colnames(sc.imp)]
#batch.covs <- ev.melt.uniqueID$lcbatch[ev.melt.uniqueID$id %in% colnames(sc.imp)]
mod<-data.frame(ev.melt.uniqueID$celltype[ev.melt.uniqueID$id %in% colnames(sc.imp)]); colnames(mod)<-"celltype"
mod<-model.matrix(~as.factor(celltype), data=mod)

# Correct
#matrix.sc.batch <- ComBat(sc.imp, batch=batch.covs, mod=mod, par.prior=F)
matrix.sc.batch <- ComBat(sc.imp, batch=batch.covs, mod=mod, par.prior=T)

t6<-matrix.sc.batch
 
####################################################################################################################################

## 10 cells

batch.covs <- ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id %in% colnames(c10.imp)]
mod<-data.frame(ev.melt.uniqueID$celltype[ev.melt.uniqueID$id %in% colnames(c10.imp)]); colnames(mod)<-"celltype"
mod<-model.matrix(~as.factor(celltype), data=mod)

matrix.c10.batch <- ComBat(c10.imp, batch=batch.covs, mod=mod)

####################################################################################################################################

## 100 cells

batch.covs <- ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id %in% colnames(c100.imp)]
mod<-data.frame(ev.melt.uniqueID$celltype[ev.melt.uniqueID$id %in% colnames(c100.imp)]); colnames(mod)<-"celltype"
mod<-model.matrix(~as.factor(celltype), data=mod)

matrix.c100.batch <- ComBat(c100.imp, batch=batch.covs, mod=mod)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# visual sanity checks post-imputation: 
hist(c(t5), breaks=b.t, xlim=xlim.t)
hist(c(t6), breaks=b.t, xlim=xlim.t)

par(mfrow=c(1,1))




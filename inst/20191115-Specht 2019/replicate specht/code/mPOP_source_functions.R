####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Check for presence of required packages for the import script
req.pkg<-c("plyr", "gridExtra", "stringr", "ggpubr", "Rcpp", "reshape2",
           "ggplot2", "gridExtra", "gridGraphics", "grid", "ggridges",
           "matrixStats", "scales","stylo")
# "patchwork"

present.pkg<-installed.packages()[,1]

if(any(!req.pkg%in%present.pkg)){

  install.packages( req.pkg[ !req.pkg%in%present.pkg ] )

}

## If installing impute or sva for the first time, uncomment the 3 lines below: 

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("impute", version = "3.8")


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

### load libraries and functions
library(reshape2)
library(scales)
library(ggplot2)
library(gridGraphics)
library(grid)
library(gridExtra)
library(stringr)
library(ggpubr)
library(ggridges)
library(matrixStats)
# library(patchwork)
library(impute)
# library("sva")
library(stylo)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Column and row normalize

cr_norm<-function(dat){
  
  for(k in 1:ncol(dat)){
    
    dat[,k]<-dat[,k]/median(dat[,k], na.rm = T)
    
    
  }
  
  
  for(k in 1:nrow(dat)){
    
    dat[k,]<-dat[k,]/mean(dat[k,], na.rm = T)
    
  }
  
  return(dat)
  
}



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Column and row normalize on log scale

cr_norm_log<-function(dat){
  
  for(k in 1:ncol(dat)){
    
    dat[,k]<-dat[,k]-median(dat[,k], na.rm = T)
    
    
  }
  
  
  for(k in 1:nrow(dat)){
    
    dat[k,]<-dat[k,]-mean(dat[k,], na.rm = T)
    
  }
  
  return(dat)
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# K - nearest neighbors imputation

hknn<-function(dat, k){
  
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      #print(length(closest.columns))
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
        
      }
      
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ] )
        
      }
      
      
    }
    
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
    
  }
  
  return(dat.imp)
  
}



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Calculate coefficient of variation for each and every protein with >= N peptides

prot.cv<-function(mat.t, meta, npeps){
  
  # Get rid of modified peptides
  modified.ind<-grep("(",row.names(mat.t), fixed=T)
  mat.t<-mat.t[-modified.ind,]
  
  # Row normalization
  for(k in 1:nrow(mat.t)){
    
    mat.t[k,]<-mat.t[k,]/mean(mat.t[k,], na.rm = T)
    #mat.t[k,]<-mat.t[k,] - mean(mat.t[k,], na.rm = T)
    
    
  }
  
  prots<-unique(ev.melt$protein)
  
  cv.mat<-matrix(data=NA, nrow=length(prots), ncol=ncol(mat.t))
  
  for(i in 1:nrow(cv.mat)){
    
    peps<-unique( ev.melt$sequence[ev.melt$protein%in%prots[i]] )
    
    if(length(peps)>=npeps){
      
      values.t<-mat.t[rownames(mat.t)%in%peps, ]
      
      # add matrix that is count of peptides going into each cv value
      cvs<-t(sqrt( rowVars(t(values.t), na.rm=T) ) / rowMeans(t(values.t), na.rm=T) )
      
      cv.mat[i,]<-cvs
      
    }
    
    
  }
  
  rownames(cv.mat)<-prots
  colnames(cv.mat)<-colnames(mat.t)
  
  return(cv.mat)
  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Calculate coefficient of variation for each and every protein with >= N peptides

prot.cv.log<-function(mat.t, meta, npeps){
  
  # Get rid of modified peptides
  modified.ind<-grep("(",row.names(mat.t), fixed=T)
  mat.t<-mat.t[-modified.ind,]
  
  # Row normalization
  for(k in 1:nrow(mat.t)){
    
    #mat.t[k,]<-mat.t[k,]/mean(mat.t[k,], na.rm = T)
    mat.t[k,]<-mat.t[k,] - mean(mat.t[k,], na.rm = T)
    
    
  }
  
  prots<-unique(ev.melt$protein)
  
  cv.mat<-matrix(data=NA, nrow=length(prots), ncol=ncol(mat.t))
  
  for(i in 1:nrow(cv.mat)){
    
    peps<-unique( ev.melt$sequence[ev.melt$protein%in%prots[i]] )
    
    if(length(peps)>=npeps){
      
      values.t<-mat.t[rownames(mat.t)%in%peps, ]
      
      # add matrix that is count of peptides going into each cv value
      cvs<-t(sqrt( rowVars(t(values.t), na.rm=T) ) / rowMeans(t(values.t), na.rm=T) )
      
      cv.mat[i,]<-cvs
      
    }
    
    
  }
  
  rownames(cv.mat)<-prots
  colnames(cv.mat)<-colnames(mat.t)
  
  return(cv.mat)
  
}



################################################################################################################
# remove.duplicates - remove duplicates across multiple columns
################################################################################################################

# Find unique entries by column1 and column2: 
# Ex remove.duplicates(TMT50,c("Sequence","Charge"))
remove.duplicates<-function(data,Cols){
  
  return(data[!duplicated(data[,Cols]),])
  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Take mean or median, removing NA values
median.na<-function(y){ median(y, na.rm = T) }
mean.na<-function(x){ mean(x, na.rm=T) }


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Color palette for cells / control well

my_cell_colors<-c(rgb(0,0,0,0.9),rgb(0,0,1,0.7),rgb(1,0,0,0.7))


# Calculate False Discovery rate from PEP values
calc_fdr <- function(pep) {
  
  return( (cumsum(pep[order(pep)]) / seq(1, length(pep)))[order(order(pep))] )
  
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Count the number of NA values in each row of a matrix
na.count<-function(data){
  
  na.v<-c()
  for(i in 1:nrow(data)){
    
    na.v<-c(na.v, length(which(is.na(data[i, ]))) )
    
  }
  
  return(na.v)
  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Matrix missing data filter function
# Filter remove COLUMNS then ROWS that have percent missing data
# (NA) greater than the specified thresholds
# Return the filtered matrix


filt.mat.cr<-function(mat, pct.r,pct.c){
  
  kc<-c()
  for(k in 1:ncol(mat)){
    
    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)
    
    
  }
  
  mat<-mat[,kc]
  
  kr<-c()
  for(k in 1:nrow(mat)){
    
    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)
    
    
  }
  
  mat<-mat[kr,]
  
  
  
  return(mat)
  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Matrix missing data filter function
# Filter remove ROWS then COLUMNS that have percent missing data
# (NA) greater than the specified thresholds
# Return the filtered matrix


filt.mat.rc<-function(mat, pct.r,pct.c){
  
  kr<-c()
  for(k in 1:nrow(mat)){
    
    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)
    
  }
  
  mat<-mat[kr,]
  
  kc<-c()
  for(k in 1:ncol(mat)){
    
    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)
    
  }
  
  mat<-mat[,kc]
  

  
  
  
  return(mat)
  
}



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Function to perform collapse, by median

collapse_to_protein<-function(ev.matrix.t, ev.melt, row_norm){
  
  ev.matrix.t.prot<-ev.matrix.t[0,]
  
  if(row_norm){
    
    # Row normalization
    for(k in 1:nrow(ev.matrix.t)){
      
      #ev.matrix.t[k,]<-ev.matrix.t[k,]/median(ev.matrix.t[k,], na.rm = T)
      ev.matrix.t[k,]<-ev.matrix.t[k,] - mean(ev.matrix.t[k,], na.rm = T)
      
    }
  }
  
  prot.v<-c()
  for(X in unique(as.character(ev.melt$protein))){
    
    rows.t<- which(rownames(ev.matrix.t) %in% ev.melt$sequence[ev.melt$protein%in%X])
    
    if(length(rows.t)>0){
      
      #print(length(rows.t))
      mat.t<-matrix(nrow = 1, ncol = ncol(ev.matrix.t))
      if(length(rows.t)>1){   mat.t<-apply(ev.matrix.t[rows.t,],2,median.na) }
      if(length(rows.t)==1){   mat.t<-ev.matrix.t[rows.t,] }
      
      ev.matrix.t.prot<-rbind(ev.matrix.t.prot, mat.t )
      
      prot.v<-c(prot.v, X)
      
    }
    
  }
  
  row.names(ev.matrix.t.prot)<-prot.v
  
  print( dim(ev.matrix.t.prot) )
  
  return( ev.matrix.t.prot )
  
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

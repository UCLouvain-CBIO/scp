path.t<-"./inst/scripts/replicate specht/"
setwd(path.t)

### Code settigngs: 
load_fresh = TRUE
only_fp94_97 = TRUE
# Remove peptides that are more abundant that 10% of the carrier channel from the single cell runs
remove_abundant_peptides<-T
# Use dart-enhanced spectra or spectra
dart_or_spectra<-"dart"
# Normalize the 10 x 10 cell and 10 x 100 cell runs to their median value
norm_10_100_median<-T
# Correct for isotopic cross-contamination from TMT11 -- this will not work properly because MQ did the correction already 
# for the data I am using
iso_cor<-F
# Minimum number of peptides observed to consider an experiment worth doing further analysis: 
thres_ids_sc <- 300
thres_ids_c10c100 <- 200
# Figure dimensions in inches
w.t<-5
h.t<-5
# Filter to less than X missing data per column or row: 
na.col<-0.99
na.row<-0.99
# Imputation
k.t<-3


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Execute all scripts in this order, some later scripts contain variables created in previous scripts

source("./code/mPOP_source_functions.R") 
source("./code/import.R") # 
source("./code/score.R") # 
source("./code/data_transformations.R") # 


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Go into the code 

### load MaxQuant data, experimental design, and meta-data (batch annotation, et cetera)...
ev<-read.delim(paste0(path.t,"ev_updated.txt"))
design<-read.csv(paste0(path.t,"annotation_fp60-97.csv"))
batch<-read.csv(paste0(path.t,"batch_fp60-97.csv"))

# Attach batch data to protein data 
ev[,colnames(batch)[-1]]<-NA
for(X in batch$set){
  
  ev$lcbatch[grep(X,ev$Raw.file)] <- as.character(batch$lcbatch[batch$set%in%X])
  ev$sortday[grep(X,ev$Raw.file)] <- as.character(batch$sortday[batch$set%in%X])
  ev$digest[grep(X,ev$Raw.file)] <- as.character(batch$digest[batch$set%in%X])
  
}

save(ev, file = "../extdata/specht2019/ev_updated_preloaded.rda")
load("../extdata/specht2019/ev_updated_preloaded.rda")

# Create unique peptide+charge column:
ev$modseq<-paste0(ev$Modified.sequence,ev$Charge)

# If somehow we have NA for raw files?! Remove 'em
ev<-ev[!is.na(ev$Raw.file),]

# Add X in front of experiment names because R doesn't like column names starting with numbers
ev$Raw.file<-paste0("X",ev$Raw.file)

# Group the experimental runs by type: single cell, 10 cell, 100 cell, 1000 cell, or ladder
all.runs<-unique(ev$Raw.file)
c10.runs<-c( all.runs[grep("col19", all.runs)] , all.runs[grep("col20", all.runs)] )
c100.runs<-c( all.runs[grep("col21", all.runs)] , all.runs[grep("col22", all.runs)] )
c1000.runs<-c( all.runs[grep("col23", all.runs)] )
ladder.runs<-c( all.runs[grep("col24", all.runs)] )
sc.runs<-all.runs[!all.runs%in%c(c10.runs,c1000.runs,c100.runs,ladder.runs)]

# Remove blank runs, if any
if(length(grep("blank",sc.runs))>0){
  sc.runs<-sc.runs[-grep("blank",sc.runs)]
}

# Consider only experiments FP94, 97
i1<-grep("FP94", sc.runs)
i2<-grep("FP97", sc.runs)
remove.sc.runs<-sc.runs[!sc.runs%in%sc.runs[c(i1,i2)]]
ev<-ev[!ev$Raw.file%in%remove.sc.runs, ]

# Which columns hold the TMT Reporter ion (RI) data
ri.index<-which(colnames(ev)%in%paste0("Reporter.intensity.",0:10))

if(iso_cor){ # Correct for isotopic cross-contamination from TMT11 -- this will not work properly because MQ did the correction already for the data I am using so iso_cor == F
  # Correct for Isotopic cross contamination
  cor.mat<-as.matrix(read.csv("te269088_lot_correction.csv")[,-1])
  corrected.ri<-t( solve(cor.mat) %*% t( ev[,ri.index] ) ) 
  # Remove 0 or negative values
  corrected.ri[corrected.ri<0.1]<-NA
  ev[,ri.index]<-corrected.ri
}

# Filtering

# Make sure all runs are described in design, if not, print and remove them:
not.described<-unique(ev$Raw.file)[ !unique(ev$Raw.file) %in% colnames(design) ]
ev<-ev[!ev$Raw.file%in%not.described,]

# Filter out reverse hits, contaminants, and contaminated spectra... 
ev<-ev[which(ev$Reverse!="+"),]
ev<-ev[-grep("REV", ev$Leading.razor.protein),]
ev<-ev[ev$Potential.contaminant!="+",]
ev<-ev[ev$PIF>0.8,]
ev<-ev[!is.na(ev$PIF),]

write.csv(ev, "ev_prefiltering.csv")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Take DART-ID or spectra peptide identifications from proteins that have been detected with < 1% FDR

if(dart_or_spectra=="spectra"){
  
  qprots<-unique(ev[ev$razor_protein_fdr<0.01,"Leading.razor.protein"])
  ev<-ev[ev$Leading.razor.protein%in%qprots, ]
  
  ev<-ev[ev$PEP<0.02,]
  
  
} else if(dart_or_spectra=="dart"){
  
  qprots<-unique(ev[ev$dart_qval<0.01,"Leading.razor.protein"])
  ev<-ev[ev$Leading.razor.protein%in%qprots, ]
  
  ev<-ev[ev$dart_PEP<0.02,]
  
}


#write.csv(ev, "ev_postfiltering.csv")

print(paste0("After filtering, there are ", length(unique(ev$Leading.razor.protein)), 
             " unique proteins identified across ", length(unique(ev$Raw.file)), " runs..."))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Remove sets with less than X peptides from 10, 100 cell runs:
ev<-ev[ !( (ev$Raw.file%in%names(table(ev$Raw.file))[table(ev$Raw.file)<thres_ids_c10c100])&(ev$Raw.file%in%c(c10.runs, c100.runs)) ) , ]

# Remove sets with less than X peptides from single cell runs:
ev<-ev[ !( (ev$Raw.file%in%names(table(ev$Raw.file))[table(ev$Raw.file)<thres_ids_sc])&(ev$Raw.file%in%sc.runs) ), ]


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Remove peptides that are more the 10% the intensity of the carrier in the single cell runs (only)

if(remove_abundant_peptides) {
  sc.index <- which(ev$Raw.file %in% sc.runs)
  
  remove.ind <- c()
  for (i in sc.index) {
    rat.t <-
      mean(as.numeric(ev[i, ri.index[4:11]] / ev[i, ri.index[1]]), na.rm = T)
    
    if (!any(is.na(rat.t))) {
      if (any(rat.t > 0.1)) {
        remove.ind <- c(remove.ind, i)
        
      }
    }
  }
  
  print(
    paste0(
      "Removing ",
      round(length(remove.ind) / length(sc.index), 2) * 100,
      "% of peptides, these are more abundant than expected for single cells."
    )
  )
  
  ev <- ev[-remove.ind,]
  
}




####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Numeric normalizations: 

# Normalize single cell runs to normalization channel
ev[ev$Raw.file%in%sc.runs, ri.index] <- ev[ev$Raw.file%in%sc.runs, ri.index] / ev[ev$Raw.file%in%sc.runs, ri.index[2]]

# Normalize single cell runs to carrier
evc<-ev
evc[evc$Raw.file%in%sc.runs, ri.index] <- evc[evc$Raw.file%in%sc.runs, ri.index] / evc[evc$Raw.file%in%sc.runs, ri.index[1]]


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Organize data into a more convenient data structure:

# Create empty data frame
ev.melt<-melt(ev[0, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev)[ri.index]) ], 
              id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest")) 

colnames(ev.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")


# Record mapping of cell type to Channel:
ct.v<-c()
qt.v<-c() 

# Create a unique ID string
unique.id.numeric<-1:11
unique.id<-paste0("i",unique.id.numeric)

# Reorganize data into a more useable format

for(X in unique(ev$Raw.file)){
  
  # Subset data by X'th experiment
  ev.t<-ev[ev$Raw.file%in%X, ]
  
  # Which experimental design used
  design.index<-which(colnames(design)%in%X)
  
  if(is.na(X)){next}
  
  # Name the RI columns by what sample type they are: carrier, single cell, unused, etc... 
  colnames(ev.t)[ri.index]<-paste0(as.character(design[, design.index]),"-", unique.id)
  
  RI_keep<-ri.index
  
  if(length(RI_keep)>0){
    
    if( X%in%c(c10.runs, c100.runs) ){
      
      ev.t[,RI_keep]<-ev.t[,RI_keep] / apply(ev.t[, RI_keep], 1, median.na)
      
    }
    
    # Melt it! and combine with other experimental sets
    ev.t.melt<-melt(ev.t[, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev.t)[RI_keep]) ],
                    id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"));
    
    # Record mapping of cell type to Channel:
    ct.v<-c(ct.v, unique.id[which(ri.index%in%RI_keep)] )
    qt.v<-c(qt.v, colnames(ev)[RI_keep] ) 
    
    colnames(ev.t.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")
    
    ev.melt<-rbind(ev.melt, ev.t.melt)
    
  }
  
  # Update unique ID string
  unique.id.numeric<-unique.id.numeric + 11
  unique.id<-paste0("i", unique.id.numeric)
  
}

c2q<-data.frame(ct.v, qt.v)
colnames(c2q)<-c("celltype","channel")

# Grab the unique number associate to each and every cell, carrier channel, and empty channel
ev.melt$id<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(2,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$celltype<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(1,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$id<-as.factor(ev.melt$id)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Creation of data matrix

# Remove duplicate observations of peptides from a single experiment
ev.melt<-remove.duplicates(ev.melt,c("sequence","id") )

# Create data frame of peptides x cells, populated by quantitation
ev.unmelt<-dcast(ev.melt, sequence ~ id, value.var = "quantitation", fill=NA)

# Also create matrix of same shape
ev.matrix<-as.matrix(ev.unmelt[,-1]); row.names(ev.matrix)<-ev.unmelt$sequence

# Create a matrix where 0s observed will be set to a low number
ev.matrix.lowbin<-ev.matrix
ev.matrix.lowbin[ev.matrix.lowbin==0]<-0.001
ev.matrix.lowbin[ev.matrix.lowbin==Inf]<-0.001
ev.matrix.lowbin[ev.matrix.lowbin==-Inf]<-0.001

# Create a matrix where 0s observed will be set to a low number
ev.matrix.lowbin0<-ev.matrix
ev.matrix.lowbin0[ev.matrix.lowbin0==0]<-0
ev.matrix.lowbin0[ev.matrix.lowbin0==Inf]<-0
ev.matrix.lowbin0[ev.matrix.lowbin0==-Inf]<-0

# Replace all 0s with NA
ev.matrix[ev.matrix==0]<-NA
ev.matrix[ev.matrix==Inf]<-NA
ev.matrix[ev.matrix==-Inf]<-NA


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Create meta data matrices

ev.melt.uniqueID<-remove.duplicates(ev.melt,"id")
ev.melt.pep<-remove.duplicates(ev.melt, c("sequence","protein") )

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Divide and Save

# Divide matrix into single cells (including intentional blanks) and carriers
sc_cols<-unique(ev.melt$id[(ev.melt$celltype%in%c("sc_u","sc_m0", "sc_0"))&(ev.melt$Raw.file%in%sc.runs)])
ev.matrix.sc<-ev.matrix[, sc_cols]
#ev.matrix.sc.lb<-ev.matrix.lowbin[, sc_cols]
#ev.matrix.sc.lb0<-ev.matrix.lowbin0[, sc_cols]

# QC: PCA to spearate macor and mono
test <- ev.matrix.sc[rowSums(is.na(ev.matrix.sc)) != ncol(ev.matrix.sc), ]
ed_specht <- read.csv("../../../scpdata/inst/extdata/specht2019/peptides-raw-RI.csv", row.names = 1)
pca <-  nipals::nipals(test, ncomp = 5, center = TRUE, scale = TRUE)
PCs <- pca$loadings
celltype <- ev.melt[, c("celltype", "id")]
celltype <- celltype[!duplicated(celltype), ]
rownames(celltype) <- celltype$id
celltype <- celltype[colnames(test), "celltype"]
ggplot(data = data.frame(PCs, celltype)) +
  geom_point(aes(x = PC1, y = PC4, col= celltype)) +
  ggtitle("PCA on the expression data")
# No separation ! 

# # Divide matrix into single cells (including intentional blanks) and carriers
# c0_cols<-unique(ev.melt$id[(ev.melt$celltype%in%c("sc_0"))&(ev.melt$Raw.file%in%sc.runs)])
# ev.matrix.0<-ev.matrix[, c0_cols]
# 
# # Carriers
# c_cols<-unique(ev.melt$id[grep("carrier",ev.melt$celltype)])
# ev.matrix.carrier<-ev.matrix[, c_cols]

# 10 cells
b_cols<-unique(ev.melt$id[ev.melt$celltype%in%c("u_10","m0_10")&(ev.melt$Raw.file%in%c10.runs)])
ev.matrix.10<-ev.matrix[, b_cols]

# 100 cells
d_cols<-unique(ev.melt$id[ev.melt$celltype%in%c("u_100","m0_100")&(ev.melt$Raw.file%in%c100.runs)])
ev.matrix.100<-ev.matrix[, d_cols]

# # Ladder cells
# e_cols<-unique(ev.melt$id[ev.melt$Raw.file%in%c(ladder.runs)])
# ev.matrix.ladder<-ev.matrix[, e_cols]

# # Double check we grabbed the right samples for each matrix
# ev.melt.uniqueID$celltype[ev.melt.uniqueID$id%in%colnames(ev.matrix.sc)]
# ev.melt.uniqueID$celltype[ev.melt.uniqueID$id%in%colnames(ev.matrix.10)]
# ev.melt.uniqueID$celltype[ev.melt.uniqueID$id%in%colnames(ev.matrix.100)]
# 
# Double check we grabbed the right samples for each matrix
ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id%in%colnames(ev.matrix.sc)]
ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id%in%colnames(ev.matrix.10)]
ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id%in%colnames(ev.matrix.100)]


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Outlier removal

# Remove FP97_col20 RI 6 --> due to missing data causing imputation artifacts

ev.matrix.10<-ev.matrix.10[ , !colnames(ev.matrix.10)=="i1206" ]

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Repeat the above with data normalized to carrier channel--> no column or row normalization, called "evc"

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Organize data into a more convenient data structure:

# Create empty data frame
evc.melt<-melt(evc[0, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(evc)[ri.index]) ], 
               id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest")) 

colnames(evc.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")


# Record mapping of cell type to Channel:
ct.v<-c()
qt.v<-c() 

# Create a unique ID string
unique.id.numeric<-1:11
unique.id<-paste0("i",unique.id.numeric)

# Tranform data into a more useable format

for(X in unique(evc$Raw.file)){
  
  # Subset data by X'th experiment
  evc.t<-evc[evc$Raw.file%in%X, ]
  
  # Which experimental design used
  design.index<-which(colnames(design)%in%X)
  
  if(is.na(X)){next}
  
  # Name the RI columns by what sample type they are: carrier, single cell, unused, etc... 
  colnames(evc.t)[ri.index]<-paste0(as.character(design[, design.index]),"-", unique.id)
  
  RI_keep<-ri.index
  
  if(length(RI_keep)>0){
    
    if( X%in%c(c10.runs, c100.runs) ){
      
      evc.t[,RI_keep]<-evc.t[,RI_keep] / apply(evc.t[, RI_keep], 1, median.na)
      
    }
    
    # Melt it! and combine with other experimental sets
    evc.t.melt<-melt(evc.t[, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(evc.t)[RI_keep]) ],
                     id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"));
    
    # Record mapping of cell type to Channel:
    ct.v<-c(ct.v, unique.id[which(ri.index%in%RI_keep)] )
    qt.v<-c(qt.v, colnames(evc)[RI_keep] ) 
    
    colnames(evc.t.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")
    
    evc.melt<-rbind(evc.melt, evc.t.melt)
    
  }
  
  # Update unique ID string
  unique.id.numeric<-unique.id.numeric + 11
  unique.id<-paste0("i", unique.id.numeric)
  
}

c2q<-data.frame(ct.v, qt.v); colnames(c2q)<-c("celltype","channel")

# Grab the unique number associate to each and evcery cell, carrier channel, and empty channel
evc.melt$id<-unlist(strsplit(as.character(evc.melt$celltype),"-"))[seq(2,length(unlist(strsplit(as.character(evc.melt$celltype),"-"))),2)]
evc.melt$celltype<-unlist(strsplit(as.character(evc.melt$celltype),"-"))[seq(1,length(unlist(strsplit(as.character(evc.melt$celltype),"-"))),2)]
evc.melt$id<-as.factor(evc.melt$id)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Creation of data matrix

# Remove duplicate observations of peptides from a single experiment
evc.melt<-remove.duplicates(evc.melt,c("sequence","id") )

# Create data frame of peptides x cells, populated by quantitation
evc.unmelt<-dcast(evc.melt, sequence ~ id, value.var = "quantitation", fill=NA)

# Also create matrix of same shape
evc.matrix<-as.matrix(evc.unmelt[,-1]); row.names(evc.matrix)<-evc.unmelt$sequence

# Create a matrix where 0s observed will be set to a low number
evc.matrix.lowbin<-evc.matrix
evc.matrix.lowbin[evc.matrix.lowbin==0]<-0.001
evc.matrix.lowbin[evc.matrix.lowbin==Inf]<-0.001
evc.matrix.lowbin[evc.matrix.lowbin==-Inf]<-0.001

# Create a matrix where 0s observed will be set to a low number
evc.matrix.lowbin0<-evc.matrix
evc.matrix.lowbin0[evc.matrix.lowbin0==0]<-0
evc.matrix.lowbin0[evc.matrix.lowbin0==Inf]<-0
evc.matrix.lowbin0[evc.matrix.lowbin0==-Inf]<-0

# Replace all 0s with NA
evc.matrix[evc.matrix==0]<-NA
evc.matrix[evc.matrix==Inf]<-NA
evc.matrix[evc.matrix==-Inf]<-NA


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Create meta data matrix
evc.melt.uniqueID<-remove.duplicates(evc.melt,"id")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Divide and Save

# Divide matrix into single cells (including intentional blanks) and carriers
sc_cols<-unique(evc.melt$id[(evc.melt$celltype%in%c("sc_u","sc_m0", "sc_0"))&(evc.melt$Raw.file%in%sc.runs)])
evc.matrix.sc<-evc.matrix[, sc_cols]
evc.matrix.sc.lb<-evc.matrix.lowbin[, sc_cols]
evc.matrix.sc.lb0<-evc.matrix.lowbin0[, sc_cols]







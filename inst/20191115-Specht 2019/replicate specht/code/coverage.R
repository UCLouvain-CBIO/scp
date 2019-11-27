# Coverage boxplot

# Calculate old FDR from old pep values, filter by FDR:
ev2<-read.csv("ev_prefiltering.csv")

# Calculate FDR
ev2$spectra_fdr<-calc_fdr(ev2$PEP)

# Remove peptides/proteins that have missing quantitation for the single cells channels
ev3<-ev2[-which(rowSums(ev2[,ri.index[4:11]], na.rm=T)==0), ]

# Filter by PEP
ev.spectra<-ev3[ev3$PEP<0.02, ]; nrow(ev.spectra)
ev.dart<-ev3[ev3$dart_PEP<0.02, ]; nrow(ev.dart)

## For spectra: 

# Grab single cell runs and count observations of (confident) peptides per run: 
sc.runs.f<-unique(ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id%in%colnames(ev.matrix.sc.f)])
sc.s<-table(ev.spectra[ev.spectra$Raw.file%in%sc.runs.f, "Raw.file"])
sc.s<-sc.s[sc.s>0]
sc.s.f<-rep("sc_spectra", length(sc.s))
df.s<-data.frame(numcells<-c(sc.s.f), npeps<-c(sc.s)); colnames(df.s)<-c("numcells","npeps")

## For dart: 

# Grab single cell runs and count observations of (confident) peptides per run: 
sc.s<-table(ev.dart[ev.dart$Raw.file%in%sc.runs.f, "Raw.file"])
sc.s<-sc.s[sc.s>0]
sc.s.f<-rep("sc_dart", length(sc.s))
df.d<-data.frame(numcells<-c(sc.s.f), npeps<-c(sc.s)); colnames(df.d)<-c("numcells","npeps")

## Combine results for spectra-only and dart, reorganize: 
df.f<-rbind(df.s,df.d)

df.f$numcells<-factor(df.f$numcells, levels=c("sc_spectra","sc_dart","c10_spectra","c10_dart","c100_spectra","c100_dart"))

df.f$type<-"Spectra only"
df.f$type[grep("dart", df.f$numcells)]<-"DART-ID"
df.f$type<-factor(df.f$type, levels=c("Spectra only","DART-ID"))

df.f$celltype<-NA
df.f$celltype[grep("sc", df.f$numcells)]<-"1-cell"
df.f$celltype[grep("c100_", df.f$numcells)]<-"100-cell"
df.f$celltype[grep("c10_", df.f$numcells)]<-"10-cell"

df.f1<-df.f

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Do the same as above for proteins

ev.spectra<-remove.duplicates(ev3[ev3$spectra_fdr<0.01, ], c("Leading.razor.protein", "Raw.file")); nrow(ev.spectra)
ev.dart<-remove.duplicates(ev3[ev3$razor_protein_fdr<0.01, ], c("Leading.razor.protein", "Raw.file")); nrow(ev.dart)

sc.s<-table(ev.spectra[ev.spectra$Raw.file%in%sc.runs.f, "Raw.file"])
sc.s<-sc.s[sc.s>0]
sc.s.f<-rep("sc_spectra", length(sc.s))
df.s<-data.frame(numcells<-c(sc.s.f), npeps<-c(sc.s)); colnames(df.s)<-c("numcells","npeps")

sc.s<-table(ev.dart[ev.dart$Raw.file%in%sc.runs.f, "Raw.file"])
sc.s<-sc.s[sc.s>0]
sc.s.f<-rep("sc_dart", length(sc.s))
df.d<-data.frame(numcells<-c(sc.s.f), npeps<-c(sc.s)); colnames(df.d)<-c("numcells","npeps")

df.f<-rbind(df.s,df.d)

df.f$numcells<-factor(df.f$numcells, levels=levels(df.f$numcells)[c(3,6,1,4,2,5)])

df.f$type<-"Spectra only"
df.f$type[grep("dart", df.f$numcells)]<-"DART-ID"
df.f$type<-factor(df.f$type)
df.f$type<-factor(df.f$type, levels=c("Spectra only","DART-ID"))

df.f$celltype<-NA
df.f$celltype[grep("sc", df.f$numcells)]<-"1-cell"
df.f$celltype[grep("c100_", df.f$numcells)]<-"100-cell"
df.f$celltype[grep("c10_", df.f$numcells)]<-"10-cell"

df.f2<-df.f


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Plotting! 


df.f1$molecule<-"Peptides"
df.f2$molecule<-"Proteins"

dfx<-rbind(df.f1, df.f2)

dfx.1c<-dfx[dfx$celltype=="1-cell",]

ggboxplot(dfx.1c, x="molecule", y="npeps", fill="type") + 
  theme(text=element_text(size=20)) + 
  xlab("")+
  ylab("# quantified / 1h run\n") +
  font("ylab",size=30)+
  font("xlab",size=30)+
  scale_fill_manual(values = c("white","lightpink") ) +
  rremove("x.ticks") + 
  font("x.text", size=20) +
  theme(axis.text.x  = element_text(angle=30, vjust=0.5)) + 
  theme(legend.position = c(0.7, 0.9)) + 
  theme(legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", 
                                         colour ="white")) +
  font("legend.title", size= 20) + 
  font("legend.text", size= 15) + 
  guides(colour = guide_legend(override.aes = list(shape = 15))) + 
  theme(legend.key.size = unit(1.2,"line")) + 
  rremove("legend.title") +   
  rremove("x.ticks") + 
  #theme(axis.text.x  = element_text(angle=30, vjust=0.5)) + 
  #theme(legend.position = c(0.35, 0.7)) + 
  #theme(legend.background = element_rect(fill="white",
  #size=1, linetype="solid", 
  #colour ="black")) +
  #font("legend.title", size= 20) + 
  #font("legend.text", size= 20) + 
  #guides(colour = guide_legend(override.aes = list(shape = 15))) + 
  #theme(legend.key.size = unit(2,"line")) + 
  rremove("legend.title") +
  rremove("xlab") + 
  ylim(0,3000) 


summary(dfx.1c[dfx.1c$type=="DART-ID" & dfx.1c$molecule=="Proteins","npeps" ])
median(dfx.1c[dfx.1c$type=="Spectra only" & dfx.1c$molecule=="Proteins","npeps" ])


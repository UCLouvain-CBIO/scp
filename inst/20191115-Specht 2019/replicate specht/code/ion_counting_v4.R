

# Import mRNA counts for CD14 monocyte cluster from Seurat tutorial
mono.rna<-read.csv("mono.csv")
mono.rna.m<-as.matrix(mono.rna[,-1]); rownames(mono.rna.m)<-mono.rna[,1]; rv<-mono.rna.m

# Import data (Proteome Discoverer (demo) search, confident PSM-level output, 1%FDR):
tev<-read.delim("pd_ioncounts_PSMs.txt")

# Take only single cell runs
tev<-tev[-grep("col",tev$Spectrum.File), ]

# Convert to MQ format sequence
ul1<-unlist(strsplit(as.character(tev$Annotated.Sequence), ".", fixed=T))
tev$modseq<-paste0("_", toupper( ul1[seq(2,length(ul1),3)] ), "_", tev$Charge)

tev$Spectrum.File<-paste0("X",tev$Spectrum.File)
colnames(tev)[30:40]<-colnames(ev)[ri.index]


# Remove overly-abundant peptides (same as import.R)
tev<-tev[!tev$modseq%in%unique(ev$modseq[remove.ind]),]

# Map TMT channels to cell IDs
ev.melt2<-ev.melt.uniqueID
ev.melt2$ri<-c2q$channel[match(ev.melt2$id, c2q$celltype)]
k.id<-colnames(ev.matrix.sc.f.n)
tev2<-tev
tev2[,30:40]<-NA

for(X in k.id){
  
  raw.t<-paste0(ev.melt2$Raw.file[ev.melt2$id==X], ".raw")
  ri.t<-as.character(ev.melt2$ri[ev.melt2$id==X])

  
  tev2[tev2$Spectrum.File==raw.t, ri.t] <-  tev[which(tev$Spectrum.File==raw.t), ri.t]
  
}

tev<-tev2

# Mapped mRNA genes from CD14 monocyte cluster to UniProt identifier 
r2p<-read.delim("rna_2_hprot_uniprot.tab")

# Filter mRNA data to common genes 
rv2<-rv
row.names(rv2)<-r2p$Entry[match(row.names(rv2), r2p[,1] )]
rv2<-rv2[row.names(rv2)%in%tev$Master.Protein.Accessions,]

# Omit 0 values from quantitative consideration
rv2[rv2==0]<-NA

# Take mean mRNA counts / single cell
rm<-rowMeans(rv2, na.rm=T)
names(rm)<-row.names(rv2)

# Filter protein data to common genes
tev<-tev[tev$Master.Protein.Accessions%in%row.names(rv2),]

# Assign Gene names
tev$gn<-r2p[match(tev$Master.Protein.Accessions,r2p$Entry),1]

# Denote reporter ion columns: 
RI<-33:40

# Unique peptide-charge entries
tev$ms<-paste0(tev$Annotated.Sequence,tev$Charge)

# Remove duplicated peptides / experiment
tev<-remove.duplicates(tev, c("ms", "Spectrum.File"))

# Melt protein data to obtain peptide-level S/N in a more convenient form
ppm<-melt(tev[,c(RI, which(colnames(tev)%in%c("ms", "Master.Protein.Accessions")))])

# Calculate mean S/N / single cell
ppv<-aggregate(data = ppm, value~ms, FUN=mean)

# Map peptides to parent protein (protein inference done by Proteome Discoverer)
ppm.u<-remove.duplicates(ppm, c("ms", "Master.Protein.Accessions"))
ppv$prot<-ppm.u$Master.Protein.Accessions[ match(ppv$ms, ppm.u$ms ) ]

# Create protein quantitation (sum of peptide S/N)
psum<-aggregate(data = ppv, value~prot, FUN=sum)

common<-intersect(names(rm), psum$prot)

rm<-rm[names(rm)%in%common]

# Create numeric vectors of S/N 
pepv<-as.numeric(ppv$value)
psum.n<-as.numeric(psum$value)

# Merge mRNA, protein, and peptide data: 
df<-data.frame( matrix(NA, nrow= max(c(length(psum.n), length(pepv), length(rm) ) ), ncol=3 ) ); colnames(df)<-c("rm","pm","ppm")

df$pm[1:length(psum.n)]<-psum.n;
df$rm[1:length(rm)]<-rm
df$ppm[1:length(pepv)]<-pepv

head(df)

# Convert S/N to number of ions on QE: 
noise<- 3.5*sqrt(240000 / 70000)
df[,2:3]<-df[,2:3]*noise

# Rearrange data for convenience
dfm<-melt(df)

# Calculate Poisson error: 
pep.error<-round( 1 / sqrt( median ( df$ppm , na.rm=T)),2)*100
prot.error<-round( 1 / sqrt( median ( df$pm , na.rm=T)),2)*100
rna.error<-round( 1 / sqrt( median ( df$rm , na.rm=T)),2)*100

# Log10 transform data
dfm$valuelog10<-log10(dfm$value)


library(scales)
#install.packages("devtools")
#devtools::install_github("thomasp85/patchwork")
library(patchwork)

pois_error<-function(x){ 1/sqrt(x) }

x1<-log10(seq(1,10000,1))
y1<-pois_error(10^(x1))

dat<-data.frame(x1,y1); colnames(dat)<-c("x","y")

# Visualize
p1<-ggplot(dfm, aes(x = valuelog10, y = variable, fill=variable)) +
  geom_density_ridges() + theme_ridges(center_axis_labels = TRUE) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(-0.1, 0)) +
  rremove("legend") +
  xlab("\ Copy number measured\n per single cell") +
  scale_fill_manual(values = c("blue", "red", "orange"))+
  scale_x_continuous(limits = c(0, 4.6), breaks = c(0,1,2,3,4),
                     labels = c(bquote(10^0),bquote(10^1),bquote(10^2),bquote(10^3),bquote(10^4) )) +
  annotate("text", x=3.5, y =2.5, label=paste0(length(psum.n), " Proteins"), color="red", size=9) +
  annotate("text", x=3.5, y =1.5, label=paste0(length(rm), " mRNAs"), color="blue", size=9)+
  annotate("text", x=3.5, y =3.35, label=paste0(length(pepv), " Peptides"), color="orange3", size=9)+
  ylab("Density") +
  theme(text=element_text(size=15)) +
  theme(text = element_text(size = 25)) + rremove("y.text") + rremove("y.ticks")+
  #ggtitle("Sampling molecules in single cells\n") +
  font("title", size=20) +
  font("x.text", size= 18)  +
  font("ylab",size= 30)


p2<-ggplot(dat, aes(x,y)) + geom_line(size=1, color="black") +
  scale_x_continuous(limits = c(0, 4.6), breaks = c(0,1,2,3,4)) + 
  ylab("CV,  %") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix="")) + 
  theme_ridges(center_axis_labels = TRUE) +
  #theme(panel.grid.major.y = element_blank())+
    rremove("x.text") + 
    rremove("xlab") + 
  geom_point(aes(x=log10( median ( df$rm , na.rm=T) ), y=rna.error/100), colour="blue", size=7) + 
  geom_point(aes(x=log10( median ( df$pm , na.rm=T) ), y=prot.error/100), colour="red", size=7) +
  geom_point(aes(x=log10( median ( df$ppm , na.rm=T) ), y=pep.error/100), colour="orange", size=7) +
  annotate("text",x=3.5, y=0.45, label = expression(paste(CV==frac(1, sqrt(n)))), size=8) + 
  annotate("text",x=3.5, y=0.9, label = "Counting error:", size=9) + 
  font("title", size=20) + 
  font("ylab",size=28) 
  #theme(axis.line.x = element_line(color="black")) 

plot_spacer() + p2 +plot_spacer()+ p1 + plot_layout(ncol = 1, heights = c(0.1,2,0.2, 5))

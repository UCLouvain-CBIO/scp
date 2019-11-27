# export

# Peptides-raw.csv
peptides.raw<-data.frame(t3)
peptides.raw$peptide<-rownames(t3)
peptides.raw$protein<-ev.melt.pep$protein[match(rownames(t3), ev.melt.pep$sequence)]

write.csv(peptides.raw[,c(ncol(peptides.raw),ncol(peptides.raw)-1 , 1:(ncol(peptides.raw)-2))], "final/export/Peptides-raw.csv", row.names = F)

# Proteins-processed.csv
prot.pro<-data.frame(matrix.sc.batch)
prot.pro$protein<-rownames(matrix.sc.batch)
write.csv(prot.pro, "final/export/Proteins-processed.csv")


# Cells.csv
head(ev.melt.uniqueID)
cells<-ev.melt.uniqueID[ev.melt.uniqueID$id%in%colnames(ev.matrix.sc.f.n), c("id","celltype","digest","sortday","lcbatch","Raw.file")]
row.names(cells) <- NULL
colnames(cells)<-c("id","celltype","batch_digest","batch_sort","batch_chromatography","raw.file")
cells.t<-t(cells)
colnames(cells.t)<-cells$id
cells.t<-cells.t[-1,]

write.csv(cells.t, "final/export/Cells.csv", row.names = T)


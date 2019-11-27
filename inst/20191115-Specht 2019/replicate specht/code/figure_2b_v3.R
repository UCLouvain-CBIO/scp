####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Protein level matrix, non-imputed, non-batch-corrected: 
mat<-ev.matrix.sc.f.n

# Seen in greater than or equal to X cells: 
ncells.v<-c(5,10,25,50, 100)

# Total number of single cells in data set
intervals<-c(50,75,100,125,150,200,250,300,ncol(mat))

# Create empty data frame: 
df<-data.frame(size.data<-c(), ncells<-c(), nprot<-c() ) 

# Set randomization seed: 
set.seed(20)

# Count the number of NA values per column and re-order the columns in mat by those values: 
mat2<-mat[, order(na.count(t(mat)), decreasing =F) ]

# Sample: 
for(Y in intervals){
  
  # Sampling approach: 
  mat.t<-mat[, sample(1:ncol(mat), Y, replace=F) ]
  
  #mat.t<-mat2[, 1:Y]
  
  ncells<-1:ncol(mat.t)
  
  # Number of observationsof each protein for the given matrix
  prot.obs<- ncol(mat.t) - na.count(mat.t)
  
  nprot<-c()
  for(X in ncells){
    
    nprot<-c(nprot, length(which(prot.obs>=X)))
    
  }
  
  fact.t<-rep(Y, length(nprot))
  
  df.t<-data.frame(size.data<-fact.t, ncells, nprot ) 
  
  df<-rbind(df, df.t)

}

colnames(df)<-c("size","ncells","nprot")
df$size<-as.factor(df$size)

head(df)

df2<-df[df$ncells%in%ncells.v,]
df2<-df2[df2$size%in%c(50,75, 100,125, 150,200,250, 300, ncol(mat)),]
df2$ncells<-as.factor(df2$ncells)

df3<-df2
df3$ncells<-as.numeric(as.character(df3$ncells))
df3$size<-as.numeric(as.character(df3$size))
df3$size<-as.numeric(as.character(df3$size))
df3$ncellsf<-as.factor(df3$ncells)
types<-unique(df2$ncells)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

my_cp<-gg_color_hue(length(unique(ncells.v)))

ggscatter(data = df3, x ="size", y="nprot", color = "ncellsf",size=4) + geom_line(aes(x=size, y=nprot, color=ncellsf)) + 
  xlab("\nTotal # single cells analyzed") +
  ylab("Proteins, 1% FDR\n") + 
  scale_x_continuous(limits=c(25,ncol(mat)+190))+
  theme(text=element_text(size=20))+
  ylim(c(0,  max(df3$nprot)+400))+
  #annotate("text", x= ncol(mat)+70, y= 6800, label = 'bold("No missing data in:")', size=8,color = rgb(0.1,0.1,1,1), parse=T)+
  annotate("text", x= ncol(mat)+45, y= max(df3$nprot)+350, label = 'bold("No missing data in:")', size=7, parse=T)+
  rremove("legend")+
  annotate("text", x=ncol(mat)+85, y= df3$nprot[(df3$size==ncol(mat))&(df3$ncells==ncells.v[1])], label = paste0("\u2265", ncells.v[1]," cells"), size=8, color=my_cp[1], fontface=2)+
  annotate("text", x= ncol(mat)+95, y= df3$nprot[(df3$size==ncol(mat))&(df3$ncells==ncells.v[2])]+10, label = paste0("\u2265", ncells.v[2]," cells"), size=8, color=my_cp[2], fontface=2)+
  annotate("text", x= ncol(mat)+95, y= df3$nprot[(df3$size==ncol(mat))&(df3$ncells==ncells.v[3])]-10, label = paste0("\u2265", ncells.v[3]," cells"), size=8, color=my_cp[3], fontface=2)+
  annotate("text", x= ncol(mat)+95, y= df3$nprot[(df3$size==ncol(mat))&(df3$ncells==ncells.v[4])]-10, label = paste0("\u2265", ncells.v[4]," cells"), size=8, color=my_cp[4], fontface=2)+
  annotate("text", x= ncol(mat)+105, y= df3$nprot[(df3$size==ncol(mat))&(df3$ncells==ncells.v[5])]-10, label = paste0("\u2265", ncells.v[5]," cells"), size=8, color=my_cp[5], fontface=2)+
  theme(plot.title = element_text(hjust = 1.5, vjust=0)) +
  font("xylab", size=25)  +
  theme(panel.grid.minor = element_blank(), 
       panel.grid.major = element_line("gray80", size = 0.1) )


  
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


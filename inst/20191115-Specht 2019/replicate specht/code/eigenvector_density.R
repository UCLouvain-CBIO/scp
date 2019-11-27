### Plot eigenvectors of lapalacians

# Eigenvectors used to perform spectral clustering on both macrophage-like + monocyte, and macrophage-like only data matrices: 
mac_vec<-c(mac_vec, rep(NA, 97))
dfxx<-data.frame(mac_mono_vec, mac_vec)
dfxm<-melt(dfxx)

# Plot! 
ggplot(data=dfxm, aes(x=value, y=variable)) + geom_density_ridges(aes(fill=variable)) + theme_pubr() + 
  scale_fill_manual(values=my_colors[c(1,2)]) + 
  xlab("Eigenvector value") + ylab("Density") + rremove("y.ticks") + rremove("y.text") + 
  font("xylab", size=20) +
  font("x.text", size=20) + 
  xlim(c(-0.15, 0.35)) + 
  #annotate("text", x=0.25, y= 1.5, label="Monocyte and\n macrophage-like", size=6)+ 
  #annotate("text", x=0.25, y= 3, label="Macrophage-like", size=6) + 
  rremove("legend")
  


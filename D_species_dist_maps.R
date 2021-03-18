################################################################################################################
###                                                                                                          ###
###                       Comparing Multiple Diversities Across Drainages                                    ###
###                                           ZD ZBINDEN et al. 2020                                         ###
###                       *** This script is for mapping species distributions ***                           ###
################################################################################################################


# Multiple plot function
##############################################################################################
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
############################################################################################




## Reading spatial data ###
############################
setwd("C:/Users/zbind/Desktop/shelf/Projects/Drainage_Compare_SWG/R_workdir")
library(ggplot2)


# reading HydroRIVERS shapefile into R (modified this file in ArcGIS)
library(raster)
hydroRiv <- shapefile("./hydroRivers/hydrosheds_reduced_Project.shp", verbose=TRUE, warnPRJ=TRUE)
plot(hydroRiv) # check map looks OK
# reading fish collection site coordinates N=142 (this file was also projected in ArcGIS)
sites <- shapefile("sites_shapefiles/Sites_project.shp")
points(sites) # check map looks OK


# read environmental data (edited after site shapefile was made)
env <- read.csv("env_data_basin_comp.csv")

# read in fish abundance data file
taxon.div <- read.csv("fish_data_basin_comp.csv", row.names = 1)
taxon.div <- t(taxon.div)
taxon.div <- as.data.frame(taxon.div)
taxon.incidence <- taxon.div
taxon.incidence[taxon.incidence >0] <-1
### snapping points to river lines ###
library(sf)
sites<-st_as_sf(sites)
hydroRiv <-st_as_sf(hydroRiv)
join<-st_join(sites, hydroRiv, st_nearest_feature,left=TRUE) # combined in situ ENV variables with hydroATLAS variables
join <- cbind.data.frame(join, taxon.incidence)
join <- st_as_sf(join)



#####################################################################
###                                                                 #
###     ICTALURIDAE -- CATFISH                                      #
###                                                                 #
####################################################################
#### Ameiurus melas
library(sp)
pdf("./species_dist_maps/Ictaluridae.pdf", paper = "a4r", width = 0, height = 0)

### Ameiurus melas
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Ameiurus_melas), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                   high = "red",
                   aesthetics = "fill")+
  ggtitle("Ameiurus melas") +
  theme_classic() +
  theme(legend.position = "none")

### Ameiurus natalis
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Ameiurus_natalis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Ameiurus natalis") +
  theme_classic() +
  theme(legend.position = "none")

### Ictalurus punctatus
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Ictalurus_punctatus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Ictalurus punctatus") +
  theme_classic() +
  theme(legend.position = "none")

### Noturus exilus
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Noturus_exilus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Noturus exilus") +
  theme_classic() +
  theme(legend.position = "none")

### Noturus gyrinus
e<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Noturus_gyrinus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Noturus gyrinus") +
  theme_classic() +
  theme(legend.position = "none")

### Noturus nocturnus
f<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Noturus_gyrinus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Noturus nocturnus") +
  theme_classic() +
  theme(legend.position = "none")

# multiplot
multiplot(a,b,c,d,e,f, cols=2)
dev.off()






#####################################################################
###                                                                 #
###     Catostomidae -- Suckers                                      #
###                                                                 #
####################################################################
library(sp)
pdf("./species_dist_maps/Catostomidae.pdf", paper = "a4r", width = 0, height = 0)
### Erimyzon oblongus
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Erimyzon_oblongus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Erimyzon oblongus") +
  theme_classic() +
  theme(legend.position = "none")

### Minytrema melanops
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Minytrema_melanops), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Minytrema melanops") +
  theme_classic() +
  theme(legend.position = "none")

### Moxostoma duqesnei 
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Moxostoma_duquesnei), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Moxostoma duquesnei") +
  theme_classic() +
  theme(legend.position = "none")

### Moxostoma erythrurum
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Moxostoma_erythrurum), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Moxostoma erythrurum") +
  theme_classic() +
  theme(legend.position = "none")

# multiplot
multiplot(a,b,c,d, cols=2)
dev.off()





#####################################################################
###                                                                 #
###     Lepomis  -- Sunfish                                         #
###                                                                 #
####################################################################
library(sp)
pdf("./species_dist_maps/Lepomis.pdf", paper = "a4r", width = 0, height = 0)

### Lepomis cyanellus
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Lepomis_cyanellus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Lepomis cyanellus") +
  theme_classic() +
  theme(legend.position = "none")

### Lepomis gulosus
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Lepomis_gulosus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Lepomis gulosus") +
  theme_classic() +
  theme(legend.position = "none")

### Lepomis humilis
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Lepomis_humilis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Lepomis humilis") +
  theme_classic() +
  theme(legend.position = "none")

### Lepomis macrochirus
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Lepomis_macrochirus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Lepomis macrochirus") +
  theme_classic() +
  theme(legend.position = "none")

### Lepomis megalotis
e<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Lepomis_megalotis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Lepomis megalotis") +
  theme_classic() +
  theme(legend.position = "none")

### Lepomis microlophus
f<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Lepomis_microlophus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Lepomis microlophus") +
  theme_classic() +
  theme(legend.position = "none")

# multiplot
multiplot(a,b,c,d,e,f, cols=2)
dev.off()



#####################################################################
###                                                                 #
###     Micropterus & Pomoxis                                        #
###                                                                 #
####################################################################
library(sp)
pdf("./species_dist_maps/Centrarchidae.pdf", paper = "a4r", width = 0, height = 0)

### Micropterus punctulatus
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Micropterus_punctulatus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Micropterus punctulatus") +
  theme_classic() +
  theme(legend.position = "none")

### Micropterus salmoides
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Micropterus_salmoides), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Micropterus salmoides") +
  theme_classic() +
  theme(legend.position = "none")

### Pomoxis annularis
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Pomoxis_annularis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Pomoxis annularis") +
  theme_classic() +
  theme(legend.position = "none")

### Pomoxis nigromaculatus
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Pomoxis_nigromaculatus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Pomoxis nigromaculatus") +
  theme_classic() +
  theme(legend.position = "none")

# multiplot
multiplot(a,b,c,d, cols=2)
dev.off()



#####################################################################
###                                                                 #
###     Notropis -- minnows                                        #
###                                                                 #
####################################################################
library(sp)
pdf("./species_dist_maps/Notropis.pdf", paper = "a4r", width = 0, height = 0)

### Notropis athernoides
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notropis_atherinoides), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notropis atherinoides") +
  theme_classic() +
  theme(legend.position = "none")

### LNotropis atrocaudalis
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notropis_atrocaudalis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notropis atrocaudalis") +
  theme_classic() +
  theme(legend.position = "none")

### Notropis boops
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notropis_boops), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notropis boops") +
  theme_classic() +
  theme(legend.position = "none")

### Notropis buchanani
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notropis_buchanani), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notropis buchanani") +
  theme_classic() +
  theme(legend.position = "none")

### Notropis ortenburgeri
e<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notropis_ortenburgeri), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notropis ortenburgeri") +
  theme_classic() +
  theme(legend.position = "none")

### Notropis stramineus
f<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notropis_stramineus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notropis stramineus") +
  theme_classic() +
  theme(legend.position = "none")

### Notropis suttkusi
g<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notropis_suttkusi), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notropis suttkusi") +
  theme_classic() +
  theme(legend.position = "none")

### Notropis volucellus
h<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notropis_volucellus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notropis volucellus") +
  theme_classic() +
  theme(legend.position = "none")



# multiplot
multiplot(a,b,c,d,e,f,g,h, cols=2)
dev.off()





#####################################################################
###                                                                 #
###     MINNOWS I                                                   #
###                                                                 #
####################################################################
library(sp)
pdf("./species_dist_maps/Minnows1.pdf", paper = "a4r", width = 0, height = 0)

### Campostoma anomalum
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Campostoma_anomalum), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Campostoma anomalum") +
  theme_classic() +
  theme(legend.position = "none")

### Campostoma spadiceum
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Campostoma_spadiceum), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Campostoma spadiceum") +
  theme_classic() +
  theme(legend.position = "none")

### Cyprinella lutrensis
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Cyprinella_lutrensis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Cyprinella lutrensis") +
  theme_classic() +
  theme(legend.position = "none")

### Cyprinella venusta
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Cyprinella_venusta), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Cyprinella venusta") +
  theme_classic() +
  theme(legend.position = "none")

### Cyprinella whipplei
e<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Cyprinella_whipplei), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Cyprinella whipplei") +
  theme_classic() +
  theme(legend.position = "none")

### Luxilus chrysocephalus 
f<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Luxilus_chrysocephalus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Luxilus chrysocephalus") +
  theme_classic() +
  theme(legend.position = "none")

### Lythrurus snelsoni
g<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Lythrurus_snelsoni), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Lythrurus snelsoni") +
  theme_classic() +
  theme(legend.position = "none")

### Lythrurus umbratilis
h<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Lythrurus_umbratilis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Lythrurus umbratilis") +
  theme_classic() +
  theme(legend.position = "none")



# multiplot
multiplot(a,b,c,d,e,f,g,h, cols=2)
dev.off()





#####################################################################
###                                                                 #
###     MINNOWS II                                                   #
###                                                                 #
####################################################################
library(sp)
pdf("./species_dist_maps/Minnows2.pdf", paper = "a4r", width = 0, height = 0)

### Notemigonus crysoleucas
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Notemigonus_crysoleucas), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Notemigonus crysoleucas") +
  theme_classic() +
  theme(legend.position = "none")

### Phenacobius mirabilis
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Phenacobius_mirabilis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Phenacobius mirabilis") +
  theme_classic() +
  theme(legend.position = "none")

### Pimephales notatus
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Pimephales_notatus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Pimephales notatus") +
  theme_classic() +
  theme(legend.position = "none")

### Pimephales promelas
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Pimephales_promelas), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Pimephales promelas") +
  theme_classic() +
  theme(legend.position = "none")

### Pimephales vigilax
e<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Pimephales_vigilax), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Pimephales vigilax") +
  theme_classic() +
  theme(legend.position = "none")

### Semotilus atromaculatus
f<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Semotilus_atromaculatus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Semotilus atromaculatus") +
  theme_classic() +
  theme(legend.position = "none")


# multiplot
multiplot(a,b,c,d,e,f, cols=2)
dev.off()



####################################################################
###                                                                 #
###     Darters                                                    #
###                                                                 #
####################################################################
library(sp)
pdf("./species_dist_maps/Darters.pdf", paper = "a4r", width = 0, height = 0)

### Etheostoma chlorosoma
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Etheostoma_chlorosoma), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Etheostoma chlorosoma") +
  theme_classic() +
  theme(legend.position = "none")

### Etheostoma gracile
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Etheostoma_gracile), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Etheostoma gracile") +
  theme_classic() +
  theme(legend.position = "none")

### Etheostoma nigrum
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Etheostoma_nigrum), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Etheostoma nigrum") +
  theme_classic() +
  theme(legend.position = "none")

### Etheostoma radiosum
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Etheostoma_radiosum), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Etheostoma radiosum") +
  theme_classic() +
  theme(legend.position = "none")

### Etheostoma spectabile
e<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Etheostoma_spectabile), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Etheostoma spectabile") +
  theme_classic() +
  theme(legend.position = "none")

### Percina caprodes
f<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Percina_caprodes), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Percina caprodes") +
  theme_classic() +
  theme(legend.position = "none")

### Percina copelandi
g<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Percina_copelandi), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Percina copelandi") +
  theme_classic() +
  theme(legend.position = "none")

### Percina phoxocephala
h<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Percina_phoxocephala), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Percina phoxocephala") +
  theme_classic() +
  theme(legend.position = "none")

### Percina sciera
i<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Percina_sciera), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Percina sciera") +
  theme_classic() +
  theme(legend.position = "none")



# multiplot
multiplot(a,b,c,d,e,f,g,h,i, cols=3)
dev.off()



####################################################################
###                                                                 #
###     ODDBALLS                                                    #
###                                                                 #
####################################################################
library(sp)
pdf("./species_dist_maps/Oddballs.pdf", paper = "a4r", width = 0, height = 0)

### Aphredoderus sayanus
a<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Aphredoderus_sayanus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Aphredoderus sayanus") +
  theme_classic() +
  theme(legend.position = "none")

### Dorosoma cepedianum
b<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Dorosoma_cepedianum), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Dorosoma cepedianum") +
  theme_classic() +
  theme(legend.position = "none")

### Esox americanus
c<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Esox_americanus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Esox americanus") +
  theme_classic() +
  theme(legend.position = "none")

### Fundulus notatus
d<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Fundulus_notatus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Fundulus notatus") +
  theme_classic() +
  theme(legend.position = "none")

### Fundulus olivaceuous
e<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Fundulus_olivaceous), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Fundulus olivaceous") +
  theme_classic() +
  theme(legend.position = "none")

### Gambusia affins
f<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Gambusia_affinis), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Gambusia affinis") +
  theme_classic() +
  theme(legend.position = "none")

### Labidesthes sicculus
g<-ggplot() +
  geom_sf(data = hydroRiv, col="black", size=0.01) +
  geom_sf(data = join, aes(pch=21, fill=Labidesthes_sicculus), size=1) +
  scale_shape_identity() + 
  scale_fill_gradient(low = "white",
                      high = "red",
                      aesthetics = "fill")+
  ggtitle("Labidesthes sicculus") +
  theme_classic() +
  theme(legend.position = "none")


# multiplot
multiplot(a,b,c,d,e,f,g, cols=2)
dev.off()


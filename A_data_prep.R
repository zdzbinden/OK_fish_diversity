### BEGINNING ---------------------------------------------------
### Comparing Multiple Diversities Across Drainages
### ZD ZBINDEN                        
### *** This script is for preparing the data for analyses ***
### use document outline -->>>>> for easy viewing
### set directory ---------------------------------------------------
setwd("C:/Users/zbind/Desktop/shelf/Projects/Drainage_Compare_SWG/R_workdir_v2/")


### functions ---------------------------------------------------
### function: multiplot
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

### function: select.envfit
select.envfit<-function(fit, p.select){ #needs two sorts of input: fit= result of envfit, r.select= numeric, correlation minimum threshold
  for (i in 1:length(fit$vectors$pvals)) { #run for-loop through the entire length of the column r in object fit$vectors$r starting at i=1
    if (fit$vectors$pvals[i]>p.select) { #Check wether r<r.select, i.e. if the correlation is weaker than the threshold value. Change this Parameter for r-based selection
      fit$vectors$arrows[i,]=NA #If the above statement is TRUE, i.e. r is smaller than r.select, then the coordinates of the vectors are set to NA, so they cannot be displayed
      i=i+1 #increase the running parameter i from 1 to 2, i.e. check the next value in the column until every value has been checked
    } #close if-loop
  } #close for-loop
  return(fit) #return fit as the result of the function
} #close the function


### 1 | ENVIRONMENTAL DATA ---------------------------------------------------
### _1.1 | Organize data --------------------------------------------------- 
library(sf)
streams <- st_read("./data/hydroRivers/streams.shp", stringsAsFactors=TRUE)
plot(streams)
sites <- st_read("./data/sites_shapefiles/sites.shp", stringsAsFactors=TRUE)
plot(sites)
### combined in situ ENV variables with hydroATLAS variables
env.join<-st_join(sites, streams, st_nearest_feature,left=TRUE) 
### can plot any of the ENV vars
plot(env.join[,10:13])
env.join <-st_set_geometry(env.join, NULL)
### organizing data
site_info <- env.join[,c(1:6,30:35,40:42)]
env.rivs <- env.join[,-c(1:6,30:32,34:35,40:43,320:324)]
row.names(env.rivs)<- site_info$study_no
### remove invariant columns
env.rivs <-env.rivs[apply(env.rivs, 2, sd) > 0]
### divide env vars into major classes 
env.instream.vars <- env.rivs[,c(1:23)]
env.hydrophysio.vars <- env.rivs[,c(24:54)]
env.climate.vars <- env.rivs[,c(55:146)]
env.landcover.vars <- env.rivs[,c(147:195)]
env.geology.vars <- env.rivs[,c(196:222)]
env.anthropogenic.vars <- env.rivs[,c(222:235)]

### _1.2 | Principal Components Analysis --------------------------------------------------- 
### ___1.2.1 | instream ---------------------------------------------------                  
### change names to make clearer
names(env.instream.vars)[names(env.instream.vars) == "width"] <- "max.stream.width"
names(env.instream.vars)[names(env.instream.vars) == "depth"] <- "max.stream.depth"
names(env.instream.vars)[names(env.instream.vars) == "pool"] <- "pool.percent"
names(env.instream.vars)[names(env.instream.vars) == "riffle"] <- "riffle.percent"
names(env.instream.vars)[names(env.instream.vars) == "channel"] <- "channel.percent"
names(env.instream.vars)[names(env.instream.vars) == "backwater"] <- "backwater.percent"
names(env.instream.vars)[names(env.instream.vars) == "oxygen"] <- "dissolved.oxygen"
names(env.instream.vars)[names(env.instream.vars) == "current"] <- "current.speed"
names(env.instream.vars)[names(env.instream.vars) == "mud"] <- "mud.percent"
names(env.instream.vars)[names(env.instream.vars) == "sand"] <- "sand.percent"
names(env.instream.vars)[names(env.instream.vars) == "gravel"] <- "gravel.percent"
names(env.instream.vars)[names(env.instream.vars) == "cobble"] <- "cobble.percent"
names(env.instream.vars)[names(env.instream.vars) == "bedrock"] <- "bedrock.percent"
names(env.instream.vars)[names(env.instream.vars) == "canopy_cov"] <- "canopy.cover.percent"
names(env.instream.vars)[names(env.instream.vars) == "fil_algae"] <- "algae.presence"
names(env.instream.vars)[names(env.instream.vars) == "macrophytes"] <- "macrophyte.presence"
names(env.instream.vars)[names(env.instream.vars) == "woody_stru"] <- "sub.structure.presence"
names(env.instream.vars)[names(env.instream.vars) == "boulders"] <- "boulders.presence"
names(env.instream.vars)[names(env.instream.vars) == "cwd"] <- "cwd.presence"
names(env.instream.vars)[names(env.instream.vars) == "bank_incis"] <- "bank.incision"
names(env.instream.vars)[names(env.instream.vars) == "bank_stabi"] <- "stable.bank"
names(env.instream.vars)[names(env.instream.vars) == "Pasture"] <- "pasture.riparian"
names(env.instream.vars)[names(env.instream.vars) == "Woodland"] <- "woodland.riparian"


######### This argument affects all environmental PCS and determines if loadings will be plotted
loadings <- FALSE
### principal components
vif <- usdm::vifstep(env.instream.vars)
env.instream.vars<- usdm::exclude(env.instream.vars, th=10)
env.instream.vars<-as.data.frame(lapply(env.instream.vars, as.numeric))
nFactors::nBartlett(env.instream.vars, nrow(env.instream.vars))
NPC<-PCDimension::rndLambdaF(env.instream.vars, B=9999)
NPC
NPC<-NPC[2]

pca <- prcomp(env.instream.vars, center = TRUE, scale. = TRUE)
pca.cumulative <-summary(pca)
pca.cumulative <- pca.cumulative$importance[,1:NPC]
write.csv(pca.cumulative, "./outputs/env_analyses/PCA_loadings/instream.cumulative.csv")
pca.loadings <-pca$rotation[,1:NPC]
write.csv(pca.loadings, "./outputs/env_analyses/PCA_loadings/instream.loadings.csv")
env.instream.PCs <- as.data.frame(pca$x)[,1:NPC]
colnames(env.instream.PCs)<-paste("instream", colnames(env.instream.PCs), sep = "_")
write.csv(env.instream.PCs, "./outputs/env_analyses/PCA_loadings/instream.PCs.csv")
### plot
library(ggplot2)
library(ggfortify)
### Plotting and coloring by drainage
pdf("./outputs/env_analyses/PCA_figs/instreamPCA.pdf", paper = "a4r", width = 0, height = 0)
PC_instream <-autoplot(pca, data=site_info, colour ='drainage', 
                          loadings = loadings, loadings.label = loadings, loadings.colour = "black", loadings.label.colour = "black",
                          frame = TRUE)
PC_instream <- PC_instream + scale_colour_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_shape_manual(values=c(21,24,23,22))
PC_instream <- PC_instream + theme_classic() + ggtitle("Instream")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(legend.text = element_text(size=20))+
  theme(legend.title = element_text(size = 20))
PC_instream
dev.off()
PC_instream

### ___1.2.2 | hydro-physio ---------------------------------------------------              
### principal components
vif <- usdm::vifstep(env.hydrophysio.vars)
env.hydrophysio.vars<- usdm::exclude(env.hydrophysio.vars, th=10)
env.hydrophysio.vars<-as.data.frame(lapply(env.hydrophysio.vars, as.numeric))
nFactors::nBartlett(env.hydrophysio.vars, nrow(env.hydrophysio.vars))
NPC<-PCDimension::rndLambdaF(env.hydrophysio.vars, B=9999)
NPC
NPC<-NPC[2]

pca <- prcomp(env.hydrophysio.vars, center = TRUE, scale. = TRUE)
pca.cumulative <-summary(pca)
pca.cumulative <- pca.cumulative$importance[,1:NPC]
write.csv(pca.cumulative, "./outputs/env_analyses/PCA_loadings/hydrophysio.cumulative.csv")
pca.loadings <-pca$rotation[,1:NPC]
write.csv(pca.loadings, "./outputs/env_analyses/PCA_loadings/hydrophysio.loadings.csv")
env.hydrophysio.PCs <- as.data.frame(pca$x)[,1:NPC]
colnames(env.hydrophysio.PCs)<-paste("hydrophysio", colnames(env.hydrophysio.PCs), sep = "_")
write.csv(env.hydrophysio.PCs, "./outputs/env_analyses/PCA_loadings/hydrophysio.PCs.csv")
### plot
### Plotting and coloring by drainage
pdf("./outputs/env_analyses/PCA_figs/hydrophysioPCA.pdf", paper = "a4r", width = 0, height = 0)
PC_hydrophysio <-autoplot(pca, data=site_info, colour ='drainage', 
                      loadings = loadings, loadings.label = loadings, loadings.colour = "black", loadings.label.colour = "black",
                      frame = TRUE)
PC_hydrophysio <- PC_hydrophysio + scale_colour_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_shape_manual(values=c(21,24,23,22))
PC_hydrophysio <- PC_hydrophysio + theme_classic() + ggtitle("Hydro-physio")+ 
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))  
PC_hydrophysio
dev.off()
PC_hydrophysio

### ___1.2.3 | climate --------------------------------------------------- 
## principal components
vif <- usdm::vifstep(env.climate.vars)
env.climate.vars<- usdm::exclude(env.climate.vars, th=10)
env.climate.vars<-as.data.frame(lapply(env.climate.vars, as.numeric))
nFactors::nBartlett(env.climate.vars, nrow(env.climate.vars))
NPC<-PCDimension::rndLambdaF(env.climate.vars, B=9999)
NPC
NPC<-NPC[2]

pca <- prcomp(env.climate.vars, center = TRUE, scale. = TRUE)
pca.cumulative <-summary(pca)
pca.cumulative <- pca.cumulative$importance[,1:NPC]
write.csv(pca.cumulative, "./outputs/env_analyses/PCA_loadings/climate.cumulative.csv")
pca.loadings <-pca$rotation[,1:NPC]
write.csv(pca.loadings, "./outputs/env_analyses/PCA_loadings/climate.loadings.csv")
env.climate.PCs <- as.data.frame(pca$x)[,1:NPC]
colnames(env.climate.PCs)<-paste("climate", colnames(env.climate.PCs), sep = "_")
write.csv(env.climate.PCs, "./outputs/env_analyses/PCA_loadings/climate.PCs.csv")
### plot
### Plotting and coloring by drainage
pdf("./outputs/env_analyses/PCA_figs/climatePCA.pdf", paper = "a4r", width = 0, height = 0)
PC_climate <-autoplot(pca, data=site_info, colour ='drainage', 
                  loadings = loadings, loadings.label = loadings, loadings.colour = "black", loadings.label.colour = "black",
                  frame = TRUE)
PC_climate <- PC_climate + scale_colour_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_shape_manual(values=c(21,24,23,22))
PC_climate <- PC_climate + theme_classic() + ggtitle("Climate")+ 
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))  
PC_climate
dev.off()
PC_climate

### ___1.2.4 | landcover --------------------------------------------------- 
## principal components
vif <- usdm::vifstep(env.landcover.vars)
env.landcover.vars<- usdm::exclude(env.landcover.vars, th=10)
env.landcover.vars<-as.data.frame(lapply(env.landcover.vars, as.numeric))
nFactors::nBartlett(env.landcover.vars, nrow(env.landcover.vars))
NPC<-PCDimension::rndLambdaF(env.landcover.vars, B=9999)
NPC
NPC<-NPC[2]

pca <- prcomp(env.landcover.vars, center = TRUE, scale. = TRUE)
pca.cumulative <-summary(pca)
pca.cumulative <- pca.cumulative$importance[,1:NPC]
write.csv(pca.cumulative, "./outputs/env_analyses/PCA_loadings/landcover.cumulative.csv")
pca.loadings <-pca$rotation[,1:NPC]
write.csv(pca.loadings, "./outputs/env_analyses/PCA_loadings/landcover.loadings.csv")
env.landcover.PCs <- as.data.frame(pca$x)[,1:NPC]
colnames(env.landcover.PCs)<-paste("landcover", colnames(env.landcover.PCs), sep = "_")
write.csv(env.landcover.PCs, "./outputs/env_analyses/PCA_loadings/landcover.PCs.csv")
### plot
### Plotting and coloring by drainage
pdf("./outputs/env_analyses/PCA_figs/landcoverPCA.pdf", paper = "a4r", width = 0, height = 0)
PC_landcover <-autoplot(pca, data=site_info, colour ='drainage', 
                      loadings = loadings, loadings.label = loadings, loadings.colour = "black", loadings.label.colour = "black",
                      frame = TRUE)
PC_landcover <- PC_landcover + scale_colour_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_shape_manual(values=c(21,24,23,22))
PC_landcover <- PC_landcover + theme_classic() + ggtitle("Landcover")+ 
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))  
PC_landcover
dev.off()
PC_landcover

### ___1.2.5 | geology --------------------------------------------------- 
### principal components
vif <- usdm::vifstep(env.geology.vars)
env.geology.vars<- usdm::exclude(env.geology.vars, th=10)
env.geology.vars<-as.data.frame(lapply(env.geology.vars, as.numeric))
nFactors::nBartlett(env.geology.vars, nrow(env.geology.vars))
NPC<-PCDimension::rndLambdaF(env.geology.vars, B=9999)
NPC
NPC<-NPC[2]+1

pca <- prcomp(env.geology.vars, center = TRUE, scale. = TRUE)
pca.cumulative <-summary(pca)
pca.cumulative <- pca.cumulative$importance[,1:NPC]
write.csv(pca.cumulative, "./outputs/env_analyses/PCA_loadings/geology.cumulative.csv")
pca.loadings <-pca$rotation[,1:NPC]
write.csv(pca.loadings, "./outputs/env_analyses/PCA_loadings/geology.loadings.csv")
env.geology.PCs <- as.data.frame(pca$x)[,1:NPC]
colnames(env.geology.PCs)<-paste("geology", colnames(env.geology.PCs), sep = "_")
write.csv(env.geology.PCs, "./outputs/env_analyses/PCA_loadings/geology.PCs.csv")
### plot
### Plotting and coloring by drainage
pdf("./outputs/env_analyses/PCA_figs/geologyPCA.pdf", paper = "a4r", width = 0, height = 0)
PC_geology <-autoplot(pca, data=site_info, colour ='drainage', 
                        loadings = loadings, loadings.label = loadings, loadings.colour = "black", loadings.label.colour = "black",
                        frame = TRUE)
PC_geology <- PC_geology + scale_colour_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_shape_manual(values=c(21,24,23,22))
PC_geology <- PC_geology + theme_classic() + ggtitle("Geology")+ 
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))  
PC_geology
dev.off()
PC_geology

### ___1.2.6 | anthropogenic --------------------------------------------------- 
### principal components
vif <- usdm::vifstep(env.anthropogenic.vars)
env.anthropogenic.vars<- usdm::exclude(env.anthropogenic.vars, th=10)
env.anthropogenic.vars<-as.data.frame(lapply(env.anthropogenic.vars, as.numeric))
nFactors::nBartlett(env.anthropogenic.vars, nrow(env.anthropogenic.vars))
NPC<-PCDimension::rndLambdaF(env.anthropogenic.vars, B=9999)
NPC
NPC<-NPC[2]

pca <- prcomp(env.anthropogenic.vars, center = TRUE, scale. = TRUE)
pca.cumulative <-summary(pca)
pca.cumulative <- pca.cumulative$importance[,1:NPC]
write.csv(pca.cumulative, "./outputs/env_analyses/PCA_loadings/anthropogenic.cumulative.csv")
pca.loadings <-pca$rotation[,1:NPC]
write.csv(pca.loadings, "./outputs/env_analyses/PCA_loadings/anthropogenic.loadings.csv")
env.anthropogenic.PCs <- as.data.frame(pca$x)[,1:NPC]
colnames(env.anthropogenic.PCs)<-paste("anthropogenic", colnames(env.anthropogenic.PCs), sep = "_")
write.csv(env.anthropogenic.PCs, "./outputs/env_analyses/PCA_loadings/anthropogenic.PCs.csv")
### plot
### Plotting and coloring by drainage
pdf("./outputs/env_analyses/PCA_figs/anthropogenicPCA.pdf", paper = "a4r", width = 0, height = 0)
PC_anthropogenic <-autoplot(pca, data=site_info, colour ='drainage', 
                      loadings = loadings, loadings.label = loadings, loadings.colour = "black", loadings.label.colour = "black",
                      frame = TRUE)
PC_anthropogenic <- PC_anthropogenic + scale_colour_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"), name="Drainage")+
  scale_shape_manual(values=c(21,24,23,22))
PC_anthropogenic <- PC_anthropogenic + theme_classic() + ggtitle("Anthropogenic")+ 
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14)) 
PC_anthropogenic
dev.off()
PC_anthropogenic

### Make an 8 panel composite plot of the "by drainage" PCAs  ###
legend <- ggpubr::get_legend(PC_instream)
legend <- ggpubr::as_ggplot(legend)
PC_instream <- PC_instream + theme(legend.position = "none")
PC_hydrophysio <- PC_hydrophysio + theme(legend.position = "none")
PC_climate <- PC_climate + theme(legend.position = "none")
PC_landcover <- PC_landcover + theme(legend.position = "none")
PC_geology <- PC_geology + theme(legend.position = "none")
PC_anthropogenic <- PC_anthropogenic + theme(legend.position = "none")

pdf("./outputs/env_analyses/PCA_figs/composite.pdf", paper = "USr", width = 0, height = 0)
multiplot(PC_instream, PC_hydrophysio, PC_climate, PC_landcover, PC_geology, PC_anthropogenic,legend, cols=3)
dev.off()
multiplot(PC_instream, PC_hydrophysio, PC_climate, PC_landcover, PC_geology, PC_anthropogenic,legend, cols=3)


### cleaning up r object workspace
rm(legend,PC_anthropogenic,PC_climate,PC_geology,PC_hydrophysio,PC_instream,PC_landcover,pca,
   pca.cumulative,pca.loadings,loadings)

### 2 | FISH TAXON DATA --------------------------------------------------- 
### _2.1 | Organize data --------------------------------------------------- 
### read in fish abundance data file
taxon.div <- read.csv("./data/fish_data_basin_comp.csv", row.names = 1)
taxon.div <- t(taxon.div)
taxon.div <- as.data.frame(taxon.div)

### Remove any outlier sites here:
### As above...
taxon.div <- taxon.div[-c(11,137),]

### _2.2 | NMDS --------------------------------------------------- 
library(vegan)
### calculate NMDS using Bray-Curtis index on incidence data ~ Sorensen 
#taxon.div<-decostand(taxon.div, method = "hellinger")
nmds <- metaMDS(taxon.div, distance = "bray", k=3, try = 1000, trymax = 1000, 
                autotransform = FALSE, center=TRUE, expand = TRUE)
stressplot(nmds)
nmds

### _2.3 | Plotting --------------------------------------------------- 
Drainage <- site_info$drainage
data.scores <- as.data.frame(scores(nmds))  
data.scores$site <- rownames(data.scores)
data.scores$grp <- Drainage
head(data.scores)
species.scores <- as.data.frame(scores(nmds, "species")) 
species.scores$species <- rownames(species.scores)
head(species.scores)

### ___2.3.1 | axes 1  & 2 --------------------------------------------------- 
library(ggplot2)
library(dplyr)
library(ggrepel)
#### subset to only important species
### can change SP threshold 
SP <- 1.0
temp.scores <-species.scores %>% filter_at(vars(1:2), any_vars(abs(.) > SP))
# Plot nmds 1 & 2
library(ggplot2)
### calculate hull values
clear <- data.scores[data.scores$grp == "Clear", ][chull(data.scores[data.scores$grp == "Clear", c("NMDS1", "NMDS2")]), ]  
muddy <- data.scores[data.scores$grp == "Muddy", ][chull(data.scores[data.scores$grp == "Muddy", c("NMDS1", "NMDS2")]), ] 
kiamichi <- data.scores[data.scores$grp == "Kiamichi", ][chull(data.scores[data.scores$grp == "Kiamichi", c("NMDS1", "NMDS2")]), ]  
little <- data.scores[data.scores$grp == "Little", ][chull(data.scores[data.scores$grp == "Little", c("NMDS1", "NMDS2")]), ]  
hull.data <- rbind(clear,muddy,kiamichi,little)
hull.data

pdf("./outputs/div.plots/NMDS12.pdf", paper = "US", width = 0, height = 0)
gg12 <- ggplot() + 
  #geom_point(data=temp.scores, aes(x=NMDS1,y=NMDS2))+
  #geom_text_repel(data=temp.scores,inherit.aes = FALSE, aes(x=NMDS1,y=NMDS2,label=species),size=3)+  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Drainage,colour=Drainage),size=1.5) + # add the point markers
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.10, show.legend = FALSE) +
  scale_shape_manual(values=c(16,17,18,15)) +
  scale_colour_manual(values=c("red3","gold3","blue3","darkgreen")) +
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"))+
  coord_equal() +
  theme_classic()+ theme(legend.position = "none" ) + xlim(-2,2) + ylim(-2,2)+
  theme(text = element_text(size=12))+
  theme(axis.text = element_text(size=12))+
  theme(plot.title = element_text(size=12))+
  theme(legend.text = element_text(size=12))
gg12
dev.off()
gg12

### ___2.3.2 | axes 1  & 3 --------------------------------------------------- 
#### calculate hull values
#### subset to only important species
temp.scores <-species.scores %>% filter_at(vars(1,3), any_vars(abs(.) > SP))

clear <- data.scores[data.scores$grp == "Clear", ][chull(data.scores[data.scores$grp == "Clear", c("NMDS1", "NMDS3")]), ]  
muddy <- data.scores[data.scores$grp == "Muddy", ][chull(data.scores[data.scores$grp == "Muddy", c("NMDS1", "NMDS3")]), ] 
kiamichi <- data.scores[data.scores$grp == "Kiamichi", ][chull(data.scores[data.scores$grp == "Kiamichi", c("NMDS1", "NMDS3")]), ]  
little <- data.scores[data.scores$grp == "Little", ][chull(data.scores[data.scores$grp == "Little", c("NMDS1", "NMDS3")]), ]  
hull.data <- rbind(clear,muddy,kiamichi,little)
hull.data

gg13 <- ggplot() + 
  #geom_point(data=temp.scores, aes(x=NMDS1,y=NMDS3))+
  #geom_text_repel(data=temp.scores,inherit.aes = FALSE, aes(x=NMDS1,y=NMDS3,label=species),size=3)+
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS3,shape=Drainage,colour=Drainage),size=1.5) + # add the point markers
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS3,fill=grp,group=grp),alpha=0.10, show.legend = FALSE) +
  scale_shape_manual(values=c(16,17,18,15)) +
  scale_colour_manual(values=c("red3","gold3","blue3","darkgreen")) +
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"))+
  coord_equal() +
  theme_classic()+ theme(legend.position = "none" ) + xlim(-2,2) + ylim(-2,2)+
  theme(text = element_text(size=12))+
  theme(axis.text = element_text(size=12))+
  theme(plot.title = element_text(size=12))+
  theme(legend.text = element_text(size=12))
gg13
pdf("./outputs/div.plots/NMDS13.pdf", paper = "US", width = 0, height = 0)
gg13
dev.off()
### ___2.3.3 | axes 2  & 3 --------------------------------------------------- 
### calculate hull values
#### subset to only important species
temp.scores <-species.scores %>% filter_at(vars(2,3), any_vars(abs(.) > SP))
clear <- data.scores[data.scores$grp == "Clear", ][chull(data.scores[data.scores$grp == "Clear", c("NMDS2", "NMDS3")]), ]  
muddy <- data.scores[data.scores$grp == "Muddy", ][chull(data.scores[data.scores$grp == "Muddy", c("NMDS2", "NMDS3")]), ] 
kiamichi <- data.scores[data.scores$grp == "Kiamichi", ][chull(data.scores[data.scores$grp == "Kiamichi", c("NMDS2", "NMDS3")]), ]  
little <- data.scores[data.scores$grp == "Little", ][chull(data.scores[data.scores$grp == "Little", c("NMDS2", "NMDS3")]), ]  
hull.data <- rbind(clear,muddy,kiamichi,little)
hull.data

gg23 <- ggplot() + 
  #geom_point(data=temp.scores, aes(x=NMDS2,y=NMDS3))+
  #geom_text_repel(data=temp.scores,inherit.aes = FALSE, aes(x=NMDS2,y=NMDS3,label=species),size=3)+
  geom_point(data=data.scores,aes(x=NMDS2,y=NMDS3,shape=Drainage,colour=Drainage),size=1.5) + # add the point markers
  geom_polygon(data=hull.data,aes(x=NMDS2,y=NMDS3,fill=grp,group=grp),alpha=0.10, show.legend = FALSE) +
  scale_shape_manual(values=c(16,17,18,15)) +
  scale_colour_manual(values=c("red3","gold3","blue3","darkgreen")) +
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"))+
  coord_equal() +
  theme_classic()+ theme(legend.position = "right" ) + xlim(-2,2) + ylim(-2,2)+
  theme(text = element_text(size=12))+
  theme(axis.text = element_text(size=12))+
  theme(plot.title = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))
 gg23
 pdf("./outputs/div.plots/NMDS23.pdf", paper = "US", width = 0, height = 0)
 gg23
 dev.off()

### ___2.3.4 | multiplot --------------------------------------------------- 
#### multiplot **must have function installed** see section above
legend <- ggpubr::get_legend(gg23)
legend <- ggpubr::as_ggplot(legend)
gg23 <- gg23 + theme(legend.position = "none")
pdf("./outputs/div.plots/CompositeNMDS.pdf", paper = "USr", width = 0, height = 0)
multiplot(gg12,gg13,gg23,legend, cols=3)
dev.off()
multiplot(gg12,gg13,gg23,legend, cols=3)

### cleaning up r object workspace
rm(gg12,gg13,gg23,hull.data,data.scores, kiamichi,little,muddy,clear,nmds,species.scores,legend, Drainage,SP, temp.scores)

### 3 | FUNCTIONAL TRAIT DATA --------------------------------------------------- 
### _3.1 | Organize data --------------------------------------------------- 
functional <- read.csv("./data/traits_data_basin_comp.csv", row.names = 1)
### pick out family ID
FAMILY <- functional[,83]
### Remove variables not of interest
functional <- functional[,-c(14:25,27:51,54:83)]
### remove invariant columns
functional <-functional[apply(functional, 2, sd) > 0]
### center and scale variables
functional <-BBmisc::normalize(functional, method="standardize")
#summary(functional)

### rename variables to be more intuitive
names(functional)[names(functional) == "BENTHIC"] <- "benthic.feed"
names(functional)[names(functional) == "SURWCOL"] <- "surface.feed"
names(functional)[names(functional) == "ALGPHYTO"] <- "algae.feed"
names(functional)[names(functional) == "MACVASCU"] <- "macrophyte.feed"
names(functional)[names(functional) == "DETRITUS"] <- "detritus.feed"
names(functional)[names(functional) == "FSHCRCRB"] <- "piscivore.feed"
names(functional)[names(functional) == "EGGS"] <- "eggs.feed"
names(functional)[names(functional) == "MAXTL"] <- "max.length"
names(functional)[names(functional) == "MATUAGE"] <- "maturation.age"
names(functional)[names(functional) == "LONGEVITY"] <- "longevity"
names(functional)[names(functional) == "FECUNDITY"] <- "fecundity"
names(functional)[names(functional) == "SERIAL"] <- "serial.spawn"
names(functional)[names(functional) == "SEASON"] <- "length.spawn"
names(functional)[names(functional) == "MINTEMP"] <- "min.temp"
names(functional)[names(functional) == "MAXTEMP"] <- "min.temp"

### _3.2 | PCOA --------------------------------------------------- 
library(vegan)
###used Gower distance because mix of binary and continuous variables
fundist <- vegdist(functional, method = "gower", zerodist="add")
###write.csv(as.matrix(fundist),"./functional.plots/functional.distance.csv")
ade4::is.euclid(fundist)
PCOA <-capscale(fundist ~ 1, distance = NULL, data=NULL, sqrt.dist = TRUE ,add = FALSE)
summary(PCOA)
functional.pcoa <-scores(PCOA, choices = c(1:4), display = "sites")
colnames(functional.pcoa) <- c("PCoA1","PCoA2","PCoA3","PCoA4")
###write.csv(functional.eigen, "./functional.plots/functional_eigen.csv")
stressplot(PCOA)

### _3.3 | NMDS --------------------------------------------------- 
library(vegan)
### used Gower distance because mix of binary and continuous variables
functional.nmds <- metaMDS(functional, distance="gower", zerodist="add", 
                           k=4, try = 1000, trymax = 1000, autotransform = FALSE, noshare = FALSE,
                           center=TRUE, pc=TRUE, halfchange=FALSE)

stressplot(functional.nmds)
functional.nmds <- as.data.frame(functional.nmds$points)
colnames(functional.nmds) <- c("NMDS1","NMDS2","NMDS3","NMDS4")


### _3.4 | Plotting --------------------------------------------------- 
library(vegan)
library(ggplot2)
### ___3.4.1 | axes 1 & 2 --------------------------------------------------- 
#### display traits that significantly correlate with axes
traits.fit <- envfit(functional.nmds[,c(1,2)],functional, permutations = 999)
traits.fit <- select.envfit(traits.fit, 0.01)
en_coord_cont = as.data.frame(scores(traits.fit,"vectors")) #* ordiArrowMul(traits.fit)

colors <-c("aquamarine","cyan","orangered","blue4","gold","cyan4","lightgreen","orange","darkred","blueviolet")
gg12 <- ggplot() +
  geom_point(data=functional.nmds,aes(x=NMDS1,y=NMDS2,colour=FAMILY),size=2, alpha=0.8) + 
  scale_color_manual(values=colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1*0.75, yend = NMDS2*0.75), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey70") +
  geom_text(data = en_coord_cont, aes(x = NMDS1*0.75, y = NMDS2*0.75), colour = "black", size=4,
            label = row.names(en_coord_cont)) + 
  theme_classic() + xlab("NMDS1") + ylab("NMDS2")+ #xlim(-0.7,0.7) + ylim(-0.7,0.7)+
  theme(text = element_text(size=12))+
  theme(axis.text = element_text(size=12))+
  theme(plot.title = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))
gg12


### ___3.4.2 | axes 1 & 3 --------------------------------------------------- 
###display traits that significantly correlate with axes
traits.fit <- envfit(functional.nmds[,c(1,3)],functional, permutations = 999)
traits.fit <- select.envfit(traits.fit, 0.01)
en_coord_cont = as.data.frame(scores(traits.fit,"vectors")) #* ordiArrowMul(traits.fit)

gg13 <- ggplot() +
  geom_point(data=functional.nmds,aes(x=NMDS1,y=NMDS3,colour=FAMILY),size=2, alpha=0.8) + 
  scale_color_manual(values=colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1*0.75, yend = NMDS3*0.75), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey70") +
  geom_text(data = en_coord_cont, aes(x = NMDS1*0.75, y = NMDS3*0.75), colour = "black", size=4,
            label = row.names(en_coord_cont)) + 
  theme_classic() + xlab("NMDS1") + ylab("NMDS3")+ #xlim(-0.7,0.7) + ylim(-0.7,0.7)+
  theme(text = element_text(size=12))+
  theme(axis.text = element_text(size=12))+
  theme(plot.title = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))
gg13

### ___3.4.3 | axes 2 & 3 --------------------------------------------------- 
### display traits that significantly correlate with axes
traits.fit <- envfit(functional.nmds[,c(2,3)],functional, permutations = 999)
traits.fit <- select.envfit(traits.fit, 0.01)
en_coord_cont = as.data.frame(scores(traits.fit,"vectors")) #* ordiArrowMul(traits.fit)

gg23 <- ggplot() +
  geom_point(data=functional.nmds,aes(x=NMDS2,y=NMDS3,colour=FAMILY),size=2, alpha=0.8) + 
  scale_color_manual(values=colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS2*0.75, yend = NMDS3*0.75), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey70") +
  geom_text(data = en_coord_cont, aes(x = NMDS2*0.75, y = NMDS3*0.75), colour = "black", size=4,
            label = row.names(en_coord_cont)) + 
  theme_classic() + xlab("NMDS2") + ylab("NMDS3")+ #xlim(-0.7,0.7) + ylim(-0.7,0.7)+
  theme(text = element_text(size=12))+
  theme(axis.text = element_text(size=12))+
  theme(plot.title = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))
gg23

### ___3.4.4 | multiplot --------------------------------------------------- 
legend <- ggpubr::get_legend(gg12)
legend <- ggpubr::as_ggplot(legend)
gg12 <- gg12 + theme(legend.position = "none")
gg13 <- gg13 + theme(legend.position = "none")
gg23 <- gg23 + theme(legend.position = "none")

pdf("./outputs/functional.plots/FUN.CompositeNMDS.pdf", paper = "USr")
multiplot(gg12,gg13,gg23,legend, cols=2)
dev.off()
multiplot(gg12,gg13,gg23,legend, cols=2)

### cleanup workspace
rm(traits.fit,colors,FAMILY,fundist,gg12,gg13,gg23, en_coord_cont, PCOA,legend)

### 4 | PHYLOGENETIC DATA --------------------------------------------------- 
### _4.1 | Get fishtree --------------------------------------------------- 
library(fishtree)
fish.phylo <-fishtree_phylogeny(
c("Ameiurus melas",        "Ameiurus natalis",        "Aphredoderus sayanus",         "Campostoma anomalum",    
"Campostoma oligolepis",   "Cyprinella lutrensis",    "Cyprinella venusta",           "Cyprinella whipplei",    
"Dorosoma cepedianum",     "Erimyzon oblongus",       "Esox americanus vermiculatus", "Etheostoma chlorosomum",  
"Etheostoma gracile",      "Etheostoma nigrum",       "Etheostoma radiosum",          "Etheostoma spectabile",  
"Fundulus notatus",        "Fundulus olivaceus",      "Gambusia affinis",             "Ictalurus punctatus",    
"Labidesthes sicculus",    "Lepomis cyanellus",       "Lepomis gulosus",              "Lepomis humilis",        
"Lepomis macrochirus",     "Lepomis megalotis",       "Lepomis microlophus",          "Luxilus chrysocephalus", 
"Lythrurus snelsoni",      "Lythrurus umbratilis",    "Micropterus punctulatus",      "Micropterus salmoides",  
"Minytrema melanops",      "Moxostoma valenciennesi", "Moxostoma erythrurum",         "Notemigonus crysoleucas",
"Notropis atherinoides",   "Notropis heterolepis",    "Notropis boops",               "Notropis buchanani",     
"Notropis nubilus",        "Notropis stramineus",     "Notropis suttkusi",            "Notropis volucellus",    
"Noturus exilis",          "Noturus gyrinus",         "Noturus nocturnus",            "Percina caprodes",       
"Percina copelandi",       "Percina phoxocephala",    "Percina sciera",               "Phenacobius mirabilis",  
"Pimephales notatus",      "Pimephales promelas",     "Pimephales vigilax",           "Pomoxis annularis",      
"Pomoxis nigromaculatus",  "Semotilus atromaculatus"),
type = "phylogram")
plot(fish.phylo)

### _4.2 | deal with tips ---------------------------------------------------
### OK 9 species not recognized in the FishTree data... 2 were due to nameing mismatches .. 7 just simply are
### missing .. Two mismatches were found and will be made consistent .. For the other 7 species, sister species
###Will be used instead, except for N. ortenburgeri for which systematics is not known, will make sister to 
### N. boops. 
tips.sub <-fish.phylo$tip.label
tips.replace <-fish.phylo$tip.label
tips.change <- cbind(tips.sub,tips.replace)
tips.change <- as.data.frame(tips.change, row.names = tips.replace)
### naming discrepencies but systematics OK
row.names(tips.change)[row.names(tips.change) == "Etheostoma_chlorosomum"] <- "Etheostoma_chlorosoma"
row.names(tips.change)[row.names(tips.change) == "Esox_americanus_vermiculatus"] <- "Esox_americanus"
### original species not included in FishTree so chose siter species to represent them
row.names(tips.change)[row.names(tips.change) == "Moxostoma_valenciennesi"] <- "Moxostoma_duquesnei"
row.names(tips.change)[row.names(tips.change) == "Campostoma_oligolepis"] <- "Campostoma_spadiceum"
row.names(tips.change)[row.names(tips.change) == "Notropis_heterolepis"] <- "Notropis_atrocaudalis"
row.names(tips.change)[row.names(tips.change) == "Notropis_nubilus"] <- "Notropis_ortenburgeri"
row.names(tips.change)[row.names(tips.change) == "Campostoma_anomalum"] <- "Campostoma_plumbeum"
row.names(tips.change)[row.names(tips.change) == "Erimyzon_oblongus"] <- "Erimyzon_claviformis"
row.names(tips.change)[row.names(tips.change) == "Etheostoma_spectabile"] <- "Etheostoma_pulchellum"
tips.change[,2] <- row.names(tips.change)
row.names(tips.change) <- c(1:58)
library(phylotools)
library(phytools)
fish.phylo <-sub.taxa.label(fish.phylo, tips.change)

### _4.3 | plot tree ---------------------------------------------------
pdf("./outputs/trees/species_phylogeny.pdf", paper = "USr", width = 0, height = 0)
plotTree(fish.phylo)
dev.off()
plotTree(fish.phylo)
### cleanup workspace
rm(tips.replace, tips.change, tips.sub)
detach(package:phylotools)
detach(package:phytools)


### 5 | SPATIAL DATA ---------------------------------------------------
library(raster)
sites.shp <- shapefile("data/sites_shapefiles/sites.shp")
library(adespatial)
library(ade4)
library(adegraphics)
### _5.1 | MEM overland dist ---------------------------------------------------
land_dist <- vegdist(sites.shp@coords, method = "euclidean")
SP_MEM_LAND <- dbmem(land_dist, thresh = NULL, MEM.autocor = "positive", store.listw = TRUE, silent = FALSE)
colnames(SP_MEM_LAND)<-paste("LAND", colnames(SP_MEM_LAND), sep = "_")
plot(SP_MEM_LAND, SpORcoords=sites.shp)
vif<- usdm::vifstep(as.data.frame(SP_MEM_LAND), th=10)
SP_MEM_LAND <- usdm::exclude(as.data.frame(SP_MEM_LAND), vif)

### _5.2 | MEM hydro network ---------------------------------------------------
######### MEM on HYDRO NETWORK DISTANCE
### reading HydroRIVERS shapefile into R (modified this file in ArcGIS)
library(riverdist)
library(rgdal)
riv.network <- line2network(path = "./data/hydroRivers/streams.shp", layer = "streams", tolerance = 250)
plot(riv.network)

### INTERACTIVE 
### Dissolve = Y; Insert Vet = Y; min dist = 500; seg mouth = 921; vet mouth = 94; Accept mouth = Y; 
### Remove add seg = N; Build Seg routes = Y
riv.network <-cleanup(riv.network)

localities <- pointshp2segvert("./data/sites_shapefiles", rivers= riv.network)
hist(localities$snapdist)
plot(riv.network)
zoomtoseg(seg = localities$seg, rivers = riv.network)
points(sites.shp, pch=16, col="red")
riverpoints(seg = localities$seg, vert = localities$vert,rivers = riv.network, pch = 15, col = "blue")
hydro_dist<-riverdistancemat(localities$seg, localities$vert, riv.network, algorithm = "segroutes")
row.names(hydro_dist) <- row.names(env.rivs)
colnames(hydro_dist) <- row.names(env.rivs)
write.csv(as.matrix(hydro_dist),"./outputs/hydro_dist.csv")
hydro_dist <- as.dist(hydro_dist)


### Spatial eigenfunction analysis of fluvial distance matrix to generate eigenvectors representing "river space"
SP_MEM_HYDRO <- dbmem(hydro_dist, thresh = NULL, MEM.autocor = "positive", store.listw = TRUE, silent = FALSE)
colnames(SP_MEM_HYDRO)<-paste("HYDRO", colnames(SP_MEM_HYDRO), sep = "_")
plot(SP_MEM_HYDRO, sites.shp)
vif<- usdm::vifstep(as.data.frame(SP_MEM_HYDRO), th=10)
SP_MEM_HYDRO <- usdm::exclude(as.data.frame(SP_MEM_HYDRO), vif)
### _5.3 | MEM upstream net ---------------------------------------------------
##################   FLOW SENSITIVE HYDRO DIST MEM (RIVERDIST PACKAGE)
######  Using NET
### distances defined as +/- (upriver/downriver)
### net = TRUE --- computes difference between up and down distances (if applicable)
### If route was down river 10 km and up a trib 20 km net upriver distance = 10 km 
ustream_dist<-upstreammat(localities$seg, localities$vert, riv.network, 
                          net=TRUE, flowconnected=FALSE, algorithm="segroutes")
write.csv(as.matrix(ustream_dist),"./outputs/ustream_dist_net.csv")
mantel(hydro_dist, ustream_dist)
row.names(ustream_dist) <- row.names(env.rivs)
colnames(ustream_dist) <- row.names(env.rivs)
ustream_dist <- as.dist(ustream_dist)
SP_MEM_UPSTR_NET <- dbmem(ustream_dist, thresh = NULL, MEM.autocor = "positive", store.listw = TRUE, silent = FALSE)
colnames(SP_MEM_UPSTR_NET)<-paste("UPNET", colnames(SP_MEM_UPSTR_NET), sep = "_")
plot(SP_MEM_UPSTR_NET, sites.shp)
vif<- usdm::vifstep(as.data.frame(SP_MEM_UPSTR_NET), th=10)
SP_MEM_UPSTR_NET <- usdm::exclude(as.data.frame(SP_MEM_UPSTR_NET), vif)

### _5.4 | MEM upstream total ---------------------------------------------------
########### Using TOTAL
### distances defined as +/- (upriver/downriver)
### net = FALSE --- computes total  up and down distances
### If route was down river 10 km and up a trib 20 km net upriver distance = 30 km
ustream_dist<-upstreammat(localities$seg, localities$vert, riv.network, 
                          net=FALSE, flowconnected=FALSE, algorithm="segroutes")
write.csv(as.matrix(ustream_dist),"./outputs/ustream_dist_total.csv")
mantel(hydro_dist, ustream_dist)
row.names(ustream_dist) <- row.names(env.rivs)
colnames(ustream_dist) <- row.names(env.rivs)
ustream_dist <- as.dist(ustream_dist)
SP_MEM_UPSTR_TOTAL<- dbmem(ustream_dist, thresh = NULL, MEM.autocor = "positive", store.listw = TRUE, silent = FALSE)
colnames(SP_MEM_UPSTR_TOTAL)<-paste("UPTOTAL", colnames(SP_MEM_UPSTR_TOTAL), sep = "_")
plot(SP_MEM_UPSTR_TOTAL, sites.shp)
vif<- usdm::vifstep(as.data.frame(SP_MEM_UPSTR_TOTAL), th=10)
SP_MEM_UPSTR_TOTAL <- usdm::exclude(as.data.frame(SP_MEM_UPSTR_TOTAL), vif)

###########################     Asymmetric MEM 
### read in data necessary
### 2-column linked segments, integer ID, coords, and original ID 
link<-read.csv("./data/links.csv", header = FALSE, row.names = NULL)
coords <-link[,c(1,6,5)]
segs <- link[,c(1,2)]
dams <-as.matrix(link[,7]) # inverse of dams present on seg (i.e. 0 = present)

### _5.5 | AEM tail-down unweighted ---------------------------------------------------
### build matrices for AEM
### rot.angle controls direction of the "asymmetric process" from bottom to top of plot area
mat <-aem.build.binary(nb.object = NULL, coords, link=segs, 
                 rm.same.y = FALSE, unit.angle = "degrees", rot.angle = 0, plot.connexions = TRUE)
### AEM analysis
SP_AEM_TDN_UW <- aem(mat, rm.link0 = FALSE)
### tidy up data and order based on original ID order in other dataframes 
SP_AEM_TDN_UW <- as.data.frame(SP_AEM_TDN_UW$vectors)
row.names(SP_AEM_TDN_UW) <- link$V3
SP_AEM_TDN_UW <- SP_AEM_TDN_UW[order(match(row.names(SP_AEM_TDN_UW), as.matrix(row.names(taxon.div))[,1])),]
colnames(SP_AEM_TDN_UW)<-paste("AEM_TD_UW", colnames(SP_AEM_TDN_UW), sep = "_")
### select AEMS with positive spatial structure
library(spdep)
listw<-mat2listw(as.matrix(mat$se.mat))
morans <-moran.randtest(SP_AEM_TDN_UW, listw)
morans.pos <-cbind(morans$names, morans$pvalue)
morans.pos<-morans.pos[morans.pos[,2] < 0.05,]
SP_AEM_TDN_UW<-subset(SP_AEM_TDN_UW, select=morans.pos[,1])
# remove "ghost linkages"
SP_AEM_TDN_UW <- SP_AEM_TDN_UW[-c(139:141),]


### _5.6 | AEM tail-down dam ---------------------------------------------------
### AEM analysis
### weight based on presence of dam on segment. DAM = 0 and NODAM = 1 (based on similarities)
SP_AEM_TDN_DAMW <- aem(mat, weight = dams[,1], rm.link0 = FALSE)
### tidy up data and order based on original ID order in other dataframes 
SP_AEM_TDN_DAMW <- as.data.frame(SP_AEM_TDN_DAMW$vectors)
row.names(SP_AEM_TDN_DAMW) <- link$V3
SP_AEM_TDN_DAMW <- SP_AEM_TDN_DAMW[order(match(row.names(SP_AEM_TDN_DAMW), as.matrix(row.names(taxon.div))[,1])),]
colnames(SP_AEM_TDN_DAMW)<-paste("AEM_TD_DAM", colnames(SP_AEM_TDN_DAMW), sep = "_")
### select AEMS with positive spatial structure
library(spdep)
listw<-mat2listw(as.matrix(mat$se.mat))
morans <-moran.randtest(SP_AEM_TDN_DAMW, listw)
morans.pos <-cbind(morans$names, morans$pvalue)
morans.pos<-morans.pos[morans.pos[,2] < 0.05,]
SP_AEM_TDN_DAMW<-subset(SP_AEM_TDN_DAMW, select=morans.pos[,1])
### remove "ghost" linkages
SP_AEM_TDN_DAMW <- SP_AEM_TDN_DAMW[-c(139:141),]


### _5.7 | AEM tail-down distance ---------------------------------------------------
### AEM analysis
### weight based inverse squared distance
dmat <-as.matrix(dist(coords[,c(2:3)]))
edges <- mat$edges[-1,]
length.edge <- vector(length = nrow(segs))
for (i in 1:nrow(edges))
{
  length.edge[i] <- dmat[edges[i,1], edges[i,2]]
}
weight.vec <- 1-(length.edge/max(length.edge))^2
SP_AEM_TDN_DISTW <- aem(mat, weight = weight.vec, rm.link0 = FALSE)
### tidy up data and order based on original ID order in other dataframes 
SP_AEM_TDN_DISTW <- as.data.frame(SP_AEM_TDN_DISTW$vectors)
row.names(SP_AEM_TDN_DISTW) <- link$V3
SP_AEM_TDN_DISTW <- SP_AEM_TDN_DISTW[order(match(row.names(SP_AEM_TDN_DISTW ), as.matrix(row.names(taxon.div))[,1])),]
colnames(SP_AEM_TDN_DISTW)<-paste("AEM_TD_DIST", colnames(SP_AEM_TDN_DISTW), sep = "_")
### select AEMS with positive spatial structure
library(spdep)
listw<-mat2listw(as.matrix(mat$se.mat))
morans <-moran.randtest(SP_AEM_TDN_DISTW, listw)
morans.pos <-cbind(morans$names, morans$pvalue)
morans.pos<-morans.pos[morans.pos[,2] < 0.05,]
SP_AEM_TDN_DISTW<-subset(SP_AEM_TDN_DISTW, select=morans.pos[,1])
### remove "ghost" linkages
SP_AEM_TDN_DISTW <- SP_AEM_TDN_DISTW[-c(139:141),]


### _5.8 | AEM tail-up unweighted ---------------------------------------------------
#### connect all terminal nodes to a "ghost" site upstream = origin of flow/direction
### without doing this the directional process is only modeled through 1 site with highest latitude
### find terminal sites/segments (i.e. most upstream site)
term.segs<-segs$V1[!segs$V1 %in% segs$V2]
term.segs.mat<- matrix(nrow = (length(term.segs)), ncol = 2)
### add connections between terminal and "ghost" segment
oS <- nrow(segs)+1
term.segs.mat[,1] <- rep(oS,(length(term.segs)))
term.segs.mat[,2] <- term.segs
### add made up location info for "ghost" segment to "coords" object
### select a central location above all sites 
segs <- rbind(segs,term.segs.mat)
coords[oS,1]<- nrow(coords)+1
coords[oS,2]<- -95.5 # just set this so long is roughly in the middle of extent
coords[oS,3]<- 35    # just set this so lat is above all other sites

### build matrices for AEM
### rot.angle controls direction of the "asymmetric process" from bottom to top of plot area
mat <-aem.build.binary(nb.object = NULL, coords, link=segs, 
                       rm.same.y = FALSE, unit.angle = "degrees", rot.angle = 180, plot.connexions = TRUE)
### AEM analysis
SP_AEM_TUP_UW <- aem(mat, rm.link0 = FALSE)
### tidy up data and order based on original ID order in other dataframes 
SP_AEM_TUP_UW <- as.data.frame(SP_AEM_TUP_UW$vectors)
link[oS,3]<- "REMOVEUP"
row.names(SP_AEM_TUP_UW) <- link$V3
SP_AEM_TUP_UW <- SP_AEM_TUP_UW[order(match(row.names(SP_AEM_TUP_UW), as.matrix(row.names(taxon.div))[,1])),]
colnames(SP_AEM_TUP_UW)<-paste("AEM_TU_UW", colnames(SP_AEM_TUP_UW), sep = "_")
### select AEMS with positive spatial structure
library(spdep)
morans <-moran.randtest(SP_AEM_TUP_UW[-142,], listw)
morans.pos <-cbind(morans$names, morans$pvalue)
morans.pos<-morans.pos[morans.pos[,2] < 0.05,]
SP_AEM_TUP_UW<-subset(SP_AEM_TUP_UW, select=morans.pos[,1])
### remove "ghost" linkages
SP_AEM_TUP_UW <- SP_AEM_TUP_UW[-c(139:142),]

### _5.9 | AEM tail-up dam ---------------------------------------------------
### AEM analysis
### add a weight to the "ghost" segment
dams <- as.matrix(append(dams, 1))
SP_AEM_TUP_DAMW <- aem(mat, weight = dams[,1], rm.link0 = FALSE)
### tidy up data and order based on original ID order in other dataframes 
SP_AEM_TUP_DAMW <- as.data.frame(SP_AEM_TUP_DAMW$vectors)
link[oS,3]<- "REMOVEUP"
row.names(SP_AEM_TUP_DAMW) <- link$V3
SP_AEM_TUP_DAMW <- SP_AEM_TUP_DAMW[order(match(row.names(SP_AEM_TUP_DAMW), as.matrix(row.names(taxon.div))[,1])),]
colnames(SP_AEM_TUP_DAMW)<-paste("AEM_TU_DAM", colnames(SP_AEM_TUP_DAMW), sep = "_")
### select AEMS with positive spatial structure
morans <-moran.randtest(SP_AEM_TUP_DAMW[-142,], listw)
morans.pos <-cbind(morans$names, morans$pvalue)
morans.pos<-morans.pos[morans.pos[,2] < 0.05,]
SP_AEM_TUP_DAMW<-subset(SP_AEM_TUP_DAMW, select=morans.pos[,1])
### remove "ghost" linkages
SP_AEM_TUP_DAMW <- SP_AEM_TUP_DAMW[-c(139:142),]

### _5.10 | AEM tail-up distance ---------------------------------------------------
### AEM analysis
### weight based inverse squared distance
dmat <-as.matrix(dist(coords[,c(2:3)]))
edges <- mat$edges[-1,]
length.edge <- vector(length = nrow(segs))
for (i in 1:nrow(edges))
{
  length.edge[i] <- dmat[edges[i,1], edges[i,2]]
}
weight.vec <- 1-(length.edge/max(length.edge))^2
SP_AEM_TUP_DISTW <- aem(mat, weight = weight.vec, rm.link0 = FALSE)
### tidy up data and order based on original ID order in other dataframes 
SP_AEM_TUP_DISTW <- as.data.frame(SP_AEM_TUP_DISTW$vectors)
row.names(SP_AEM_TUP_DISTW) <- link$V3
SP_AEM_TUP_DISTW <- SP_AEM_TUP_DISTW[order(match(row.names(SP_AEM_TUP_DISTW), as.matrix(row.names(taxon.div))[,1])),]
colnames(SP_AEM_TUP_DISTW)<-paste("AEM_TU_DIST", colnames(SP_AEM_TUP_DISTW), sep = "_")
### select AEMS with positive spatial structure
morans <-moran.randtest(SP_AEM_TUP_DISTW[-142,], listw)
morans.pos <-cbind(morans$names, morans$pvalue)
morans.pos<-morans.pos[morans.pos[,2] < 0.05,]
SP_AEM_TUP_DISTW<-subset(SP_AEM_TUP_DISTW, select=morans.pos[,1])
### remove "ghost" linkages
SP_AEM_TUP_DISTW <- SP_AEM_TUP_DISTW[-c(139:142),]


### tidy environment
rm(vif,loadings,coords,dams,dmat,edges,link,mat,segs,sites,sites.shp,streams,term.segs.mat, hydro_dist, i, land_dist, length.edge, oS, term.segs, ustream_dist, weight.vec)

##########################                    END           ########################################
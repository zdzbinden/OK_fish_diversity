# BEGINNING ---------------------------------------------------------------
### Comparing Multiple Diversities Across Drainages                        
### ZD ZBINDEN                                       
### ____ This script is for data analyses  ____                        
### run all of script "A_data_prep" and have objects in Global Environment
### use document outline -->>>>> for easy viewing
# set directory -----------------------------------------------------------
setwd("C:/Users/zbind/Desktop/shelf/Projects/Drainage_Compare_SWG/R_workdir_v2/")

### 1 | PAIRWISE BETA MATRICES ---------------------------------------------------
### _ 1.1 | taxon beta -----------------------------------------------------------
### generate beta diversity pairwise matrices from taxon incidence & abundance
library("betapart")
taxon.incidence <- taxon.div
taxon.incidence[taxon.incidence >0]<-1
### taxon beta diversity based on incidence data (presence/absence)
beta.matrices <-beta.pair(taxon.incidence, index.family = "sorensen")
taxon.beta.sor <-beta.matrices$beta.sor
taxon.beta.sim <- beta.matrices$beta.sim
taxon.beta.sne <- beta.matrices$beta.sne
### taxon beta diversity based on abundance data
beta.matrices <-beta.pair.abund(taxon.div, index.family = "bray")
taxon.beta.bray.bal <- beta.matrices$beta.bray.bal
taxon.beta.bray.gra <- beta.matrices$beta.bray.gra
taxon.beta.bray <- beta.matrices$beta.bray




### _ 1.2 | functional beta ---------------------------------------------------
### generate beta diversity pairwise matrices from functional data
library("betapart")
library("doParallel")
### this takes a while to run, if already ran and saved you can skip to read in saved object below
nc <-detectCores()
beta.functional.core <-functional.betapart.core(taxon.incidence, functional.nmds, multi = FALSE, warning.time = FALSE,
                         return.details = TRUE, fbc.step = FALSE,
                         parallel = TRUE, opt.parallel = beta.para.control(nc=nc))
saveRDS(beta.functional.core, file="./data/beta.fun.core.rds")
#### Skip to here: 
beta.functional.core <- readRDS(file="./data/beta.fun.core.rds")
beta.matrices <-functional.beta.pair(beta.functional.core,index.family = "sorensen")
functional.beta.sor <- beta.matrices$funct.beta.sor
functional.beta.sim <- beta.matrices$funct.beta.sim
functional.beta.sne <- beta.matrices$funct.beta.sne


### _ 1.3 | phylogenetic beta -------------------------------------------------
### phylogenetic beta diversity
beta.matrices <-phylo.beta.pair(taxon.incidence, fish.phylo, index.family = "sorensen")
phylo.beta.sor <- beta.matrices$phylo.beta.sor
phylo.beta.sim <- beta.matrices$phylo.beta.sim
phylo.beta.sne <- beta.matrices$phylo.beta.sne


### 2 | PLOT BETA DISTRIBUTIONS ---------------------------------------------
### _ 2.1 | Overall Turnover Plots ------------------------------------------
### ___ 2.1.1 | taxon incidence (Bsor) ---------------------------------------
### get data from distance matrix into long format for plotting 
### taxon.beta.sor (Clear Boggy)
d<-as.matrix(taxon.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
taxon.list.dist <- t(combn(colnames(dd),2))
taxon.list.dist<-data.frame(taxon.list.dist, dist=dd[taxon.list.dist])
taxon.list.dist[ ,4]<- rep("Clear", nrow(taxon.list.dist))
### taxon.beta.sor (Muddy Boggy)
d<-as.matrix(taxon.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Muddy", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sor (Kiamchi)
d<-as.matrix(taxon.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sor (Little)
d<-as.matrix(taxon.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Little", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sor (Total)
d<-as.matrix(taxon.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Total", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
############## PLOT WITH GGPLOT2
colnames(taxon.list.dist)<-c("one", "two", "Bsor", "Drainage")
ggtaxonBsor <-ggplot(taxon.list.dist, aes(x=Bsor, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("TD Bsor")+
  ylab("TD (INCIDENCE)")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.x = element_blank())  
ggtaxonBsor


### ___ 2.1.2 | taxon abundance (Bray) --------------------------------------------
### get data from distance matrix into long format for plotting 
### taxon.beta.bray (Clear Boggy)
d<-as.matrix(taxon.beta.bray)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
taxon.list.dist <- t(combn(colnames(dd),2))
taxon.list.dist<-data.frame(taxon.list.dist, dist=dd[taxon.list.dist])
taxon.list.dist[ ,4]<- rep("Clear", nrow(taxon.list.dist))
### taxon.beta.bray (Muddy Boggy)
d<-as.matrix(taxon.beta.bray)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Muddy", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray (Kiamchi)
d<-as.matrix(taxon.beta.bray)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray (Little)
d<-as.matrix(taxon.beta.bray)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Little", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray (Total)
d<-as.matrix(taxon.beta.bray)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Total", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### plot with GGPLOT2
colnames(taxon.list.dist)<-c("one", "two", "Bray", "Drainage")
ggtaxonBray <-ggplot(taxon.list.dist, aes(x=Bray, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  labs(title = 'OVERALL TURNOVER')+
  xlab("TD Bray")+
  ylab("TD (ABUNDANCE)")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14, hjust = 0.5))+
  theme(axis.text.x = element_blank())  
ggtaxonBray


### ___ 2.1.3 | functional (Bsor) ---------------------------------------
### get data from distance matrix into long format for plotting 
### functional.beta.sor (Clear Boggy)
d<-as.matrix(functional.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
functional.list.dist <- t(combn(colnames(dd),2))
functional.list.dist<-data.frame(functional.list.dist, dist=dd[functional.list.dist])
functional.list.dist[ ,4]<- rep("Clear", nrow(functional.list.dist))
### functional.beta.sor (Muddy Boggy)
d<-as.matrix(functional.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Muddy", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sor (Kiamchi)
d<-as.matrix(functional.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sor (Little)
d<-as.matrix(functional.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Little", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sor (Total)
d<-as.matrix(functional.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Total", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### plot with GGPLOT2
colnames(functional.list.dist)<-c("one", "two", "Bsor", "Drainage")
ggfunctionalBsor <-ggplot(functional.list.dist, aes(x=Bsor, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("FD Bsor")+
  ylab("FUNCTIONAL")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.x = element_blank())  
ggfunctionalBsor


### ___ 2.1.4 | phylogenetic (Bsor) ---------------------------------------
### get data from distance matrix into long format for plotting 
### phylo.beta.sor (Clear Boggy)
d<-as.matrix(phylo.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
phylo.list.dist <- t(combn(colnames(dd),2))
phylo.list.dist<-data.frame(phylo.list.dist, dist=dd[phylo.list.dist])
phylo.list.dist[ ,4]<- rep("Clear", nrow(phylo.list.dist))
### phylo.beta.sor (Muddy Boggy)
d<-as.matrix(phylo.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Muddy", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sor (Kiamchi)
d<-as.matrix(phylo.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sor (Little)
d<-as.matrix(phylo.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Little", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sor (Total)
d<-as.matrix(phylo.beta.sor)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Total", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### plot with GGPLOT2
colnames(phylo.list.dist)<-c("one", "two", "Bsor", "Drainage")
ggphyloBsor <-ggplot(phylo.list.dist, aes(x=Bsor, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("PD Bsor")+
  ylab("PHYLOGENETIC")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14)) 
  ggphyloBsor

### _ 2.2 | Replacement Plots ---------------------------------------
### ___ 2.2.1 | taxon incidence (Bsim) ---------------------------------------
### get data from distance matrix into long format for plotting 
### taxon.beta.sim (Clear Boggy)
d<-as.matrix(taxon.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
taxon.list.dist <- t(combn(colnames(dd),2))
taxon.list.dist<-data.frame(taxon.list.dist, dist=dd[taxon.list.dist])
taxon.list.dist[ ,4]<- rep("Clear", nrow(taxon.list.dist))
### taxon.beta.sim (Muddy Boggy)
d<-as.matrix(taxon.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Muddy", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sim (Kiamchi)
d<-as.matrix(taxon.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sim (Little)
d<-as.matrix(taxon.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Little", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sim (Total)
d<-as.matrix(taxon.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Total", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
############## PLOT WITH GGPLOT2
colnames(taxon.list.dist)<-c("one", "two", "Bsim", "Drainage")
ggtaxonBsim <-ggplot(taxon.list.dist, aes(x=Bsim, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("TD Bsim")+
  ylab("")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                        quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                        jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_blank()) 
ggtaxonBsim
  
  
  
### ___ 2.2.2 | taxon abundance (Bray bal) ---------------------------------------
### get data from distance matrix into long format for plotting 
### taxon.beta.bray.bal (Clear Boggy)
d<-as.matrix(taxon.beta.bray.bal)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
taxon.list.dist <- t(combn(colnames(dd),2))
taxon.list.dist<-data.frame(taxon.list.dist, dist=dd[taxon.list.dist])
taxon.list.dist[ ,4]<- rep("Clear", nrow(taxon.list.dist))
### taxon.beta.bray.bal (Muddy Boggy)
d<-as.matrix(taxon.beta.bray.bal)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Muddy", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray.bal (Kiamchi)
d<-as.matrix(taxon.beta.bray.bal)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray.bal (Little)
d<-as.matrix(taxon.beta.bray.bal)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Little", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray.bal (Total)
d<-as.matrix(taxon.beta.bray.bal)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Total", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### plot with GGPLOT2
colnames(taxon.list.dist)<-c("one", "two", "BrayBal", "Drainage")
ggtaxonBrayBal <-ggplot(taxon.list.dist, aes(x=BrayBal, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  labs(title = 'REPLACEMENT')+
  xlab("TD Bray Balanced")+
  ylab("")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14, hjust = 0.5))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_blank()) 
ggtaxonBrayBal


### ___ 2.2.3 | functional (Bsim) ---------------------------------------
### get data from distance matrix into long format for plotting 
### functional.beta.sim (Clear Boggy)
d<-as.matrix(functional.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
functional.list.dist <- t(combn(colnames(dd),2))
functional.list.dist<-data.frame(functional.list.dist, dist=dd[functional.list.dist])
functional.list.dist[ ,4]<- rep("Clear", nrow(functional.list.dist))
### functional.beta.sim (Muddy Boggy)
d<-as.matrix(functional.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Muddy", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sim (Kiamchi)
d<-as.matrix(functional.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sim (Little)
d<-as.matrix(functional.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Little", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sim (Total)
d<-as.matrix(functional.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Total", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### plot with GGPLOT2
colnames(functional.list.dist)<-c("one", "two", "Bsim", "Drainage")
ggfunctionalBsim <-ggplot(functional.list.dist, aes(x=Bsim, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("FD Bsim")+
  ylab("")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_blank())  
ggfunctionalBsim


### ___ 2.2.4 | phylogenetic (Bsim) ---------------------------------------
### get data from distance matrix into long format for plotting 
### phylo.beta.sim (Clear Boggy)
d<-as.matrix(phylo.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
phylo.list.dist <- t(combn(colnames(dd),2))
phylo.list.dist<-data.frame(phylo.list.dist, dist=dd[phylo.list.dist])
phylo.list.dist[ ,4]<- rep("Clear", nrow(phylo.list.dist))
### phylo.beta.sim (Muddy Boggy)
d<-as.matrix(phylo.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Muddy", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sim (Kiamchi)
d<-as.matrix(phylo.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sim (Little)
d<-as.matrix(phylo.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Little", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sim (Total)
d<-as.matrix(phylo.beta.sim)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Total", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### plot with GGPLOT2
colnames(phylo.list.dist)<-c("one", "two", "Bsim", "Drainage")
ggphyloBsim <-ggplot(phylo.list.dist, aes(x=Bsim, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("PD Bsim")+
  ylab("")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.y = element_blank()) 
ggphyloBsim
  
  
### _ 2.3 | Nestedness Plots ---------------------------------------
### ___ 2.3.1 | taxon incidence (Bsne) ---------------------------------------
### get data from distance matrix into long format for plotting 
### taxon.beta.sne (Clear Boggy)
d<-as.matrix(taxon.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
taxon.list.dist <- t(combn(colnames(dd),2))
taxon.list.dist<-data.frame(taxon.list.dist, dist=dd[taxon.list.dist])
taxon.list.dist[ ,4]<- rep("Clear", nrow(taxon.list.dist))
### taxon.beta.sne (Muddy Boggy)
d<-as.matrix(taxon.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Muddy", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sne (Kiamchi)
d<-as.matrix(taxon.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sne (Little)
d<-as.matrix(taxon.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Little", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.sne (Total)
d<-as.matrix(taxon.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Total", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
############## PLOT WITH GGPLOT2
colnames(taxon.list.dist)<-c("one", "two", "Bsne", "Drainage")
ggtaxonBsne <-ggplot(taxon.list.dist, aes(x=Bsne, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("TD Bsne")+
  ylab("")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_blank())  
ggtaxonBsne

### ___ 2.3.2 | taxon abundance (Bray gra) ---------------------------------------
### get data from distance matrix into long format for plotting 
### taxon.beta.bray.gra (Clear Boggy)
d<-as.matrix(taxon.beta.bray.gra)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
taxon.list.dist <- t(combn(colnames(dd),2))
taxon.list.dist<-data.frame(taxon.list.dist, dist=dd[taxon.list.dist])
taxon.list.dist[ ,4]<- rep("Clear", nrow(taxon.list.dist))
### taxon.beta.bray.gra (Muddy Boggy)
d<-as.matrix(taxon.beta.bray.gra)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Muddy", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray.gra (Kiamchi)
d<-as.matrix(taxon.beta.bray.gra)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray.gra (Little)
d<-as.matrix(taxon.beta.bray.gra)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Little", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### taxon.beta.bray.gra (Total)
d<-as.matrix(taxon.beta.bray.gra)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
taxon.list.dist.tmp <- t(combn(colnames(dd),2))
taxon.list.dist.tmp <-data.frame(taxon.list.dist.tmp, dist=dd[taxon.list.dist.tmp])
taxon.list.dist.tmp[ ,4]<- rep("Total", nrow(taxon.list.dist.tmp))
taxon.list.dist <-rbind(taxon.list.dist, taxon.list.dist.tmp)
### plot with GGPLOT2
colnames(taxon.list.dist)<-c("one", "two", "BrayGra", "Drainage")
ggtaxonBrayGra <-ggplot(taxon.list.dist, aes(x=BrayGra, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  labs(title = 'NESTEDNESS')+
  xlab("TD Bray Gradient")+
  ylab("")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14, hjust = 0.5))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_blank())  
ggtaxonBrayGra


### ___ 2.3.3 | functional (Bsne) ---------------------------------------
### get data from distance matrix into long format for plotting 
### functional.beta.sne (Clear Boggy)
d<-as.matrix(functional.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
functional.list.dist <- t(combn(colnames(dd),2))
functional.list.dist<-data.frame(functional.list.dist, dist=dd[functional.list.dist])
functional.list.dist[ ,4]<- rep("Clear", nrow(functional.list.dist))
### functional.beta.sne (Muddy Boggy)
d<-as.matrix(functional.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Muddy", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sne (Kiamchi)
d<-as.matrix(functional.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sne (Little)
d<-as.matrix(functional.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Little", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### functional.beta.sne (Total)
d<-as.matrix(functional.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
functional.list.dist.tmp <- t(combn(colnames(dd),2))
functional.list.dist.tmp <-data.frame(functional.list.dist.tmp, dist=dd[functional.list.dist.tmp])
functional.list.dist.tmp[ ,4]<- rep("Total", nrow(functional.list.dist.tmp))
functional.list.dist <-rbind(functional.list.dist, functional.list.dist.tmp)
### plot with GGPLOT2
colnames(functional.list.dist)<-c("one", "two", "Bsne", "Drainage")
ggfunctionalBsne <-ggplot(functional.list.dist, aes(x=Bsne, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("FD Bsne")+
  ylab("")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_blank())  
ggfunctionalBsne

### ___ 2.3.4 | phylogenetic (Bsne) ---------------------------------------
### get data from distance matrix into long format for plotting 
### phylo.beta.sne (Clear Boggy)
d<-as.matrix(phylo.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[1:C]))
phylo.list.dist <- t(combn(colnames(dd),2))
phylo.list.dist<-data.frame(phylo.list.dist, dist=dd[phylo.list.dist])
phylo.list.dist[ ,4]<- rep("Clear", nrow(phylo.list.dist))
### phylo.beta.sne (Muddy Boggy)
d<-as.matrix(phylo.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+1):(C+M)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Muddy", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sne (Kiamchi)
d<-as.matrix(phylo.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+1):(C+M+K)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Kiamchi", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sne (Little)
d<-as.matrix(phylo.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)[(C+M+K+1):(C+M+K+L)]))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Little", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### phylo.beta.sne (Total)
d<-as.matrix(phylo.beta.sne)
dd<-as.matrix(usedist::dist_subset(d, row.names(env.rivs)))
phylo.list.dist.tmp <- t(combn(colnames(dd),2))
phylo.list.dist.tmp <-data.frame(phylo.list.dist.tmp, dist=dd[phylo.list.dist.tmp])
phylo.list.dist.tmp[ ,4]<- rep("Total", nrow(phylo.list.dist.tmp))
phylo.list.dist <-rbind(phylo.list.dist, phylo.list.dist.tmp)
### plot with GGPLOT2
colnames(phylo.list.dist)<-c("one", "two", "Bsne", "Drainage")
ggphyloBsne <-ggplot(phylo.list.dist, aes(x=Bsne, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlab("PD Bsne")+
  ylab("")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size=1, alpha=0.3, point_alpha=0.8,
                       jittered_points=FALSE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen","grey"))+
  xlim(0,1)+
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.y = element_blank()) 
ggphyloBsne

### ___ 2.4 | Multiplot ---------------------------------------
pdf("./outputs/div.plots/betadiv.pdf", paper = "USr", width = 0, height = 0)
multiplot(ggtaxonBray, ggtaxonBsor, ggfunctionalBsor, ggphyloBsor, 
          ggtaxonBrayBal, ggtaxonBsim, ggfunctionalBsim, ggphyloBsim,
          ggtaxonBrayGra, ggtaxonBsne, ggfunctionalBsne, ggphyloBsne,
          cols = 3)
dev.off()
multiplot(ggtaxonBray, ggtaxonBsor, ggfunctionalBsor, ggphyloBsor, 
          ggtaxonBrayBal, ggtaxonBsim, ggfunctionalBsim, ggphyloBsim,
          ggtaxonBrayGra, ggtaxonBsne, ggfunctionalBsne, ggphyloBsne,
          cols = 3)
  


### 3 | MODEL SELECTION  ---------------------------------------
### _ 3.1 | Taxon Abundance ---------------------------------------
### __  3.1.1 | Overall Turnover  ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(taxon.beta.bray ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
mod0 <- capscale(taxon.beta.bray ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
mod1 <- capscale(taxon.beta.bray ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                   R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
model$anova
terms <- attr(model$terminfo$terms,"term.labels")
red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(taxon.beta.bray ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(taxon.beta.bray ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(taxon.beta.bray ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(taxon.beta.bray ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(taxon.beta.bray ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(taxon.beta.bray ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(taxon.beta.bray ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(taxon.beta.bray ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(taxon.beta.bray ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(taxon.beta.bray ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(taxon.beta.bray ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(taxon.beta.bray ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(taxon.beta.bray ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(taxon.beta.bray ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(taxon.beta.bray ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(taxon.beta.bray ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(taxon.beta.bray ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(taxon.beta.bray, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(taxon.beta.bray ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(taxon.beta.bray ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(taxon.beta.bray ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(taxon.beta.bray ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(taxon.beta.bray ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### ___ 3.1.2 | Replacement ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(taxon.beta.bray.bal ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(taxon.beta.bray.bal ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(taxon.beta.bray.bal ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(taxon.beta.bray.bal ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(taxon.beta.bray.bal ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(taxon.beta.bray.bal ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(taxon.beta.bray.bal ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(taxon.beta.bray.bal ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.bal ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.bal ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(taxon.beta.bray.bal ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(taxon.beta.bray.bal, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(taxon.beta.bray.bal ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(taxon.beta.bray.bal ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(taxon.beta.bray.bal ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(taxon.beta.bray.bal ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(taxon.beta.bray.bal ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### ___ 3.1.3 | Nestedness ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-capscale(taxon.beta.bray.gra ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(taxon.beta.bray.gra ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(taxon.beta.bray.gra ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(taxon.beta.bray.gra ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(taxon.beta.bray.gra ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(taxon.beta.bray.gra ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(taxon.beta.bray.gra ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(taxon.beta.bray.gra ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.bray.gra ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.bray.gra ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(taxon.beta.bray.gra ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(taxon.beta.bray.gra, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(taxon.beta.bray.gra ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(taxon.beta.bray.gra ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(taxon.beta.bray.gra ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(taxon.beta.bray.gra ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(taxon.beta.bray.gra ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### _3.2 | Taxon Incidence ---------------------------------------
### __  3.2.1 | Overall Turnover  ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(taxon.beta.sor ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(taxon.beta.sor ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(taxon.beta.sor ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(taxon.beta.sor ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(taxon.beta.sor ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(taxon.beta.sor ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(taxon.beta.sor ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(taxon.beta.sor ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(taxon.beta.sor ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(taxon.beta.sor ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(taxon.beta.sor ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(taxon.beta.sor ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(taxon.beta.sor ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(taxon.beta.sor ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(taxon.beta.sor ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(taxon.beta.sor ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(taxon.beta.sor ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sor ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sor ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(taxon.beta.sor ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(taxon.beta.sor, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(taxon.beta.sor ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(taxon.beta.sor ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(taxon.beta.sor ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(taxon.beta.sor ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(taxon.beta.sor ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### ___ 3.2.2 | Replacement ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(taxon.beta.sim ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(taxon.beta.sim ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(taxon.beta.sim ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(taxon.beta.sim ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(taxon.beta.sim ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(taxon.beta.sim ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(taxon.beta.sim ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(taxon.beta.sim ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(taxon.beta.sim ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(taxon.beta.sim ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(taxon.beta.sim ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(taxon.beta.sim ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(taxon.beta.sim ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(taxon.beta.sim ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(taxon.beta.sim ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(taxon.beta.sim ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(taxon.beta.sim ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sim ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sim ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(taxon.beta.sim ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(taxon.beta.sim, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(taxon.beta.sim ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(taxon.beta.sim ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(taxon.beta.sim ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(taxon.beta.sim ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(taxon.beta.sim ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### ___ 3.2.3 | Nestedness ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(taxon.beta.sne ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(taxon.beta.sne ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(taxon.beta.sne ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(taxon.beta.sne ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(taxon.beta.sne ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(taxon.beta.sne ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(taxon.beta.sne ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(taxon.beta.sne ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(taxon.beta.sne ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(taxon.beta.sne ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(taxon.beta.sne ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(taxon.beta.sne ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(taxon.beta.sne ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(taxon.beta.sne ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(taxon.beta.sne ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(taxon.beta.sne ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(taxon.beta.sne ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(taxon.beta.sne ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(taxon.beta.sne ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(taxon.beta.sne ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(taxon.beta.sne, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(taxon.beta.sne ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(taxon.beta.sne ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(taxon.beta.sne ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(taxon.beta.sne ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(taxon.beta.sne ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)




### _3.3 | Functional ---------------------------------------
### ___ 3.3.1 | overall turnover ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(functional.beta.sor ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(functional.beta.sor ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(functional.beta.sor ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(functional.beta.sor ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(functional.beta.sor ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(functional.beta.sor ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(functional.beta.sor ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(functional.beta.sor ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(functional.beta.sor ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(functional.beta.sor ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(functional.beta.sor ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(functional.beta.sor ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(functional.beta.sor ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(functional.beta.sor ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(functional.beta.sor ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(functional.beta.sor ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(functional.beta.sor ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sor ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sor ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(functional.beta.sor ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(functional.beta.sor, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(functional.beta.sor ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(functional.beta.sor ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(functional.beta.sor ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(functional.beta.sor ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(functional.beta.sor ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### ___ 3.3.2 | replacement ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(functional.beta.sim ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(functional.beta.sim ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(functional.beta.sim ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(functional.beta.sim ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(functional.beta.sim ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(functional.beta.sim ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(functional.beta.sim ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(functional.beta.sim ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(functional.beta.sim ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(functional.beta.sim ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(functional.beta.sim ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(functional.beta.sim ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(functional.beta.sim ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(functional.beta.sim ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(functional.beta.sim ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(functional.beta.sim ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(functional.beta.sim ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sim ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sim ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(functional.beta.sim ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(functional.beta.sim, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(functional.beta.sim ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(functional.beta.sim ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(functional.beta.sim ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(functional.beta.sim ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(functional.beta.sim ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### ___ 3.3.2 | nestedness ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(functional.beta.sne ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(functional.beta.sne ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(functional.beta.sne ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(functional.beta.sne ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(functional.beta.sne ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(functional.beta.sne ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(functional.beta.sne ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(functional.beta.sne ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(functional.beta.sne ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(functional.beta.sne ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(functional.beta.sne ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(functional.beta.sne ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(functional.beta.sne ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(functional.beta.sne ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(functional.beta.sne ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(functional.beta.sne ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(functional.beta.sne ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(functional.beta.sne ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(functional.beta.sne ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(functional.beta.sne ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(functional.beta.sne, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(functional.beta.sne ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(functional.beta.sne ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(functional.beta.sne ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(functional.beta.sne ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(functional.beta.sne ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### _3.4 | Phylogenetic ---------------------------------------
### ___ 3.4.1 | overall turnover ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(phylo.beta.sor ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(phylo.beta.sor ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(phylo.beta.sor ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(phylo.beta.sor ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(phylo.beta.sor ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(phylo.beta.sor ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(phylo.beta.sor ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(phylo.beta.sor ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(phylo.beta.sor ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(phylo.beta.sor ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(phylo.beta.sor ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(phylo.beta.sor ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(phylo.beta.sor ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(phylo.beta.sor ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(phylo.beta.sor ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(phylo.beta.sor ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(phylo.beta.sor ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sor ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sor ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(phylo.beta.sor ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(phylo.beta.sor, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(phylo.beta.sor ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(phylo.beta.sor ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(phylo.beta.sor ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(phylo.beta.sor ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(phylo.beta.sor ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### ___ 3.4.2 | replacement ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(phylo.beta.sim ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(phylo.beta.sim ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(phylo.beta.sim ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(phylo.beta.sim ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(phylo.beta.sim ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(phylo.beta.sim ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(phylo.beta.sim ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(phylo.beta.sim ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(phylo.beta.sim ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(phylo.beta.sim ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(phylo.beta.sim ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(phylo.beta.sim ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(phylo.beta.sim ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(phylo.beta.sim ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(phylo.beta.sim ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(phylo.beta.sim ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(phylo.beta.sim ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sim ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sim ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(phylo.beta.sim ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(phylo.beta.sim, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(phylo.beta.sim ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(phylo.beta.sim ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(phylo.beta.sim ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(phylo.beta.sim ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(phylo.beta.sim ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)


### ___ 3.4.2 | nestedness ---------------------------------------
### ______ environmental ---------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-dbrda(phylo.beta.sne ~ ., data= env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=env.instream.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-dbrda(phylo.beta.sne ~ ., data= env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=env.hydrophysio.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-dbrda(phylo.beta.sne ~ ., data= env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=env.climate.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-dbrda(phylo.beta.sne ~ ., data= env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=env.landcover.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-dbrda(phylo.beta.sne ~ ., data= env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=env.geology.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-dbrda(phylo.beta.sne ~ ., data= env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=env.anthropogenic.PCs, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.anthropogenic.PCs<- subset(env.anthropogenic.PCs, select=terms)
} else {
  red.env.anthropogenic.PCs<-data.frame(row.names = rownames(env.anthropogenic.PCs)) 
}
model$anova
remove(mod0,mod1,model)

#### combine sig factors from model selection into one dataframe
ENV_select <- cbind.data.frame(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)

rda<-dbrda(phylo.beta.sne ~ ., data= ENV_select, sqrt.dist = FALSE, add = FALSE)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(ENV_select, th=10)
ENV_select<-usdm::exclude(ENV_select, vif)
remove(rda,sig.terms,terms,vif)


### ______ spatial ---------------------------------------
####### SP_MEM_LAND
rda<-dbrda(phylo.beta.sne ~ ., data= SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_MEM_LAND, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-dbrda(phylo.beta.sne ~ ., data= SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_MEM_HYDRO, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-dbrda(phylo.beta.sne ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_MEM_UPSTR_TOTAL, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-dbrda(phylo.beta.sne ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_MEM_UPSTR_NET, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-dbrda(phylo.beta.sne ~ ., data= SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_AEM_TDN_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-dbrda(phylo.beta.sne ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_AEM_TDN_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-dbrda(phylo.beta.sne ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_AEM_TDN_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-dbrda(phylo.beta.sne ~ ., data= SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_AEM_TUP_UW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW
rda<-dbrda(phylo.beta.sne ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_AEM_TUP_DAMW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-dbrda(phylo.beta.sne ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- capscale(phylo.beta.sne ~1, data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  mod1 <- capscale(phylo.beta.sne ~., data=SP_AEM_TUP_DISTW, sqrt.dist=FALSE, add=FALSE)
  model <-ordiR2step(mod0,mod1, Pin=0.05, permutations=how(nperm=999), 
                     R2permutations=999, direction="forward", R2scope=TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########## Combine select terms from different sets
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
### multiple regression test to confirm all terms are significantly related to response
rda<-dbrda(phylo.beta.sne ~., data= SP_select, sqrt.dist=FALSE, add = FALSE)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)

### ______  variation partitioning ---------------------------------------
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(phylo.beta.sne, SP_select, ENV_select, sqrt.dist = F)
partial.reg
plot(partial.reg, cex=1.2, bg=c("royalblue1","olivedrab"), Xnames=c("Spatial","Environmental"))
#### TESTING SIGNIFICANCE OF FRACTIONS
#### factions [a+b+c]:
comb <- cbind(SP_select, ENV_select)
rda.all<-dbrda(phylo.beta.sne ~ ., data = comb, sqrt.dist = F)
RsquareAdj(rda.all)
anova(rda.all)
#### fraction [a]:
rda.pure.spatial<-dbrda(phylo.beta.sne ~ . + Condition(as.matrix(ENV_select)),data=SP_select, sqrt.dist = F)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=999)
#### fraction [c]:
rda.pure.environ<-dbrda(phylo.beta.sne ~ . + Condition(as.matrix(SP_select)),data=ENV_select, sqrt.dist = F)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=999)
#### fraction [a+b]:
rda.spatial<-dbrda(phylo.beta.sne ~ ., data= SP_select, sqrt.dist = F)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-dbrda(phylo.beta.sne ~ ., data= ENV_select, sqrt.dist = F)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)



### 4 | MATRIX CORRELATION ---------------------------------------
protest(taxon.beta.bray, taxon.beta.bray.bal)
protest(taxon.beta.bray, taxon.beta.bray.gra)
protest(taxon.beta.bray, taxon.beta.sor)
protest(taxon.beta.bray, taxon.beta.sim)
protest(taxon.beta.bray, taxon.beta.sne)
protest(taxon.beta.bray, functional.beta.sor)
protest(taxon.beta.bray, functional.beta.sim)
protest(taxon.beta.bray, functional.beta.sne)
protest(taxon.beta.bray, phylo.beta.sor)
protest(taxon.beta.bray, phylo.beta.sim)
protest(taxon.beta.bray, phylo.beta.sne)

protest(taxon.beta.bray.bal, taxon.beta.bray.gra)
protest(taxon.beta.bray.bal, taxon.beta.sor)
protest(taxon.beta.bray.bal, taxon.beta.sim)
protest(taxon.beta.bray.bal, taxon.beta.sne)
protest(taxon.beta.bray.bal, functional.beta.sor)
protest(taxon.beta.bray.bal, functional.beta.sim)
protest(taxon.beta.bray.bal, functional.beta.sne)
protest(taxon.beta.bray.bal, phylo.beta.sor)
protest(taxon.beta.bray.bal, phylo.beta.sim)
protest(taxon.beta.bray.bal, phylo.beta.sne)

protest(taxon.beta.bray.gra, taxon.beta.sor)
protest(taxon.beta.bray.gra, taxon.beta.sim)
protest(taxon.beta.bray.gra, taxon.beta.sne)
protest(taxon.beta.bray.gra, functional.beta.sor)
protest(taxon.beta.bray.gra, functional.beta.sim)
protest(taxon.beta.bray.gra, functional.beta.sne)
protest(taxon.beta.bray.gra, phylo.beta.sor)
protest(taxon.beta.bray.gra, phylo.beta.sim)
protest(taxon.beta.bray.gra, phylo.beta.sne)

protest(taxon.beta.sor, taxon.beta.sim)
protest(taxon.beta.sor, taxon.beta.sne)
protest(taxon.beta.sor, functional.beta.sor)
protest(taxon.beta.sor, functional.beta.sim)
protest(taxon.beta.sor, functional.beta.sne)
protest(taxon.beta.sor, phylo.beta.sor)
protest(taxon.beta.sor, phylo.beta.sim)
protest(taxon.beta.sor, phylo.beta.sne)

protest(taxon.beta.sim, taxon.beta.sne)
protest(taxon.beta.sim, functional.beta.sor)
protest(taxon.beta.sim, functional.beta.sim)
protest(taxon.beta.sim, functional.beta.sne)
protest(taxon.beta.sim, phylo.beta.sor)
protest(taxon.beta.sim, phylo.beta.sim)
protest(taxon.beta.sim, phylo.beta.sne)

protest(taxon.beta.sne, functional.beta.sor)
protest(taxon.beta.sne, functional.beta.sim)
protest(taxon.beta.sne, functional.beta.sne)
protest(taxon.beta.sne, phylo.beta.sor)
protest(taxon.beta.sne, phylo.beta.sim)
protest(taxon.beta.sne, phylo.beta.sne)

protest(functional.beta.sor, functional.beta.sim)
protest(functional.beta.sor, functional.beta.sne)
protest(functional.beta.sor, phylo.beta.sor)
protest(functional.beta.sor, phylo.beta.sim)
protest(functional.beta.sor, phylo.beta.sne)

protest(functional.beta.sim, functional.beta.sne)
protest(functional.beta.sim, phylo.beta.sor)
protest(functional.beta.sim, phylo.beta.sim)
protest(functional.beta.sim, phylo.beta.sne)

protest(functional.beta.sne, phylo.beta.sor)
protest(functional.beta.sne, phylo.beta.sim)
protest(functional.beta.sne, phylo.beta.sne)

protest(phylo.beta.sor, phylo.beta.sim)
protest(phylo.beta.sor, phylo.beta.sne)

protest(phylo.beta.sim, phylo.beta.sne)
### END ---------------------------------------
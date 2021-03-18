### BEGINNING ---------------------------------------------------
### Comparing Multiple Diversities Across Drainages                                
### ZD ZBINDEN                                    
### *** This script is for data analyses  ***                                 
### run all of script "A_data_prep" and have resulting objects in Global Environment
### use document outline -->>>>> for easy viewing
### packages ---------------------------------------------------


### set directory ---------------------------------------------------
### specify project directory
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


### 1 | MULTI-SITE COMPARISON ---------------------------------------------------
### _1.1 | Diversity Summary Table ---------------------------------------------------
### create presence absence matrix and calculate species richness
taxon.incidence <- taxon.div
taxon.incidence[taxon.incidence >0]<-1
write.csv(taxon.incidence, "./outputs/taxon.incidence.csv")
### Create a summary table for diversity analyses by drainage
diversity.summary <- as.data.frame(matrix(nrow = 5, ncol = 13))
row.names(diversity.summary) <-c("Clear","Muddy","Kiamichi","Little","Total")
colnames(diversity.summary) <-c("N sites", "Taxon Gamma", "Function Gamma","Phylo Gamma","Taxon Alpha", "Function Alpha", "Phylo Alpha","Taxon Bsor", "Taxon Bsim", "Taxon Bsne","Phylo Bsor", "Phylo Bsim", "Phylo Bsne")
### Add number of sites to matrix
diversity.summary[1,1] <- nrow(site_info[site_info$drainage == "Clear",])
diversity.summary[2,1] <- nrow(site_info[site_info$drainage == "Muddy",])
diversity.summary[3,1] <- nrow(site_info[site_info$drainage == "Kiamichi",])
diversity.summary[4,1] <- nrow(site_info[site_info$drainage == "Little",])
diversity.summary[5,1] <- nrow(site_info)
C <-diversity.summary[1,1]
M <-diversity.summary[2,1]
K <-diversity.summary[3,1]
L <-diversity.summary[4,1]
T <-diversity.summary[5,1]
### ___1.1.1 | Taxon Summaries ---------------------------------------------------
######## compute gamma diversity
library("betapart")
### gamma of total study
taxon.gamma <-betapart.core(taxon.incidence) # core function
diversity.summary[5,2] <-taxon.gamma$St # total species in matrix
### gamma clear boggy
taxon.gamma<-betapart.core(taxon.incidence[1:C,])
diversity.summary[1,2] <-taxon.gamma$St
### gamma muddy boggy
taxon.gamma <-betapart.core(taxon.incidence[(C+1):(C+M),])
diversity.summary[2,2] <-taxon.gamma$St
### gamma kiamichi
taxon.gamma <-betapart.core(taxon.incidence[(C+M+1):(C+M+K),])
diversity.summary[3,2] <-taxon.gamma$St
### gamma little river
taxon.gamma <-betapart.core(taxon.incidence[(C+M+K+1):(C+M+K+L),])
diversity.summary[4,2] <-taxon.gamma$St
######## compute alpha diversity (average species richness)
diversity.summary[5,5] <- round(mean(rowSums(taxon.incidence)), digits = 1) # total
diversity.summary[1,5] <- round(mean(rowSums(taxon.incidence[1:C,])), digits = 1) # clear
diversity.summary[2,5] <- round(mean(rowSums(taxon.incidence[(C+1):(C+M),])), digits = 1) # muddy
diversity.summary[3,5] <- round(mean(rowSums(taxon.incidence[(C+M+1):(C+M+K),])), digits = 1) # kiamichi
diversity.summary[4,5] <- round(mean(rowSums(taxon.incidence[(C+M+K+1):(C+M+K+L),])), digits = 1) # little
taxon.alpha <- rowSums(taxon.incidence)
taxon.alpha <- as.data.frame(taxon.alpha, row.names = rownames(taxon.incidence))
######## compute beta diversity
### use beta sample rather than multi because of different sample sizes across drainages
### total
beta.sum <- beta.sample.abund(taxon.div, index.family = "bray", sites = 25, samples = 1000)
diversity.summary[5,8] <- round(beta.sum$mean.values[3], digits = 3)
diversity.summary[5,9] <- round(beta.sum$mean.values[1], digits = 3)
diversity.summary[5,10] <- round(beta.sum$mean.values[2], digits = 3)
##### clear
beta.sum <- beta.sample.abund(taxon.div[1:C,], index.family = "bray", sites=25, samples = 1000)
diversity.summary[1,8] <- round(beta.sum$mean.values[3], digits = 3)
diversity.summary[1,9] <- round(beta.sum$mean.values[1], digits = 3)
diversity.summary[1,10] <- round(beta.sum$mean.values[2], digits = 3)
#### muddy
beta.sum <- beta.sample.abund(taxon.div[(C+1):(C+M),], index.family = "bray", sites=25, samples = 1000)
diversity.summary[2,8] <- round(beta.sum$mean.values[3], digits = 3)
diversity.summary[2,9] <- round(beta.sum$mean.values[1], digits = 3)
diversity.summary[2,10] <- round(beta.sum$mean.values[2], digits = 3)
#### kiamichi
beta.sum <- beta.sample.abund(taxon.div[(C+M+1):(C+M+K),], index.family = "bray", sites=25, samples = 1000)
diversity.summary[3,8] <- round(beta.sum$mean.values[3], digits = 3)
diversity.summary[3,9] <- round(beta.sum$mean.values[1], digits = 3)
diversity.summary[3,10] <- round(beta.sum$mean.values[2], digits = 3)
#### little
beta.sum <- beta.sample.abund(taxon.div[(C+M+K+1):(C+M+K+L),], index.family = "bray", sites=25, samples = 1000)
diversity.summary[4,8] <- round(beta.sum$mean.values[3], digits = 3)
diversity.summary[4,9] <- round(beta.sum$mean.values[1], digits = 3)
diversity.summary[4,10] <- round(beta.sum$mean.values[2], digits = 3)


### ___1.1.2 | Phylogenetic Summaries ---------------------------------------------------
######## Total phylo gamma
phylo.core <- phylo.betapart.core(taxon.incidence, fish.phylo)
diversity.summary[5,4]<-round(phylo.core$St, digits = 1)
#### Clear boggy phylo gamma
phylo.core <- phylo.betapart.core(taxon.incidence[1:C,], fish.phylo)
diversity.summary[1,4]<-round(phylo.core$St, digits = 1)
#### Muddy boggy phylo gamma
phylo.core <- phylo.betapart.core(taxon.incidence[(C+1):(C+M),], fish.phylo)
diversity.summary[2,4]<-round(phylo.core$St, digits = 1)
#### Kiamichi phylo gamma
phylo.core <- phylo.betapart.core(taxon.incidence[(C+M+1):(C+M+K),], fish.phylo)
diversity.summary[3,4]<-round(phylo.core$St, digits = 1)
### Kiamichi phylo gamma
phylo.core <- phylo.betapart.core(taxon.incidence[(C+M+K+1):(C+M+K+L),], fish.phylo)
diversity.summary[4,4]<-round(phylo.core$St, digits = 1)
######## Mean phylo alpha
library("picante")
######## Total phylo alpha mean
phylogenetic.alpha <-picante::pd(taxon.incidence, fish.phylo)
diversity.summary[5,7]<-round(mean(phylogenetic.alpha$PD), digits = 1)
### Clear phylo alpha mean
phylo.alpha <-picante::pd(taxon.incidence[1:C,], fish.phylo)
diversity.summary[1,7]<-round(mean(phylo.alpha$PD), digits = 1)
### Muddy phylo alpha mean
phylo.alpha <-picante::pd(taxon.incidence[(C+1):(C+M),], fish.phylo)
diversity.summary[2,7]<-round(mean(phylo.alpha$PD), digits = 1)
### Kiamichi phylo alpha mean
phylo.alpha <-picante::pd(taxon.incidence[(C+M+1):(C+M+K),], fish.phylo)
diversity.summary[3,7]<-round(mean(phylo.alpha$PD), digits = 1)
### Little phylo alpha mean
phylo.alpha <-picante::pd(taxon.incidence[(C+M+K+1):(C+M+K+L),], fish.phylo)
diversity.summary[4,7]<-round(mean(phylo.alpha$PD), digits = 1)
########## compute beta diversity
### No beta.sample.abund exists for phylo diversity so this uses a loop to sample 25 sites 100 times and take the average
### total
### number of sites to sample
n <- 25
### number of replications
N <- 100
p.b.sor <- rep(NA, N)
p.b.sim <- rep(NA, N)
p.b.sne <- rep(NA, N)
for (i in 1:N){
  X <- sample(c(1:T), size = n, replace = FALSE)
  beta.sum <- phylo.beta.multi(taxon.incidence[c(X),], fish.phylo, index.family = "sorensen")
  p.b.sor[i] <- beta.sum$phylo.beta.SOR
  p.b.sim[i] <- beta.sum$phylo.beta.SIM
  p.b.sne[i] <- beta.sum$phylo.beta.SNE
}
diversity.summary[5,11] <- round(mean(p.b.sor), digits = 3)
diversity.summary[5,12] <- round(mean(p.b.sim), digits = 3)
diversity.summary[5,13] <- round(mean(p.b.sne), digits = 3)

### clear
### number of sites to sample
n <- 25
### number of replications
N <- 100
p.b.sor <- rep(NA, N)
p.b.sim <- rep(NA, N)
p.b.sne <- rep(NA, N)
for (i in 1:N){
  X <- sample(c(1:C), size = n, replace = FALSE)
  beta.sum <- phylo.beta.multi(taxon.incidence[c(X),], fish.phylo, index.family = "sorensen")
  p.b.sor[i] <- beta.sum$phylo.beta.SOR
  p.b.sim[i] <- beta.sum$phylo.beta.SIM
  p.b.sne[i] <- beta.sum$phylo.beta.SNE
}
diversity.summary[1,11] <- round(mean(p.b.sor), digits = 3)
diversity.summary[1,12] <- round(mean(p.b.sim), digits = 3)
diversity.summary[1,13] <- round(mean(p.b.sne), digits = 3)

### muddy 
### number of sites to sample
n <- 25
### number of replications
N <- 100
p.b.sor <- rep(NA, N)
p.b.sim <- rep(NA, N)
p.b.sne <- rep(NA, N)
for (i in 1:N){
  X <- sample(c((C+1):(C+M)), size = n, replace = FALSE)
  beta.sum <- phylo.beta.multi(taxon.incidence[c(X),], fish.phylo, index.family = "sorensen")
  p.b.sor[i] <- beta.sum$phylo.beta.SOR
  p.b.sim[i] <- beta.sum$phylo.beta.SIM
  p.b.sne[i] <- beta.sum$phylo.beta.SNE
}
diversity.summary[2,11] <- round(mean(p.b.sor), digits = 3)
diversity.summary[2,12] <- round(mean(p.b.sim), digits = 3)
diversity.summary[2,13] <- round(mean(p.b.sne), digits = 3)

### kiamichi
### number of sites to sample
n <- 25
### number of replications
N <- 100
p.b.sor <- rep(NA, N)
p.b.sim <- rep(NA, N)
p.b.sne <- rep(NA, N)
for (i in 1:N){
  X <- sample(c((C+M+1):(C+M+K)), size = n, replace = FALSE)
  beta.sum <- phylo.beta.multi(taxon.incidence[c(X),], fish.phylo, index.family = "sorensen")
  p.b.sor[i] <- beta.sum$phylo.beta.SOR
  p.b.sim[i] <- beta.sum$phylo.beta.SIM
  p.b.sne[i] <- beta.sum$phylo.beta.SNE
}
diversity.summary[3,11] <- round(mean(p.b.sor), digits = 3)
diversity.summary[3,12] <- round(mean(p.b.sim), digits = 3)
diversity.summary[3,13] <- round(mean(p.b.sne), digits = 3)

### Little 
### number of sites to sample
n <- 25
### number of replications
N <- 100
p.b.sor <- rep(NA, N)
p.b.sim <- rep(NA, N)
p.b.sne <- rep(NA, N)
for (i in 1:N){
  X <- sample(c((C+M+K+1):(C+M+K+L)), size = n, replace = FALSE)
  beta.sum <- phylo.beta.multi(taxon.incidence[c(X),], fish.phylo, index.family = "sorensen")
  p.b.sor[i] <- beta.sum$phylo.beta.SOR
  p.b.sim[i] <- beta.sum$phylo.beta.SIM
  p.b.sne[i] <- beta.sum$phylo.beta.SNE
}
diversity.summary[4,11] <- round(mean(p.b.sor), digits = 3)
diversity.summary[4,12] <- round(mean(p.b.sim), digits = 3)
diversity.summary[4,13] <- round(mean(p.b.sne), digits = 3)



### ___1.1.3 | Functional Summaries ---------------------------------------------------
### this function from betapart provides functional richness but not as a proportion of the whole volume
#fun.core <-functional.betapart.core(taxon.incidence, functional.pcoa, multi = FALSE, warning.time = FALSE,
#                         return.details = TRUE, fbc.step = FALSE,
#                         parallel = TRUE, opt.parallel = beta.para.control())
### how do results compare to multidimFD function? 
#fun.core$details$CH$FRi

### multi-site functional calculation does not finish -> illustrate with distributions instead
#fun.mutli <- functional.beta.multi(taxon.incidence, functional.pcoa, 
#                                   index.family = "sorensen", warning.time = FALSE)
#fun.mutli$funct.beta.SOR
#fun.mutli$funct.beta.SIM
#fun.mutli$funct.beta.SNE


######## Functional richness using function multidimFD
source("./FD/functions/multidimFD.R")
### calculate functional richness for each site
functional.richness <-as.data.frame(multidimFD(coord = as.matrix(functional.nmds), weight = as.matrix(taxon.div)))
functional.richness <- as.data.frame(functional.richness$FRic, rownames(taxon.div))
colnames(functional.richness) <- "FRi"
### compute functional alpha diversity (average functional richness)
diversity.summary[5,6] <- round(mean(functional.richness$FRi), digits =2) # total
diversity.summary[1,6] <- round(mean(functional.richness[1:C,]), digits =2) # clear
diversity.summary[2,6] <- round(mean(functional.richness[(C+1):(C+M),]), digits =2) # muddy
diversity.summary[3,6] <- round(mean(functional.richness[(C+M+1):(C+M+K),]), digits =2) # kiamichi
diversity.summary[4,6] <- round(mean(functional.richness[(C+M+K+1):(C+M+K+L),]), digits =2) # little
### Create a summary table for diversity analyses by drainage
functional.gamma <- as.data.frame(matrix(nrow = 4, ncol = 58))
row.names(functional.gamma) <- c("clearBoggy","MuddyBoggy","Kiamichi","LittleRiver")
colnames(functional.gamma) <- colnames(taxon.incidence)
### calculate overall functional richness proportions for drainages separately
functional.gamma[1,] <- colSums(taxon.incidence[1:C,])
functional.gamma[2,] <- colSums(taxon.incidence[(C+1):(C+M),])
functional.gamma[3,] <- colSums(taxon.incidence[(C+M+1):(C+M+K),])
functional.gamma[4,] <- colSums(taxon.incidence[(C+M+K+1):(C+M+K+L),])
functional.gamma.richness <-as.data.frame(multidimFD(coord = as.matrix(functional.nmds), weight = as.matrix(functional.gamma)))
functional.gamma.richness <- as.data.frame(functional.gamma.richness$FRic, rownames(functional.gamma))
colnames(functional.gamma.richness) <- "FRi"
diversity.summary[1,3] <- round(functional.gamma.richness[1,], digits = 2)
diversity.summary[2,3] <- round(functional.gamma.richness[2,], digits = 2)
diversity.summary[3,3] <- round(functional.gamma.richness[3,], digits = 2)
diversity.summary[4,3] <- round(functional.gamma.richness[4,], digits = 2)
diversity.summary[5,3] <- 1.0
write.csv(diversity.summary, "./outputs/divesity_summary.csv")


### make table of correlations among alpha diversity facets
data <- cbind(taxon.alpha, functional.richness, phylogenetic.alpha$PD)
alpha.cor <-cor(data, method = "spearman")
row.names(alpha.cor) <- c("Taxon Alpha", "Functional Alpha", "Phylo Alpha")
colnames(alpha.cor) <- c("Taxon Alpha", "Functional Alpha", "Phylo Alpha")
write.csv(alpha.cor, "./outputs/alpha_correlation_table.csv")


### _1.2 | Multifacet Plots ---------------------------------------------------
library(ggplot2)
library(ggridges)
### taxon
taxon.alpha[1:C,2] <- rep("Clear", C)
taxon.alpha[(C+1):(C+M),2] <- rep("Muddy", M)
taxon.alpha[(C+M+1):(C+M+K),2] <- rep("Kiamichi", K)
taxon.alpha[(C+M+K+1):(C+M+K+L),2] <- rep("Little", L)
colnames(taxon.alpha)<-c("Species_Richness", "Drainage")
ggTD <-ggplot(taxon.alpha, aes(x=Species_Richness, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"))+
  labs(title = 'Alpha Diversity')+
  xlab("Taxon Richness")+
  ylab("TAXON")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size = 1, alpha=0.2, point_alpha=0.8,
                       jittered_points=TRUE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen"))+
  scale_x_continuous(breaks = c(5,10,15,20,25)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))
ggTD

### functional
functional.richness[1:C,2] <- rep("Clear", C)
functional.richness[(C+1):(C+M),2] <- rep("Muddy", M)
functional.richness[(C+M+1):(C+M+K),2] <- rep("Kiamichi", K)
functional.richness[(C+M+K+1):(C+M+K+L),2] <- rep("Little", L)
colnames(functional.richness)<-c("Functional_Richness", "Drainage")
ggFD <-ggplot(functional.richness, aes(x=Functional_Richness, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"))+
  labs(title = '')+
  xlab("Functional Richness")+
  ylab("FUNCTIONAL")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size = 1, alpha=0.2, point_alpha=0.8,
                       jittered_points=TRUE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen"))+
  scale_x_continuous(breaks = c(0,0.25,0.5)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))
ggFD

### phylogenetic
phylogenetic.alpha[1:C,2] <- rep("Clear", C)
phylogenetic.alpha[(C+1):(C+M),2] <- rep("Muddy", M)
phylogenetic.alpha[(C+M+1):(C+M+K),2] <- rep("Kiamichi", K)
phylogenetic.alpha[(C+M+K+1):(C+M+K+L),2] <- rep("Little", L)
colnames(phylogenetic.alpha)<-c("Phylogenetic_Diversity", "Drainage")
ggPD <-ggplot(phylogenetic.alpha, aes(x=Phylogenetic_Diversity, y=Drainage, group=Drainage, fill=Drainage))+
  scale_fill_manual(values=c("red3","gold3","blue3","darkgreen"))+
  labs(title = '')+
  xlab("Phylogenetic Diversity")+
  ylab("PHYLOGENETIC")+
  geom_density_ridges2(aes(point_color = Drainage, point_fill = Drainage),
                       quantile_lines = TRUE, quantiles = 2,vline_size = 1, alpha=0.2, point_alpha=0.8,
                       jittered_points=TRUE, rel_min_height=0)+
  scale_discrete_manual(aesthetics = "point_color",values=c("red3","gold3","blue3","darkgreen"))+
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(legend.position = "none")+
  theme(text = element_text(size=14))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(size=14))
ggPD

### combine plots
pdf("./outputs/div.plots/alphadiv.pdf", paper = "US", width = 0, height = 0)
multiplot(ggTD, ggFD, ggPD, cols=1)
dev.off()
multiplot(ggTD, ggFD, ggPD, cols=1)

### 2 | MODEL SELECTION ---------------------------------------------------
### _2.1 | Taxon Alpha ---------------------------------------------------
### ___ environmental selection ---------------------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-rda(taxon.alpha$Species_Richness ~ ., data= env.instream.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
mod0 <- rda(taxon.alpha$Species_Richness ~1, data=env.instream.PCs)
mod1 <- rda(taxon.alpha$Species_Richness ~., data=env.instream.PCs)
model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                   R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
model$anova
terms <- attr(model$terminfo$terms,"term.labels")
red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-rda(taxon.alpha$Species_Richness ~ ., data= env.hydrophysio.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=env.hydrophysio.PCs)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=env.hydrophysio.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-rda(taxon.alpha$Species_Richness ~ ., data= env.climate.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=env.climate.PCs)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=env.climate.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-rda(taxon.alpha$Species_Richness ~ ., data= env.landcover.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=env.landcover.PCs)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=env.landcover.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-rda(taxon.alpha$Species_Richness ~ ., data= env.geology.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=env.geology.PCs)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=env.geology.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-rda(taxon.alpha$Species_Richness ~ ., data= env.anthropogenic.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=env.anthropogenic.PCs)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=env.anthropogenic.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
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
usdm::vif(ENV_select)
rda<-rda(taxon.alpha$Species_Richness ~., data= ENV_select)
RsquareAdj(rda)
anova(rda)
#sig.terms<-anova(rda, by="term", permutations=99999)
#sig.terms
### remove non-sig vars
#terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
#terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
#terms<-terms[terms$V2<0.1,]
#ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
#vif<-usdm::vifstep(ENV_select, th=10)
#ENV_select<-usdm::exclude(ENV_select, vif)
#remove(rda,sig.terms,terms,vif)


### ___ spatial selection ---------------------------------------------------
####### forward model selection for each variable set
######## SP_MEM_LAND
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_MEM_LAND, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_MEM_LAND)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_MEM_LAND)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_MEM_HYDRO, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_MEM_HYDRO)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_MEM_HYDRO)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_MEM_UPSTR_TOTAL)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_MEM_UPSTR_TOTAL)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_MEM_UPSTR_NET)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_MEM_UPSTR_NET)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_AEM_TDN_UW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_AEM_TDN_UW)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_AEM_TDN_UW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_AEM_TDN_DAMW)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_AEM_TDN_DAMW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_AEM_TDN_DISTW)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_AEM_TDN_DISTW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_AEM_TUP_UW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_AEM_TUP_UW)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_AEM_TUP_UW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW 
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_AEM_TUP_DAMW)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_AEM_TUP_DAMW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-rda(taxon.alpha$Species_Richness ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(taxon.alpha$Species_Richness ~1, data=SP_AEM_TUP_DISTW)
  mod1 <- rda(taxon.alpha$Species_Richness ~., data=SP_AEM_TUP_DISTW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### GET SELECTED TERMS
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
### multiple regression test to confirm all terms are significantly related to response
rda<-rda(taxon.alpha$Species_Richness~., data= SP_select)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=99999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
#vif<-usdm::vifstep(SP_select, th=10)
#SP_select<-usdm::exclude(SP_select, vif)
SP_select<-SP_select$AEM_TU_DIST_V1
SP_select<- as.data.frame(SP_select)
colnames(SP_select)<-c("AEM_TU_DIST_V1")
#remove(rda, sig.terms, terms, vif)

### ___ variation partitioning ---------------------------------------------------
### Partial redundancy analysis of models (taxon alpha)
### partial RDA of full environmental set
comb <- cbind(SP_select, ENV_select)
rda.all<-rda(taxon.alpha$Species_Richness ~ ., data=comb)
RsquareAdj(rda.all)
anova(rda.all)
#anova(rda.all, by="term", permutations=99999)
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(taxon.alpha$Species_Richness, SP_select, ENV_select)
partial.reg

#### TESTING SIGNIFICANCE OF FRACTIONS
#### fraction [a]:
rda.pure.spatial<-rda(taxon.alpha$Species_Richness ~ . + Condition(as.matrix(ENV_select)), data=SP_select)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=99999)
#### fraction [c]:
rda.pure.environ<-rda(taxon.alpha$Species_Richness ~ . + Condition(as.matrix(SP_select)), data=ENV_select)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=99999)
#### fraction [a+b]:
rda.spatial<-rda(taxon.alpha$Species_Richness ~., data= SP_select)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-rda(taxon.alpha$Species_Richness ~ ., data= ENV_select)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=99999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)

### _2.2 | Functional Alpha ---------------------------------------------------
### ___ environmental selection ---------------------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-rda(functional.richness$Functional_Richness ~ ., data= env.instream.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=env.instream.PCs)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=env.instream.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-rda(functional.richness$Functional_Richness ~ ., data= env.hydrophysio.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=env.hydrophysio.PCs)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=env.hydrophysio.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-rda(functional.richness$Functional_Richness ~ ., data= env.climate.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=env.climate.PCs)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=env.climate.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-rda(functional.richness$Functional_Richness ~ ., data= env.landcover.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=env.landcover.PCs)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=env.landcover.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-rda(functional.richness$Functional_Richness ~ ., data= env.geology.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=env.geology.PCs)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=env.geology.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-rda(functional.richness$Functional_Richness ~ ., data= env.anthropogenic.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=env.anthropogenic.PCs)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=env.anthropogenic.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
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
usdm::vif(ENV_select)
rda<-rda(functional.richness$Functional_Richness ~., data= ENV_select)
RsquareAdj(rda)
anova(rda)
#sig.terms<-anova(rda, by="term", permutations=99999)
#sig.terms
### remove non-sig vars
#terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
#terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
#terms<-terms[terms$V2<0.1,]
#ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
#vif<-usdm::vifstep(ENV_select, th=10)
#ENV_select<-usdm::exclude(ENV_select, vif)
#remove(rda,sig.terms,terms,vif)


### ___ spatial selection ---------------------------------------------------
####### forward model selection for each variable set
######## SP_MEM_LAND
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_MEM_LAND, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_MEM_LAND)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_MEM_LAND)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_MEM_HYDRO, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_MEM_HYDRO)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_MEM_HYDRO)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_MEM_UPSTR_TOTAL)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_MEM_UPSTR_TOTAL)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_MEM_UPSTR_NET)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_MEM_UPSTR_NET)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_AEM_TDN_UW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_AEM_TDN_UW)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_AEM_TDN_UW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_AEM_TDN_DAMW)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_AEM_TDN_DAMW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_AEM_TDN_DISTW)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_AEM_TDN_DISTW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_AEM_TUP_UW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_AEM_TUP_UW)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_AEM_TUP_UW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW 
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_AEM_TUP_DAMW)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_AEM_TUP_DAMW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-rda(functional.richness$Functional_Richness ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(functional.richness$Functional_Richness ~1, data=SP_AEM_TUP_DISTW)
  mod1 <- rda(functional.richness$Functional_Richness ~., data=SP_AEM_TUP_DISTW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### GET SELECTED TERMS
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
### multiple regression test to confirm all terms are significantly related to response
rda<-rda(functional.richness$Functional_Richness~., data= SP_select)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=99999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
remove(rda, sig.terms, terms, vif)


### ___ variation partitioning ---------------------------------------------------
### Partial redundancy analysis of models (taxon alpha)
### partial RDA of full environmental set
comb <- cbind(SP_select, ENV_select)
rda.all<-rda(functional.richness$Functional_Richness ~ ., data=comb)
RsquareAdj(rda.all)
anova(rda.all)
#anova(rda.all, by="term", permutations=99999)
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(functional.richness$Functional_Richness, SP_select, ENV_select)
partial.reg

#### TESTING SIGNIFICANCE OF FRACTIONS
#### fraction [a]:
rda.pure.spatial<-rda(functional.richness$Functional_Richness ~ . + Condition(as.matrix(ENV_select)), data=SP_select)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=99999)
#### fraction [c]:
rda.pure.environ<-rda(functional.richness$Functional_Richness ~ . + Condition(as.matrix(SP_select)), data=ENV_select)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=99999)
#### fraction [a+b]:
rda.spatial<-rda(functional.richness$Functional_Richness ~., data= SP_select)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-rda(functional.richness$Functional_Richness ~ ., data= ENV_select)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=99999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)
### _2.3 | Phylogenetic Alpha ---------------------------------------------------
### ___ environmental selection ---------------------------------------------------
library("vegan")
####### forward model selection for each variable set
### instream
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= env.instream.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] < 0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=env.instream.PCs)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=env.instream.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.instream.PCs<- subset(env.instream.PCs, select=terms)
} else {
  red.env.instream.PCs<-data.frame(row.names = rownames(env.instream.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### hydrophysio
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= env.hydrophysio.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=env.hydrophysio.PCs)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=env.hydrophysio.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.hydrophysio.PCs<- subset(env.hydrophysio.PCs, select=terms)
} else {
  red.env.hydrophysio.PCs<-data.frame(row.names = rownames(env.hydrophysio.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### climate
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= env.climate.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=env.climate.PCs)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=env.climate.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.climate.PCs<- subset(env.climate.PCs, select=terms)
} else {
  red.env.climate.PCs<-data.frame(row.names = rownames(env.climate.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### landcover
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= env.landcover.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=env.landcover.PCs)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=env.landcover.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.landcover.PCs<- subset(env.landcover.PCs, select=terms)
} else {
  red.env.landcover.PCs<-data.frame(row.names = rownames(env.landcover.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### geology
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= env.geology.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=env.geology.PCs)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=env.geology.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.env.geology.PCs<- subset(env.geology.PCs, select=terms)
} else {
  red.env.geology.PCs<-data.frame(row.names = rownames(env.geology.PCs)) 
}
model$anova
remove(mod0,mod1,model)

### anthropogenic
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= env.anthropogenic.PCs, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=env.anthropogenic.PCs)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=env.anthropogenic.PCs)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
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
usdm::vif(ENV_select)
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data= ENV_select)
RsquareAdj(rda)
anova(rda)
#sig.terms<-anova(rda, by="term", permutations=99999)
#sig.terms
### remove non-sig vars
#terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
#terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
#terms<-terms[terms$V2<0.1,]
#ENV_select<-subset(ENV_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
#vif<-usdm::vifstep(ENV_select, th=10)
#ENV_select<-usdm::exclude(ENV_select, vif)
#remove(rda,sig.terms,terms,vif)


### ___ spatial selection ---------------------------------------------------
####### forward model selection for each variable set
######## SP_MEM_LAND
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_MEM_LAND, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_MEM_LAND)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_MEM_LAND)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_LAND<- subset(SP_MEM_LAND, select=terms)
} else {
  red.SP_MEM_LAND<-data.frame(row.names = rownames(SP_MEM_LAND)) 
}
model$anova
remove(mod0,mod1,model)

######## SP_MEM_HYDRO
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_MEM_HYDRO, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_MEM_HYDRO)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_MEM_HYDRO)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_HYDRO<- subset(SP_MEM_HYDRO, select=terms)
} else {
  red.SP_MEM_HYDRO<-data.frame(row.names = rownames(SP_MEM_HYDRO)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_TOTAL
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_MEM_UPSTR_TOTAL, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_MEM_UPSTR_TOTAL)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_MEM_UPSTR_TOTAL)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_TOTAL<- subset(SP_MEM_UPSTR_TOTAL, select=terms)
} else {
  red.SP_MEM_UPSTR_TOTAL<-data.frame(row.names = rownames(SP_MEM_UPSTR_TOTAL)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_MEM_UPSTR_NET
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_MEM_UPSTR_NET, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_MEM_UPSTR_NET)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_MEM_UPSTR_NET)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_MEM_UPSTR_NET<- subset(SP_MEM_UPSTR_NET, select=terms)
} else {
  red.SP_MEM_UPSTR_NET<-data.frame(row.names = rownames(SP_MEM_UPSTR_NET)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_UW
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_AEM_TDN_UW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_AEM_TDN_UW)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_AEM_TDN_UW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_UW<- subset(SP_AEM_TDN_UW, select=terms)
} else {
  red.SP_AEM_TDN_UW<-data.frame(row.names = rownames(SP_AEM_TDN_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DAMW
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_AEM_TDN_DAMW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_AEM_TDN_DAMW)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_AEM_TDN_DAMW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DAMW<- subset(SP_AEM_TDN_DAMW, select=terms)
} else {
  red.SP_AEM_TDN_DAMW<-data.frame(row.names = rownames(SP_AEM_TDN_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TDN_DISTW
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_AEM_TDN_DISTW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_AEM_TDN_DISTW)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_AEM_TDN_DISTW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TDN_DISTW<- subset(SP_AEM_TDN_DISTW, select=terms)
} else {
  red.SP_AEM_TDN_DISTW<-data.frame(row.names = rownames(SP_AEM_TDN_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_UW
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_AEM_TUP_UW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_AEM_TUP_UW)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_AEM_TUP_UW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_UW<- subset(SP_AEM_TUP_UW, select=terms)
} else {
  red.SP_AEM_TUP_UW<-data.frame(row.names = rownames(SP_AEM_TUP_UW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DAMW 
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_AEM_TUP_DAMW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_AEM_TUP_DAMW)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_AEM_TUP_DAMW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DAMW<- subset(SP_AEM_TUP_DAMW, select=terms)
} else {
  red.SP_AEM_TUP_DAMW<-data.frame(row.names = rownames(SP_AEM_TUP_DAMW)) 
}
model$anova
remove(mod0,mod1,model)

########### SP_AEM_TUP_DISTW
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= SP_AEM_TUP_DISTW, sqrt.dist = FALSE, add = FALSE)
sig<-anova(rda)
sig
if (sig$`Pr(>F)`[1] <0.1){
  mod0 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~1, data=SP_AEM_TUP_DISTW)
  mod1 <- rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data=SP_AEM_TUP_DISTW)
  model <-ordiR2step(mod0,mod1, Pin=0.1, permutations=how(nperm=99999), 
                     R2permutations=99999, direction="forward", R2scope = TRUE, trace=FALSE)
  model$anova
  terms <- attr(model$terminfo$terms,"term.labels")
  red.SP_AEM_TUP_DISTW<- subset(SP_AEM_TUP_DISTW, select=terms)
} else {
  red.SP_AEM_TUP_DISTW<-data.frame(row.names = rownames(SP_AEM_TUP_DISTW)) 
}
model$anova
remove(mod0,mod1,model)

########### GET SELECTED TERMS
SP_select<- cbind.data.frame(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, 
                             red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, 
                             red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
usdm::vif(SP_select)
### multiple regression test to confirm all terms are significantly related to response
rda<-rda(phylogenetic.alpha$Phylogenetic_Diversity~., data= SP_select)
RsquareAdj(rda)
anova(rda)
sig.terms<-anova(rda, by="term", permutations=99999)
sig.terms
### remove non-sig vars
terms <- as.data.frame(attr(rda$terminfo$terms,"term.labels"))
terms[,2]<-sig.terms$`Pr(>F)`[-(length(sig.terms$`Pr(>F)`-1))]
terms<-terms[terms$V2<0.05,]
SP_select<-subset(SP_select, select=terms$`attr(rda$terminfo$terms, "term.labels")`)
vif<-usdm::vifstep(SP_select, th=10)
SP_select<-usdm::exclude(SP_select, vif)
SP_select<-as.data.frame(SP_select[,-2])
remove(rda, sig.terms, terms, vif)
### ___ variation partitioning ---------------------------------------------------
### Partial redundancy analysis of models (taxon alpha)
### partial RDA of full environmental set
comb <- cbind(SP_select, ENV_select)
rda.all<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data=comb)
RsquareAdj(rda.all)
anova(rda.all)
#anova(rda.all, by="term", permutations=99999)
#### Partial Redundancy Analysis of selected model from procedure above
partial.reg <- varpart(phylogenetic.alpha$Phylogenetic_Diversity, SP_select, ENV_select)
partial.reg

#### TESTING SIGNIFICANCE OF FRACTIONS
#### fraction [a]:
rda.pure.spatial<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ . + Condition(as.matrix(ENV_select)), data=SP_select)
RsquareAdj(rda.pure.spatial)
anova(rda.pure.spatial)
anova(rda.pure.spatial, by="term", permutations=99999)
#### fraction [c]:
rda.pure.environ<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ . + Condition(as.matrix(SP_select)), data=ENV_select)
RsquareAdj(rda.pure.environ)
anova(rda.pure.environ)
anova(rda.pure.environ, by="term", permutations=99999)
#### fraction [a+b]:
rda.spatial<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~., data= SP_select)
RsquareAdj(rda.spatial)
anova(rda.spatial)
#### fraction [b+c]:
rda.envrion<-rda(phylogenetic.alpha$Phylogenetic_Diversity ~ ., data= ENV_select)
RsquareAdj(rda.envrion)
anova(rda.envrion)
anova(rda.envrion, by="term", permutations=99999)

#### clean global environment
remove(red.env.instream.PCs,red.env.hydrophysio.PCs,red.env.climate.PCs,red.env.landcover.PCs,red.env.geology.PCs,red.env.anthropogenic.PCs)
remove(red.SP_MEM_LAND, red.SP_MEM_HYDRO, red.SP_MEM_UPSTR_TOTAL, red.SP_MEM_UPSTR_NET, red.SP_AEM_TDN_UW, red.SP_AEM_TDN_DAMW, red.SP_AEM_TDN_DISTW, red.SP_AEM_TUP_UW, red.SP_AEM_TUP_DAMW, red.SP_AEM_TUP_DISTW)
remove(ENV_select, SP_select)
remove(rda,sig.terms,terms,mod0,mod1,model)
##########################          END           #####################################################
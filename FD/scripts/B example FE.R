rm(list=ls()) # cleaning memory

##############################################################################################
# sourcing the "homemade" R functions

# setting working directory
my_path<-"/Volumes/Data/tutorial_FD" #  <= PUT direction of the folder where you saved the zipfile containing functions and data

# loading functions
setwd( paste(my_path,"/functions", sep="") ) # folder where R functions have been saved
source("species_to_FE.R")
source("FE_metrics.R")
source("quality_funct_space.R")
source("plot_funct_space.R")
source("multidimFD.R")

##############################################################################################
# IMPORTING DATA that have already been imported in "script example fruits.R"

setwd( paste(my_path,"/results", sep="") ) # folder where outputs from script "A script example FD" have been saved
load("traits_fruits")
load("weight_fruits_baskets")
load("species_codes")

# folder where outputs will be saved
setwd( paste(my_path,"/results/FE", sep="") ) 

##############################################################################################
# raw trait matrix
summary(traits_fruits) # Sugar_content is the only continuous trait, not considered hereafter for convenience

# new trait matrix with no continuous traits
traits_fruits_Q<-traits_fruits[,1:4]
summary(traits_fruits_Q) # all traits are qualitative or semi-quantitative


# decreasing the number of modalities per trait for convenience (to have less unique combinations of trait values)

# size grouped into only 3 categories
traits_fruits_Q[,"Size"]<-as.character(traits_fruits_Q[,"Size"])
traits_fruits_Q[which(traits_fruits_Q[,"Size"]=="very_small"),"Size"]<-"small"
traits_fruits_Q[which(traits_fruits_Q[,"Size"]=="very_large"),"Size"]<-"large"
traits_fruits_Q[,"Size"]<-factor(traits_fruits_Q[,"Size"], levels=c("small", "medium", "large" ), ordered = TRUE )

# Plant type grouped into only 2 categories
traits_fruits_Q[,"Plant"]<-as.character(traits_fruits_Q[,"Plant"])
traits_fruits_Q[which(traits_fruits_Q[,"Plant"]!="tree"),"Plant"]<-"Not_tree"
traits_fruits_Q[,"Plant"]<-factor(traits_fruits_Q[,"Plant"], levels=c("Not_tree", "tree"), ordered = TRUE )

# Plant Origin grouped into only 2 categories
traits_fruits_Q[,"Origin"]<-as.character(traits_fruits_Q[,"Origin"])
traits_fruits_Q[which(traits_fruits_Q[,"Origin"]!="temperate"),"Origin"]<-"tropical"
traits_fruits_Q[,"Origin"]<-factor(traits_fruits_Q[,"Origin"], levels=c("temperate", "tropical"), ordered = TRUE )


# names of two seed type are similar=> replacing them for next steps
traits_fruits_Q[,"Seed"]<-factor(traits_fruits_Q[,"Seed"], levels=c("none","pip", "pit" ), labels = c("no", "pp", "pt" ), ordered = FALSE )

# checking the new trait database
summary(traits_fruits_Q)

# saving it
save(traits_fruits_Q, file="traits_fruits_Q")

##############################################################################################
# computing FUNCTIONAL ENTITIES = unique combinations of trait values

# grouping species into FE
species_to_FE_fruits<-species_to_FE(traits_fruits_Q)

# codes of FE
species_to_FE_fruits$FE_codes # the 25 species were grouped into 13 Funct entities

# number of species per FE
apply(species_to_FE_fruits$FE_sp_01,1,sum  ) # from 1 to 6 species per FE

# FE to which each fruit species belongs
FE<-species_to_FE_fruits$FE
FE_fruits_01<-species_to_FE_fruits$FE_sp_01

# trait values of the FE
FE_traits<-species_to_FE_fruits$FE_traits

# matrix of fruit occurence (0/1) in baskets
basket_fruits_occ<-weight_fruits_baskets
basket_fruits_occ[which(basket_fruits_occ>0)]<-1

##################################
# matrix of FE biomass in baskets
basket_FE_weight<-matrix(0, nrow(weight_fruits_baskets), nrow(FE_fruits_01), 
                            dimnames=list( row.names(weight_fruits_baskets), row.names(FE_fruits_01) ) ) # empty matrix
for (k in row.names(FE_fruits_01) ) # loop on FE
{
     sp_k<-names(which(FE_fruits_01[k,]==1))   
     if( length(sp_k)==1 ) {basket_FE_weight[,k]<-weight_fruits_baskets[,sp_k] } # if only one species in FE k
     if (length(sp_k)>1 ) {basket_FE_weight[,k]<-apply(weight_fruits_baskets[,sp_k],1,sum)  } # if more than 1 species in FE k
}# end of k

sum(weight_fruits_baskets)==sum(basket_FE_weight) # check total biomass kept constant

# matrix of FE occurence (0/1) in baskets
basket_FE_occ<-basket_FE_weight
basket_FE_occ[which(basket_FE_occ>0)]<-1

##################################
# computing diversity metrics based on Funct ent for the set of fruits baskets studied
baskets_FE_metrics<-FE_metrics(species_to_FE_fruits, basket_fruits_occ, check_species_pool=TRUE, 
  folder_plot="/Volumes/Data/tutorial_FD/results/FE/plot_FE_metrics", nm_asb_plot=row.names( basket_fruits_occ) )
round(baskets_FE_metrics,3)
# plots illustrating distribution of species into Funct entities are in the subfolder ".../results/plot_FE" in your working directory.

##############################################################################################

# computing MULTIDIMENSIONAL FUNCTIONAL DIVERISTY INDICES based on FE position in a functional space to assess 
# how weight is distributed in the functional space independently from packing of species into FE (which is assessed by metrics presented above)

# computing Functional space based on trait value of Functional entities (not based on species trait values 
# because we want to represent distance between combinations of trait values independently from their frequency among species, i.e. to give same weight to each FE whatever its number of species)
qual_funct_space_FE<-quality_funct_space(FE_traits, traits_weights=NULL, nbdim=10, metric="Gower", dendro=TRUE, plot="quality_funct_space_FE")
qual_funct_space_FE$meanSD # => best space has 3 dimensions + have a look to the "quality_funct_space_FE.jpeg" file in the ".../results/FE" folder

#  FE coordinates in the best space
coord_FE_3D<-qual_funct_space_FE$details_funct_space$mat_coord[,1:3]

# species coordinates in the space according to those of FE
coord_sp_3D<-coord_FE_3D[FE,]
row.names(coord_sp_3D)<-names(FE)

##################################
# FD indices

# computing FD indices according to species weights in the 5 first baskets
FD_baskets_sp<-multidimFD(coord_sp_3D, weight_fruits_baskets[c(1:5),], check_species_pool=FALSE, verb=TRUE  )

# computing FD indices according to FE weights in the 5 first baskets
FD_baskets_FE<-multidimFD(coord_FE_3D, basket_FE_weight[c(1:5),], check_species_pool=FALSE, verb=TRUE  )

##################################
# comparing values between species- and FE-based FD indices 

# for FRic that does not account for weight
round( rbind( FRic_sp=FD_baskets_sp[,"FRic"], FRic_FE= FD_baskets_FE[,"FRic"]) ,3)
# => identical by defintion since FRic accounts only for coordinates of the most extreme points 
# and most extreme FEs are made by most extreme species


# for metrics accounting for weight

# Average position along the first axis 
round( rbind( FIde_PC1_sp=FD_baskets_sp[,"FIde_PC1"], FIde_PC1_FE= FD_baskets_FE[,"FIde_PC1"]) ,3) 
# => no difference since summing weights of the species that belong to the same FE (and hence have the same coordinates in the FE-based space) do not change weighted-average coordinates for each axis

# Functional evenness
round( rbind( FEve_sp=FD_baskets_sp[,"FEve"], FEve_FE= FD_baskets_FE[,"FEve"]) ,3)
# => grouping species into FE changes the regularity in the distribution of weights among the nodes of the Minimum Spanning Tree 
# linking all points (species or FE). Here FEve computed with species is lower than FEve computed with FE for all baskets but 'basket 1'

# Functional specialization
round( rbind( FSpe_sp=FD_baskets_sp[,"FSpe"], FSpe_FE= FD_baskets_FE[,"FSpe"]) ,3)
# => species coordinates came from coordinates computed for FEs trait values, hence average position of FEs on each axis is 0 but since FEs group different number of species, average position of species on each axis is not 0. Eventually FSpe differs between the 2 approaches

# GENERAL COMMENT:
# computing FD indices based on species or FE weights are both correct, choice should be done according to question addressed :
#   - species-based indices are relevant to detect assembly rules since processes such as dispersion and competition act on species
#  - FE-based indices are relevant to assess the links between FD and ecosystem processes since redundant species are expected to have same ecological roles

##############################################################################################
# END of script
##############################################################################################


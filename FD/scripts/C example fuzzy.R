rm(list=ls()) # cleaning memory

##############################################################################################
# sourcing "homemade" R functions

# setting working directory
my_path<-"/Volumes/Data/tutorial_FD" #  <= PUT direction of the folder where you saved the zipfile containing functions and data


setwd( paste(my_path,"/functions", sep="") ) # folder where R functions have been saved
source("quality_funct_space_fromdist.R")


##############################################################################################
# IMPORTING DATASETS

# importing matrix with 5 functional traits of fruits (from script "A example FD"), renaming it
setwd( paste(my_path,"/results", sep="") ) # folder where R functions have been saved
load( "traits_fruits" )
traits5_fruits<-traits_fruits
summary(traits5_fruits)

# importing an additional trait "fruit use", i.e. proportion of fruit eaten raw, in cake or as jam, that was coded using 3 variables (% of consumption)
use_fruits<-read.table( paste(my_path,"/data/use_fruits.txt", sep="") , header=T, row.names = "Species")
summary(use_fruits) # 100% for each row

# checking that species names are the same in the two matrices
sum( row.names(traits5_fruits)== row.names(use_fruits) ) == nrow(traits5_fruits)

# preparing a single trait matrix using  library "ade4"
library(ade4)

# fuzzy coding of proportion of each use (sum of each row equals 1)
use_fruits_fuzz<-prep.fuzzy( as.data.frame(use_fruits), col.blocks=3, label="Use")
summary(use_fruits_fuzz)

# list of the 6 traits as an object of type "ktab"
traits6_fruits_ktab<- ktab.list.df(list( as.data.frame(traits5_fruits$Size), as.data.frame(traits5_fruits$Plant), as.data.frame(traits5_fruits$Origin), 
                                          as.data.frame(traits5_fruits$Seed), as.data.frame(traits5_fruits$Sugar_content), use_fruits_fuzz)  )
traits6_fruits_ktab

# computing Gower's distance between all species using function from Pavoine et al (2009)
# need to set type of each trait: Size (ordinal), 
Gower_fruits<- dist.ktab(traits6_fruits_ktab, type=c("O","N","O", "O", "Q","F"), option=c("scaledBYrange") )

Gower_fruits # a distance matrix without species names

# adding species names
Gower_fruits<-as.matrix(Gower_fruits)
row.names(Gower_fruits)<-row.names(traits5_fruits)
colnames(Gower_fruits)<-row.names(traits5_fruits)
Gower_fruits<-as.dist(Gower_fruits)
Gower_fruits # OK
##############################################################################################

# FUNCTIONAL SPACE

setwd( paste(my_path,"/results/plot_funct_space", sep="") )  # setting working directory for results

# computing all functional spaces based on dendrogram or PCoA (up to 10 axes)
# NOTE: you need to be connected to internet for a correct functioning of quality_funct_space function

qual_funct_space_6traits<-quality_funct_space_fromdist(Gower_fruits, nbdim=10, plot=NA)
qual_funct_space_6traits$meanSD # => 10 dimensional space is the best but 4D-space hs a meanSD<0.01

# species coordinates in the 4D space
coord_fruits_4D<-qual_funct_space_6traits$details_funct_space$mat_coord[,1:4]

##############################################################################################

# than see script "A example FD" to compute FD indices based on species coordinates in the functional space and species dominance in communities

##############################################################################################
# END
##############################################################################################
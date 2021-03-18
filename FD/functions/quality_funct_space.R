#################################################################################################################################
## quality_funct_space: R function for computing the quality of functional dendrogramm and multidimensional functional spaces  ##
##                                                                                                                             ##
##    This function is an updated version of the Appendix S1 associated to Maire et al. 2015 (Global Ecology and Biogeography) ##
##                                                                                                                             ##
##                  Given a functional matrix, the function computes the quality (i.e. mean squared-deviation between          ##
##                  initial functional distance and standardized distance in the functional space) for the best functional     ##
##                  dendrogram and all the multidimensional functional spaces from 2 to N dimensions (N selected by the user). ## 
##                  A graphical output illustrating the quality of each functional space is also provided.                     ##
##                                                                                                                             ##
##  Code by Eva Maire & Sébastien Villéger (sebastien.villeger@cnrs.fr). Last update on 2016-04-29                           ##
##  Best functional dendrogram is computed using the GFD function written by F. Guilhaumon (see Mouchet et al. 2008, Oikos)    ##
##                                                                                                                             ##
##                                                                                                                             ##
##     INPUTS:                                                                                                                 ##
##                                                                                                                             ##
##      - "mat_funct" : a species x functional traits matrix (NA are not allowed, at least 3 species and 3 traits)             ##
##                        Traits could be of different types (e.g. numeric, ordinal, nominal)                                  ##
##                                                                                                                             ##
##      - "traits_weights" : a numeric vector (NA are not allowed) with weights for all traits in "mat_funct"                  ##
##                           Applied only with Gower's distance. Default is same weight for all traits                         ##
##                                                                                                                             ##
##      - "nbdim" : maximum number of dimensions for multidimensional functional spaces. By default, nbdim=7                   ##
##      		Final number of dimensions depends on the number of positive eigenvalues obtained with PCoA                        ##
##                                                                                                                             ##
##      - "metric" : metric to be used to compute functional distance, "Euclidean" or "Gower" (=default)                       ##
##                                                                                                                             ##
##      - "dendro" : a logical value indicating whether the best functional dendrgoam shoudl be looked for (default is TRUE)   ##
##      		        Setting value to FALSE will save computation time, especially when working with >100 species               ##
##                                                                                                                             ##
##      - "plot" : character string to set the name of the jpeg file for plots illustrating the quality of functional spaces   ##
##                  NA means no plot                                                                                           ##
##                                                                                                                             ##
##     NB: 1/  high value for 'nbdim' can increase computation time                                                            ##
##         2/ if at least one trait is not numeric, 'metric' should be set as 'Gower'                                          ##
##         3/ if metric=Euclidean, functional traits are scaled (mean=0, sd=1) before computing functional distances           ##
##         4/ R libraries 'ape', 'clue', 'cluster', 'geometry', 'gtools' are required                                          ##
##                                                                                                                             ##
##                                                                                                                             ##
##      OUTPUTS: a list with                                                                                                   ##
##                                                                                                                             ##
##       - $ meanSD : a vector with mean squared deviation values for all the functional spaces tested							           ##
##            names are "t_clustering algorithm" for the best tree and 'm_kD' (with k=2:nbdim) for multidimensional spaces     ##
##                                                                                                                             ##
##       - $ details_funct_space : a list with details about functional spaces                                                 ##
##                                                                                                                             ##
##          - $ mat_dissim : matrix of functional dissimilarities between species based on their trait values                 ##
##                                                                                                                             ##
##          - $ mat_coord : coordinates of species in the nbdim multidimensional functional space (PCoA)                       ##
##                                                                                                                             ##
##       	  - $ alg_best_tree : name of the clustering method that produced the best tree                                      ##
##                                                                                                                             ##
##          - $ best_tree : the best tree obtained with the GFD function.                                                      ##
##                                                                                                                             ##
##          - $dist_raw and $dist_st : lists with raw and standardized distances between species in each functional space      ##
##                       names of these distance matrices are as names in meanSD (e.g. $dist_raw$t_UPGMA or $dist_raw$m_3D)    ##
##                                                                                                                             ##
##       - a jpeg file in the current working directory with :                                                                 ##
##                      - a barplot showing the meanSD for all functional spaces                                               ##
##                      - 'nbdim' panels (or only 15 if nbdim>15) illustrating the quality of each functional space            ##
##                        Points represent species pairs. Mean squared deviation (mSD) is provided at the top of each panel.	 ##
##                                                                                                                             ##
#################################################################################################################################

quality_funct_space <- function(mat_funct, traits_weights=NULL, nbdim=7, 
                                metric="Gower", dendro=TRUE, plot="quality_funct_space") 
  {

#loading required libraries
require(ape)
require(clue)
require(cluster)
require(geometry)
require(gtools)

# sourcing GFD function to compute the best dendrogram
source("http://villeger.sebastien.free.fr/R%20scripts/GFD_matcomm.R")  ; GFD<-GFD_matcomm
################################################################################################################################# 

# checking data
if (sum(is.na(mat_funct))!=0)   {  stop(" NA are not allowed in 'funct_mat' ")     }
if (nbdim<2)   {  stop(" 'nbdim' must be higher than 1")     }
if (ncol(mat_funct)<3)   {  stop(" there must be at least 3 traits in 'funct_mat' ")     }
if (nrow(mat_funct)<3)   {  stop(" there must be at least 3 species in 'funct_mat' ")     }

if (metric=="Euclidean" & sum(apply(mat_funct,2,is.numeric))!=ncol(mat_funct) )  { stop("using Euclidean distance requires that all traits are continuous")   }
if (metric=="Euclidean" & nbdim>ncol(mat_funct) )  { stop("when using Euclidean distance, number of dimensions tested should be lower than number of traits")   }
if (metric=="Gower" & nbdim>nrow(mat_funct) )  { stop("when using Gower's distance, number of dimensions tested should be lower than number of species")   }
if (metric=="Euclidean" & is.null(traits_weights)==FALSE )  { warning("traits weights were ignored when computing Euclidean distance")   }

if(is.null(traits_weights) ) {traits_weights<-rep(1,ncol(mat_funct) ) }  # if no traits weighting provided, setting same weight for all traits
if (is.numeric(traits_weights)==FALSE)  {  stop("'traits_weight' should be numeric ")     }
if (sum(is.na(traits_weights))!=0)   {  stop(" NA are not allowed in 'traits_weights' ")     }
if (length(traits_weights)!=ncol(mat_funct))   {  stop("length of 'traits_weight' should match number of columns of 'mat_funct' ")     }

################################################################
# computing functional dissimilarity between species given their trait values
if (metric=="Euclidean") mat_funct <-scale(mat_funct) # scaling if continuous traits
mat_dissim<-daisy(mat_funct, metric=tolower(metric), weights=traits_weights )

# lists to store distances matrices
dist_raw<-list()
dist_st<-list()
  
################################################################
# computing PCoA
mat_pcoa<-pcoa(mat_dissim)

# changing number of dimensions given number of positive eigenvalues
nbdim<-min(nbdim,ncol(mat_pcoa$vectors) )

# keeping species coordoinates on the 'nbdim' axes
mat_coord<-mat_pcoa$vectors[,1:nbdim]
row.names(mat_coord)<-row.names(mat_funct)
colnames(mat_coord)<-paste("PC",1:nbdim,sep="")

# computing Euclidean distances between species in the (nbdim-1) multidimensionnal functional spaces 
for (k in 2:nbdim) 
    {
    eval(parse(text=paste("dist_",k,"D<-dist(mat_coord[,1:",k,"],method='euclidean')", sep="")))
    eval(parse(text=paste("dist_raw$m_",k,"D<-dist_",k,"D", sep="")))
    } # end of k
  
################################################################
# computing the best functional dendrogram using the GFD function and keeping the cophenetic distances between species on this tree
alg_best_tree<-NA
best_tree<-NA
dist_best_tree<-NA
  if (dendro==TRUE)
  {
  best_tree_pool <- GFD( mat_funct, NA, distances=tolower(metric), dissimMethod="spectral")
  alg_best_tree<-best_tree_pool$methods
  best_tree<-best_tree_pool$tree
  dist_best_tree<-best_tree_pool$ultram
  }# end of if compute dendrograms
eval(parse(text=paste("dist_raw$t_",alg_best_tree,"<-dist_best_tree", sep="")))

################################################################
# computing mean squared deviation between initial distance and standardized final distance in the functional space
meanSD<-rep(NA,nbdim) ; names(meanSD)<-c(paste("t_",alg_best_tree,sep="") ,paste("m_",2:nbdim,"D",sep=""))

x<-mat_dissim # initial distance
S<-nrow(mat_funct) # species richness

if (dendro==TRUE)
{
  # for tree
  y<-dist_best_tree
  yst<- y/max(y) * max(x)
  eval(parse(text=paste("dist_st$t_",alg_best_tree,"<-yst", sep="")))
  meanSD[paste("t_",alg_best_tree,sep="")]<-round( ( (sum((x-yst)^2)) / (S*(S-1)/2) ) ,6)
}# end of if compute dendrograms


# for muldimensionnal spaces
for (k in 2:nbdim)  
  {
  eval(parse(text=paste("y<-dist_",k,"D",sep="")))
  yst<- y/max(y) * max(x)
  eval(parse(text=paste("dist_st$m_",k,"D<-dist_",k,"D", sep="")))
  meanSD[paste("m_",k,"D",sep="")]<-round( ( (sum((x-yst)^2)) / (S*(S-1)/2) ) ,6)
  }  # end of k

# list of outputs
res<-list(meanSD=meanSD, details_funct_space=list(mat_dissim=mat_dissim, mat_coord=mat_coord, 
              alg_best_tree=alg_best_tree, best_tree=best_tree, dist_raw=dist_raw, dist_st=dist_st )  )
  
################################################################################################################################
# GRAPHICS if plot has a name
################################################################################################################################

if (is.na(plot)==FALSE)
{
# window
if (nbdim<=3) {jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=600)
                layout(matrix(c(1:(nbdim+1),rep(0,3-nbdim)),1,4,T)) ; layout.show(nbdim+1)  }
if (nbdim>3 & nbdim<=7) {jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=1200)
                          layout(matrix(c(1:(nbdim+1),rep(0,7-nbdim)),2,4,T)) ; layout.show(nbdim+1) }
if (nbdim>7 & nbdim<=11) {jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=1800)
                          layout(matrix(c(1:(nbdim+1),rep(0,11-nbdim)),3,4,T)) ; layout.show(nbdim+1) }
if (nbdim>11 & nbdim<=15) {jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=2400)
                          layout(matrix(c(1:(nbdim+1),rep(0,15-nbdim)),4,4,T)) ; layout.show(nbdim+1) }
if (nbdim>15) { jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=2400)
  layout(matrix(1:16,4,4,T)) ; layout.show(16) ; meanSD_plot<-meanSD[1:15]}	

par(mar=c(4,4,3,3))
  
# change in meanSD with increasing number of dimensions  
barplot(height=meanSD,names.arg=names(meanSD), xlab="Functional space", ylab= "Quality (Mean SD)", 
        space=0, cex.names=0.7, col=c("red", rep("blue",nbdim-1) ) )

# quality of each functional space
x<-mat_dissim # functional distances

# dendrogram quality
if (dendro==TRUE)
  {
  eval(parse(text=paste("yst<-dist_st$t_",alg_best_tree, sep="")))
  plot(x,yst, xlab=paste(metric, "distance"), ylab= "Cophenetic distance", xlim=c(0,max(x)), ylim=c(0,max(yst)), 
        pch=21, col="red", bg="red", cex=0.3, cex.axis=0.8, cex.lab=0.9 )	
  abline(a=0,b=1)
  title(main=paste(alg_best_tree, "   mSD=",round(meanSD[paste("t_",alg_best_tree,sep="")],4),sep=""), cex.main=1.1, col.main="red", line=0.5   )
  }# end of if compute dendrograms

# multidimensional spaces quality
for (k in 2:min(nbdim,15) )  
  {
  eval(parse(text=paste("yst<-dist_st$m_",k,"D",sep="")))
  plot(x,yst, xlab=paste(metric, "distance"), ylab= "Euclidean distance", xlim=c(0,max(x)), ylim=c(0,max(yst)), 
          pch=21, col="blue", bg="blue", cex=0.3, cex.axis=0.8, cex.lab=0.9  )	
  abline(a=0,b=1)
  title(main=paste(paste(k,"D",sep=""),"   mSD=",round(meanSD[paste("m_",k,"D",sep="")],4),sep=""), cex.main=1.1, col.main="blue", line=0.5  )
  }  # end of k	


graphics.off()
  
} # end of of plot
 
################################################################################################################################
################################################################################################################################
   
invisible(res)

} # end of function quality_funct_space

################################################################################################################################
################################################################################################################################



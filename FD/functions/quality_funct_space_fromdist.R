#################################################################################################################################
## R function for computing the quality of functional dendrogramm and multidimensional functional spaces                       ##
##                                                                                                                             ##
##    This function is a simplified version of the Appendix S1 associated to Maire et al. 2015 (Global Ecol. and Biogeogr.)    ##
##        i.e. functional distance matrix is provided by user and only UPGMA tree is computed                                  ##
##                                                                                                                             ##
##                  Given a functional distance matrix, the function computes the quality                                      ## 
##                  (i.e. mean squared-deviation between initial functional distance and standardized distance in the          ##
##                  functional space) for UPGMA dendrogram and all the multidimensional functional spaces                      ##
##                  from 2 to N dimensions (N selected by the user).                                                           ## 
##                  A graphical output illustrating the quality of each functional space is also provided.                     ##
##                                                                                                                             ##
##                                                                                                                             ##
##     INPUTS:                                                                                                                 ##
##                                                                                                                             ##
##      - "dist_funct" : a species x species functional distance matrix. Object should be of class 'dist', for instance output ##
##                       of function 'daisy' (package 'cluster') or of function 'dist.ktab' (package ade4) that computes       ##
##                        generalized Gower's distance accounting for all types of traits (e.g. fuzzy-coded)                   ##
##                                                                                                                             ##
##      - "nbdim" : maximum number of dimensions for multidimensional functional spaces. By default, nbdim=7                   ##
##      		Final number of dimensions depends on the number of positive eigenvalues obtained with PCoA                        ##
##                                                                                                                             ##
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
##          - $ mat_coord : coordinates of species in the nbdim multidimensional functional space (PCoA)                       ##
##                                                                                                                             ##
##          - $ upgma_tree : dendrogram built with the UPGMA clustering algorithm                                              ##
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

quality_funct_space_fromdist <- function( dist_funct,  nbdim=7,   plot="quality_funct_space") 
{

#loading required libraries
require(ape)
require(clue)
require(cluster)
require(geometry)
require(gtools)

################################################################################################################################# 

  # checking data
  if ( ("dist" %in% class(dist_funct) ) ==FALSE )   {  stop(" 'dist_funct' must be of class 'dist'")     }
  if (length(dist_funct)<3)   {  stop(" there must be at least 3 species in 'dist_funct' ")     }
  if (sum(is.na(dist_funct))!=0)   {  stop(" NA are not allowed in 'dist_funct' ")     }
  if (nbdim<2)   {  stop(" 'nbdim' must be higher than 1")     }
  
  # functional distance 
  mat_dissim<-dist_funct
  
  # species names
  nm_sp<-row.names(as.matrix(dist_funct))
  
################################################################
# computing PCoA
mat_pcoa<-pcoa(mat_dissim)

# changing number of dimensions given number of positive eigenvalues
nbdim<-min(nbdim,ncol(mat_pcoa$vectors) )

# keeping species coordoinates on the 'nbdim' axes
mat_coord<-mat_pcoa$vectors[,1:nbdim]
row.names(mat_coord)<-nm_sp
colnames(mat_coord)<-paste("PC",1:nbdim,sep="")

# lists to store distance matrices
dist_raw<-list()
dist_st<-list()

# computing Euclidean distances between species in the (nbdim-1) multidimensionnal functional spaces 
for (k in 2:nbdim) 
    {
    eval(parse(text=paste("dist_",k,"D<-dist(mat_coord[,1:",k,"],method='euclidean')", sep="")))
    eval(parse(text=paste("dist_raw$m_",k,"D<-dist_",k,"D", sep="")))
    } # end of k
  
################################################################
# computing an UPGMA-dendrogram then  cophenetic distances between species on this tree
alg_best_tree<-"UPGMA"
best_tree<-NA
dist_best_tree<-NA

    best_tree <-hclust( mat_dissim , method="average")
    
    dist_best_tree<-cl_ultrametric(best_tree)
  
eval(parse(text=paste("dist_raw$t_",alg_best_tree,"<-dist_best_tree", sep="")))

################################################################ 
# computing mean squared deviation between initial distance and standardized final distance in the functional space
meanSD<-rep(NA,nbdim) ; names(meanSD)<-c(paste("t_",alg_best_tree,sep="") ,paste("m_",2:nbdim,"D",sep=""))

x<-mat_dissim # initial distance
S<-length(nm_sp) # species richness

  # for tree
  y<-dist_best_tree
  yst<- y/max(y) * max(x)
  eval(parse(text=paste("dist_st$t_",alg_best_tree,"<-yst", sep="")))
  meanSD[paste("t_",alg_best_tree,sep="")]<-round( ( (sum((x-yst)^2)) / (S*(S-1)/2) ) ,6)


# for muldimensionnal spaces
for (k in 2:nbdim)  
  {
  eval(parse(text=paste("y<-dist_",k,"D",sep="")))
  yst<- y/max(y) * max(x)
  eval(parse(text=paste("dist_st$m_",k,"D<-dist_",k,"D", sep="")))
  meanSD[paste("m_",k,"D",sep="")]<-round( ( (sum((x-yst)^2)) / (S*(S-1)/2) ) ,6)
  }  # end of k

# list of outputs
res<-list(meanSD=meanSD, details_funct_space=list( mat_coord=mat_coord, upgma_tree=best_tree, dist_raw=dist_raw, dist_st=dist_st )  )
  
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
  
# plotting change in meanSD with increasing number of dimensions  
barplot(height=meanSD,names.arg=names(meanSD), xlab="Functional space", ylab= "Quality (Mean SD)", 
        space=0, cex.names=0.7, col=c("red", rep("blue",nbdim-1) ) )

# plotting quality of each functional space

# functional distances
x<-as.dist( mat_dissim) 

# dendrogram quality
  eval(parse(text=paste("yst<-dist_st$t_",alg_best_tree, sep="")))
  plot(x,yst , xlab="Initial distance", ylab= "Cophenetic distance", xlim=c(0,max(x)), ylim=c(0,max(yst)), 
        pch=21, col="red", bg="red", cex=0.3, cex.axis=0.8, cex.lab=0.9 )	
  abline(a=0,b=1)
  title(main=paste(alg_best_tree, "   mSD=",round(meanSD[paste("t_",alg_best_tree,sep="")],4),sep=""), cex.main=1.1, col.main="red", line=0.5   )

# multidimensional spaces quality
for (k in 2:min(nbdim,15) )  
  {
  eval(parse(text=paste("yst<-dist_st$m_",k,"D",sep="")))
  plot(x,yst, xlab="Initial distance", ylab= "Euclidean distance", xlim=c(0,max(x)), ylim=c(0,max(yst)), 
          pch=21, col="blue", bg="blue", cex=0.3, cex.axis=0.8, cex.lab=0.9  )	
  abline(a=0,b=1)
  title(main=paste(paste(k,"D",sep=""),"   mSD=",round(meanSD[paste("m_",k,"D",sep="")],4),sep=""), cex.main=1.1, col.main="blue", line=0.5  )
  }  # end of k	


graphics.off()
  
} # end of of plot
 
################################################################################################################################
################################################################################################################################
   
invisible(res)

} # end of function quality_funct_space_fromdist

################################################################################################################################
################################################################################################################################



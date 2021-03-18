###########################################################################################################################################
#
# 'multidimFbetaD': function to compute and illustrate multdimensional functional beta-diversity indices for pairs of species assemblages
#
# For details about indices formulas see Villéger et al. 2013, Global Ecology and biogeography (22:671-681).
# This function uses R library 'betapart' for indices computation, see its help for more details.
#                 
# INPUTS: 
#       - 'coord': a matrix with coordinates of at least 3 species (rows) along 2,3 or 4 axes (columns) of a functional space
#		    - 'occ': a matrix with presence/absence (coded as 0/1) of species (columns) in at least 2 assembalges (rows). 
#       - 'check_species_pool" : a logical value indicating wether two tests are performed to check that all columns of 'occ' have a sum stricly positive
#                               (i.e. all species should be present in at least one assemblage). Default is TRUE.
#		    - 'verb': a logical value (default=TRUE) indicating whether printing progress of computation process
#
#     NB :  
#           Column names of 'occ' should be the same than row names of 'coord' (i.e. same codes and same order).
#           If 'coord' has no names, default names will be set ('Axis_1', 'Axis_2',...,'Axis_n').
#           'occ' should have at least 2 rows and row names
#           NA are not allowed in 'coord' or in 'occ'.
#           Only 0/1 are allowed in 'occ'.
#           All assemblages should have more species than number of axes so that convex-hull-based indices could be computed.
#           Computation time when there are 4 axes and many species could be high. This is also why it is not possible to work with 5-dimensional spaces.
#
# 		  - 'nm_asb_plot': a vector with the names of at least 2 assemblages for which FD indices need to be illustrated. 
#                       Default is NULL (i.e. no plot). To plot FD of all assemblages, set to 'row.names(occ)'.
#       - 'folder_plot': a character string for setting the working directory where plots will be saved, default is current working directory.
#       - 'Faxes_plot': a vector with names of axes to plot. Should be of length from 2 to 4. 
#               Default is NULL which means that the four first axes will be plotted.
#               Value of 'Faxes_plot' does not affect the way FD indices are computed (i.e. acccording to all columns of 'coord').
#       - 'Faxes_nm_plot': a vector with titles of axes (default is titles as 'Faxes_plot') 
# 		  - 'plot_pool': a logical value indicating whether all species present in 'coord' need to be illustrated on all plots.
#                      If yes (default value), space filled by the pool of species is in white while background is filled 
#                         with color specified in 'col_bg' (default is light grey), and absent species are plotted with a different symbol (see below).
#       - 'col_bg': a R color name. See above.
#       - "pch_sp_pool", "cex_sp_pool", 'col_sp_pool': Shape, size and color of symbol to plot absent species. See above.
#       - 'pch_sp': two numeric value coding shape of symbol, as in function 'points' (default is square and dot), to plot species of each assemblage.
#       - 'col_sp': two character string with hexadecimal code (e.g. from www.colorpicker.com) for color used for symbols and convex hull of each assemblage.
#       - 'cex_sp': two numeric value coding size of symbol for each assemblage, default is dot slighlty smaller than square (see 'pch_sp') to see both symbols
# 		  - 'transp': a single numeric value indicating the percentage of transparency for convex hull filling. Default is 50%.
#
#
# OUTPUTS: 
#     => a list of 4 dist objects with values of functional beta-diversity indices (based on convex hulls overlap) for all pairs of assemblages
#       - 'F_beta': Jaccard-like dissimilarity index 
#       - 'F_turn': turnover component of the Jaccard-like dissimilarity index 
#       - 'F_nest_res': nestedness-resultant component of the Jaccard-like dissimilarity index 
#       - 'F_pturn': contribution of turonver ot fucntinola beta-diversity (F_pturn=F_turn/F_beta)
#
#   => for each of the assemblage listed in 'nm_asb_plot': a jpeg file named 'beta_AssemblageA_AseemblageB.jpeg' 
#         illustrating Functional beta-Diversity indices for assemblages 'A' and 'B'
#         File has 1, 3 or 6 panels illustrating all combinations of 2, 3 or 4 axes, respectively.
#
#       All plots (i.e. all assemblages and all pairs of axes) have the same axis scale to faithfully represent FD.
#       Species present in the assemblage are shown with 'pch_sp' symbol. 
#			  The colored convex polygons are projection of the multidimensional convex hulls of the two assemblages in 2D. 
#       Filled symbols are species being vertices in the multidimensional space. 
#      
##########################################################################################################################################

multidimFbetaD<-function( coord,  occ, check_species_pool=TRUE, verb=TRUE,
                          folder_plot=NULL, nm_asb_plot=NULL, Faxes_plot=NULL, Faxes_nm_plot=NULL, 
                          plot_pool=TRUE, col_bg="grey90", col_sp_pool="grey30", pch_sp_pool="+", cex_sp_pool=1,
                          pch_sp=c(22,21), col_sp=c("#504AE8","#FA1900"), cex_sp=c(1.5,1.2), transp=50  ) 
  {

  # library required for indices computation
  require (geometry)
  require(betapart)
  
  # saving name of current working directory
  current_wd<-getwd()
  
  ##############################################################################
  # checking inputs
  
  # coordinates in the functional space
  if( nrow(coord)<3 ) stop(paste(" error: there must be at least 2 species in the dataset"))
  if( ncol(coord)<2 ) stop(paste(" error: there must be at least 2 functional axes"))
  if( ncol(coord)>4 ) stop(paste(" error: there must be less than 4 functional axes "))
  if ( is.null(colnames(coord)) ) { colnames(coord)<-paste("Axis", 1:ncol(coord),sep="_") }   # if no column names in 'coord' default value
  if( is.numeric(coord)==FALSE ) stop(paste(" error: 'coord' is not numeric"))
  if( is.na(sum(coord)) ) stop(paste(" error: NA in 'coord'"))
  
  # occurences of species
  if( is.matrix(occ)==FALSE ) stop( " 'occ' should be an object of type matrix")
  if( nrow(occ)<2 ) stop(paste(" error: there must be at least 2 assemblages"))
  if( is.numeric(occ)==FALSE ) stop(paste(" error: 'occ' is not numeric"))
  if( is.na(sum(occ)) ) stop(paste(" error: NA in 'occ'"))
  if( length(which(occ==0 | occ==1) )!=( nrow(occ)*ncol(occ) )  ) stop(paste(" error: only 0 or 1 allowed in 'occ'"))
  if(min(apply(occ,1,sum))==0 ) 
    stop(paste(" error: all rows of 'occ' should have a sum stricly positive, i.e. all assemblage must host at least one species"))
  
  # match between datasets
  if( sum(colnames(occ) == row.names(coord))!= nrow(coord) ) stop(paste(" error: 'occ' does not have the same column names than row names of 'coord'"))
  if (length(which( apply(occ,1,sum)<=ncol(coord) )>0 ) ) stop(paste(" error: all assemblages should have more species than number of functional axes"))
  
  # checking graphical parameters
  if( sum( nm_asb_plot %in% row.names(occ) ) != length(nm_asb_plot) ) stop( paste(" error: 'nm_asb_plot' should be subset of 'coord' row names") )

  # checking species pool
  if (check_species_pool==TRUE)
  {
    if(min(apply(occ,2,sum))==0 ) 
      stop(paste(" error: all columns of 'occ' should have a sum stricly positive, i.e. all species must occur in at least one assemblage"))
  }# end of check species pool
  
  # checking graphical parameters
  if( length(pch_sp)!=2 ) stop(paste(" error:'pch_sp' should contain 2 values"))
  if( length(col_sp)!=2 ) stop(paste(" error:'col_sp' should contain 2 values"))
  if( length(cex_sp)!=2 ) stop(paste(" error:'cex_sp' should contain 2 values"))
  if( length(col_bg)!=1 ) stop(paste(" error:'col_bg' should contain only one value"))
  if( length(col_sp_pool)!=1 ) stop(paste(" error:'col_sp_pool' should contain only one value"))
  if( length(pch_sp_pool)!=1 ) stop(paste(" error:'pch_sp_pool' should contain only one value"))
  if( length(cex_sp_pool)!=1 ) stop(paste(" error:'cex_sp_pool' should contain only one value"))
  if( length(transp)!=1 ) stop(paste(" error:'transp' should contain only one value"))
  
  ##############################################################################
  # preliminary computation at the species pool level
  
  # number and names of axes
  nm_axes<-colnames(coord)
  nb_axes<-length(nm_axes)
  
  # number and names of assemblages
  nm_asb<-row.names(occ)
  nb_asb<-length(nm_asb)
  
  # end of preliminary computation
  
  ######################################################################################

  
  # computing convex hulls and their intersections for all pairs of assemblages
  F_betapart_core<-functional.betapart.core(x=occ , traits=coord, multi=FALSE, return.details=TRUE)
  names(F_betapart_core$details$CH$coord_vertices)<-names(F_betapart_core$details$CH$FRi) # names of assemblages in sublist of matrices of vertices
  
  # computing functional beta diversity indices for 
  F_beta_pair<-functional.beta.pair( F_betapart_core, index.family="sorensen")
  F_beta_pair
  
  # storing results
  F_beta<-F_beta_pair$funct.beta.sor
  F_turn<-F_beta_pair$funct.beta.sim
  F_nest_res<-F_beta_pair$funct.beta.sne
  F_pturn<-F_turn/F_beta
  FbetaD<-list( F_beta=F_beta, F_turn=F_turn, F_nest_res=F_nest_res, F_pturn=F_pturn)
  
  # printing step achieved
  if (verb==TRUE) print( "Funct. beta-diversity indices have been computed for all pairs of assemblages " )
  ######################################################################################

  # graphical outputs if needed
  if ( is.null(nm_asb_plot)==FALSE  )
  {
    ##############################################################################
    # setting same graphical parameters for all assemblages
    
    # setting working directory to store jpeg files
    if (is.null(folder_plot) ) {  folder_plot<-current_wd }
    
    # setting folder for saving jpeg files
    test_folder_plot<-try( setwd(folder_plot) , silent=TRUE)
    if ( class(test_folder_plot) =="try-error") {
      folder_plot<-current_wd 
      print(paste(" /!\    WARNING: '",folder_plot," does not exist', jpeg files have been saved in '",current_wd, "'",sep="" ))
    } # end of if
    setwd(folder_plot)
    
    # range of 'coord' extended by 5%
    rge_coord<-range(coord)
    extrge_coord<-c( rge_coord[1]-0.05*(rge_coord[2]-rge_coord[1]) , rge_coord[2]+0.05*(rge_coord[2]-rge_coord[1]))
    
    # pretty axes labels within observed range
    lab<-pretty(extrge_coord, n=4) # labels
    lab<-lab[which(lab>=rge_coord[1] & lab<=rge_coord[2]) ]# filtering
    
    # axes limits a,d range
    Faxes_lim<-extrge_coord
    rge_Faxes_lim<-Faxes_lim[2]-Faxes_lim[1]
    
    # setting margins size of plot
    par( pty="s", mar=c(4,4,2,2) ) 
    
    ###########################################
    # function to draw functional space given axes limits and background color, option: lengend for species occs
    functional_space<-function(axes_xy, nm_axes_xy, col_bg, plot_pool=plot_pool) 
    {
      # setting margins size of plot
      par( pty="s", mar=c(3,3,3,3) ) 
      
      plot(Faxes_lim,Faxes_lim,type="n",axes=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=Faxes_lim,ylim=Faxes_lim) # window
      rect(Faxes_lim[1],Faxes_lim[1],Faxes_lim[2],Faxes_lim[2], col="white")   # background
      
      # customized X and Y axes
      for (k in 1:2)
      {
        axis(side=k, at=lab, labels=F, tcl=-0.3, pos=Faxes_lim[1])  # ticks
        mtext(side=k, lab, at=lab, line=-0.2, cex=0.9, las=1) # labels
        mtext(side=k, nm_axes_xy[k], cex=1,line=1.5, font=2) # title  
      } # end of k
      
      
      if (plot_pool==TRUE)
      {
        # grey background
        rect(Faxes_lim[1],Faxes_lim[1],Faxes_lim[2],Faxes_lim[2], col=col_bg)   
        
        # convex hull of species pool
        vert0<-convhulln( coord[,axes_xy] ,"Fx TO 'vert.txt'")
        vert1<-scan("vert.txt",quiet=T)
        vert_ij<-(vert1+1)[-1]
        polygon(coord[vert_ij,axes_xy], border=NA,col="white")   
        
        # all species
        points( coord[,axes_xy[1] ], coord[,axes_xy[2] ] , pch=pch_sp_pool, col=col_sp_pool, cex=cex_sp_pool)
        
      }# end of if plot_pool
      
      
    }# end of function functional_space
    ###########################################
    # axes to plot if not specified: up to first 4 axes 
    if ( is.null(Faxes_plot) ) { Faxes_plot<-colnames(coord)[ 1:(min(c(ncol(coord),4)))]  }
  
    # if not specififed default axes names
    if ( is.null(Faxes_nm_plot) ) { Faxes_nm_plot<-Faxes_plot }
  
    # checking inputs

    if( length(nm_asb_plot)<2 ) stop(paste(" error: 'nm_asb_plot' should contain names of at least 2 assemblages"))
    
    if( sum( nm_asb_plot %in% row.names(occ) ) != length(nm_asb_plot) ) stop(paste(" error: 'nm_asb_plot' should be subset of 'coord' row names"))
    
    if( length(Faxes_plot) <2 | length(Faxes_plot) >4) stop(paste("length of 'Faxes_plot' should be 2, 3 or 4 "))
    if( sum( Faxes_plot %in% colnames(coord) ) != length(Faxes_plot) ) stop(paste(" error: 'Faxes_plot' should be subset of 'coord' column names"))
    if( length(Faxes_plot) != length(Faxes_nm_plot) ) stop(paste("length of 'Faxes_plot' should match length of 'Faxes_nm_plot' "))
  
    # checking graphical parameters
    if( length(pch_sp)!=2 ) stop(paste(" error:'pch_sp' should contain two values"))
    if( length(col_sp)!=2 ) stop(paste(" error:'col_sp' should contain two values"))
    if( length(cex_sp)!=2 ) stop(paste(" error:'cex_sp' should contain two values"))
    if( length(col_bg)!=1 ) stop(paste(" error:'col_bg' should contain two values"))
    if( length(col_sp_pool)!=1 ) stop(paste(" error:'col_sp_pool' should contain only one value"))
    if( length(pch_sp_pool)!=1 ) stop(paste(" error:'pch_sp_pool' should contain only one value"))
    if( length(cex_sp_pool)!=1 ) stop(paste(" error:'cex_sp_pool' should contain only one value"))
    
    
    # shortening object names
    Faxes<-Faxes_plot
    Faxes_nm<-Faxes_nm_plot
    
    # number of axes
    nb_Faxes<-length(Faxes_plot)

    # number of assemblages
    nb_asb_plot<-length(nm_asb_plot)
    
    ################################
    # loop on pairs of assemblages
    for (i in 1:(nb_asb_plot-1) )
    for (j in (i+1):nb_asb_plot )
    {
        # species names and coordinates for each assemblages
        occ_i<-occ[nm_asb_plot[i],]
        nm_sp_i<-row.names(coord)[which(occ_i>0)]
        coord_sp_i<-coord[nm_sp_i,]
        
        occ_j<-occ[nm_asb_plot[j],]
        nm_sp_j<-row.names(coord)[which(occ_j>0)]
        coord_sp_j<-coord[nm_sp_j,]
        
        
        # creating jpeg file with 
        nm_jpeg<-paste("beta_",nm_asb_plot[i] ,"_",nm_asb_plot[j],".jpeg",sep="")
        
        # size of jpeg file and organisation of panels according to number of functional axes to plot
        if(nb_Faxes==2 ) { jpeg(file=nm_jpeg, res=300, width=2400, height=1200)
          layout(matrix(c(1:2), 1,2,F)) ; layout.show(2)  ; cext<-1.2 } # end of if 2 axes
        
        if(nb_Faxes==3 ){ jpeg(file=nm_jpeg, res=300, width=2400, height=2400)
          layout(matrix(c(1:2,4,3),2,2,F)) ; layout.show(4) ; cext<-1.5 } # end of if 3 axes
        
        if(nb_Faxes==4 ) { jpeg(file=nm_jpeg, res=300, width=3600, height=3600)
          layout(matrix(c(1:3, 0,4:5 ,7,0,6),3,3,F)) ; layout.show(7) ; cext<-2   } # end of if 4 axes
        
        # setting margins size of plot
        par( pty="s", mar=c(3,3,3,3) ) 
        
            # loop on pairs of axes
            for (x in 1:(nb_Faxes-1) )
            for (y in (x+1):nb_Faxes )
            {
            functional_space( c(Faxes[x] ,Faxes[y]), c(Faxes_nm[x] ,Faxes_nm[y]), col_bg=col_bg, plot_pool=plot_pool )
              
            # projected convex hulls in 2D
            vert0<-convhulln(coord_sp_i[,Faxes[c(x,y)]],"Fx TO 'vert.txt'")
            vert1<-scan("vert.txt",quiet=T) ; vertices2D_i<-(vert1+1)[-1]
            polygon(coord_sp_i[vertices2D_i,Faxes[c(x,y)]],border=NA, col=paste(col_sp[1],transp,sep=""))
          
            vert0<-convhulln(coord_sp_j[,Faxes[c(x,y)]],"Fx TO 'vert.txt'")
            vert1<-scan("vert.txt",quiet=T) ; vertices2D_j<-(vert1+1)[-1]
            polygon(coord_sp_j[vertices2D_j,Faxes[c(x,y)]],border=NA, col=paste(col_sp[2],transp,sep=""))
          
          
            # all points (empty) then filling points being vertices in nD
            points(coord_sp_i[,Faxes[x]], coord_sp_i[,Faxes[y]], pch=pch_sp[1], bg="white", col=col_sp[1],cex=cex_sp[1])
            nm_vertices_i<-row.names( F_betapart_core$details$CH$coord_vertices[[ nm_asb_plot[i] ]] )
            points(coord_sp_i[nm_vertices_i,Faxes[x] ], coord_sp_i[nm_vertices_i,Faxes[y]], pch=pch_sp[1], bg=col_sp[1], col=col_sp[1],cex=cex_sp[1])
          
            points(coord_sp_j[,Faxes[x]], coord_sp_j[,Faxes[y]], pch=pch_sp[2], bg="white", col=col_sp[2],cex=cex_sp[2])
            nm_vertices_j<-row.names( F_betapart_core$details$CH$coord_vertices[[ nm_asb_plot[j] ]] )
            points(coord_sp_j[nm_vertices_j,Faxes[x] ], coord_sp_j[nm_vertices_j,Faxes[y]], pch=pch_sp[2], bg=col_sp[2], col=col_sp[2],cex=cex_sp[2])
          
            }# end of x,y (pair of axes)
    
    # indices values  in the last panel  
    plot(1:5,1:5,type="n", axes=F,xlab="", ylab="") # empty plot
    text( 1, 4,  paste("Funct. beta-div.=", round( as.matrix(F_beta)[ nm_asb_plot[i],  nm_asb_plot[j] ],3) , sep=""), cex=cext, adj=0 )     
    text( 1, 3,  paste("Funct. turn.=", round( as.matrix(F_turn)[ nm_asb_plot[i],  nm_asb_plot[j] ],3) , sep=""), cex=cext, adj=0 )     
    text( 1, 2,  paste("Funct. nest.-res.=", round( as.matrix(F_nest_res)[ nm_asb_plot[i],  nm_asb_plot[j] ],3) , sep=""), cex=cext, adj=0 )     

    # closing jpeg
    graphics.off()
    
    # printing step achieved
    if (verb==TRUE) print( paste("plotting Funct. beta-diversity for assemblages: ", i , " and ", j, sep="") )
    
    }# end of working on pair of assemblages i,j
    ################################
    
  } # end of if plot of Functional beta-diversity
  ######################################################################################  
  # returning to current working directory
  setwd(current_wd)
  
  # returning results	
  return(FbetaD)	
  
} # end of function multidimFbetaD
########################################################################################################################################
########################################################################################################################################

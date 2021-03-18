##################################################################################################################################################################
# 'plot_funct_space': R function to illustrate species positions in a functional space
#						
# INPUTS:
#     coord: a numeric matrix with species (as rows) coordinates on functional axes, i.e. traits or PC axes from PCA or PCoA (as columns).
#     Faxes: a vector with names of axes to plot. Should be of length from 2 to 4. 
#               Default is NULL which means that the four first axes will be plotted.
#     Faxes_nm: a vector with labels of axes (default means as Faxes) 
#     Faxes_lim: a vector with minimum and maximum for values for axes. 
#                Default is c(NA,NA) which means that range is according to range of values among all axes.
#                 To have a fair representation of position of species in all plots, all axes have the same range.
#
#     convex_pool: a logical value (default=FALSE) indicating whether the convex hull filled by all species should be shown. 
#                 If yes, space filled by the pool of species is in white and background is filled with color specified in 'col_bg' (default is light grey) 
#
#     col_bg: a R color name. See above
#     pch_sp: symbol to use to show position of species. Numeric value means symbol shape as in function 'points' (default is 21, i.e. points).
#             A vector of characters of same length than 'length of 'coord' columns means species are illustrated with text (e.g. species names)
#     cex_sp: size of symbol of text
#     col_sp: color of symbol or text (as in 'points' and 'text' functions)
#
#
#     nm_jpeg: a character to set the name of the jpeg file, without the '.jpeg' extension.
#     close_jpeg: a logical value (default is TRUE) indicating whether jpeg file should be finalized. 
#                   If set to FLASE it is then possible to add a legend on the last plot (top-right corner )
#
# OUTPUTS: a jpeg file with plots showing the position of species in all 2D spaces made by pairs of axes. 
#           An empty plot is available on top-right corner to plot information, e.g. legend of points colors.
#
# NB: library "geometry" is required for plotting convex hulls
#
##################################################################################################################################################################


plot_funct_space<-function( coord , Faxes=NULL, Faxes_nm=Faxes, Faxes_lim=c(NA, NA), convex_pool=FALSE,
                            col_bg="grey90", pch_sp=21, cex_sp=0.8, col_sp="blue",  nm_jpeg="Fspace_.jpeg", close_jpeg=TRUE )

  {
  
  # if no column names in 'coord' default value
  if ( is.null(colnames(coord)) ) { colnames(coord)<-paste("Axis", 1:ncol(coord),sep="_") } 

  # if no axes selected, default is to keep the first 4 axes
  if ( is.null(Faxes)  ) { Faxes<- colnames(coord)[ 1:(min(c(ncol(coord),4)))] ; Faxes_nm<-Faxes }

  # checking inputs
  if( length(Faxes) <2 | length(Faxes) >4 ) 
    stop( "Check number of axes to plot: it should be from 2 to 4 axes")
  if ( sum( Faxes %in% colnames(coord) ) != length(Faxes) )  
    stop("Check names of axes to plot: they should be present in 'coord' column names")
  if ( length(Faxes) != length(Faxes_nm) )  
    stop("Check 'Faxes_nm': it should have the same length than 'Faxes'")
  if ( (length(pch_sp) != 1) & ( length(pch_sp) !=nrow(coord) ) )  
    stop("Check 'pch_sp': it should have a single value or the same length than number of rows in 'coord' ")
  
  
  # size of jpeg file and organisation of panels according to number of functional axes to plot
  nb_Faxes<-length(Faxes)
  if(nb_Faxes==2 ) { jpeg(file=paste(nm_jpeg,".jpeg",sep=""), res=300, width=2400, height=1200)
                      layout(matrix(c(1:2), 1,2,F)) ; layout.show(2)    } # end of if 2 axes
    
  if(nb_Faxes==3 ){ jpeg(file=paste(nm_jpeg,".jpeg",sep=""), res=300, width=2400, height=2400)
                      layout(matrix(c(1:2,4,3),2,2,F)) ; layout.show(4)    } # end of if 3 axes
    
  if(nb_Faxes==4 ) { jpeg(file=paste(nm_jpeg,".jpeg",sep=""), res=300, width=3600, height=3600)
                      layout(matrix(c(1:3, 0,4:5 ,7,0,6),3,3,F)) ; layout.show(7)    } # end of if 4 axes

  # range of values on focal axes
  rge_coord<-range(coord[,Faxes])
  
  # if axes limits not provided, values set according to range of coordinates for all axes 
  if ( sum( is.na(Faxes_lim) ) !=0 ) { 
          Faxes_lim<-c( rge_coord[1]-0.05*(rge_coord[2]-rge_coord[1]) , rge_coord[2]+0.05*(rge_coord[2]-rge_coord[1])) # range extended by 5%
          } # end of if 
  
  # axes labels
  lab<-pretty(Faxes_lim, n=4) # pretty labels
  lab<-lab[which(lab>=Faxes_lim[1] & lab<=Faxes_lim[2]) ]# filtering
  
  
  # background in col_bg if convex of all species to be drawn, else set to white
  if( convex_pool==FALSE)  { col_bg<-NA } # white
  if( convex_pool==TRUE)  { require(geometry) } # loading library for further use

    # loop on pairs of axes
    for (i in 1:(nb_Faxes-1) )
      for (j in (i+1):nb_Faxes )
      {
        # setting margins size of plot
        par( pty="s", mar=c(4,4,1,1) ) 
        
        # customized empty plot (with only axes) to be sure axes limits fit chosen limits
        plot(Faxes_lim,Faxes_lim,type="n",axes=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=Faxes_lim,ylim=Faxes_lim) # window
        rect(Faxes_lim[1],Faxes_lim[1],Faxes_lim[2],Faxes_lim[2], col=col_bg)   # border and background

          # customized X and Y axes
          for (k in 1:2)
          {
          axis(side=k, at=lab, labels=F, tcl=-0.3, pos=Faxes_lim[1])  # ticks
          mtext(side=k, lab, at=lab, line=-0.2, cex=0.9, las=1) # labels
          mtext(side=k, c(Faxes_nm[i] ,Faxes_nm[j])[k], cex=1,line=1.5, font=2) # title  
          } # end of k
        
        
        # projected  convex hull of all species filled in white
        if( convex_pool==TRUE)  { 
          vert0<-convhulln( coord[,c(Faxes[i],Faxes[j])] ,"Fx TO 'vert.txt'")
          vert1<-scan("vert.txt",quiet=T)
          vert_ij<-(vert1+1)[-1]
          polygon(coord[vert_ij,c(Faxes[i],Faxes[j])], border=NA,col="white")   
        } # end of if convex
        

        # plotting species with points or text
        if ( is.numeric(pch_sp) )   { points( coord[,c(Faxes[i],Faxes[j])], pch=pch_sp, cex=cex_sp, col=col_sp, bg=col_sp ) } # points
        if ( is.character(pch_sp) ) { text( coord[,c(Faxes[i],Faxes[j])], pch_sp, cex=cex_sp, col=col_sp ) } # text

      } # end of i,j

        # closing jpeg file. IF NOT extra panel on top right coudl be used for plotting legend
        if( close_jpeg==TRUE) { graphics.off() } 
  
  }# end of function
    
##################################################################################################################################################################
##################################################################################################################################################################

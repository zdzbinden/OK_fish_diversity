################################################################################################################################################
# R function to compute redundancy and vulnerability metrics based on distribution of species into Functional Entities (FE) for a set of assemblages
#
#  INPUTS: 
#   - species_and_FE: a list of objects provided by the function "species_to_FE" (see corresponding R script)
#   - asb_sp_01: a matrix of species (columns) occurences in a set of assemblages (rows)
#                 All species in "asb_sp_01" should be present in the objects of "species_and_FE"
#   - check_species_pool: a logical value indicating wether two tests are performed to check that all columns of 'asb_sp_01' have a sum stricly positive 
#                         (i.e. all species should be present in at least one assemblage). Default is TRUE. 
#                         
#   - 'folder_plot': a character string for setting the working directory where plots will be saved. NULL means current working directory.
#   - 'nm_asb_plot': a vector with the names of assemblages for which FE metrics need to be illustrated. 
#                       Default is NULL (i.e. no plot). To plot metrics for all assemblages, set 'nm_asb_plot' to 'row.names(asb_sp_01)'.
#   - 'col_FE': a R color name for the filling of bars representing number of species per FE. Default is light grey.
#   - 'border_FE': a R color name for the border of bars representing number of species per FE. Default is black. NA means no border.
#   - 'col_overR': a R color name for illustrating functional over-redundancy, as the proportion of species in excess in species-rich FE
#                  i.e. FE with mroe species than average functional redundancy (dashed horizontal line). Default is green.
#   - 'col_vuln': a R color name for the horizontal arrow illustrating Funct. Vulnerability as the proportion of FE with a single species. 
#                   Default is red.
#   - 'print_FE': a logical value indicating whether names of FE are shown on plot. Default is TRUE.
#
#  OUTPUTS: 
#     - a matrix of FE-based metrics (columns) for each assemblage (rows)
#
#     - for each of the assemblages listed in 'nm_asb_plot': a graphical file illustrating distribution of species in FE and associated metrics. 
#       Jpeg File has a resolution of 150dpi and dimensions of 1800x1200 pixels which means a size of around 150ko.
#       Funct ent. are sorted according to their number of species (decreasing from left to right). They are represented as colored bars. 
#       Functional redundancy is illustrated as horizontal dashed line, functional over-redundancy as color-filled area over this line 
#            and functional vulnerability as an horizontal double-head arrow.
#
################################################################################################################################################


 FE_metrics<-function(species_and_FE, asb_sp_01, check_species_pool=TRUE,
                        folder_plot=NULL, nm_asb_plot=NULL, 
                        col_FE="grey60", border_FE="black", col_overR="green4", col_vuln="red2", print_FE=TRUE) {
 

  ##############################################################################
  # checking inputs
    
   # species and FE
    if ( sum(names(species_and_FE)%in%c("FE", "FE_codes",  "FE_sp_01", "FE_traits"))!=4 ) 
        stop(paste(" error: 'species_and_FE' should be a list of objects built with the R function 'species_to_FE' "))
  
   # assemblages
    if( is.matrix(asb_sp_01)==FALSE ) stop( " 'asb_sp_01' should be an object of type matrix")
    if( sum(colnames(asb_sp_01) %in% names(species_and_FE$FE))!= ncol(asb_sp_01) ) stop(paste(" error: all species in 'asb_sp_01' should be in 'species_and_FE'"))
    if( is.numeric(asb_sp_01)==FALSE ) stop(paste(" error: 'asb_sp_01' is not numeric"))
    if( length( which(asb_sp_01==0 | asb_sp_01==1) )!=length(asb_sp_01) ) stop(paste(" error: 'asb_sp_01' must contain only 0 and 1"))
    if(min(apply(asb_sp_01,1,sum))==0 ) stop(paste(" error: all assemblages must host at least 1 species"))

    # checking species pool
    if (check_species_pool==TRUE)
    {
      if(min(apply(asb_sp_01,2,sum))==0 ) 
        stop(paste(" error: all columns of 'asb_sp_01' should have a sum stricly positive, i.e. all species must occur in at least one assemblage"))
    }# end of check species pool

   # checking graphical parameters
   if( length(col_FE)!=1 ) stop(paste(" error:'col_FE' should contain only one value"))
   if( length(border_FE)!=1 ) stop(paste(" error:'border_FE' should contain only one value"))
   if( length(col_overR)!=1 ) stop(paste(" error:'col_overR' should contain only one value"))
   if( length(col_vuln)!=1 ) stop(paste(" error:'col_vuln' should contain only one value"))
   
  ##############################################################################
    
   # names and number of assemblages
   nm_asb<-row.names(asb_sp_01)
   nb_asb<-length(nm_asb)
   
    # matrix to store results
    indices<-c( "Nb_sp", "Nb_FE", "F_Redundancy", "F_OverRedundancy", "F_Vulnerability" )
    asb_indices<-matrix(NA, nb_asb, length(indices), dimnames=list(nm_asb,indices))
    
    ##############################################################################
    # setting same graphical parameters for all assemblages
    
    # saving name of current working directory
    current_wd<-getwd()
    
    # setting working directory to store jpeg files
    if (is.null(folder_plot) ) {  folder_plot<-current_wd }
    
    # setting folder for saving jpeg files
    test_folder_plot<-try( setwd(folder_plot) , silent=TRUE)
    if ( class(test_folder_plot) =="try-error") {
      folder_plot<-current_wd 
      print(paste(" /!\    WARNING: '",folder_plot," does not exist', jpeg files have been saved in '",current_wd, "'",sep="" ))
    } # end of if
    setwd(folder_plot)
    
    ##############################################################################
    # loop on assemblages for computing indices and plotting distribution of species in FE
  
    for (k in nm_asb)
    { 
      # species present
      sp_k<-names(which(asb_sp_01[k,]==1))
      nb_sp_k<-length(sp_k)
      
      # FE present
      FE_k<-unique(species_and_FE$FE[sp_k])
      nb_FE_k<-length(FE_k)
      
      # number of species per FE
      if( length(FE_k)==1) { FE_nbsp_k<-sum(species_and_FE$FE_sp_01[FE_k,sp_k]) }# end of if >1 FE
      if( length(FE_k)>1) { FE_nbsp_k<-apply(species_and_FE$FE_sp_01[FE_k,sp_k],1,sum) }# end of if >1 FE
      
      # Functional redundancy = average number of species per FE
      F_Redundancy_k<-nb_sp_k/nb_FE_k
      
      # Functional over-redundancy = proportion of species in excess in species-rich FE
      F_OverRedundancy_k<-sum( sapply(FE_nbsp_k, function(x) { max( c(x, F_Redundancy_k ) - F_Redundancy_k ) } )   ) / nb_sp_k
        
      # Functional vulnerability = proportion of FE with only 1 species
      F_Vulnerability_k<-length(which(FE_nbsp_k==1))/nb_FE_k
      
      # indices values
      asb_indices[k,]<-c( nb_sp_k, nb_FE_k, F_Redundancy_k, F_OverRedundancy_k, F_Vulnerability_k)

      
    ##############################################################################
    # plot of distribution of species into FE
    
    if(k %in% nm_asb_plot ) {
      
      # creating a jpeg file
      jpeg(file=paste(k,"_nbsp_per_FE.jpeg",sep=""), res=150, width=1800, height=1200)

      # setting margins size of plot
      par( mar=c(2,2,2,1) ) 
      
      # setting X axis limit and labels
      lim_x<-c(0,nb_FE_k)+0.5
      lab_x<-round(pretty(c(1,nb_FE_k), n=4)) # labels by default given number of FE
      lab_x<-unique(c(1,lab_x[which(lab_x>=lim_x[1] & lab_x<=lim_x[2]) ], nb_FE_k)) # customizing to ensure 1 and nb of FE present
      
      # setting Y axis limit and labels
      lim_y<-c(0,max(FE_nbsp_k)+0.5)
      lab_y<-round( pretty(FE_nbsp_k, n=4) ) # labels by default given number of species per FE
      lab_y<-unique(c(0,1,lab_y[which(lab_y>=lim_y[1] & lab_y<=lim_y[2]) ], max(FE_nbsp_k)) ) # customizing to ensure 0 and max nb of species present
      
    
      # window
      plot(lim_x,lim_y,type="n",axes=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=lim_x,ylim=lim_y) # default
      rect(lim_x[1],lim_y[1],lim_x[2],lim_y[2], col="white") # customized window
      
      # customized X axis
      axis(side=1, at=lab_x, labels=F, tcl=-0.4, pos=lim_y[1])  # ticks
      mtext(side=1, lab_x, at=lab_x, line=-0.8, cex=1.2, las=1) # labels
      mtext(side=1, "Rank of Functional Entities", cex=1.5,line=0.7, font=2) # title  

      # customized Y axis
      axis(side=2, at=lab_y, labels=F, tcl=-0.4, pos=lim_x[1])  # ticks
      mtext(side=2, lab_y, at=lab_y, line=-1.5, cex=1.2, las=1) # labels
      mtext(side=2, "Number of species", cex=1.5 ,line=0.4, font=2) # title  
      
      # sorting FE according to decreasing number of species
      FE_nbsp_k_D<-sort(FE_nbsp_k,T)
      
      # bar to illustrate species richness in each FE
      rect( (1:nb_FE_k)-0.5, 0,(1:nb_FE_k)+0.5, FE_nbsp_k_D, col=col_FE , border=NA )
      
      # filling some bars to illustrate over-redundancy
      rect( (1:nb_FE_k)[which(FE_nbsp_k_D>=F_Redundancy_k)]-0.5, F_Redundancy_k, (1:nb_FE_k)[which(FE_nbsp_k_D>=F_Redundancy_k)]+0.5,
        FE_nbsp_k_D[which(FE_nbsp_k_D>=F_Redundancy_k)] , col=col_overR , border=NA )
      
      # horizontal arrow for illustrating vulnerability
      if (min(FE_nbsp_k_D)==1)
        { arrows( min(which(FE_nbsp_k_D==1))-0.5, 1+0.03*(lim_y[2]-1), max(which(FE_nbsp_k_D==1))+0.5, 1+0.03*(lim_y[2]-1) , code=3, col=col_vuln , lwd=3 ) }
      
      # horizontal dashed line for redundancy
      segments( lim_x[1], F_Redundancy_k, lim_x[2], F_Redundancy_k, lwd=4, lty=2)
      
      # border of bars (if needed)
      if (is.na(border_FE)==FALSE) { rect( (1:nb_FE_k)-0.5, 0,(1:nb_FE_k)+0.5, FE_nbsp_k_D, col=NA , border=border_FE ) }
      
      # names of FE if needed
      if (print_FE==TRUE) text( 1:nb_FE_k, 0.2, names(FE_nbsp_k_D), srt=90, adj=c(0,0.5) , cex=1 )
      
      # title = indices values
      mtext( side=3, paste(" Nb. sp=", nb_sp_k, "  Nb. FE=", nb_FE_k, 
            "  Funct. Redund.=", round(F_Redundancy_k,1),"  FOver-Redund=", round(F_OverRedundancy_k,1), "  Funct. Vuln.=", round(F_Vulnerability_k,1) , sep=""), 
            cex=1.5 ,line=0, font=2)

      # closing jpeg
      graphics.off()
      
    }# end of if plot
    
    
      
    }# end of k (loop on assemblages)

    ##############################################################################
    # returning to current working directory
    setwd(current_wd)
    
    # returning results	
     return(asb_indices)	

 } # end of function FE_metrics
########################################################################################################################################
########################################################################################################################################




################################################################################################################################################
# R function to classify species into Functional Entities based on non-continuous traits
#
#  INPUT: 
#     - 'traits' a matrix with values of functional traits (columns) for a set of species (rows). 
#         All traits should be objects of type 'factor' or 'ordered' ('numeric' objects are not allowed). NA are not allowed.
#         2 first letters of trait names should be unique, and for each trait 2 first letters of levels (i.e. categories) should be unique 
#         NB: Trait names and trait levels could be coded with a single character.
#
#  OUTPUT: a list of objects
#             - FE: Funct Ent to which each species belongs (names of vector elements are row names of 'traits')
#             - FE_codes: a vector with Funct. Ent. names
#             - FE_sp_01: a presence-absence (0/1) like matrix of species (rows) in Funct Ent (columns)
#             - FE_traits: a matrix with trait values of the Funct. Ent.
#
#  Names of Funct Ent are made as a character chain of up to 2 letters for trait name in upper case font then up to 2 letters for level in lower case font,
#   separated by "_" between traits; examples: ("TAc2_TBxx_TCyy" & "TAc3_TBff_TCyy") or ("A2_Bx_Cy" & "A3_Bf_Cy")
#   NB: to ensure that names of Funct. entities are as short as possible (while being unique), trait names are abbreviated to a single letter whenever possible
#
################################################################################################################################################


species_to_FE<-function(traits) {

########################################
# checking input
  
# size of trait matrix  
if ( nrow(traits)<2 ) stop("Error: 'traits' should have at least 2 rows")
if ( ncol(traits)<2 ) stop("Error: 'traits' should have at least 2 columns")
  
# checking coding of trait and absence of NA  
for ( t in names(traits) )
{
  if ( is.numeric(traits[,t]) ) stop( paste("Error: trait '",t,"' is coded as 'numeric' ", sep="") )
  if( sum(is.na((traits[,t])))>0 ) stop(paste("Error: NA in trait '",t,"' ", sep=""))
  
}# end of t  
  
  
# checking that trait names have different first letter or at least unique 2 first letters
traits_codes<-substr(names(traits),1,1)  # default is one letter

if ( length(unique(traits_codes))!=ncol(traits) ) {
  traits_codes<-substr(names(traits),1,2) # if some traits have the same first letter, then up to 2 first letters kept
  if ( length(unique(traits_codes))!=ncol(traits) ) stop("Error: 2 first letters of trait names should be unique ")
}# end of t
names(traits_codes)<-names(traits)

# checking for each trait that levels have unique two first letters
for ( t in names(traits) )
{
  # levels of trait t
  mod_t<-unique(traits[,t])
  
  # checking the uniqueness of their first 2 letters
  if ( length(unique(substr(mod_t,1,2)))!=length(mod_t) )  stop( paste("Error: some levels of trait '",t,"' have the same 2 first letters", sep="") )
}# end of t

########################################

# defining Functional Entities as unique combinations of trait values
FE<-paste( toupper(traits_codes[1]), tolower(substr(traits[,1], 1,2) ) ,sep="")
for (t in names(traits)[-1] )
{
  FE<-paste( FE, paste( toupper(traits_codes[t]), tolower(substr(traits[,t], 1,2) ) ,sep="") , sep="_")
} # end of t
names(FE)<-row.names(traits)

########################################

# codes of FE
FE_codes<-unique(FE)

########################################
# matrix of species occurence in FE
FE_sp_01<-matrix(0, length(FE_codes), nrow(traits), dimnames=list( FE_codes, row.names(traits) ) )
for (f in FE_codes)
{ 
  FE_sp_01[f,names(which(FE==f))]<-1
}# end of f

########################################

# trait values for FE
FE_traits<-traits[apply(FE_sp_01, 1, function(x) {names(which(x==1))[1]} ), ]
row.names(FE_traits)<-row.names(FE_sp_01)

########################################
# results in a single list
res<-list( FE=FE, FE_codes=FE_codes, FE_sp_01=FE_sp_01, FE_traits=FE_traits)
return(res)

########################################

 } # end of function


##############################################################################################
# END of function
##############################################################################################
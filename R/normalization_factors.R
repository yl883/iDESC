#' normalization_factors: Calculating different types of normalization factors
#'
#' @param mat Count matrix
#' @param norm_opt  Option for normalizing factors
#' @param user_sf Option for user-specific normalizing factors
#'
#' @return A vector of normalizing factors
#' @export

normalization_factors<-function(mat,norm_opt,user_sf){
  mat<-as.matrix(mat)
  if(norm_opt!="User"&(!is.null(user_sf))){
    print('Disable user-specific normalization, please set the norm_opt="User".')
    norm_opt="None"
  }
  if(norm_opt=="User"&(is.null(user_sf))){
    print('No user-specific normalization provided, disable normalization.')
    norm_opt="None"
  }
  if(norm_opt=="SeqDepth"){
    sf=colSums(mat)
  }else if(norm_opt=="SizeFactor"){
    GEOmean_cell<-apply(mat,1,function(gene){
      exp(sum(log(gene[which(gene>0)]), na.rm = TRUE)/length(gene[which(gene>0)]))
    })
    sf<-apply(mat,2,function(cell){
      stats::median((cell/GEOmean_cell)[which(cell/GEOmean_cell != 0)])
    })
    rm(GEOmean_cell)
  }else if(norm_opt=="User"&(!is.null(user_sf))){
    if(ncol(mat)!=length(user_sf)){
      print("User-specific normalization cannot be applied to all cells due to data filtering")
      user_sf<-user_sf[colnames(mat)]
    }
    sf=user_sf
  }else if(norm_opt=="None"){sf<-rep(1,ncol(mat))}
  return(sf)
}

#zp_diagnosis is desinged for 1/2 groups
#' zp_diagnosis: Diagnosis about group differences on inflated zeros
#'
#' @param vec A vector of data (one gene) that need diagnosis
#' @param zp_group_var The name of group/disease information for differences of inflated zeros in meta
#' @param zp.thresh Threshold of zero proportions
#' @param meta Data frame including information for cells
#'
#' @return A string with the format "Two Pi Model" or "One Pi Model"
#' @export
#'
#' @examples
zp_diagnosis<-function(vec,meta,zp_group_var,zp.thresh=zp.thresh){
  if(length(levels(factor(meta[,zp_group_var])))<=1){
    return("One Pi Model")
  }
  if(length(levels(factor(meta[,zp_group_var])))>2){
    return("Two Pi Model")
  }
  zp<-sapply(levels(factor(meta[,zp_group_var])),function(group){
    c(mean(vec[which(meta[,zp_group_var]==group)]==0,na.rm=T),sum(meta[,zp_group_var]==group,na.rm=T))
  })
  z_zp<-(zp[1,1]-zp[1,2])/sqrt((zp[1,1]*(1-zp[1,1])/zp[2,1])+(zp[1,2]*(1-zp[1,2])/zp[2,2]))
  return(ifelse(2*stats::pnorm(z_zp,lower.tail = F)<=zp.thresh,"Two Pi Model","One Pi Model"))
}

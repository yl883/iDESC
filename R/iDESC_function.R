#' iDESC: main function for DE analysis using zero-inflated negative binomial mixed model
#'
#' @param mat Count matrix
#' @param meta Data frame including information for cells
#' @param subject_var The name of subject information in meta
#' @param group_var The name of group/disease information for DE analysis in meta
#' @param norm_opt Option for normalizing factors
#' @param user_sf Option for user-specific normalizing factors
#' @param sub_cell_filtering Filtering on cells within each subject
#' @param gene_sub_filtering Filtering on genes based on expression in each subject
#' @param gene_cell_filtering Filtering on genes based on expression across all cells
#' @param ncell_filtering Filtering on cells based on the number of genes expressed
#' @param span smoothing parameter for LOESS curve
#'
#' @return A result table contains parameter estimation
#' @export
#'
#' @examples
#'
#' library(iDESC)
#' data(IPF_example)
#' mat=IPF_example$mat
#' meta=IPF_example$meta
#' sequencing_depth=IPF_example$meta$sequencing_depth
#' result=iDESC(mat,meta,subject_var="subject",group_var="disease",norm_opt="User",user_sf = sequencing_depth,span = 0.7)
iDESC<-function(mat,meta,subject_var,group_var,norm_opt=c("SeqDepth","SizeFactor","User","None"),user_sf=NULL,
                sub_cell_filtering=5,gene_sub_filtering=0,gene_cell_filtering=0.05,ncell_filtering=1,span=0.05){
  if(!is.null(user_sf)){
    if(is.null(names(user_sf))){names(user_sf)<-colnames(mat)}else{user_sf<-user_sf[colnames(mat)]}
  }
  dat<-preprocessing(mat,meta,subject_var,group_var,sub_cell_filtering=sub_cell_filtering,gene_sub_filtering=gene_sub_filtering,gene_cell_filtering=gene_cell_filtering,ncell_filtering=ncell_filtering)
  mat<-dat$mat
  meta<-dat$meta
  norm_factor<-normalization_factors(mat,norm_opt,user_sf)
  group<-meta[,group_var]
  subject<-meta[,subject_var]

  predict_pi<-zp_prediction(mat,norm_factor,span)
  stopifnot("LOESS: zero-width neighborhood. make span bigger"=sum(is.na(predict_pi))==0)

  res_tb<-Reduce(plyr::rbind.fill,lapply(1:nrow(mat),function(g){
    gene<-mat[g,]
    predict_pi_offset<-predict_pi[g]
    tmp.df<-data.frame(y=gene,norm_sf=norm_factor,predict_pi_offset=predict_pi_offset,ind_zero=1*(gene==0),group=group,sub=as.numeric(factor(subject)))

    f1<-try(glmmTMB::glmmTMB(y ~ group + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = glmmTMB::nbinom2, zi = ~ offset(predict_pi_offset)-1+ind_zero+(1 | sub)))

    res1<-try(c(summary(f1)$coefficients$cond[,"Estimate"],1/summary(f1)$sigma,exp(f1$fit$par["theta"])^2,exp(f1$fit$par["thetazi"])^2,
                summary(f1)$coefficients$zi[1,"Estimate"],summary(f1)$coefficients$cond[-1,"Pr(>|z|)"],summary(f1)$AICtab["deviance"]))
    if('try-error' %in% class(res1)){
      res1<-c(Alpha=NA)
    }else{
      names(res1)<-c("Alpha",paste0("Beta_",names(summary(f1)$coefficients$cond[,"Estimate"])[-1]), "Dispersion","Sigma2","Sigma2_ZI","Theta",
                     paste0("Pval_Beta_",names(summary(f1)$coefficients$cond[,"Estimate"])[-1]),"Deviance")

    }

    data.frame(t(res1))
  }))
  rownames(res_tb)<-rownames(mat)
  res_tb<-res_tb[which(!is.na(res_tb[,1])&!is.na(res_tb[,2])&!is.na(res_tb[,7])),]
  return(res_tb)
}





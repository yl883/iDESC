#' iDESC: main function for DE analysis using zero-inflated negative binomial mixed model
#'
#' @param mat Count matrix
#' @param meta Data frame including information for cells
#' @param subject_var The name of subject information in meta
#' @param group_var The name of group/disease information for DE analysis in meta
#' @param zp_group_var The name of group/disease information for differences of inflated zeros in meta
#' @param norm_opt Option for normalizing factors
#' @param user_sf Option for user-specific normalizing factors
#' @param sub_cell_filtering Filtering on cells within each subject
#' @param gene_sub_filtering Filtering on genes based on expression in each subject
#' @param gene_cell_filtering Filtering on genes based on expression across all cells
#' @param zp.thresh Threshold of zero proportions
#' @param diagnosis Option for doing diagnosis on zero proportions
#'
#' @return A result table contains parameter estimation
#' @export
#'
#' @examples

iDESC<-function(mat,meta,subject_var,group_var,zp_group_var,norm_opt=c("SeqDepth","SizeFactor","User","None"),user_sf=NULL,
                sub_cell_filtering=5,gene_sub_filtering=1,gene_cell_filtering=0.05,zp.thresh=0.05,diagnosis=F){
  dat<-preprocessing(mat,meta,subject_var,group_var,sub_cell_filtering=sub_cell_filtering,gene_sub_filtering=gene_sub_filtering,gene_cell_filtering=gene_cell_filtering)
  mat<-dat$mat
  meta<-dat$meta
  norm_factor<-normalization_factors(mat,norm_opt,user_sf)
  group<-meta[,group_var]
  subject<-meta[,subject_var]
  if(diagnosis){
    model_select_diagnosis<-apply(mat,1,function(g){
      zp_diagnosis(g,meta,zp_group_var,zp.thresh)
    })
  }
  res_tb<-Reduce(plyr::rbind.fill,lapply(1:nrow(mat),function(g){
    gene<-mat[g,]
    model_select<-model_select_diagnosis[g]
    tmp.df<-data.frame(y=gene,norm_sf=norm_factor,group=group,zp_group=meta[,zp_group_var],sub=as.numeric(factor(subject)))
    if(model_select=="One Pi Model"){
      f1<-try(glmmTMB::glmmTMB(y ~ group + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = glmmTMB::nbinom2, zi = ~ 1))
      res1<-try(c(summary(f1)$coefficients$cond[,"Estimate"],1/summary(f1)$sigma,exp(f1$fit$par["theta"]),
                  boot::inv.logit(c(summary(f1)$coefficients$zi[1,"Estimate"],summary(f1)$coefficients$zi[-1,"Estimate"]+summary(f1)$coefficients$zi[1,"Estimate"])),summary(f1)$coefficients$cond[-1,"Pr(>|z|)"],
                  summary(f1)$coefficients$zi[-1,"Pr(>|z|)"],f1$fit$convergence))
      if('try-error' %in% class(res1)){
        res1<-c(Convergence=NA)
      }else{
        names(res1)<-c("Alpha",paste0("Beta_",names(summary(f1)$coefficients$cond[,"Estimate"])[-1]), "Dispersion","Sigma",
                       "Pi_.Intercept.",
                       paste0("Pval-Beta_",names(summary(f1)$coefficients$cond[,"Estimate"])[-1]),"Convergence")
        if(sum(is.na(res1))){res1["Convergence"]<-1}
      }
    }else{
      f1<-try(glmmTMB::glmmTMB(y ~ group + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = glmmTMB::nbinom2, zi = ~ zp_group))
      res1<-try(c(summary(f1)$coefficients$cond[,"Estimate"],1/summary(f1)$sigma,exp(f1$fit$par["theta"]),
                  boot::inv.logit(c(summary(f1)$coefficients$zi[1,"Estimate"],summary(f1)$coefficients$zi[-1,"Estimate"]+summary(f1)$coefficients$zi[1,"Estimate"])),summary(f1)$coefficients$cond[-1,"Pr(>|z|)"],
                  summary(f1)$coefficients$zi[-1,"Pr(>|z|)"],f1$fit$convergence))
      if('try-error' %in% class(res1)){
        res1<-c(Convergence=NA)
      }else{
        names(res1)<-c("Alpha",paste0("Beta_",names(summary(f1)$coefficients$cond[,"Estimate"])[-1]), "Dispersion","Sigma",
                       paste0("Pi_",names(summary(f1)$coefficients$zi[,"Estimate"])),
                       paste0("Pval-Beta_",names(summary(f1)$coefficients$cond[,"Estimate"])[-1]),
                       paste0("Pval-Pi_",names(summary(f1)$coefficients$zi[,"Estimate"])[-1]),"Convergence")
        if(sum(is.na(res1))){res1["Convergence"]<-1}
      }
    }
    data.frame(t(res1))
  }))
  rownames(res_tb)<-rownames(mat)
  res_tb<-cbind(res_tb,Model=model_select_diagnosis)
  res_tb<-res_tb[res_tb$Convergence==0,-which(colnames(res_tb)=="Convergence")]
  return(res_tb)
}





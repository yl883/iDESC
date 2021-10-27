
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
      median((cell/GEOmean_cell)[which(sample_j != 0)])
    })
    rm(GEOmean_cell)
  }else if(norm_opt=="User"&(!is.null(user_sf))){
    stopifnot("User-specific normalization cannot be applied to all cells"=ncol(mat)==length(user_sf))
    sf=user_sf
  }else if(norm_opt=="None"){sf<-rep(1,ncol(mat))}
  return(sf)
}

preprocessing<-function(mat,meta,subject_var,group_var,sub_cell_filtering,gene_sub_filtering,gene_cell_filtering){
  #meta is a df with rownames=cell_barcode, columns includes:subject,group
  #sub_cell_filtering is an integer
  stopifnot("Gene filtering should be based on an indicated proportion"=(gene_sub_filtering>=0&gene_sub_filtering<=1))
  stopifnot("Gene filtering should be based on an indicated proportion"=(gene_cell_filtering>=0&gene_cell_filtering<=1))
  stopifnot("Cell barcodes in meta data and matrix should be the same"=(colnames(mat)==rownames(meta)))

  mat<-as.matrix(mat)

  #Remove subject with <"sub_cell_filtering" cells
  barcode_orig<-rownames(meta)
  subject_name<-names(table(meta[,subject_var]))
  remain_sub<-subject_name[which(table(meta[,subject_var])>=sub_cell_filtering)]
  remain_cell<-barcode_orig[which(meta[,subject_var]%in%remain_sub)]
  meta<-meta[remain_cell,]

  #Keep genes expressed at least a% cells in b% subject in either group
  if(length(levels(factor(meta[,group_var])))>3){
    print('"Group" has more than 3 levels and was treated as continuous.')
    group_sub_list<-list(subject_name)
  }else{
    group_sub_list<-lapply(1:length(levels(factor(meta[,group_var]))),function(i){
      unique(meta[which(meta[,group_var]==levels(factor(meta[,group_var]))[i]),subject_var])
    })
  }

  #sub_prop is a gene*sub matrix
  sub_prop<-sapply(remain_sub,function(sub){
    rowMeans(mat[,barcode_orig[which(meta[,subject_var]==sub)]]!=0,na.rm=T)
  })

  remain_gene<-rownames(mat)[which(rowSums(sapply(group_sub_list,function(sub_list){
    rowSums(sub_prop[,as.character(sub_list)]>=gene_cell_filtering,na.rm=T)>=(gene_sub_filtering*length(sub_list))
  }))>0)]
  return(list(mat=mat[remain_gene,remain_cell],meta=meta))
}

#zp_diagnosis is desinged for 1/2 groups
zp_diagnosis<-function(vec,zp_group_var,zp.thresh=zp.thresh){
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
  return(ifelse(2*pnorm(z_zp,lower.tail = F)<=zp.thresh,"Two Pi Model","One Pi Model"))
}



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
      zp_diagnosis(g,zp_group_var,zp.thresh)
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





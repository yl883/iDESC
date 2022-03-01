#' preprocessing: filtering step on genes/cells
#'
#' @param mat Count matrix
#' @param meta Data frame including information for cells
#' @param subject_var The name of subject information in meta
#' @param group_var The name of group/disease information for DE analysis in meta
#' @param sub_cell_filtering Filtering on cells within each subject
#' @param gene_sub_filtering Filtering on genes based on expression in each subject
#' @param gene_cell_filtering Filtering on genes based on expression across all cells
#'
#' @return A list of processed count matrix and data frame
#' @export
#'
#' @examples
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

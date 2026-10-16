#' preprocessing: filtering step on genes/cells
#'
#' @param mat  Count matrix
#' @param meta Data frame including information for cells
#' @param subject_var The name of subject information in meta
#' @param group_var The name of group/disease information for DE analysis in meta
#' @param sub_cell_filtering Filtering on cells within each subject
#' @param gene_sub_filtering Filtering on genes based on expression in each subject: a gene is kept
#' if, in at least one group, it is expressed in at least `gene_cell_filtering` of the cells in at
#' least `gene_sub_filtering` x (number of subjects of that group) subjects. With 0 every gene is kept.
#' @param gene_cell_filtering Minimum proportion of cells of a subject expressing the gene (used
#' together with `gene_sub_filtering`)
#' @param ncell_filtering Filtering on cells based on the number of genes expressed
#' @param cores Number of cores for parallel processing (default is 1)
#'
#' @return A list with the processed count matrix (`mat`) and the matching data frame of cells (`meta`)
#' @importFrom Matrix colSums
#' @export

preprocessing<-function(mat,meta,subject_var,group_var,sub_cell_filtering,gene_sub_filtering,gene_cell_filtering,ncell_filtering,cores=1){
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
  sub_cell<-barcode_orig[which(meta[,subject_var]%in%remain_sub)]
  expressed_cell<-barcode_orig[which(Matrix::colSums(mat!=0)>0)]
  remain_cell<-intersect(sub_cell,expressed_cell)
  meta<-meta[remain_cell,]
  #barcodes of the cells that remain (positions refer to the filtered meta)
  barcode_kept<-rownames(meta)

  #Keep genes expressed at least a% cells in b% subject in either group
  if(length(levels(factor(meta[,group_var])))>3){
    print('"Group" has more than 3 levels and was treated as continuous.')
    group_sub_list<-list(unique(as.character(meta[,subject_var])))
  }else{
    group_sub_list<-parallel::mclapply(1:length(levels(factor(meta[,group_var]))),function(i){
      unique(meta[which(meta[,group_var]==levels(factor(meta[,group_var]))[i]),subject_var])
    },mc.cores=cores)
  }

  #sub_prop is a gene*sub matrix
  sub_prop_list<-parallel::mclapply(remain_sub,function(sub){
    rowMeans(mat[,barcode_kept[which(meta[,subject_var]==sub)],drop=FALSE]!=0,na.rm=T)
  },mc.cores=cores)

  # list to matrix
  sub_prop <- do.call(cbind, sub_prop_list)
  colnames(sub_prop) <- remain_sub

  group_pass <- do.call(cbind, lapply(group_sub_list,function(sub_list){
    rowSums(sub_prop[,as.character(sub_list),drop=FALSE]>=gene_cell_filtering,na.rm=T)>=(gene_sub_filtering*length(sub_list))
  }))
  remain_gene<-rownames(mat)[which(rowSums(group_pass)>0)]

  # filter cells by number of expressed genes (ncell_filtering)
  mat_filtered <- mat[remain_gene, remain_cell, drop = FALSE]
  n_genes_per_cell <- Matrix::colSums(mat_filtered != 0)
  final_cells <- names(n_genes_per_cell[n_genes_per_cell >= ncell_filtering])

  # Subset matrix and metadata to keep only those cells
  mat_filtered <- mat_filtered[, final_cells, drop = FALSE]
  meta_filtered <- meta[final_cells, ]

  return(list(mat = mat_filtered, meta = meta_filtered))
}

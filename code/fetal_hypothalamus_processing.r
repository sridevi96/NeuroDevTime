library(Seurat)
library(dplyr)
library(caret)
library(readr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(tidyr)

#data is downloaded from
#'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoSamples_AllCells.rda'
#metadata is downloaded from S3: https://www.science.org/doi/suppl/10.1126/sciadv.adf6251/suppl_file/sciadv.adf6251_sm.pdf


load("data/Herb_hypothalamus_2023/Herb_hypothalamus_2023.Rdata") #LIST OF SEURAT OBJECTS CALLED HYPODAT
load("data/Herb_hypothalamus_2023/Sample_metadata.Rdata")#METADATA DOWNLOADED FROM PAPER SUPPLEMENTARY; includes cells from another study as well
final_common_genes<-readRDS("models/commongenes.rds")#list of common genes


Herb_metadata<-Herb_hypothalamus_celltypelabels%>%filter(study=="Herb")
split_Herb_metadata<-split(Herb_metadata,f = Herb_metadata$sample)
names(split_Herb_metadata)<-lapply(split_Herb_metadata,function(x){unique(x$dev_age)})%>%unlist()

#Add metadata to each Herb object and filter to keep cells with metadata assigned cell types only
Herb_filtered_data<-lapply(1:length(HypoDat),function(x){
    sample_name<-names(HypoDat)[[x]]
    data=HypoDat[[sample_name]]
    metadata=split_Herb_metadata[[sample_name]]
    filtered_data<-data[,metadata$Barcode]
    filtered_data<-AddMetaData(filtered_data,metadata = metadata[match(Cells(filtered_data),metadata$Barcode),])
    return(filtered_data)
})

#Aggregate expression at the celltype level
names(Herb_filtered_data)<-lapply(Herb_filtered_data,function(x){unique(x$dev_age)})
merged_Herb_filtered<-merge(Herb_filtered_data[[1]],Herb_filtered_data[-1])

Herb_pseudobulk_celltypes_rawcounts<- lapply(Herb_filtered_data,function(x){
    AggregateExpression(x,assays = "RNA",group.by = "CellTypes",return.seurat = TRUE)})

age_carnegiestage<-c(6,7,7,10,10)#CS 13,14,15,22 respectively
names(Herb_pseudobulk_celltypes_rawcounts)[1:5]<-paste0(names(Herb_pseudobulk_celltypes_rawcounts)[1:5],"_GW",age_carnegiestage[1:5])


#Function to create a list of matrices, one for each celltype with columns corresponding to different samples
create_matrices_by_strings <- function(list_of_matrices, strings) {
  new_list <- lapply(strings, function(str) {
    lapply(list_of_matrices, function(mat) {
      selected_cols <- which(colnames(mat)==str)
      if (any(selected_cols)) {
        mat[, selected_cols, drop = FALSE]
      } else {
        matrix(NA, nrow = nrow(mat), ncol = 0)
      }
    })
  })
  new_list
}

Herb_celltypes<-unique(unlist(lapply(Herb_pseudobulk_celltypes_rawcounts,colnames)))
Herb_pb_celltypesbyage<-create_matrices_by_strings(Herb_pseudobulk_celltypes_rawcounts,Herb_celltypes)
names(Herb_pb_celltypesbyage)<-Herb_celltypes


#get a merged matrix of all cell type transcriptomes per sample
Herb_pb_combinedcelltypes<-lapply(Herb_pb_celltypesbyage,function(seuratobjlist){
    #each is a list of seurat objects with different times, some maybe empty
    indextoskip<-lapply(seuratobjlist,is.logical)%>%unlist()
    filteredlist<-seuratobjlist[!indextoskip]
    temp<-names(filteredlist)
    tempdata<-lapply(filteredlist,function(x){x@assays[["RNA"]]@layers[["data"]]})%>%do.call(cbind,.)
    rownames(tempdata)<-Features(Herb_pb_celltypesbyage$Astrocytes$CS13_GW6)#set rownames as genes from one sample
    return(tempdata)
    })%>%do.call(cbind,.)%>%.[final_common_genes,]#restrict to only the common genes

#get list of all cell types in the combined matrix
subtract_Herb<-lapply(Herb_pb_celltypesbyage,function(x){lapply(x,is.logical)%>%unlist%>%sum}) #some of the list elements are empty
celltypes_Herb<-list()
for (i in 1:length(Herb_pb_celltypesbyage)){
  celltypes_Herb[[i]]<-rep(names(Herb_pb_celltypesbyage)[i],11-subtract_Herb[[i]])
}
celltypes_Herb<-unlist(celltypes_Herb)


Herb_metadata<-data.frame(study=rep("Herb",ncol(Herb_pb_combinedcelltypes)),
                          obs_age=str_extract(colnames(Herb_pb_combinedcelltypes), "(?<=GW)\\d+"),
                          celltypes=celltypes_Herb)
save(Herb_pb_combinedcelltypes,Herb_metadata,file="processed_data/Herb_processed.Rdata")


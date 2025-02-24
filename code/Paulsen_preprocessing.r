library(Seurat)
library(dplyr)
library(Matrix)
library(readr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(tidyr)

setwd("/home/sridevi/inkwell03_sridevi//metadevorganoid/werneranalysis/DevTime_gitub_repo/") #change to the working directory

final_common_genes<-readRDS("models/commongenes.rds")#list of common genes

#files were downloaded from https://www.synapse.org/Synapse:syn26346581 
filepath="/home/sridevi/inkwell03_sridevi/metadevorganoid/celltypeagepredictor/data/disease_organoid_datasets/Paulsen_autism_arlottaorganoids_2022/"

#files in Paulsen et al.2022
filenames=list.files(filepath) 
datasets=filenames[str_detect(filenames,".rds")]%>%str_split(.,".rds")%>%lapply(.,function(x) x[[1]])
                                                                                
#filter the names to get gene, cell line, age, replicate
temp<-str_split(datasets,"_")

metadata <- do.call(rbind, lapply(temp, function(x) {
  length(x) <- 4  # Ensure all vectors have 4 elements, filling with NA if necessary
  return(x)
}))
metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
colnames(metadata)<-c("gene","line","age","replicate")
metadata

files=paste0(filepath,filenames[str_detect(filenames,".rds")])
basename(files)

#load each file, pseudobulk by organoid,celltype and disease state
Paulsen_pb_celltypes<-lapply(files,function(file){
    data<-readRDS(file)
    pb_data<-AggregateExpression(data,group.by = c("org","CellType","treat"),return.seurat = TRUE)
    rm(data)
    return(pb_data)
})
  
names(Paulsen_pb_celltypes)<-datasets

save(Paulsen_pb_celltypes,file="processed_data/Paulsen2022_processed.Rdata")
save(metadata,file="processed_data/Paulsen2022_metadata.Rdata")
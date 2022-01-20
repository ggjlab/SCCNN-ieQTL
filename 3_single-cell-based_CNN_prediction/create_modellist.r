library(stringi)
library(argparser, quietly=TRUE)
p <- add_argument(p, "--model_save_directory", help="", default="")
p <- add_argument(p, "--model_list_directory", help="", default="")
argv <- parse_args(p)

model_save_directory <- p$model_save_directory
model_list_directory <- p$model_list_directory

file_list <- list.files(path = model_save_directory , pattern = 'save')
 
all <- unlist(strsplit(file_list,'\\.'))
a <- grep(all, pattern='save')
all_celltypes <- all[a-1]

modellist <- data.frame(ModelName=paste0(model_save_directory,file_list) , Tissue=all_celltypes )

write.table(modellist, file=paste0(model_list_directory,'modellist'), sep='\t', quote=F, row.names=F, col.names=T)
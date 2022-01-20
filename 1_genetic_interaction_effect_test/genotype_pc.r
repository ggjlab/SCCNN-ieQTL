library(argparser, quietly=TRUE)
p <- arg_parser("")
p <- add_argument(p, "--deconv_result_dir", help="")
p <- add_argument(p, "--peer_dir", help="", default="")
p <- add_argument(p, "--tissue", short="-c", help="tissue name")
p <- add_argument(p, "--out_dir", short="-o", help="out_dir")

argv <- parse_args(p)

deconvolution_result_directory = argv$deconv_result_dir
PEER_directory = argv$peer_dir
out_dir = argv$out_dir

tissue = as.character(argv$tissue)

setwd(paste0(PEER_directory, tissue))

pca <- read.table("./.eigenvec",sep = "\t",header = F)
rownames(pca) <- pca$V2
pca <- pca[,3:ncol(pca)]
colnames(pca) <- paste0('pca','',c(1:5))

cell_fraction <- read.table(paste0(deconvolution_result_directory, tissue,'_fraction.txt') ,header = 1, row.names=1, sep='\t')
choosesample <- rownames(cell_fraction)
   
all <- unlist(strsplit(choosesample,'-'))
a <- grep(all, pattern='GTEX')
b <- a+1
choosesample_participant <- paste0(all[a],'-',all[b])

pca_sample <- pca[choosesample_participant,]
rownames(pca_sample) <- choosesample

write.table(t(pca_sample),file=paste0(out_dir,'genotype_pcs.txt'),sep='\t',quote=F,col.names=T,row.names=T)

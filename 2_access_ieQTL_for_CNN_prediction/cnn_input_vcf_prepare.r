
#### The script is used to get all top ieQTL (after genetic interaction effect test)

# define paths
# ieqtl_output_directory='/media/ggj/Guo-4T-MS2/XiaoYanyu/qtl_result/' # ieqtl output directory
# vcf_directory='/media/ggj/My Passport1/Xiaoyanyu/CNN/vcf/OtherTissue/' # vcf output directory
# bfile_directory='/tmpfile/79471/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/GenotypeFiles/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze'
# coordinate_file='/home/ggj/hg38ToHg19.over.chain.gz'

library(argparser, quietly=TRUE)
p <- arg_parser("cnn input vcf prepare")
p <- add_argument(p, "--ieqtl_output_dir", help="")
p <- add_argument(p, "--vcf_dir", help="", default="")
p <- add_argument(p, "--bfile_dir", help="", default="")
p <- add_argument(p, "--coord_file", help="", default="")
p <- add_argument(p, "--run_cell_type", help="", default="")

argv <- parse_args(p)

ieqtl_output_directory = argv$ieqtl_output_dir
vcf_directory = argv$vcf_dir
bfile_directory = argv$bfile_dir
coordinate_file = argv$coord_file
run_cell_type = argv$run_cell_type
run_cell_type <- unlist(strsplit(run_cell_type,','))

#########the next is main_1.r##############

all_top_assoc <- data.frame()
for (i in run_cell_type){
  if (file.exists(paste0(ieqtl_output_directory, i, "/SCCNN-ieQTL.cis_qtl_top_assoc.txt.gz")) == TRUE){
    top_assoc <- read.csv(paste0(ieqtl_output_directory, i, '/SCCNN-ieQTL.cis_qtl_top_assoc.txt.gz'), sep = '\t')
    top_assoc$celltype <- rep(i, nrow(top_assoc))
    all_top_assoc <- rbind(all_top_assoc, top_assoc)
    print(i)
  }
}
all_top_assoc <- all_top_assoc[abs(all_top_assoc$tss_distance) < 20000,] # define the range of CNN predictions

all_top_assoc$gene <- gsub('\\..*', '', all_top_assoc$phenotype_id)
all_top_assoc$index <- paste0(all_top_assoc$celltype, '_', all_top_assoc$gene)
rownames(all_top_assoc) <- all_top_assoc$index

write.table(unique(all_top_assoc$variant_id), file=paste0(vcf_directory,'topsnp.txt'), row.names=F, col.names=F, quote=F)
write.table(all_top_assoc, file=paste0(vcf_directory, 'all_top_assoc.txt'), row.names=F, col.names=T, quote=F, sep='\t')

#########the next is main_2.sh##############
system(paste0("Plink2 --bfile ", bfile_directory, "--r2 --ld-snp-list ", vcf_directory ,"topsnp.txt", "--ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.6", 
               "--out ", vcf_directory),
               intern=FALSE, wait=TRUE)

#########the next is main_3.r##############
LD <- read.table(paste0(vcf_directory, 'plink.ld'), header=T)
save(LD,file=paste0(vcf_directory, 'LD.RData'))

LD1 <- LD[LD$R2>0.6, ]
LD_tmp <- data.frame(top_ieqtl=LD1$SNP_A, variant_id=LD1$SNP_B)

all_top_assoc <- read.csv(paste0(vcf_directory, 'all_top_assoc.txt'), sep='\t')
all_top_assoc_tmp <- data.frame(gene=all_top_assoc$gene, top_ieqtl=all_top_assoc$variant_id)
all_top_assoc_tmp$index <- paste0(all_top_assoc_tmp$gene, '-', all_top_assoc_tmp$top_ieqtl)
all_top_assoc_tmp <- all_top_assoc_tmp[!duplicated(all_top_assoc_tmp$index),]

hb <- merge(LD_tmp, all_top_assoc_tmp, all=T, by="top_ieqtl")
hb <- na.omit(hb)

hb$index <- paste0(hb$variant_id,'-', hb$gene)
hb <- hb[!duplicated(hb$index),]

aa<-matrix(unlist(strsplit(as.character(hb$variant_id),"_")), ncol = 5, byrow = T)
all_predict_vcf <- data.frame(chr = aa[,1], start=aa[,2], end = as.numeric(aa[,2]) + 1,
                              ref = aa[,3], mut = aa[,4], index = hb$index)

options(scipen = 200)

all_predict_vcf_bed <- all_predict_vcf
all_predict_vcf_bed$ref <- NULL
all_predict_vcf_bed$mut <- NULL

write.table(all_predict_vcf_bed, file=paste0(vcf_directory,'unique_pos_bed.txt'), sep='\t', quote=F, row.names=F, col.names=F)

#########the next is main_4.sh##############
dir.create(paste0(vcf_directory, "h19"))
system(paste0("liftOver ", vcf_directory, "unique_pos_bed.txt ", coordinate_file, " ",
              vcf_directory, "h19/unique_pos_bed.txt ", vcf_directory, "h19/'NotMap'"),
       intern=FALSE, wait=TRUE)

#########the next is main_5.r##############
h19_bed <- read.csv(paste0(vcf_directory,'h19/unique_pos_bed.txt'),sep='\t',header=F)
aa<-matrix(unlist(strsplit(as.character(h19_bed$V4),"_")),ncol = 5,byrow = T)
h19_bed$ref <- aa[,3]
h19_bed$mut <- aa[,4]
h19_bed <- h19_bed[,c('V1','V2','V3','ref','mut','V4')]
write.table(h19_bed,file=paste0(vcf_directory,'h19/unique_pos.txt'),sep='\t',quote=F,row.names=F,col.names=F)

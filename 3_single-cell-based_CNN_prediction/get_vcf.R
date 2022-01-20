library(stringi)
library(argparser, quietly=TRUE)
p <- arg_parser("get vcf for cnn input")
p <- add_argument(p, "--vcf_dir", help="", default="")
p <- add_argument(p, "--resource_dir", help="", default="")
p <- add_argument(p, "--output_dir", help="", default="")
argv <- parse_args(p)

vcf_directory = argv$vcf_dir
resource_directory = argv$resource_dir
output_directory = argv$output_dir

temp <- read.table(paste0(vcf_directory, "h19/unique_pos.txt"))
temp$id <- gsub("[-]","\t",temp$V6)
gene <- read.csv(paste0(resource_directory, "geneanno.csv"))
temp$id <- temp$V8
temp_use <- merge(temp,gene,by="id")
vcf_bed_use <- temp_use[,c(1,2,3,5,6,11,12,14)]
vcf_bed_use$V1_2 <- vcf_bed_use$V2 -1
vcf_bed_use$CAGE_representative_TSS_1 <- vcf_bed_use$CAGE_representative_TSS -1
vcf_bed_use <- vcf_bed_use[,c(2,9,3,4,5,6,10,8,7,1)] #4429218
vcf_bed_use$d <- vcf_bed_use$CAGE_representative_TSS - vcf_bed_use$V2
colnames(vcf_bed_use) <- c("chr","variant_position_1","variant_position","reference","allele","chr_2","TSS_1","TSS","strand","gene_id","distance")
vcf_bed_use$reference_length<-stri_length(vcf_bed_use$reference)
vcf_bed_use$allele_length<-stri_length(vcf_bed_use$allele)
vcf_bed_use <- vcf_bed_use[vcf_bed_use$reference_length == 1,]
vcf_bed_use <- vcf_bed_use[vcf_bed_use$allele_length == 1,]
colnames(vcf_bed_use)
#[1] "chr"                "variant_position_1" "variant_position"   "reference"          "allele"            
#[6] "chr_2"              "TSS_1"              "TSS"                "strand"             "gene_id"           
#[11] "distance"           "reference_length"   "allele_length"     
vcf_bed_use_1 <- vcf_bed_use


#-------loop start---strand plus----(from 0kb to 200kb)----------------------------------
for (i in 0:19){
  k= i*10000 #(0kb)
  h= (i+1)*10000   #(200kb)
  outname=paste("From_",k/1000,"kb_To_",h/1000,"kb",sep="")
  File_1 = paste(output_directory, "all_celltype.vcf.bed.sorted.bed_", outname,".closestgene",sep="")
  File_2 = paste(output_directory, "all_celltype_", outname, ".vcf", sep="")
  vcf_bed_use_2 <- vcf_bed_use_1[which(vcf_bed_use_1$distance >= k),]
  vcf_bed_use_2 <- vcf_bed_use_2[which(vcf_bed_use_2$distance < h),]
  dim(vcf_bed_use_2) #159722
  
  vcf_bed_use_2$chr <- gsub("chr","",vcf_bed_use_2$chr)
  vcf_bed_use_2$chr_2 <- gsub("chr","",vcf_bed_use_2$chr_2)
  vcf_bed_use_2[1:2,]
  vcf_bed_use_2 <- vcf_bed_use_2[,-c(12,13)]
  vcf_bed_use_2[1:2,]
  write.table(vcf_bed_use_2,file = File_1,row.names = F,col.names = F,sep = "\t",quote = F)
  
  vcf_use <- vcf_bed_use_2[,c(1,3,4,5)]
  vcf_use$chr <- paste("chr",vcf_use$chr,sep = "")
  vcf_use$tag <- "-"
  vcf_use<- vcf_use[,c(1,2,5,3,4)]
  vcf_use[1:2,]
  write.table(vcf_use,file = File_2,row.names = F,col.names = F,sep = "\t",quote = F)
}
#-------loop start---strand plus----(from -200kb to 0kb)----------------------------------
for(i in 0:19) {
  k= -(i+1)*10000 #(-200kb)
  h= -i*10000  #(0kb)
  outname=paste("From_",k/1000,"kb_To_",h/1000,"kb",sep="")
  File_1 = paste(output_directory, "all_celltype.vcf.bed.sorted.bed_",outname,".closestgene",sep="")
  File_2 = paste(output_directory, "all_celltype_",outname,".vcf",sep="")
  vcf_bed_use_2 <- vcf_bed_use_1[which(vcf_bed_use_1$distance >= k),]
  vcf_bed_use_2 <- vcf_bed_use_2[which(vcf_bed_use_2$distance < h),]
  dim(vcf_bed_use_2) #159722
  
  vcf_bed_use_2$chr <- gsub("chr","",vcf_bed_use_2$chr)
  vcf_bed_use_2$chr_2 <- gsub("chr","",vcf_bed_use_2$chr_2)
  vcf_bed_use_2[1:2,]
  vcf_bed_use_2 <- vcf_bed_use_2[,-c(12,13)]
  vcf_bed_use_2[1:2,]
  write.table(vcf_bed_use_2,file = File_1,row.names = F,col.names = F,sep = "\t",quote = F)
  
  vcf_use <- vcf_bed_use_2[,c(1,3,4,5)]
  vcf_use$chr <- paste("chr",vcf_use$chr,sep = "")
  vcf_use$tag <- "-"
  vcf_use<- vcf_use[,c(1,2,5,3,4)]
  vcf_use[1:2,]
  write.table(vcf_use,file = File_2,row.names = F,col.names = F,sep = "\t",quote = F)
}
##--------------------------------END-----------------------------------

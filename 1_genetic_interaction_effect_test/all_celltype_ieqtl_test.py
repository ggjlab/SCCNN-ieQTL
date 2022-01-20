
### The script is used to compute genetic interaction effect

import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import os
import argparse

# define paths
# PEER_directory = '/home/ggj/QTL/PEER/' # PEER directory
# bfile_directory = '/tmpfile/79471/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/GenotypeFiles/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze'
# expression_bed = '/home/ggj/QTL/GTEx_V8.expression.bed.gz'
# prefix = 'GTEx_V8_samples'
# celltype_list = ['AdultTransverseColon0', 'AdultTransverseColon10', 'AdultTransverseColon11'] # cell types used for ieQTL test
# celltype_tissue_anno = '/home/ggj/QTL/Celltype-GTExSampleNumber-Tissue-Annotation.txt' # cell types and their matched tissue names
# sample_participant_anno = '/home/ggj/QTL/sample_participant_lookup.txt' # bulk samples and their matched participants
# deconvolution_result_directory = '/home/ggj/QTL/deconvolution_result/absolute_fraction/Transformed/' # transformed deconvolution result directory(文件以'组织_fraction.txt'命名)
# ieqtl_output_directory='/media/ggj/Guo-4T-MS2/XiaoYanyu/qtl_result/' # ieqtl output directory

parser = argparse.ArgumentParser(description='all_celltype_ieqtl_test.py')
parser.add_argument('--peer_dir', help='')
parser.add_argument('--bfile_dir', help='')
parser.add_argument('--expr_bed_file', help='')
parser.add_argument('--deconv_result_dir', help='')
parser.add_argument('--celltype_tissue_anno', help='')
parser.add_argument('--sample_participant_anno', help='')
parser.add_argument('--qtl_output_dir', default='.', help='Output directory')
parser.add_argument('--run_cell_type', help='test celltype')

args = parser.parse_args()

PEER_directory = args.peer_dir
bfile_directory = args.bfile_dir
expression_bed = args.expr_bed_file
deconvolution_result_directory = args.deconv_result_dir
celltype_tissue_anno = args.celltype_tissue_anno # cell types and their matched tissue names
sample_participant_anno = args.sample_participant_anno # bulk samples and their matched participants
ieqtl_output_directory = args.qtl_output_dir
celltype_list = str(args.run_cell_type)
celltype_list = celltype_list.split(',')
print(celltype_list)

prefix = 'SCCNN-ieQTL'

# load phenotypes
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(bfile_directory)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

celltype_tissue_mapping = pd.read_table(celltype_tissue_anno,sep='\t',header=None)
celltype_tissue_mapping

for i in celltype_list : 
 print(i+' start')
 tissue = celltype_tissue_mapping[celltype_tissue_mapping[0]==i][1].to_string(index=False).strip()
 #covariates
 covariates_file = PEER_directory + tissue +'/GTEx_V8.combined_covariates.txt'
 covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T  
 #interaction
 ##load interaction_all 
 interaction_file = deconvolution_result_directory+tissue+'_fraction.txt'
 interaction_s_all = pd.read_csv(interaction_file, sep='\t', index_col=0) 
 interaction_s_all.index = [x.replace('.','-') for x in interaction_s_all.index]
 interaction_s = pd.Series(interaction_s_all[i].values, index=interaction_s_all.index)
 #####get genotype sample#####
 mapping = pd.read_table(sample_participant_anno,sep='\t',header=None)
 allsample = phenotype_df.columns.tolist()
 usepaticiant = genotype_df.columns.tolist()
 usesample = [x for x in allsample if mapping[mapping[0]==x][1].to_string(index=False).strip() in usepaticiant]
 ##phenotype
 phenotype_df1 = phenotype_df[usesample]
 ##covariate
 covariates_df1 = covariates_df[covariates_df.index.isin(usesample)]
 ##interaction
 interaction_s1 = interaction_s[interaction_s.index.isin(usesample)]
 #####tissue sample
 choosesample = interaction_s_all.index
 ##phenotype
 phenotype_df2 = phenotype_df1[[x for x in usesample if x in choosesample]]
 ##covariate
 covariates_df2 = covariates_df1[covariates_df1.index.isin(choosesample)]
 ##interaction
 interaction_s2 = interaction_s1[interaction_s1.index.isin(choosesample)]
 covariates_df2 = covariates_df2.reindex(phenotype_df2.columns.tolist()) 
 os.chdir(ieqtl_output_directory + i)
 cis.map_nominal(genotype_df=genotype_df, variant_df=variant_df, phenotype_df=phenotype_df2, 
                 phenotype_pos_df=phenotype_pos_df,
                 covariates_df=covariates_df2, prefix=prefix,  
                 sample_participant_lookup_file=sample_participant_anno,
                 interaction_s=interaction_s2, maf_threshold_interaction=0.05,
                 group_s=None, run_eigenmt=True, output_dir='./')
 print(i+' has been finished') 

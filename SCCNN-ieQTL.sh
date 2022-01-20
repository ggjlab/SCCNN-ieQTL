####### define user input parameter ####
run_script_folder=$(pwd)
expr_bed_file=${run_script_folder}"/data/GTEx_V8.expression.bed.gz"
cell_type_expression_profile=${run_script_folder}"/data/Tissue20_imputation_celltype_group_exp.csv"
deconv_result_dir=${run_script_folder}"/data/cell_type_compositions/"
output_dir=${run_script_folder}"/result/"
celltype_tissue_anno=${run_script_folder}'/data/Celltype-GTExSampleNumber-Tissue-Annotation.txt'
sample_participant_anno=${run_script_folder}'/data/sample_participant_lookup.txt'

function print_help_info() {
	echo "**********this page showes how to use SCCNN-ieQTL*********"
	echo "	  -b, --bfile_dir, genotype data of bulk tissue samples in PLINK format"
	echo "          -bed_file, --expr_bed_file,	phenotype data of bulk tissue samples in bed format, default=${expr_bed_file}"
	echo "          -expr_profile, --expression_profile, cell-type expression profiles from scRNA-seq, default=${cell_type_expression_profile}"
	echo "          -deconv_res, --deconv_result_dir, cell-type compositions across bulk tissue samples, default=${deconv_result_dir}"
	echo "          -output, --output_dir, SCCNN-ieQTL result directory, default=${output_dir}"
	echo "          -celltype_tissue_lookup, --celltype_tissue_anno, the first column is your cell type names, and the second column is the matched tissue names, default=${celltype_tissue_anno}"
	echo "          -sample_participant_lookup, --sample_participant_anno, the first column is your bulk tissue sample names, and the second column is the matched participant names, default=${sample_participant_anno}"
	echo "          -run_tissue, --run_tissue, tissue names you choose to run SCCNN-ieQTL"
	echo "          -h, --help,	help info"
}

if [ $# -ge 1 ]; then
	if [ $1 = "-h" ] || [ $1 = "--help" ]; then
		print_help_info; exit 0;
	fi
fi

while [ $# -ge 2 ]; do
	case "$1" in
		-b|--bfile_dir)
			bfile_directory=$2; shift 2;;
		-bed_file|--expr_bed_file)
			expr_bed_file=$2; shift2;;
		-expr_profile|--expression_profile)
			cell_type_expression_profile=$2; shift2;;
		-deconv_res|--deconv_result_dir)
			deconv_result_dir=$2; shift 2;;
    -output|--output_dir)
	    output_dir=$2; shift 2;;
    -celltype_tissue_lookup|--celltype_tissue_anno)
	    celltype_tissue_anno=$2; shift 2;;
    -sample_participant_lookup|--sample_participant_anno)
	    celltype_tissue_anno=$2; shift 2;;
		-run_tissue|--run_tissue)
	    run_tissue=$2; shift 2;;
		*)
			echo "unknown input parameter!"; exit 1;;
	esac
done

########## define internal path ###############

resource_directory=${run_script_folder}"/resource/"
coordinate_file=${resource_directory}"hg38ToHg19.over.chain.gz"


alias Rscript='/usr/bin/Rscript'
alias Plink2='/home/ggj/plink2'
alias Python3='/home/ggj/anaconda2/envs/py3.7/bin/python'

##### part I: Program_genetic_interaction_effect_test #######

echo 'extracting tissue samples'

mkdir -p ${output_dir}tmpfile/PEER
PEER_directory=${output_dir}"tmpfile/PEER/"

Rscript ${run_script_folder}/1_genetic_interaction_effect_test/extract_matched_sample.r --deconv_result_dir ${deconv_result_dir}\
                               --peer_dir ${PEER_directory}\
                               --expr_bed_file ${expr_bed_file}\
                               --run_tissue ${run_tissue}

echo 'estimating covariates'

tissuename=$(ls ${PEER_directory})
for i in ${tissuename}; do

    echo "${i}"
	
	# peer run
	Rscript ${run_script_folder}/1_genetic_interaction_effect_test/run_PEER.r ${PEER_directory}${i}/choosesample_exp.txt SCCNN-ieQTL 60 \
            --max_iter 200 --output_dir ${PEER_directory}${i};

	# get genotype
	Plink2 --allow-extra-chr --threads 20 -bfile ${bfile_directory} --pca 5 --out ${PEER_directory}${i}/ --keep \
                              ${PEER_directory}${i}/genotype_pc_analyse_sample.txt

	Rscript ${run_script_folder}/1_genetic_interaction_effect_test/genotype_pc.r --deconv_result_dir ${deconv_result_dir} \
                                    --peer_dir ${PEER_directory} \
                                    --out_dir ${PEER_directory}${i}/ \
                                    --tissue ${i};
	# combine covariates
	prefix="SCCNN-ieQTL";
	genotype_pcs="genotype_pcs.txt";

	Python3 ${run_script_folder}/1_genetic_interaction_effect_test/combine_covariates.py ${PEER_directory}${i}/${prefix}.PEER_covariates.txt ${prefix} \
                             --genotype_pcs ${PEER_directory}${i}/${genotype_pcs} --output_dir ${PEER_directory}${i}/

done


###

echo 'genetic interaction calculation start'

for tissue in $(echo ${run_tissue} | sed "s/,/ /g")
do
  all_celltype_list=${all_celltype_list}" "$(head -n 1 $(find ${deconv_result_dir} -type f -name "${tissue}*"))
done

for i in ${all_celltype_list}; do
  mkdir -p ${output_dir}ieQTL/${i}
done
qtl_output_directory=${output_dir}"ieQTL/"

all_celltype_list2=$(echo ${all_celltype_list} | sed 's/[ ][ ]*/,/g')
echo $all_celltype_list2

Python3 ${run_script_folder}/1_genetic_interaction_effect_test/all_celltype_ieqtl_test.py --peer_dir ${PEER_directory} --bfile_dir ${bfile_directory} \
                  --expr_bed_file ${expr_bed_file} --deconv_result_dir ${deconv_result_dir} \
                  --celltype_tissue_anno ${celltype_tissue_anno} \
                  --sample_participant_anno ${sample_participant_anno} \
                  --qtl_output_dir ${qtl_output_directory} \
                  --run_cell_type ${all_celltype_list2} 


##### part II: Program_CNN_input_vcf_prepare #######

echo 'preparing vcf files for CNN prediction'

mkdir -p ${output_dir}tmpfile/CNN_input
vcf_directory=${output_dir}"tmpfile/CNN_input/"

Rscript ${run_script_folder}/2_access_ieQTL_for_CNN_prediction/cnn_input_vcf_prepare.r --ieqtl_output_dir ${qtl_output_directory} \
                                --vcf_dir ${vcf_directory} \
                                --bfile_dir ${bfile_directory} \
                                --coord_file ${coordinate_file} \
                                --run_cell_type ${all_celltype_list2} 
								
##### part III: Program_CNN_predict #######
Rscript ./3_single-cell-based_CNN_prediction/get_vcf.R --vcf_dir ${vcf_directory} \
                --resource_dir ${resource_directory} \
                --output_dir ${vcf_directory}

echo 'downloading CNN resources'

if [ -d ${resource_directory}CNN/resources ]; then
  echo "resources cnn resources has been downloaded"
else
  wget -P ${resource_directory}CNN http://deepsea.princeton.edu/media/code/expecto/resources_20190807.tar.gz
  tar -xf ${resource_directory}CNN/resources_20190807.tar.gz -C ${resource_directory}CNN/
fi

mkdir -p ${output_dir}tmpfile/CNN_model/models
mkdir -p ${output_dir}tmpfile/CNN_model/pred_target
mkdir -p ${output_dir}CNN_prediction
model_save_directory=${output_dir}tmpfile/CNN_model/models/
model_list_directory=${output_dir}tmpfile/CNN_model/
pred_target_directory=${output_dir}tmpfile/CNN_model/pred_target/
cnn_predict_directory=${output_dir}CNN_prediction/


## model traning

echo 'traning single-cell based models'

INDEX=1
for i in ${all_celltype_list}
do
  Python3 ${resource_directory}CNN/train.py --expFile ${cell_type_expression_profile} \
          --targetIndex $INDEX --output ${model_save_directory}model.$i --evalFile ${pred_target_directory}model.$i.csv \
          --inputFile ${resource_directory}CNN/resources/Xreducedall.2002.npy \
          --annoFile ${resource_directory}CNN/resources/geneanno.csv
  INDEX=$((INDEX + 1))
done

##create used modellist file 

Rscript ${run_script_folder}/3_single-cell-based_CNN_prediction/create_modellist.r \
       --model_save_directory ${model_save_directory} \
       --model_list_directory ${model_list_directory}

## computes the chromatin effects of the variants

echo 'computing the chromatin effects of the variants'

cd ${vcf_directory}
for file_name in ${vcf_directory}*vcf
do
Python3 ${resource_directory}CNN/chromatin.py $file_name --cuda \
       --hg19_fa ${resource_directory}CNN/resources/hg19.fa \
       --deepsea ${resource_directory}CNN/resources/deepsea.beluga.pth
done
cd ${run_script_folder}

## model prediction

echo 'computing variant effects for each cell type'

for ((i=0;i<=19;i++))
do
	Python3 ${resource_directory}CNN/predict.py --coorFile ${vcf_directory}all_celltype_From_$[$i *10]kb_To_$[$[$i *10] +10]kb.vcf \
                       --geneFile ${vcf_directory}all_celltype.vcf.bed.sorted.bed_From_$[$i *10]kb_To_$[$[$i *10] +10]kb.closestgene \
                       --snpEffectFilePattern ${vcf_directory}all_celltype_From_$[$i *10]kb_To_$[$[$i *10] +10]kb.vcf.shift_SHIFT.diff.h5 \
                       --modelList ${model_list_directory}modellist \
                       --output ${cnn_predict_directory}all_celltype_From_$[$i *10]kb_To_$[$[$i *10] +10]kb_output.csv
done

###### part IV #################

echo 'inferring confident causal genetic regulation'

Rscript ./4_confident_causal_regulation_inferring/ --vcf_dir ${vcf_directory} \
                                                   --cnn_predict_dir ${cnn_predict_directory} \
                                                   --output_dir ${output_dir} \
                                                   --resource_dir ${resource_directory} \
                                                   --qtl_output_dir ${qtl_output_directory} \
                                                   --run_tissue ${run_tissue}


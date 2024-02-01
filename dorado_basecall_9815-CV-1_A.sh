#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=A_dorado
#SBATCH --mail-user=weishwu@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --gpus-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --time=240:00:00
#SBATCH --account=bis0
#SBATCH --partition=gpu
#SBATCH --output=/nfs/turbo/lsa-bis/DNAmet/analysis/wgbs/ont_UW/pod5/%x-%j.log

run_id=9815-CV-1_A

####################################

module load singularity

out_dir=/nfs/turbo/lsa-bis/DNAmet/analysis/wgbs/ont_UW/analysis/
dorado_bin=/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/dorado-0.5.0-linux-x64/bin/dorado
pod5_sif=/nfs/turbo/lsa-bis/envs/pod5.sif
dorado_model_bc=/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/dorado-0.5.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0
dorado_model_bc_mod=/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/dorado-0.5.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v1
pod5_dir=/nfs/turbo/lsa-bis/DNAmet/analysis/wgbs/ont_UW/pod5

#ont_dir=`ls -d /nfs/turbo/lsa-bis/DNAmet/analysis/wgbs/ont_UW/${run_id}/*/*`

#sample_id=`echo ${ont_dir} |sed 's:/$::g'| awk -F'/' '{print $(NF-1)"-"$(NF-2)}'`
sample_id=${run_id}

cd ${out_dir}

(#singularity run -B /nfs/ ${pod5_sif} pod5 \
#	convert fast5 \
#	${ont_dir}/fast5_pass/*.fast5 \
#	--output ${pod5_dir}/${sample_id}/

${dorado_bin} basecaller \
     ${dorado_model_bc} \
     --modified-bases-models ${dorado_model_bc_mod} \
     ${pod5_dir}/${sample_id}/ \
     > ${out_dir}/Sample_${sample_id}.doradoSup.basecall.bam
) 2>&1|tee >${out_dir}/logs/Sample_${sample_id}.doradoSup.basecall.log


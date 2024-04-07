#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=barcode${f}_sup
#SBATCH --mail-user=weishwu@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --gpus-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G
#SBATCH --time=240:00:00
#SBATCH --account=cgates1
#SBATCH --partition=gpu
#SBATCH --output=/nfs/mm-isilon/bioinfcore/ActiveProjects/Freeman_freemanz_ONT7-Amplicon_weishwu_9946-ZF/outputs/logs/%x-%j.log

smp_id=barcode${f}

####################################

module load singularity

pjdir=/nfs/mm-isilon/bioinfcore/ActiveProjects/Freeman_freemanz_ONT7-Amplicon_weishwu_9946-ZF/
nf_conf=${pjdir}/scripts/nextflow_resource.cfg
pod5_dir=${pjdir}/inputs/10494-ZF/pod5s_10494-ZF/pod5_pass/
samplesheet=${pjdir}/inputs/10494-ZF/samplesheet_10494-ZF.csv
smp_name=`awk -v p=${smp_id} -F ',' '{if ($2 == p) {print $3}}' ${samplesheet}`

source /nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/Anaconda3/bin/activate nextflow

bcdir=${pjdir}/outputs/basecalls_sup/
mkdir -p ${bcdir}
mkdir -p ${bcdir}/${smp_name}

cd ${bcdir}/${smp_name}
export NXF_SINGULARITY_CACHEDIR=${pjdir}/env/nextflow_singularity_images/

(nextflow run epi2me-labs/wf-basecalling -revision prerelease \
--input ${pod5_dir}/${smp_id}/ \
--fastq_only \
--out_dir ${bcdir}/${smp_name} \
--sample_name ${smp_name} \
--basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
--dorado_ext pod5 \
--qscore_filter 10 \
--merge_threads 20 \
--stats_threads 20 \
-profile singularity \
-c ${nf_conf} ) 2>&1|tee > ${pjdir}/outputs/logs/basecall_sup.${smp_name}.log



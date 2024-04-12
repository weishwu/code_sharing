#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=barcode75_sup
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
#SBATCH --output=/nfs/mm-isilon/bioinfcore/ActiveProjects/Bachman_mikebach_ONT1_weishwu/outputs/logs/%x-%j.log

smp_ids="barcode75"

####################################

module load singularity

pjdir=/nfs/mm-isilon/bioinfcore/ActiveProjects/Bachman_mikebach_ONT1_weishwu
nf_conf=${pjdir}/scripts/nextflow_resource.cfg
pod5_dir=${pjdir}/inputs/10548-SM/pod5_files_10548-SM/pod5_pass
samplesheet=${pjdir}/inputs/10548-SM/samplesheet_10548-SM.csv

source /nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/Anaconda3/bin/activate nextflow

bcdir=${pjdir}/outputs/basecalls_sup_test
mkdir -p ${bcdir}

export NXF_SINGULARITY_CACHEDIR=${pjdir}/env/nextflow_singularity_images/

for smp_id in ${smp_ids}
do
smp_name=`awk -v p=${smp_id} -F ',' '{if ($2 == p) {print $3}}' ${samplesheet}`

mkdir -p ${bcdir}/${smp_name}
cd ${bcdir}/${smp_name}
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
-c ${nf_conf} ) 2>&1|tee > ${pjdir}/outputs/logs/basecall_sup.test.${smp_name}.log
done


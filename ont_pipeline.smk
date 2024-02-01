samples = ['5555C-A-9815-CV-1.doradoSup','5555C-B-9815-CV-1.doradoSup']

out_dir = '/nfs/turbo/lsa-bis/DNAmet/analysis/wgbs/ont_UW/analysis/'
log_dir = out_dir + 'logs/'
bmk_dir = out_dir + 'benchmarks/'
tmp_dir = '/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/temp/'
ref_dir = '/nfs/turbo/lsa-bis/reference/'

nanofilt_len = 100,
nanofilt_headcrop = 50,
nanofilt_tailcrop = 20,
nanofilt_baseq = 10,

clair3_models = {}
clair3_models['5555C-A-9815-CV-1.doradoSup'] = '/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/rerio/clair3_models/r1041_e82_400bps_sup_v430'
clair3_models['5555C-B-9815-CV-1.doradoSup'] = '/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/rerio/clair3_models/r1041_e82_400bps_sup_v430'


dorado_bin = '/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/dorado-0.5.0-linux-x64/bin/dorado'
modkit_bin = '/nfs/turbo/lsa-bis/envs/modkit/modkit'
dss_code = '/nfs/turbo/lsa-bis/DNAmet/analysis/scripts/dss_ont_asm.R'
haps_allele_count_comp_code = '/nfs/turbo/lsa-bis/DNAmet/analysis/scripts/haps_allele_count_comp.py'
genome_fasta = '/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/GATKbundles/hg38_noChr11Random_primaryContigs/Homo_sapiens_assembly38_noChr11random_primaryContigs.fasta'
genome_table = '/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/GATKbundles/hg38_noChr11Random_primaryContigs/Homo_sapiens_assembly38_noChr11random_primaryContigs_genometable'

import os
import glob

def clair_model(x):
  return(clair3_models[x])


###############################################

rule all:
   input:
      expand(out_dir + "Sample_{sample}.snikt.{zoom}.html", zoom = ['50_20','80_30','120_50', '120_100'], sample = samples),
      expand(out_dir + "Sample_{sample}.mosdepth.summary.txt", sample = samples),
      expand(out_dir + "Sample_{sample}.GQ20.test_callDMR.bedGraph", sample = samples),
      expand(out_dir + "Sample_{sample}.GQ20.test_DMLtest.bw", sample = samples),
      expand(out_dir + 'Sample_{sample}.dorado_aligned.GQ20.hap{hap}.CpG.meth.bw', sample = samples, hap = ['1','2']),
 



rule bam_to_fastq:
   input: out_dir + 'Sample_{sample}.basecall.bam',
   output: out_dir + 'Sample_{sample}.basecall.fastq',
   singularity: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_2',
   threads: 20
   log: log_dir + 'Sample_{sample}.bam2fq.log',
   benchmark: bmk_dir + 'Sample_{sample}.bam2fq.txt',
   shell:
      '''(samtools fastq -@ {threads} {input} > {output}) 2>&1|tee >{log}'''


rule snikt_qc:
   input: out_dir + 'Sample_{sample}.basecall.fastq',
   output:
      expand(out_dir + "Sample_{sample}.snikt.{zoom}.html", zoom = ['50_20','80_30','120_50', '120_100'], sample = ['{sample}']),
   singularity: "docker://weishwu/snikt:0.5.0",
   log: log_dir + "Sample_{sample}.snikt.log",
   benchmark: bmk_dir + "Sample_{sample}.snikt.txt",
   params:
      work_dir = out_dir
   shell:
      '''(cd {params.work_dir}
      snikt.R --zoom5 80 --zoom3 30 --notrim -o Sample_{wildcards.sample}.snikt.80_30 {input}
      snikt.R --zoom5 50 --zoom3 20 --notrim -o Sample_{wildcards.sample}.snikt.50_20 {input}
      snikt.R --zoom5 120 --zoom3 50 --notrim -o Sample_{wildcards.sample}.snikt.120_50 {input}
      snikt.R --zoom5 120 --zoom3 100 --notrim -o Sample_{wildcards.sample}.snikt.120_100 {input}
) 2>&1|tee >{log}'''


rule nanofilt_trim:
   input: out_dir + 'Sample_{sample}.basecall.fastq',
   output: out_dir + 'Sample_{sample}.basecall.filt.fastq',
   singularity: "docker://quay.io/biocontainers/nanofilt:2.8.0--py_0",
   log: log_dir + "Sample_{sample}.nanofilt.log",
   benchmark: bmk_dir + "Sample_{sample}.nanofilt.txt"
   params:
      nanofilt_len = nanofilt_len,
      nanofilt_headcrop = nanofilt_headcrop,
      nanofilt_tailcrop = nanofilt_tailcrop,
      nanofilt_baseq = nanofilt_baseq,
   shell:
      '''(NanoFilt -l {params.nanofilt_len} --headcrop {params.nanofilt_headcrop} --tailcrop {params.nanofilt_tailcrop} -q {params.nanofilt_baseq} {input} >{output}) 2>&1|tee >{log}'''


rule fastq_to_bam:
   input: out_dir + 'Sample_{sample}.basecall.filt.fastq',
   output: out_dir + 'Sample_{sample}.basecall.filt.bam',
   singularity: 'docker://quay.io/nf-core/gatk:4.4.0.0',
   log: log_dir + 'Sample_{sample}.fq2bam.log',
   benchmark: bmk_dir + 'Sample_{sample}.fq2bam.txt',
   resources:
      mem_gb=120,
   params:
      tmp_dir = tmp_dir,
   shell:
      '''(gatk FastqToSam --java-options '-Xmx{resources.mem_gb}G' --FASTQ {input} --OUTPUT {output} --SAMPLE_NAME {wildcards.sample} --TMP_DIR {params.tmp_dir}) 2>&1|tee >{log}'''


rule donor_bam_sort_n:
   input: out_dir + 'Sample_{sample}.basecall.bam',
   output: out_dir + 'Sample_{sample}.basecall.nsort.bam',
   singularity: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_2',
   threads: 20
   log: log_dir + 'Sample_{sample}.basecall.bam.sort.log',
   benchmark: bmk_dir + 'Sample_{sample}.basecall.bam.sort.txt',
   shell:
      '''(samtools sort -@ {threads} -n -O BAM -o {output} {input}) 2>&1|tee >{log}'''


rule acceptor_bam_sort_n:
   input: out_dir + 'Sample_{sample}.basecall.filt.bam',
   output: out_dir + 'Sample_{sample}.basecall.filt.nsort.bam',
   singularity: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_2',
   threads: 20
   log: log_dir + 'Sample_{sample}.basecall.filt.bam.nsort.log',
   benchmark: bmk_dir + 'Sample_{sample}.basecall.filt.bam.nsort.txt',
   shell:
      '''(samtools sort -@ {threads} -n -O BAM -o {output} {input}) 2>&1|tee >{log}'''


rule repair_bam:
   input:
      donor = out_dir + 'Sample_{sample}.basecall.nsort.bam',
      acceptor = out_dir + 'Sample_{sample}.basecall.filt.nsort.bam',
   output: out_dir + 'Sample_{sample}.basecall.filt.nsort.repaired.bam',
   params:
      code = modkit_bin,
      logfile = log_dir + 'Sample_{sample}.bam_repair.logfile',
   log: log_dir + 'Sample_{sample}.bam_repair.log',
   benchmark: bmk_dir + 'Sample_{sample}.bam_repair.txt',
   shell:
      '''({params.code} repair --donor-bam {input.donor} --acceptor-bam {input.acceptor} --log-filepath {params.logfile} --output-bam {output}) 2>&1|tee >{log}'''


rule read_align:
   input: out_dir + 'Sample_{sample}.basecall.filt.nsort.repaired.bam',
   output: out_dir + 'Sample_{sample}.dorado_aligned.bam',
   params:
      code = dorado_bin,
      genome_fasta = genome_fasta,
   threads: 40
   log: log_dir + 'Sample_{sample}.dorado_align.log',
   benchmark: bmk_dir + 'Sample_{sample}.dorado_align.txt',
   shell:
      '''({params.code} aligner --max-reads 100000000 -t {threads} {params.genome_fasta} {input} > {output}) 2>&1|tee >{log}'''


rule bam_sort_c:
   input: out_dir + 'Sample_{sample}.dorado_aligned.bam',
   output:
      bam = out_dir + 'Sample_{sample}.dorado_aligned.sort.bam',
      index = out_dir + 'Sample_{sample}.dorado_aligned.sort.bam.bai',
   singularity: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_2',
   threads: 20
   log: log_dir + 'Sample_{sample}.dorado_align.bam.sort.log',
   benchmark: bmk_dir + 'Sample_{sample}.dorado_align.bam.sort.txt',
   shell:
      '''(samtools sort -@ {threads} -O BAM -o {output.bam} {input}
      samtools index {output.bam}) 2>&1|tee >{log}'''


rule aln_cov:
   input:
      bam = out_dir + "Sample_{sample}.dorado_aligned.sort.bam",
      index = out_dir + "Sample_{sample}.dorado_aligned.sort.bam.bai",
   output:
      out_dir + "Sample_{sample}.mosdepth.summary.txt",
   singularity: "docker://quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2",
   threads: 16
   log: log_dir + "Sample_{sample}.dorado_aligned.sort.bam.coverage.log",
   benchmark: bmk_dir + "Sample_{sample}.dorado_aligned.sort.bam.coverage.txt",
   params:
      prefix = out_dir + "Sample_{sample}",
   shell:
      '''(mosdepth -t {threads} {params.prefix} {input.bam}) 2>&1|tee >{log}'''


rule call_variants:
   input: out_dir + 'Sample_{sample}.dorado_aligned.sort.bam',
   output: out_dir + "Sample_{sample}.clair3/merge_output.vcf.gz",
   singularity: "docker://quay.io/biocontainers/clair3:1.0.5--py39hf5e1c6e_0",
   log: log_dir + "Sample_{sample}.clair3.call_variants.log",
   benchmark: bmk_dir + "Sample_{sample}.call_variants.txt",
   threads: 20
   params:
      genome_fasta = genome_fasta,
      model = lambda wildcards: clair_model(wildcards.sample),
      outdir = out_dir + 'Sample_{sample}.clair3/',
      tmp_dir = tmp_dir,
   shell:
      '''(export TMPDIR={params.tmp_dir}
  run_clair3.sh --bam_fn={input} \
  --ref_fn={params.genome_fasta} \
  --output={params.outdir} \
  --threads={threads} --platform=ont \
  --model_path={params.model}) 2>&1|tee >{log}'''


rule pass_variants:
   input: out_dir + "Sample_{sample}.clair3/merge_output.vcf.gz",
   output: out_dir + "Sample_{sample}.clair3/merge_output.pass.vcf",
   log: log_dir + "Sample_{sample}.clair3.variants.pass_filter.log",
   benchmark: bmk_dir + "Sample_{sample}.variants_filter.txt",
   shell:
      '''(gunzip -c {input} | awk '{{if (($1~"^#" || $7=="PASS")) {{print $0}}}}' > {output}) 2>&1|tee >{log}'''


rule gzip_variants:
   input: out_dir + "Sample_{sample}.clair3/merge_output.pass.vcf",
   output: out_dir + "Sample_{sample}.clair3/merge_output.pass.vcf.gz",
   singularity: 'docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0',
   shell:
      '''(bgzip {input}
      tabix {output})'''


rule filter_variants_GQ:
   input: out_dir + "Sample_{sample}.clair3/merge_output.pass.vcf.gz",
   output: out_dir + "Sample_{sample}.clair3/merge_output.pass.GQ20.vcf",
   singularity: 'docker://quay.io/biocontainers/bcftools:1.19--h8b25389_0',
   log: log_dir + "Sample_{sample}.pass.GQ20.log",
   shell:
      '''(
      bcftools view -i "FORMAT/GQ>=20" -Ov -o {output} {input}) 2>&1|tee >{log}'''


rule phase_variants:
   input: 
      vcf = out_dir + "Sample_{sample}.clair3/merge_output.pass.GQ20.vcf",
      bam = out_dir + 'Sample_{sample}.dorado_aligned.sort.bam',
   output: 
      vcf = out_dir + "Sample_{sample}.pass.phased.GQ20.vcf",
   singularity: "docker://quay.io/biocontainers/whatshap:2.1--py39h1f90b4d_0",
   log: log_dir + "Sample_{sample}.variants.whatshap.log",
   benchmark: bmk_dir + "Sample_{sample}.variants_phase.txt",
   params:
      genome_fasta = genome_fasta,
   shell:
      '''(whatshap phase --ignore-read-groups --indels --reference {params.genome_fasta} -o {output.vcf} {input.vcf} {input.bam}) 2>&1|tee >{log}'''


rule modkit_extract:
   input: out_dir + 'Sample_{sample}.dorado_aligned.sort.bam',
   output: 
      out_dir + 'Sample_{sample}.dorado_aligned.modkit_extract.CpG.tsv',  
   params:
      code = modkit_bin,
      genome_fasta = genome_fasta,
   threads: 20
   log: log_dir + 'Sample_{sample}.dorado_aligned.modkit_extract.CpG.log',
   shell:
      '''({params.code} extract --cpg --mapped-only --threads {threads} --reference {params.genome_fasta} {input} {output}) 2>&1|tee >{log}'''



rule index_methcall:
   input: out_dir + 'Sample_{sample}.dorado_aligned.modkit_extract.CpG.tsv',
   output:
      bedgz = out_dir + "Sample_{sample}.methylation_calls.bed.gz",
      index = out_dir + "Sample_{sample}.methylation_calls.bed.gz.tbi", 
   singularity: "docker://weishwu/nanomethphase:1.2.0",
   log: log_dir + "Sample_{sample}.index_methcall.log",
   benchmark: bmk_dir + "Sample_{sample}.index_methcall.txt",
   threads: 20,
   shell:
      '''(set +eu
      source /opt/conda/bin/activate nanomethphase

      echo -e "chrom\tpos\tstrand\tpos_in_strand\treadname\tread_strand\tprob_0\tprob_1\tcalled_label\tk_mer" >{input}.reformatted.tsv

      awk '{{OFS="\t";if (($12=="m")) {{print $4,$3,$6,"NA",$1,"NA",1-$11,$11,"NA",$14}}}}' {input} >> {input}.reformatted.tsv

      python /usr/share/NanoMethPhase/nanomethphase.py methyl_call_processor --tool_and_callthresh deepsignal:0.6 -mc {input}.reformatted.tsv -t {threads} | sort -k1,1 -k2,2n -k3,3n | bgzip > {output.bedgz}
      tabix -p bed {output.bedgz}) 2>&1|tee >{log}'''


rule halplotype_meth:
   input:
      aln = out_dir + 'Sample_{sample}.dorado_aligned.sort.bam',
      meth = out_dir + "Sample_{sample}.methylation_calls.bed.gz",
      vcf = out_dir + "Sample_{sample}.pass.phased.GQ20.vcf",
   output: expand(out_dir + "Sample_{sample}.GQ20_NanoMethPhase_HP{hap}_MethylFrequency.tsv", hap = ['1','2'], sample = ['{sample}']),
   singularity: "docker://weishwu/nanomethphase:1.2.0",
   params:
      genome_fasta = genome_fasta,
      out_prefix = out_dir + "Sample_{sample}.GQ20",
   threads: 20
   log: log_dir + "Sample_{sample}.meth_haplotype.GQ20.log",
   benchmark: bmk_dir + "Sample_{sample}.meth_haplotype.txt",
   shell:
      '''(set +eu
      source /opt/conda/bin/activate nanomethphase
      python /usr/share/NanoMethPhase/nanomethphase.py phase --overwrite --include_indels -b {input.aln} -v {input.vcf} -mc {input.meth} -r {params.genome_fasta} -o {params.out_prefix} -of methylcall -t {threads}) 2>&1|tee >{log}'''




# dmr calling using adjusted parameters: p0.01, smoothing pan 1000
rule diff_meth_relax:
   input:
      expand(out_dir + "Sample_{sample}.GQ20_NanoMethPhase_HP{hap}_MethylFrequency.tsv", hap = ['1','2'], sample = ['{sample}']),
   output: 
      dmr = out_dir + "Sample_{sample}.GQ20.test_callDMR.txt",
      test = out_dir + "Sample_{sample}.GQ20.test_DMLtest.txt"
   singularity: "docker://weishwu/nanomethphase:1.2.0",
   log: log_dir + "Sample_{sample}.GQ20.diff_meth.test.log",
   benchmark: bmk_dir + "Sample_{sample}.diff_meth.test.txt",
   params:
      out_dir = out_dir,
      prefix = "Sample_{sample}.GQ20.test",
   threads: 10
   shell:
      '''(set +eu
      source /opt/conda/bin/activate nanomethphase
      python /usr/share/NanoMethPhase/nanomethphase.py dma -c 1,2,4,5,7 -ca {input[0]} -co {input[1]} -o {params.out_dir} -op {params.prefix} --overwrite --pval_cutoff 0.01 --smoothing_span 1000) 2>&1|tee >{log}'''


rule prep_igv_tracks:
   input:
     dmr = out_dir + "Sample_{sample}.GQ20.test_callDMR.txt",
     test = out_dir + "Sample_{sample}.GQ20.test_DMLtest.txt",
   output:
     dmr = out_dir + "Sample_{sample}.GQ20.test_callDMR.bedGraph",
     beta = out_dir + "Sample_{sample}.GQ20.test_DMLtest.bw",
   singularity: 'docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0',
   log: log_dir + '{sample}_igv_tracks_prep.log',
   params:
      genome_table = genome_table
   shell:
      '''(awk '{{OFS="\t"; if (($1 == "chr")) {{print "#"$1,$2,$3,$8,$4,$5,$6,$7,$9}} else {{print $1,$2,$3,$8,$4,$5,$6,$7,$9}}}}' {input.dmr} | awk '{{OFS="\t";if (($4 < 0)) {{$4 = 0 - $4; print $0}} else {{print $0}}}}' > {output.dmr}
      sed '1d' {input.test} | awk '{{OFS="\t"; if (($5 < 0)) {{print $1,$2,$2+1,0-$5}} else {{print $1,$2,$2+1,$5}}}}' > {output.beta}.tmp
      bedGraphToBigWig {output.beta}.tmp {params.genome_table} {output.beta}
      rm {output.beta}.tmp) 2>&1|tee >{log}'''


rule haplotag_bam:
   input: 
      vcf = out_dir + "Sample_{sample}.pass.phased.GQ20.vcf",
      bam = out_dir + "Sample_{sample}.dorado_aligned.sort.bam",
   output: 
      bam = out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.haplotagged.bam",
      list = out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.haplolist.tsv",
   singularity: "docker://quay.io/biocontainers/whatshap:2.1--py39h1f90b4d_0",
   log: log_dir + "{sample}.sorted.GQ20.haplotagged.log",
   benchmark: bmk_dir + "{sample}.sorted.haplotagged.txt",
   params:
      genome_fasta = genome_fasta,
   threads: 20
   shell:
      '''(if [ ! -f {input.vcf}.gz ]; then
      bgzip -c {input.vcf} > {input.vcf}.gz
      tabix {input.vcf}.gz
      fi
      whatshap haplotag -o {output.bam} --reference {params.genome_fasta} --ignore-read-groups --output-haplotag-list {output.list} --ploidy 2 --output-threads {threads} {input.vcf}.gz {input.bam}) 2>&1|tee >{log}'''


rule haplosplit_bam:
   input: 
      bam = out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.haplotagged.bam",
      list = out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.haplolist.tsv",
   output: 
      bam1 = out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.hap1.bam",
      bam2 = out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.hap2.bam",
   singularity: "docker://quay.io/biocontainers/whatshap:2.1--py39h1f90b4d_0",
   log: log_dir + "{sample}.sorted.GQ20.haplosplit.log",
   benchmark: bmk_dir + "{sample}.sorted.haplosplit.txt",
   threads: 20
   shell:
      '''(whatshap split --output-h1 {output.bam1} --output-h2 {output.bam2} {input.bam} {input.list}) 2>&1|tee >{log}'''


rule index_haplobam:
   input: out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.hap{hap}.bam",
   output: out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.hap{hap}.bam.bai",
   singularity: 'docker://quay.io/biocontainers/samtools:1.17--hd87286a_2',
   threads: 20
   shell:
      '''(samtools index -@ {threads} {input})'''


rule modkit_pileup_haps:
   input:
      bam = out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.hap{hap}.bam",
      bai = out_dir + "Sample_{sample}.dorado_aligned.sort.GQ20.hap{hap}.bam.bai",
   output: 
      out_dir + 'Sample_{sample}.dorado_aligned.GQ20.hap{hap}.modkit_pileup.bed',  
   params:
      code = modkit_bin,
      genome_fasta = genome_fasta,
   threads: 20
   log: log_dir + 'Sample_{sample}.dorado_aligned.GQ20.hap{hap}.modkit_extract.log',
   shell:
      '''({params.code} pileup -t {threads} --cpg --ref {params.genome_fasta} --ignore h --combine-strands {input.bam} {output}) 2>&1|tee >{log}'''


rule meth_haps_bigwig:
   input: out_dir + 'Sample_{sample}.dorado_aligned.GQ20.hap{hap}.modkit_pileup.bed', 
   output: out_dir + 'Sample_{sample}.dorado_aligned.GQ20.hap{hap}.CpG.meth.bw',
   singularity: 'docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0'
   params:
      genome_table = genome_table, 
   shell:
      '''(awk '{{OFS="\t"; print $1,$2,$3,$11}}' {input} | sort -k1,1 -k2,2n > {output}.tmp.bed
      bedGraphToBigWig {output}.tmp.bed {params.genome_table} {output}
      rm {output}.tmp.bed)''' 



rule diff_meth_relax_p1e5:
   input:
      expand(out_dir + "Sample_{sample}.GQ20_NanoMethPhase_HP{hap}_MethylFrequency.tsv", hap = ['1','2'], sample = ['{sample}']),
   output:
      dmr = out_dir + "Sample_{sample}.GQ20_p1e5_callDMR.txt",
      test = out_dir + "Sample_{sample}.GQ20_p1e5_DMLtest.txt"
   singularity: "docker://weishwu/nanomethphase:1.2.0",
   log: log_dir + "Sample_{sample}.diff_meth.GQ20.p1e5.log",
   benchmark: bmk_dir + "Sample_{sample}.diff_meth.GQ20.p1e5.txt",
   params:
      out_dir = out_dir,
      prefix = "Sample_{sample}.GQ20_p1e5",
   threads: 10
   shell:
      '''(set +eu
      source /opt/conda/bin/activate nanomethphase
      python /usr/share/NanoMethPhase/nanomethphase.py dma -c 1,2,4,5,7 -ca {input[0]} -co {input[1]} -o {params.out_dir} -op {params.prefix} --overwrite --pval_cutoff 0.00001 --smoothing_span 1000) 2>&1|tee >{log}'''



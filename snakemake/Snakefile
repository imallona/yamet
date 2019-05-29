# todo
# ensure the bam file is not empty after subset by chrom


# import os.path as op
import os
import glob
import re
# configfile: "config.yaml"
# validate(config, "config.schema.yaml")

## yaml tests start ##########

import yaml

print('note the bamfiles are retrieved using 01_encode_bulk_run, and that the hmm segmentations as well')

# roadmap ids
# https://egg2.wustl.edu/roadmap/web_portal/meta.html
# https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
# 15 states HMM

config = """
ENCFF112TXF:
 exp: https://www.encodeproject.org/experiments/ENCSR888FON/'
 bam: https://www.encodeproject.org/files/ENCFF112TXF/@@download/ENCFF112TXF.bam
 assembly: GRCh38
 genotype: IMR90
 replicate: 1
 sequencing: single
 roadmap_id: E017
 roadmap_hmm: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E017_15_coreMarks_hg38lift_mnemonics.bed.gz

ENCFF957OIM:
 exp: https://www.encodeproject.org/experiments/ENCSR881XOU/
 bam: https://www.encodeproject.org/files/ENCFF957OIM/@@download/ENCFF957OIM.bam
 assembly: GRCh38
 genotype: HepG2
 replicate : 1
 sequencing: paired
 roadmap_id : E118
 roadmap_hmm: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E118_15_coreMarks_hg38lift_mnemonics.bed.gz

ENCFF572KNK :
 exp: https://www.encodeproject.org/experiments/ENCSR881XOU/
 bam: https://www.encodeproject.org/files/ENCFF572KNK/@@download/ENCFF572KNK.bam
 assembly: GRCh38
 genotype: HepG2
 replicate : 2
 sequencing: paired
 roadmap_id : E118
 roadmap_hmm: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E118_15_coreMarks_hg38lift_mnemonics.bed.gz

ENCFF193RVP :
 exp: https://www.encodeproject.org/experiments/ENCSR550RTN/
 bam: https://www.encodeproject.org/files/ENCFF193RVP/@@download/ENCFF193RVP.bam
 assembly: GRCh38
 genotype: HeLa-S3
 replicate : 1
 sequencing: paired
 roadmap_id: E117
 roadmap_hmm: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E117_15_coreMarks_hg38lift_mnemonics.bed.gz

ENCFF845VFH :
 exp: https://www.encodeproject.org/experiments/ENCSR550RTN/
 bam: https://www.encodeproject.org/files/ENCFF845VFH/@@download/ENCFF845VFH.bam
 assembly: GRCh38
 genotype: HeLa-S3
 replicate: 2
 sequencing: paired
 roadmap_id: E117
 roadmap_hmm: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E117_15_coreMarks_hg38lift_mnemonics.bed.gz

ENCFF079RGH :
 exp: https://www.encodeproject.org/experiments/ENCSR440MLE/
 bam: https://www.encodeproject.org/files/ENCFF079RGH/@@download/ENCFF079RGH.bam
 assembly : GRCh38
 genotype: GM23248
 replicate: 1
 sequencing: paired
 roadmap_id : E126
 roadmap_hmm : https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E126_15_coreMarks_hg38lift_mnemonics.bed.gz

ENCFF119ELB :
 exp : https://www.encodeproject.org/experiments/ENCSR440MLE/
 bam : https://www.encodeproject.org/files/ENCFF119ELB/@@download/ENCFF119ELB.bam
 assembly : GRCh38
 genotype: GM23248
 replicate: 2
 sequencing: paired
 roadmap_id : E126
 roadmap_hmm : https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E126_15_coreMarks_hg38lift_mnemonics.bed.gz
"""

config = yaml.load(config)
# print(config)


## yaml tests end ###########

# Rbin = "/home/imallona/soft/R/R-3.5.1/bin/R"
Rbin = "/usr/local/bin/R"

# CHROMS = ['chr' + str(c) for c in range(17,19)]
CHROMS = ['chr19']
MINCOVERAGE = 20
NTHREADS = 4

SAMPLES = glob.glob(os.path.join('test_data', '*.bam'))

BAMS = {os.path.splitext(os.path.basename(sample))[0] : sample for sample in SAMPLES}
# HMMS = [os.path.splitext(os.path.basename(hmm))[0] for hmm in glob.glob(os.path.join('hmm', '*.gz'))]

# HMMS = {os.path.splitext(os.path.basename(hmm))[0] : hmm for hmm in glob.glob(os.path.join('hmm', '*.gz'))}
HMMS = {os.path.basename(hmm).split(os.extsep,1)[0] : hmm for hmm in glob.glob(os.path.join('hmm', '*.gz'))}



METHTUPLE = expand("{sample}/qnsorted_{sample}_{chrom}.CG.2.tsv",
                    sample =  BAMS.keys(),
                    chrom = CHROMS)

AWK = expand("{sample}/qnsorted_{sample}_{chrom}.CG.2_cov_{cov}.tsv",
                    sample =  BAMS.keys(),
                    chrom = CHROMS,
                    cov = MINCOVERAGE)

METHYLATIONS =  {os.path.splitext(os.path.basename(sample))[0] : expand("{sample}/qnsorted_{sample}_{chrom}_cov_{cov}_methylation.bed",
                    sample =  sample,
                    chrom = CHROMS,
                    cov = MINCOVERAGE) for sample in BAMS}

ENTROPIES =  {os.path.splitext(os.path.basename(sample))[0] : expand("{sample}/qnsorted_{sample}_{chrom}_cov_{cov}_entropy.bed",
                    sample =  sample,
                    chrom = CHROMS,
                    cov = MINCOVERAGE) for sample in BAMS}

MERGED_ENTROPY = expand("{sample}/{sample}_cov_{cov}_entropy.bed.gz",
                        sample = BAMS.keys(),
                           cov = MINCOVERAGE)
MERGED_METHYLATION = expand( "{sample}/{sample}_cov_{cov}_methylation.bed.gz",
                             sample = BAMS.keys(),
                             cov = MINCOVERAGE)

# KNITTED_REPORT = expand("{sample}/01_data_visualization_{sample}_cov_{cov}_hmm_{hmm}.html",
#                         sample = BAMS.keys(),
#                         cov = MINCOVERAGE,
#                         hmm = HMMS.keys())

COLORED = expand( "{sample}/{sample}_cov_{cov}_entropy_and_meth_hmm_{hmm}.bed.gz",
                  sample = BAMS.keys(),
                  cov = MINCOVERAGE,
                  hmm = HMMS.keys())

R_OUTS = expand( "{sample}/{sample}_cov_{cov}_hmm_{hmm}.test",
                  sample = BAMS.keys(),
                  cov = MINCOVERAGE,
                  hmm = HMMS.keys())

# R_OUTS_VARIANT =  ','.join({os.path.splitext(os.path.basename(sample))[0] : expand( "{sample}/{sample}_cov_{cov}_hmm_{hmm}_to_integrate.tsv.gz",
#                   sample = BAMS.keys(),
#                   cov = MINCOVERAGE,
#                   hmm = HMMS.keys()) for sample in BAMS})

R_OUTS_VARIANT =  {os.path.splitext(os.path.basename(sample))[0] : expand( "{sample}/{sample}_cov_{cov}_hmm_{hmm}_to_integrate.tsv.gz",
                  sample = BAMS.keys(),
                  cov = MINCOVERAGE,
                  hmm = HMMS.keys()) for sample in BAMS}

INTEGRATIVE_OUTPUT = expand( "across_cov_{cov}_hmm.test",
                        cov = MINCOVERAGE)

# print(ENTROPIES)
# print(METHYLATIONS)
# print(R_OUTS_VARIANT)
# print(','.join(','.join(item) for item in R_OUTS_VARIANT.values()))

rule all:
    input:
        # BAMS, AWK, METHTUPLE, ENTROPIES, METHYLATIONS, MERGED_ENTROPY, MERGED_METHYLATION, COLORED#, R_OUTS
        AWK, METHTUPLE, MERGED_METHYLATION, MERGED_ENTROPY, COLORED#, INTEGRATIVE_OUTPUT



# onsuccess:
#     print("Workflow finished, no error")
#     params:
#         script = "scripts/multisample_report.R",
#         samples = ','.join(','.join(item) for item in R_OUTS_VARIANT.values())    
#     log:
#         "multi_sample_report.log"
#     output:
#         "across_cov_{cov}_hmm.test"
#     shell: '''{Rbin} CMD BATCH --no-restore --no-save \
#     "--args annotated='{params.samples}' \
#     output='{output}'" {params.script} {log}'''
    
# ## does not work Not all output, log and benchmark files of rule run_across_samples_report
# ##  contain the same wildcards.
# rule run_across_samples_report:
#     input:
#         # lambda wildcards: R_OUTS_VARIANT[wildcards.sample]
#         # ','.join(R_OUTS_VARIANT.values())
#         stuff =   lambda wildcards: R_OUTS_VARIANT[wildcards.sample],
#         script = "scripts/multisample_report.R"
#     params:
#         ','.join(','.join(item) for item in R_OUTS_VARIANT.values())    
#     log:
#         "multi_sample_report.log"
#     output:
#          "across_cov_{cov}_hmm.test"
#     shell: '''{Rbin} CMD BATCH --no-restore --no-save \
#     "--args annotated='{params}' \
#     output='{output}'" {input.script} {log}'''


rule run_sample_report:
    input:
        colored_entropy =  "{sample}/{sample}_cov_{cov}_entropy_and_meth_hmm_{hmm}.bed.gz",
        script = 'scripts/report.R'
    log:
         "{sample}/{sample}_cov_{cov}_hmm_{hmm}_report.log"
    output:
        "{sample}/{sample}_cov_{cov}_hmm_{hmm}.test"
    # script:
    #     "scripts/report.R"
    shell:
        '''{Rbin} CMD BATCH --no-restore --no-save \
        "--args colored_entropy='{input.colored_entropy}' \
        output='{output}'" {input.script} {log}'''

## for this proper hmm must be downloaded @todo check the shell step
rule annotate_entropies_to_hmm:
    input:
        # "{sample}/{sample}_cov_{MINCOVERAGE}_entropy"
        entropy =  "{sample}/{sample}_cov_{MINCOVERAGE}_entropy_and_meth.bed.gz",
        hmm = "hmm/{hmm}.bed.gz"
    output:
        "{sample}/{sample}_cov_{MINCOVERAGE}_entropy_and_meth_hmm_{hmm}.bed.gz"
    shell:
        """ bedtools intersect -wa -wb \
              -a {input.entropy} \
              -b {input.hmm} | \
        bedtools groupby -g 1-7 -c 11 -o concat | gzip > \
        {output}"""
        # """ touch  {output}"""


        
rule add_methylation_to_entropy_file:
    input:
        methylation = "{sample}/{sample}_cov_{MINCOVERAGE}_methylation.bed.gz",
        entropy = "{sample}/{sample}_cov_{MINCOVERAGE}_entropy.bed.gz"
    output:
        "{sample}/{sample}_cov_{MINCOVERAGE}_entropy_and_meth.bed.gz"
    shell:
        """
        # merge
        command=paste
        for i in {input.entropy} {input.methylation}
        do
          command="$command <(gzip -cd $i)"
        done
  
        eval $command | cut -f1-6,11 | dos2unix | gzip > {output}
        """

rule merge_methylations:
    input:
        lambda wildcards: METHYLATIONS[wildcards.sample]
    output:
        temp("{sample}/{sample}_cov_{cov}_methylation.bed.gz")
    shell:
        "cat {input} | bedtools sort | gzip > {output}"

rule merge_entropies:
    priority: 70
    input:
         lambda wildcards: ENTROPIES[wildcards.sample]
    output:
        temp("{sample}/{sample}_cov_{MINCOVERAGE}_entropy.bed.gz")
    shell:
        "cat {input} | bedtools sort | gzip > {output}"
        
rule get_methylation:
    input:
        "{sample}/qnsorted_{sample}_{chrom}.CG.2_cov_{MINCOVERAGE}.tsv"
    output:
        temp("{sample}/qnsorted_{sample}_{chrom}_cov_{MINCOVERAGE}_methylation.bed")
    shell:
        """
        awk '
BEGIN {{OFS=FS="\t"}}
NR == 1 {{next}}
{{
  print $1,$3,$4,"MM"$5";MU"$6";UM"$7";UU"$8,(($5+(0.5*$6)+(0.5*$7))/($5+$6+$7+$8)),$2
}}' {input} > {output}
"""
   
rule get_entropy:
    input:
        "{sample}/qnsorted_{sample}_{chrom}.CG.2_cov_{MINCOVERAGE}.tsv"
    output:
        temp("{sample}/qnsorted_{sample}_{chrom}_cov_{MINCOVERAGE}_entropy.bed")
    shell:
        """
        awk '
BEGIN {{OFS=FS="\t"}}
NR == 1 {{next}}
{{
  H["MM"] = $5
  H["UU"] = $8
  H["UM"] = $7
  H["MU"] = $6
  N  = ($5+$6+$7+$8)
  E = 0
  p = ""
  for (i in H) {{
    if (H[i] != 0) {{
      p = H[i]/N;
      E -=  p * log(p);
    }}
  }}

  print $1,$3,$4,"MM"$5";MU"$6";UM"$7";UU"$8,E,$2

}}' {input} > {output}
"""

rule filter_by_coverage:
       input:
           "{sample}/qnsorted_{sample}_{chrom}.CG.2.tsv"
       output:
           temp("{sample}/qnsorted_{sample}_{chrom}.CG.2_cov_{MINCOVERAGE}.tsv")
       shell:
           """
           awk -v mincov="{MINCOVERAGE}" '
           NR == 1 {{next}}
           {{
            if (($5+$6+$7+$8) >= mincov) print $0
           }}' {input} > {output}
           """
                   
rule run_meththuple:
    priority: 80
    input:
        "{sample}/qnsorted_{sample}_{chrom}.bam"
    output:
        "{sample}/qnsorted_{sample}_{chrom}.CG.2.tsv"
    shell: """
        set +u;
        source ~/virtenvs/methtuple/bin/activate
        methtuple -m 2 --methylation-type CG {input}    
        deactivate
        set -u;
        """

rule sort_by_queryname:
    priority: 90
    input:
        "{sample}/{sample}_{chrom}.bam"
    output:
        temp("{sample}/qnsorted_{sample}_{chrom}.bam")
    threads: NTHREADS
    shell:
        """
        samtools sort -@ {threads} -n {input} -o {output}
        # mv {output}.bam {output}"""

rule split_by_chrom:
    input:
        "{sample}/{sample}.bam"
        # lambda wildcards, attempt: attempt * 100
    output:
        temp("{sample}/{sample}_{chrom}.bam")
    threads: NTHREADS
    shell:
        "samtools view -@ {threads} -b {input} {wildcards.chrom} > {output}"
 
rule samtools_sort:
    priority: 100
              
    input:
        "test_data/{sample}.bam"
    output:
        temp("{sample}/{sample}.bam")
    threads : NTHREADS
    shell:
        """
        samtools sort -@ {threads} {input} -o {output}
        # mv {output}.bam {output}
        samtools index {output}
        """

# rule download_hmm:
#     # input:
#     #     [config[sample]['roadmap_hmm'] for sample in config.keys()]
#     output:
#         "{encode}"
#     shell:
#         "wget --quiet {output} -O hmm/{output}"
        

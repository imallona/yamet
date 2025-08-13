"""
Handles CRC data

Single-cell multiomics sequencing and analyses of human colorectal cancer

Shuhui Bian https://orcid.org/0000-0002-9662-113X, Yu Hou https://orcid.org/0000-0001-7875-7087, Xin Zhou https://orcid.org/0000-0002-4048-4017, Xianlong Li https://orcid.org/0000-0001-7745-2695, Jun Yong https://orcid.org/0000-0002-3770-2108, Yicheng Wang https://orcid.org/0000-0002-3901-1482, Wendong Wang https://orcid.org/0000-0002-8631-3109, Jia Yan https://orcid.org/0000-0002-6203-0016, Boqiang Hu https://orcid.org/0000-0002-2045-3619, Hongshan Guo https://orcid.org/0000-0003-2799-2989, Jilian Wang https://orcid.org/0000-0002-3951-6931, Shuai Gao https://orcid.org/0000-0002-8208-9197, Yunuo Mao https://orcid.org/0000-0003-1428-270X, Ji Dong, Ping Zhu, Dianrong Xiu, Liying Yan https://orcid.org/0000-0001-9572-9440, Lu Wen, Jie Qiao https://orcid.org/0000-0003-2126-1376 , Fuchou Tang https://orcid.org/0000-0002-8625-7717 , and Wei Fu https://orcid.org/0000-0001-5248-7891

https://www.science.org/doi/10.1126/science.aao3791?url_ver=Z39.88-2003

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693

This is hg19

Izaskun Mallona
GPLv3
"""

import os.path as op
from glob import glob

CRC = op.join("data", "crc")
CRC_RAW = op.join(CRC, "raw")                       ## raw, as in GEO and from GEO
CRC_HARMONIZED = op.join(CRC, "harmonized")         ## ingestable by yamet
CRC_OUTPUT = op.join(CRC, "output")                 ## output for features (genes, promoters etc)
CRC_WINDOWS_OUTPUT = op.join(CRC, 'windows_output') ## output for tiles/genomic windows

## bedfiles
ANNOTATIONS = {
    "pmd": ["pmds", "hmds"],
    "hmm": [
        "0_Enhancer",
        "2_Enhancer",
        "11_Promoter",
        "12_Promoter",
        "1_Transcribed",
        "4_Transcribed",
        "5_RegPermissive",
        "7_RegPermissive",
        "6_LowConfidence",
        "3_Quiescent",
        "8_Quiescent",
        "10_Quiescent",
        "9_ConstitutiveHet",
        "13_ConstitutiveHet",
    ],
    "chip": ["H3K27me3", "H3K9me3", "H3K4me3"],
    "lad": ["laminb1"],
}

## patient -> biopsy sites mapping
SAMPLES = {"CRC01" : ['NC', 'PT', 'LN', 'ML', 'MP'],
           "CRC02" : ['NC', 'PT', 'ML', 'PT'],
           "CRC04" : ['NC', 'PT'],
           "CRC10" : ['NC', 'PT', 'LN'],
           "CRC11" : ['NC', 'PT', 'LN'],
           "CRC13" : ['NC', 'PT', 'LN'],
           "CRC15" : ['NC', 'PT', 'LN', 'ML', 'MO']} # What is MO, a typo?


def list_raw_files_from_metadata():
    input_list_path = op.join("src", "sctrioseq_bismark_files.txt")
    with open(input_list_path, "r") as f:
        files = [line.strip() for line in f if line.strip()]
    return [op.join(CRC_RAW, file) for file in files]


rule download_crc_accessors:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        gsm=op.join(CRC, "bisulfites_gsm.txt"),
    message:
        """
        CRC Accessors download
        """
    params:
        raw=CRC_RAW
    shell:
        """
        mkdir -p {params.raw}
        
        esearch -db sra -query PRJNA382695 |
          efetch -format runinfo |
          cut -f11,13,14,15,30 -d"," |
          grep Bis |
          cut -f5 -d"," > {output.gsm}
        """


rule download_crc:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        gsm=op.join(CRC, "bisulfites_gsm.txt")
    output:
        download_flag = protected(op.join(CRC, "download.flag")),
        files = protected(list_raw_files_from_metadata())
    log:
        op.join("logs", "crc_download.log")
    params:
        raw=CRC_RAW
    message:
        """
        CRC download
        """
    shell:
        """
        ## this is extremely dirty, better control to avoid re-download needed
        ##  not robust neither, we should keep a list of the files being included into the
        ##  analysis and not glob them afterwards
        echo {output.download_flag}
        cat {output.download_flag}
        
        if [ ! -f {output.download_flag} ]
        then
          echo "downloading data "
        
          while read -r gsm
          do
            short="$(echo $gsm | cut -c1-7)"
            url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"$short"nnn/"$gsm"/suppl/
            if [ ! -e "{params.raw}"/$(basename "$url") ]
            then
              # removed a -r because it was making -nc not work
              # but might be -r is necessary, not deleting files to debug
              # so if not working add `-r`
              wget --no-directories --directory-prefix={params.raw} \
                --no-clobber --execute robots=off -k -A gz $url
              fi
          done < {input.gsm} > {log}
        
        date >> {output.download_flag}
        fi
        """

# ruleorder: harmonize_cell_report_for_yamet > run_yamet_on_separate_features > run_yamet_on_windows


checkpoint harmonize_cell_report_for_yamet:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        list_raw_files_from_metadata(),
        op.join(CRC, "download.flag")
    output:
        outdir = directory(op.join(CRC_HARMONIZED, "{patient}", "{location}"))
    params:
        raw = CRC_RAW
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.outdir}
        shopt -s nullglob
        for f in {params.raw}/G*_{wildcards.patient}_{wildcards.location}*.txt.gz; do
          base=$(basename "$f")
          gunzip -c "$f" |
            grep "CpG$" |
            awk '{{OFS="\t"; if ($4 == "-") $2 = $2 - 2; else if ($4 == "+") $2 = $2 - 1; print $1, $2, $2+1, "", "", $4, $6, $5;}}' |
            bedtools merge -c 7,8 -o sum |
            awk '{{OFS=FS="\t"; bin = ($4/$5 > 0.1) ? 1 : 0; print $1,$2,$4,$5,bin}}' |
            sort -k1,1 -k2,2n |
            gzip -c > {output.outdir}/"$base"
        done
        """

def get_harmonized_files(wildcards):
    ckpt = checkpoints.harmonize_cell_report_for_yamet.get(
        patient=wildcards.patient, location=wildcards.location
    )
    outdir = ckpt.output.outdir
    return sorted(glob(op.join(outdir, "*.txt.gz")))

        
# print(get_harmonized_files(patient = 'CRC01', location = "NC"))

rule run_yamet_on_separate_features:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        # cells=lambda wildcards: expand(
        #     "{file}",
        #     file=get_harmonized_files(wildcards.patient, wildcards.location),
        # ),
        cells = get_harmonized_files,
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        bed=op.join(HG19_BASE, "{subcat}.{cat}.bed")
    output:
        simple=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}.{location}.out"),
        det=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.det.out"),
        meth=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.meth.out")
    log:
        op.join('logs', 'yamet_{subcat}_{cat}_{patient}_{location}.log')
    group:
        "yamet"
    params:
        path = CRC_OUTPUT
    threads: 16
    shell:
        """
        mkdir -p {params.path}
        yamet \
         --cell {input.cells} \
         --reference {input.ref} \
         --intervals {input.bed} \
         --cores {threads} \
         --print-sampens F \
         --out {output.simple} \
         --det-out {output.det} \
         --meth-out {output.meth} &> {log}
        """

def list_relevant_yamet_outputs():
    res = []
    for annot in ANNOTATIONS:
        for subannot in ANNOTATIONS[annot]:
            for patient in SAMPLES.keys():
                for sampling in SAMPLES[patient]:
                    res.append(f"{subannot}_{annot}_{patient}_{sampling}.det.out")
    return [op.join(CRC_OUTPUT, item) for item in res]

rule render_crc_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        list_relevant_yamet_outputs()
    params:
        output_path=CRC_OUTPUT
    output:
        op.join(CRC, "results", "crc.html")
    log:
        log = op.join("logs", "render_crc.log")
    script:
        "src/crc.Rmd"

        
#####################################################################################################
        
## windows-based stuff ##############################################################################

#####################################################################################################
 
""" the idea is to get a common set of features, genomic tiles, to run
    some data mining techniques afterwards
    largey las in prototyping/crc_prototype.sh
"""

rule get_sizes_hg19:
    conda:
        op.join("..", "envs", "yamet.yml")
    output:
        sizes = op.join(HG19_BASE, "hg19.sizes"),
    shell:
        """
        ## the `hg19.smk` code is so convoluted we need to download it again, cannot reuse
        mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -s -e \
          'SELECT chrom,size FROM hg19.chromInfo' > {output.sizes}       
        """

rule make_windows_hg19:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        sizes = op.join(HG19_BASE, "hg19.sizes"),
    output:
        windows =  op.join(HG19_BASE, "windows_{win_size,\d+}_nt.bed")
    shell:
        """
        bedtools makewindows -g {input.sizes} \
           -w {wildcards.win_size} | sort -k1,1 -k2,2n -k3,3n > {output.windows}
        """

rule get_single_annotion_coverage_per_window:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        windows =  op.join(HG19_BASE, "windows_{win_size}_nt.bed"),
        annotation =  op.join(HG19_BASE, "{subcat}.{cat}.bed")
    output:
        header = temp(op.join(HG19_BASE, "{win_size}_nt_{subcat}.{cat}.header")),
        body = temp(op.join(HG19_BASE, "{win_size}_nt_{subcat}.{cat}.body")),
        annotated_windows =  op.join(HG19_BASE, "windows_{win_size}_nt_{subcat}_{cat}_annotation.frac")
    shell:
        """
        bedtools coverage -a {input.windows} \
            -b {input.annotation} | cut -f7 > {output.body}
        echo "{wildcards.subcat}_{wildcards.cat}" > {output.header}
        cat {output.header} {output.body} > {output.annotated_windows}
        """

def list_annotated_windows():
    res = []
    for cat in ANNOTATIONS:
        for subcat in ANNOTATIONS[cat]:
            res.append(f"windows_{{win_size}}_nt_{subcat}_{cat}_annotation.frac")
    return [op.join("hg19", item) for item in res]

print(list_annotated_windows())

rule combine_annotated_windows:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        annotated_windows =  list_annotated_windows()
    output:
        op.join(HG19_BASE, "windows_{win_size}_nt_annotation.gz")
    shell:
        """
        paste {input.annotated_windows} | gzip -c> {output}
        """

"""
issue is, we need to permute backgrounds... and/or correct analytically
favouring the analytical aproach now, skip permuting
"""

rule run_yamet_on_windows:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        # cells=lambda wildcards: expand(
        #     "{file}",
        #     file=get_harmonized_files(wildcards.patient, wildcards.location),
        # ),
        cells = get_harmonized_files,
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        windows=op.join(HG19_BASE, "windows_{win_size}_nt.bed")
    output:
        simple=op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.out"),
        det=op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.det.out"),
        meth=op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.meth.out")
    log:
        op.join('logs', 'yamet_{win_size}_{patient}_{location}.log')
    group:
        "yamet"
    params:
        path = CRC_WINDOWS_OUTPUT
    threads: 16
    shell:
        """
        mkdir -p {params.path}
        yamet \
         --cell {input.cells} \
         --reference {input.ref} \
         --intervals {input.windows} \
         --cores {threads} \
         --print-sampens F \
         --out {output.simple} \
         --det-out {output.det} \
         --meth-out {output.meth} &> {log}
        """

def list_relevant_yamet_windows_outputs():
    res = []
    for patient in SAMPLES.keys():
        for location in SAMPLES[patient]:
            res.append(f"{{win_size}}_{patient}_{location}.det.out")
    return [op.join(CRC_WINDOWS_OUTPUT, item) for item in res]

rule render_crc_windows_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        yamet_dets = list_relevant_yamet_windows_outputs(),
        annotations = op.join(HG19_BASE, "windows_{win_size,\d+}_nt_annotation.gz")
    params:
        output_path=CRC_WINDOWS_OUTPUT
    output:
        op.join(CRC, "results", "crc_windows_{win_size}_nt.html")
    log:
        log = op.join("logs", "render_crc_windows_{win_size}.log")
    script:
        "src/crc_windows.Rmd"
    # shell:
    #     """
    #     echo 'nothing done' > {output}
    #     """

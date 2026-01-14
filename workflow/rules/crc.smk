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
    "genes": ["genes"],
    "cpgIslandExt": ['cpgIslandExt'],
    "scna": ["crc01_nc_scna", "crc01_gain_scna", "crc01_lost_scna"]
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
    """
    Reports a hardcoded list of bismark meth reports, to make sure download is complete
    """
    
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

checkpoint download_crc_bismarks:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        gsm=op.join(CRC, "bisulfites_gsm.txt")
    output:
        urls = temp(op.join(CRC_RAW, 'tmp.urls')),
        download_flag = op.join(CRC, "download.flag")#,
        # files = list_raw_files_from_metadata() ## commented out so snmk does not remove all these files before re-running the rule, so the finding/skipping algo makes sense
    log:
        op.join("logs", "crc_download.log")
    params:
        raw=CRC_RAW
    threads:
        1
    message:
        """
        CRC download
        """
    shell: ## this is stupid, because it wipes out the folder before exec
       """
        exec &> {log}
        
        set -euo pipefail
        mkdir -p "{params.raw}"
        
        while read -r gsm
          do
            echo $gsm
            short="$(echo $gsm | cut -c1-7)"
            url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"$short"nnn/"$gsm"/suppl/
        
            # wget -q -e robots=off -r -k -A gz -nd -P {params.raw} $url
            matches=$(find "{params.raw}" \
                   -maxdepth 1 \
                   -type f \
                   -name "${{gsm}}*" \
                   -print)
            echo "$matches"
 
            if [[ -n "$matches" ]]; then
               echo "Found existing report for $gsm:"
               printf "  %s Skipping download\\n\\n" "${{matches[@]}}"
            else
               printf "Downloading $gsm from $url\\n\\n"

               wget \
                  --quiet \
                  --execute=robots=off \
                  --recursive \
                  --convert-links \
                  --accept=gz \
                  --no-directories \
                  --directory-prefix={params.raw} \
                  --timestamping \
                  "$url"
            
                echo "$url" >> {output.urls}
           fi
          
        done < {input.gsm}

        touch {output.urls} # in case files were already there
        touch {output.download_flag}
        """


rule harmonize_cell_report_for_yamet:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(CRC, "download.flag"),
        # list_raw_files_from_metadata(),
        op.join(CRC_RAW, "{file}")        
        # op.join(CRC_RAW, "{gsm}_{patient}_{location}_{cellid}.singleC.txt.gz")
    output:
        op.join(CRC_HARMONIZED, "{file}") ## these should tmp
    params:
        raw=CRC_RAW,
        harmonized=CRC_HARMONIZED
    threads: max(8, workflow.cores/8)
    log:
        op.join("logs", "{file}_harmonization.log")
    shell:
        """
        exec &> {log}

        zcat {params.raw}/{wildcards.file} |
          grep "CpG$" |
        awk '
            {{OFS="\\t";}} {{
                if ($4 == "-") $2 = $2 - 2;
                else if ($4 == "+") $2 = $2 - 1;
                print $1, $2, $2+1, "", "", $4, $6, $5;
            }}
        ' |
        bedtools merge -c 7,8 -o sum |
        awk '
             {{OFS=FS="\\t";
                bin = ($4/$5 > 0.1) ? 1 : 0;
                print $1,$2,$4,$5,bin}}
        ' |
        sort -k1,1 -k2,2n |
        gzip -c > {output}
        """

## location is whether is a NC, PT, LN and so on
## patient is whether CRC01, CRC01 and so on
def get_harmonized_files(patient, location):
    ckpt = checkpoints.download_crc_bismarks.get()
    
    filepaths = glob(op.join(CRC_RAW, f"G*_{patient}_{location}*.txt.gz"))
    if len(filepaths) == 0:
        raise Exception('No cytosine reports matching the specs.')
    return [op.join(CRC_HARMONIZED, op.basename(file)) for file in filepaths]



rule run_yamet_on_separate_features:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        # cells = get_harmonized_files,
        cells=lambda wildcards: expand(
            "{file}",
            file=get_harmonized_files(wildcards.patient, wildcards.location),
        ),
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        bed=op.join(HG19_BASE, "{subcat}.{cat}.bed")
    output:
        simple_uncomp = temp(op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.out")),
        det_uncomp = temp(op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.det.out")),
        meth_uncomp = temp(op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.meth.out")),
        simple=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.out.gz"),
        det=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.det.out.gz"),
        meth=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.meth.out.gz")
    log:
        op.join('logs', 'yamet_{subcat}_{cat}_{patient}_{location}.log')
    group:
        "yamet"
    params:
        path = CRC_OUTPUT
    threads: max(8, workflow.cores/8)
    shell:
        """
        mkdir -p {params.path}
        yamet \
         --cell {input.cells} \
         --reference {input.ref} \
         --intervals {input.bed} \
         --cores {threads} \
         --print-sampens F \
         --out {output.simple_uncomp} \
         --det-out {output.det_uncomp} \
         --meth-out {output.meth_uncomp} &> {log}
        
        gzip --keep -f {params.path}/{wildcards.subcat}_{wildcards.cat}_{wildcards.patient}_{wildcards.location}*out  &>> {log}
        """

def list_relevant_yamet_outputs():
    res = []
    for annot in ANNOTATIONS:
        for subannot in ANNOTATIONS[annot]:
            for patient in SAMPLES.keys():
                for sampling in SAMPLES[patient]:
                    res.append(f"{subannot}_{annot}_{patient}_{sampling}.det.out.gz")
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
        windows =  op.join(HG19_BASE, r"windows_{win_size,\d+}_nt.bed")
    shell:
        """
        bedtools makewindows -g {input.sizes} \
           -w {wildcards.win_size} | sort -k1,1 -k2,2n -k3,3n > {output.windows}
        """

## @todo make sure this ingests genes, CpGis and SCNAs!
rule get_aggregated_annotation_coverage_per_window:
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
        annotated_windows =  list_annotated_windows(),
        windows =  op.join(HG19_BASE, r"windows_{win_size,\d+}_nt.bed")
    output:
        annotation = op.join(HG19_BASE, "windows_{win_size}_nt_annotation.gz"),
        header = temp(op.join(HG19_BASE, 'header_{win_size}_chr.txt')),
        windows_with_header = temp(op.join(HG19_BASE, 'window_{win_size}_with_header.txt'))
    shell:
        """
        echo -e "chr\\tstart\\tend" > {output.header}
        cat {output.header} {input.windows} > {output.windows_with_header}
        paste {output.windows_with_header} {input.annotated_windows} | gzip -c > {output.annotation}
        """

"""
issue is, we need to permute backgrounds... and/or correct analytically
favouring the analytical aproach now, skip permuting
"""

rule run_yamet_on_windows:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        cells=lambda wildcards: expand(
            "{file}",
            file=get_harmonized_files(wildcards.patient, wildcards.location),
        ),
        # cells = get_harmonized_files,
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        windows=op.join(HG19_BASE, "windows_{win_size}_nt.bed")
    output:
        simple_uncomp=temp(op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.out")),
        det_uncomp=temp(op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.det.out")),
        meth_uncomp=temp(op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.meth.out")),
        simple=op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.out.gz"),
        det=op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.det.out.gz"),
        meth=op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.meth.out.gz")
    log:
        op.join('logs', 'yamet_{win_size}_{patient}_{location}.log')
    group:
        "yamet"
    params:
        path = CRC_WINDOWS_OUTPUT
    threads: max(8, workflow.cores/8)
    shell:
        """
        mkdir -p {params.path}
        yamet \
         --cell {input.cells} \
         --reference {input.ref} \
         --intervals {input.windows} \
         --cores {threads} \
         --print-sampens F \
         --out {output.simple_uncomp} \
         --det-out {output.det_uncomp} \
         --meth-out {output.meth_uncomp} &> {log}

        gzip --keep -f {params.path}/{wildcards.win_size}_{wildcards.patient}_{wildcards.location}*out  &>> {log}
        """

def list_relevant_yamet_windows_outputs():
    res = []
    for patient in SAMPLES.keys():
        for location in SAMPLES[patient]:
            res.append(f"{{win_size}}_{patient}_{location}.det.out.gz")
    return [op.join(CRC_WINDOWS_OUTPUT, item) for item in res]


## other reports from Atreya to be re-categorized / placed somewhere - caution not sure where the SCNA bedfile is used during the annotation phase @todo

rule get_scna_patient1_from_supplementary_data:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        op.join(".", "src", "scna_hg19.xlsx"),
    output:
        scna_bed = op.join(HG19_BASE, "patient_crc01_scna.scna.bed.gz.pre")
    script:
        "src/parse_scna.R"

rule split_patient1_crc_in_kept_lost_gained:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        scna_bed = op.join(HG19_BASE, "patient_crc01_scna.scna.bed.gz.pre"),
        genome_sizes = op.join(HG19_BASE, "genome.sizes"),
    output:
        uncomp = temp(op.join(HG19_BASE, "crc01_scna.bed")),
        gained = op.join(HG19_BASE, "crc01_gain_scna.scna.bed"),
        lost = op.join(HG19_BASE, "crc01_lost_scna.scna.bed"),
        kept = op.join(HG19_BASE, "crc01_nc_scna.scna.bed"),        
    shell:
        """
           gzip -d -c {input.scna_bed} > {output.uncomp}
           grep "deleted" {output.uncomp} | bedtools sort > {output.lost}
           grep "amplified" {output.uncomp} | bedtools sort > {output.gained}
           bedtools complement -i {output.uncomp} -g {input.genome_sizes} > {output.kept}
        """
        
rule render_crc_windows_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        yamet_dets = list_relevant_yamet_windows_outputs(),
        annotations = op.join(HG19_BASE, r"windows_{win_size,\d+}_nt_annotation.gz")
    params:
        output_path=CRC_WINDOWS_OUTPUT
    threads: workflow.cores        
    output:
        op.join(CRC, "results", "crc_windows_{win_size}_nt.html")
    log:
        log = op.join("logs", "render_crc_windows_{win_size}.log")
    script:
        "src/crc_windows.Rmd"

## from crc_windows.Rmd, to crc_windows_sce.Rmd
rule lazy_move_rds_objects:
    conda:
        op.join("..", "envs", "r.yml")
    params:
        output_path=CRC_WINDOWS_OUTPUT
    input:
        op.join(CRC, "results", "crc_windows_{win_size}_nt.html")
    output:
        sce = op.join(CRC, 'results', 'sce_windows_{win_size}_colon.rds'),
        de =  op.join(CRC, 'results', 'de_list_{win_size}.rds')
    threads: 1
    shell:
        """
        cp sce_windows_colon.rds {output.sce}
        cp de_list.rds {output.de}
        """
        
## params$corrected_sce is an output, but Rmd scripts only allow one
rule render_crc_sce_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        sce = op.join(CRC, 'results', 'sce_windows_{win_size}_colon.rds'),
        de =  op.join(CRC, 'results', 'de_list_{win_size}.rds'),
        windows_annotation = op.join(HG19_BASE, "windows_{win_size}_nt_annotation.gz")
    threads: workflow.cores
    params:
        output_path=CRC_WINDOWS_OUTPUT,
        corrected_sce = op.join(CRC, 'results', 'sce_windows_colon_corrected_{win_size}.rds')
    output:
        op.join(CRC, "results", "crc_windows_sce_{win_size}.html")
    log:
        log = op.join("logs", "render_crc_windows_sce_{win_size}.log")
    script:
        "src/crc_windows_sce.Rmd"

# rule run_crc_stats_report:
#     conda:
#         op.join("..", "envs", "r.yml")
#     input:
#         list_relevant_yamet_windows_outputs(),
#         annotation=op.join(HG19_BASE, "windows_{win_size}_nt_annotation.gz")
#     output:
#         op.join("results", "crc_stats_{win_size}.html"),
#     params:
#         output_path=CRC_WINDOWS_OUTPUT,
#     log:
#         log=op.join("logs", "render_crc_stats_{win_size}.log"),
#     threads: 16
#     script:
#         "src/crc_stats.Rmd"

# rule run_crc_deletions_report:
#     conda:
#         op.join("..", "envs", "r.yml")
#     input:
#         scna_bed =  op.join(HG19_BASE, "patient_crc01_scna.scna.bed.gz"),
#         yamet_dets=expand(
#             op.join(CRC_WINDOWS_OUTPUT, "{{win_size}}_CRC01_{location}.det.out.gz"),
#             location=SAMPLES["CRC01"],
#         ),
#         annotation=op.join(HG19_BASE, "windows_{win_size}_nt_annotation.gz"), ## but we want this on scnas themselves, not only on windows
#     output:
#         op.join("results", "crc_deletions_{win_size}.html"),
#     log:
#         log=op.join("logs", "render_crc_deletions_{win_size}.log"),
#     threads: 16
#     script:
#         "src/crc_deletions.Rmd"

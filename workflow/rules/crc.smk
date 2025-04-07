"""
Handles CRC data

Single-cell multiomics sequencing and analyses of human colorectal cancer

Shuhui Bian https://orcid.org/0000-0002-9662-113X, Yu Hou https://orcid.org/0000-0001-7875-7087, Xin Zhou https://orcid.org/0000-0002-4048-4017, Xianlong Li https://orcid.org/0000-0001-7745-2695, Jun Yong https://orcid.org/0000-0002-3770-2108, Yicheng Wang https://orcid.org/0000-0002-3901-1482, Wendong Wang https://orcid.org/0000-0002-8631-3109, Jia Yan https://orcid.org/0000-0002-6203-0016, Boqiang Hu https://orcid.org/0000-0002-2045-3619, Hongshan Guo https://orcid.org/0000-0003-2799-2989, Jilian Wang https://orcid.org/0000-0002-3951-6931, Shuai Gao https://orcid.org/0000-0002-8208-9197, Yunuo Mao https://orcid.org/0000-0003-1428-270X, Ji Dong, Ping Zhu, Dianrong Xiu, Liying Yan https://orcid.org/0000-0001-9572-9440, Lu Wen, Jie Qiao https://orcid.org/0000-0003-2126-1376 , Fuchou Tang https://orcid.org/0000-0002-8625-7717 , and Wei Fu https://orcid.org/0000-0001-5248-7891

https://www.science.org/doi/10.1126/science.aao3791?url_ver=Z39.88-2003

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693

This is hg19

"""

from glob import glob

CRC = op.join("data", "crc")
CRC_RAW = op.join(CRC, "raw")          ## from geo
CRC_HARMONIZED = op.join(CRC, "raw")   ## ingestable by yamet
CRC_OUTPUT = op.join(CRC, "output")    ## output

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
        download_flag = op.join(CRC, "download.flag")
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
        while read -r gsm
        do
          short="$(echo $gsm | cut -c1-7)"
          url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"$short"nnn/"$gsm"/suppl/
          wget --no-directories --directory-prefix={params.raw} \
            --no-clobber --execute robots=off -r -k -A gz $url
        done < {input.gsm} > {log}
        
        date >> {output.download_flag}
        """


rule parse_single_crc:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(CRC, "download.flag")
    output:
        op.join(CRC, "{file}")
    params:
        raw=CRC_RAW
    shell:
        """
        gunzip -c {params.raw}/{wildcards.file} |
          grep "CpG$" |
        awk '
            {OFS="\t";} {
                if ($4 == "-") $2 = $2 - 2;
                else if ($4 == "+") $2 = $2 - 1;
                print $1, $2, $2+1, "", "", $4, $6, $5;
            }
        ' |
        bedtools merge -c 7,8 -o sum |
        awk '
             {OFS=FS="\t";
                bin = ($4/$5 > 0.1) ? 1 : 0;
                print $1,$2,$4,$5,bin}
        ' |
        sort -k1,1 -k2,2n |
        gzip -c > {output}
        """

def get_raw_files(patient, location):
    filepaths = glob(op.join(CRC_RAW, f"G*_{patient}_{location}*.txt.gz"))
    return [op.basename(file) for file in filepaths]


rule yamet_crc_cg:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        download_flag = op.join("crc", "download.flag"),
        cells=lambda wildcards: expand(
            op.join(CRC, "{file}"),
            file=get_raw_files(wildcards.patient, wildcards.location),
        ),
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        bed=op.join(HG19_BASE, "{subcat}.{cat}.bed")
    output:
        simple=op.join(CRC_OUTPUT, "{subcat}.{cat}.{patient}.{location}.out"),
        det=op.join(CRC_OUTPUT, "{subcat}.{cat}.{patient}.{location}.det.out")
    log:
        op.join('logs', 'yamet_{subcat}_{cat}_{patient}_{location}.log')
    group:
        "yamet"
    params:
        path = CRC_OUTPUT
    threads: 16
    shell:
        """
        mkdir -p {params.output_path}
        yamet \
         --cell {input.cells} \
         --reference {input.ref} \
         --intervals {input.bed} \
         --cores {threads} \
         --print-sampens F \
         --out {output.simple} \
         --det-out {output.det} &> {log}
        """


CAT_MAP = {
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
LOCATIONS = ["NC", "PT", "LN"]
PATIENTS = ["CRC01", "CRC02", "CRC04", "CRC10", "CRC11", "CRC13", "CRC15"]


rule crc_doc:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        expand(
            op.join(CRC_OUTPUT, "{jnt}.{patient}.{location}.out"),
            jnt=[f"{subcat}.{cat}" for cat in CAT_MAP for subcat in CAT_MAP[cat]],
            location=LOCATIONS,
            patient=PATIENTS,
        ),
    params:
        yamet=CRC_OUTPUT,
    output:
        op.join("crc", "crc.html"),
    script:
        "src/crc.Rmd"

#!/bin/bash

##
## Processes the CRC data quick and dirty
##
## Started 25th Feb 2025
##
## Izaskun Mallona
## GPLv3

DATAPATH=$HOME/src/yamet/crc_data
WD=$HOME/tmp/yamet_crc

normals=$(find $DATAPATH -name "*CRC01_NC*gz" | sort | paste -sd " ")
primaries=$(find $DATAPATH -name "*CRC01_PT*gz" |  sort| paste -sd " ")

# download pmds

mkdir -p $WD
cd $WD

wget https://zhouserver.research.chop.edu/GenomeAnnotation/hg19/PMD_coordinates_hg19.bed.gz

# download / installed yamet

# wget https://github.com/imallona/yamet/releases/download/v1.1.0-rc.2/yamet-Linux-x86_64-1.1.0-rc.2.deb
# dpkg-deb -e yamet-Linux-x86_64-1.1.0-rc.2.deb yamet

# well, compiling it, caution cmake 3.31.6, from the same branch

cd method
bash build.sh
export PATH=/home/imallona/src/yamet/method/build:$PATH
# generate custom reference, this time cpgs present in the longest normal file

zcat /home/imallona/src/yamet/crc_data/ftp.ncbi.nlm.nih.gov/geo/samples/GSM2697nnn/GSM2697701/suppl/GSM2697701_scTrioSeq2Met_CRC01_NC_529.singleC.txt.gz | grep CpG | awk '{OFS=FS="\t"; print $1,$2,$2+1}'| sort |  gzip -c > custom.ref.gz

# run yamet, with process substitutions

# (yamet) imallona@mlshainan:~/tmp/yamet_crc$ zcat  /home/imallona/src/yamet/crc_data/ftp.ncbi.nlm.nih.gov/geo/samples/GSM2697nnn/GSM2697695/suppl/GSM2697695_scTrioSeq2Met_CRC01_NC_507.singleC.txt.gz | head
# #Chr    Pos     Ref     Chain   Total   Met     UnMet   MetRate Ref_context     Type
# chr7    16411   c       +       1       1       0       1       CCT     CHH
# chr7    16412   c       +       1       1       0       1       CTA     CHH
# chr7    16416   c       +       1       1       0       1       CCC     CHH
# chr7    16417   c       +       1       1       0       1       CCT     CHH
# chr7    16418   c       +       1       1       0       1       CTA     CHH
# chr7    16422   c       +       1       1       0       1       CCC     CHH
# chr7    16423   c       +       1       1       0       1       CCT     CHH

# and files are expected to be     chromosome, position, number of methylated reads, total number of reads and the rate 



yamet --cell <(zcat /home/imallona/src/yamet/crc_data/ftp.ncbi.nlm.nih.gov/geo/samples/GSM2697nnn/GSM2697695/suppl/GSM2697695_scTrioSeq2Met_CRC01_NC_507.singleC.txt.gz | grep CpG | awk '{print $1,$2,$6,$5,$6/$5}' | sort) \
      --reference custom.ref.gz \
      --intervals PMD_coordinates_hg19.bed.gz \
      --out normals.yamet.out.gz \
      --cores 10 \
      --print-sampens F

# but breaks with
# Error: in line

zcat /home/imallona/src/yamet/crc_data/ftp.ncbi.nlm.nih.gov/geo/samples/GSM2697nnn/GSM2697695/suppl/GSM2697695_scTrioSeq2Met_CRC01_NC_507.singleC.txt.gz | grep CpG | head -10000 | awk '{OFS=FS="\t"; print $1,$2,$6,$5,$6/$5}' | gzip -c > test.cpg.gz

yamet --cell test.cpg.gz \
      --reference custom.ref.gz \
      --intervals PMD_coordinates_hg19.bed.gz \
      --out normals.yamet.out.gz \
      --cores 10 \
      --print-sampens F

# pfff, again

# (yamet) imallona@mlshainan:~/tmp/yamet_crc$ yamet --cell test.cpg.gz       --reference custom.ref.gz       --intervals PMD_coordinates_hg19.bed.gz       --out normals.yamet.out.gz       --cores 10       --print-sampens F
# Error: in line


# (yamet) imallona@mlshainan:~/tmp/yamet_crc$ zcat test.cpg.gz | head
# chr7    16486   1       1       1
# chr7    24406   1       1       1
# chr7    24480   1       1       1
# chr7    30601   2       2       1
# chr7    30626   2       2       1
# chr7    30638   2       2       1
# chr7    30702   2       2       1
# chr7    30801   1       1       1
# chr7    30810   1       1       1
# chr7    30812   1       1       1
# (yamet) imallona@mlshainan:~/tmp/yamet_crc$ zcat custom.ref.gz | head
# chr10   100000004       100000005
# chr10   100000012       100000013
# chr10   10000002        10000003
# chr10   10000011        10000012
# chr10   100000172       100000173
# chr10   100001868       100001869
# chr10   10000191        10000192
# chr10   100001923       100001924
# chr10   100002038       100002039
# chr10   10000205        10000206
# (yamet) imallona@mlshainan:~/tmp/yamet_crc$ zcat PMD_coordinates_hg19.bed.gz | head
# chr1    0       99999   0.131812279885387       PMD     neither
# chr1    100000  199999  0.158459023408353       PMD     commonPMD
# chr1    200000  299999  0.153218066385386       PMD     commonPMD
# chr1    300000  399999  0.150108771667091       PMD     commonPMD
# chr1    400000  499999  0.166654936685486       PMD     commonPMD
# chr1    500000  599999  0.154852098520618       PMD     commonPMD
# chr1    600000  699999  0.15312085101117        PMD     commonPMD
# chr1    700000  799999  0.0665392757807886      HMD     commonHMD
# chr1    800000  899999  0.127529527615858       PMD     neither
# chr1    900000  999999  0.131609015581668       PMD     neither

# shouldn't be because of the sorting, I guess

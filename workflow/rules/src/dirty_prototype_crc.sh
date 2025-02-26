#!/bin/bash

##
## Processes the CRC data quick and dirty
##

# 1. automate the process substitutions for all normals
# 2. run yamet per feature
# 3. sort features per avg methylation
# 4. repeat for primary tumors
# 5. compare whether ranks are inverted
# 6. evaluate whether there is a within-cell-type or across-cell-types change in meth variability

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

gunzip PMD*gz

# download / installed yamet

# wget https://github.com/imallona/yamet/releases/download/v1.1.0-rc.2/yamet-Linux-x86_64-1.1.0-rc.2.deb
# dpkg-deb -e yamet-Linux-x86_64-1.1.0-rc.2.deb yamet

# well, compiling it, caution cmake 3.31.6, from the same branch

cd method
bash build.sh
export PATH=/home/imallona/src/yamet/method/build:$PATH
# generate custom reference, this time cpgs present in the longest normal file

zcat /home/imallona/src/yamet/crc_data/ftp.ncbi.nlm.nih.gov/geo/samples/GSM2697nnn/GSM2697701/suppl/GSM2697701_scTrioSeq2Met_CRC01_NC_529.singleC.txt.gz | \
    grep CpG | awk '{OFS=FS="\t"; print $1,$2,$2+1}'| \
    sort -k1,1 -k2,2n -k3,3n  |  gzip -c > custom.ref.gz

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

zcat /home/imallona/src/yamet/crc_data/ftp.ncbi.nlm.nih.gov/geo/samples/GSM2697nnn/GSM2697695/suppl/GSM2697695_scTrioSeq2Met_CRC01_NC_507.singleC.txt.gz | \
    grep CpG | awk '{OFS=FS="\t"; print $1,$2,$6,$5,$6/$5}' | \
    sort -k1,1 -k2,2n -k3,3n |  gzip -c > test.cpg.gz

sort -k1,1 -k2,2n -k3,3n PMD_coordinates_hg19.bed > PMD_coordinates_hg19.bed.sorted

yamet --cell test.cpg.gz \
      --reference custom.ref.gz \
      --intervals PMD_coordinates_hg19.bed.sorted \
      --out normals.yamet.out \
      --cores 10 \
      --det-out normals.yamet.det.out \
      --print-sampens F 


# instead of using process substitutions, automate yamet-zation
# these should be temporary files within the snmk

mkdir -p nc pt

for normal in $normals
do
    bn=$(basename $normal .txt.gz)
    zcat $normal |\
      grep CpG | awk '{OFS=FS="\t"; print $1,$2,$6,$5,$6/$5}' | \
      sort -k1,1 -k2,2n -k3,3n |  gzip -c > nc/"$bn".yametized.gz
done


yamet --cell nc/*NC*yametized.gz \
      --reference custom.ref.gz \
      --intervals PMD_coordinates_hg19.bed.sorted \
      --out normals.yamet.out \
      --cores 30 \
      --det-out normals.yamet.det.out \
      --print-sampens F 
      
for primary in $primaries
do
    bn=$(basename $primary .txt.gz)
    zcat $primary |\
      grep CpG | awk '{OFS=FS="\t"; print $1,$2,$6,$5,$6/$5}' | \
      sort -k1,1 -k2,2n -k3,3n |  gzip -c > pt/"$bn".yametized.gz
done

yamet --cell pt/*PT*yametized.gz
      --reference custom.ref.gz \
      --intervals PMD_coordinates_hg19.bed.sorted \
      --out primaries.yamet.out \
      --cores 30 \
      --det-out primaries.yamet.det.out \
      --print-sampens F 

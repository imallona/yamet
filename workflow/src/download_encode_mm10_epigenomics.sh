#!/bin/bash
##
## mm10 neuron-related, beds, marks
##
## https://www.encodeproject.org/search/?type=Experiment&control_type!=*&related_series.@type=ReferenceEpigenome&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&status=released&biosample_ontology.classification=tissue&biosample_ontology.term_name=forebrain&assembly=mm10&assay_title=Histone+ChIP-seq&replicates.library.biosample.life_stage=postnatal&files.file_type=bed+narrowPeak&target.label=H3K27ac&target.label=H3K4me3&target.label=H3K9me3&target.label=H3K27me3

# forebrain, postnatal

# "https://www.encodeproject.org/metadata/?control_type%21=%2A&related_series.%40type=ReferenceEpigenome&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&status=released&biosample_ontology.classification=tissue&biosample_ontology.term_name=forebrain&assembly=mm10&assay_title=Histone+ChIP-seq&replicates.library.biosample.life_stage=postnatal&files.file_type=bed+narrowPeak&target.label=H3K27ac&target.label=H3K4me3&target.label=H3K9me3&target.label=H3K27me3&target.label=H3K4me1&type=Experiment&files.analyses.status=released&files.preferred_default=true"

mkdir -p annotation/mm10/encode/
cd $_

#H3K4me3
curl https://www.encodeproject.org/files/ENCFF160SCR/@@download/ENCFF160SCR.bed.gz -o h3k4me3.bed.gz

#H3K9me3
curl https://www.encodeproject.org/files/ENCFF658QTP/@@download/ENCFF658QTP.bed.gz -o h3k9me3.bed.gz

#H3K4me1
curl https://www.encodeproject.org/files/ENCFF937JHP/@@download/ENCFF937JHP.bed.gz -o h3k4me1.bed.gz

#H3K27me3
curl https://www.encodeproject.org/files/ENCFF827BBC/@@download/ENCFF827BBC.bed.gz -o h3k27me3.bed.gz

#H3K27ac 
curl https://www.encodeproject.org/files/ENCFF442GIT/@@download/ENCFF442GIT.bed.gz -o h3k27ac.bed.gz

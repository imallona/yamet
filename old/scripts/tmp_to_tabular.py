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

for key in config.keys():
    curr = config[key].values().reverse()
    print(key + '\t' + '\t'.join(map(str, curr)))

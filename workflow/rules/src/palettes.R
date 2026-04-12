## Shared color palettes and symbol mappings for all yamet reports.
##
## Source this file in every Rmd alongside ggtheme.R. All palettes are named
## character vectors of hex colors keyed by the display-level name. Using fixed
## hex codes (rather than generating them at runtime with turbo/viridis) ensures
## that colors stay stable even when a subset of levels is plotted.
##
## Palette families:
##   - Cell/biological groups use qualitative palettes (hand-picked or
##     RColorBrewer-derived) for maximum perceptual distinctness.
##   - Genomic annotations use semantically grouped colors: histone marks
##     share a hue family, repeats another, regulatory elements another, etc.
##   - Driver categories use a colorblind-safe four-color scheme.
##
## Naming convention:
##   <dataset>_<variable>_pal    for color palettes
##   <dataset>_<variable>_shapes for point shape (pch) vectors
##   driver_pal / driver_shapes  for the shared driver categorization
##
## Usage in ggplot:
##   scale_color_manual(values = argelaguet_stage_pal)
##
## Usage in ComplexHeatmap:
##   rowAnnotation(stage = df$stage, col = list(stage = argelaguet_stage_pal))
##
## When only a subset of levels is present, subsetting works directly:
##   scale_color_manual(values = argelaguet_lineage_pal[levels_in_data])

# ---------------------------------------------------------------------------
# Argelaguet gastrulation
# ---------------------------------------------------------------------------

argelaguet_stage_pal <- c(
  "E4.5" = "#7570B3",
  "E5.5" = "#1B9E77",
  "E6.5" = "#D95F02",
  "E7.5" = "#E7298A"
)

## 22 lineages -- Paired + extra hand-picked for the tail end
argelaguet_lineage_pal <- c(
  "Anterior_Primitive_Streak" = "#A6CEE3",
  "Caudal_epiblast"           = "#1F78B4",
  "Embryonic_endoderm"        = "#B2DF8A",
  "Epiblast"                  = "#33A02C",
  "ExE_ectoderm"              = "#FB9A99",
  "ExE_mesoderm"              = "#E31A1C",
  "Gut"                       = "#FDBF6F",
  "Mature_mesoderm"           = "#FF7F00",
  "Mesenchyme"                = "#CAB2D6",
  "Mixed_mesoderm"            = "#6A3D9A",
  "Nascent_mesoderm"          = "#FFFF99",
  "Notochord"                 = "#B15928",
  "Paraxial_mesoderm"         = "#8DD3C7",
  "Parietal_endoderm"         = "#BEBADA",
  "Pharyngeal_mesoderm"       = "#FB8072",
  "Primitive_endoderm"        = "#80B1D3",
  "Primitive_Streak"          = "#FDB462",
  "Rostral_neurectoderm"      = "#B3DE69",
  "Rostral_neuroectoderm"     = "#FCCDE5",
  "Somitic_mesoderm"          = "#D9D9D9",
  "Surface_ectoderm"          = "#BC80BD",
  "Visceral_endoderm"         = "#CCEBC5"
)

## 8 broad lineage classes -- Dark2
argelaguet_lineage_class_pal <- c(
  "Ectoderm"           = "#1B9E77",
  "Endoderm"           = "#D95F02",
  "Epiblast"           = "#7570B3",
  "ExE_ectoderm"       = "#E7298A",
  "Mesoderm"           = "#66A61E",
  "Primitive_endoderm" = "#E6AB02",
  "Primitive_Streak"   = "#A6761D",
  "Visceral_endoderm"  = "#666666"
)

## Argelaguet genomic annotations: grouped by biological category
## - Gene/promoter: blues
## - Repeats (LINEs, SINEs): greens
## - ENCODE histone marks: oranges/reds
## - E7.5 enhancers: purples
## - E7.5 H3K4me3: pinks
## - ESC marks: browns
argelaguet_annotation_pal <- c(
  "Genes"              = "#2166AC",
  "Promoters"          = "#67A9CF",
  "LINEs"              = "#1B7837",
  "SINEs"              = "#7FBC41",
  "H3K4me3 (ENCODE)"   = "#E08214",
  "H3K9me3 (ENCODE)"   = "#B2182B",
  "H3K27me3 (ENCODE)"  = "#D6604D",
  "H3K4me1 (ENCODE)"   = "#F4A582",
  "H3K27ac (ENCODE)"   = "#FDDBC7",
  "Enh E7.5 Ect"      = "#762A83",
  "Enh E7.5 End"      = "#9970AB",
  "Enh E7.5 Mes"      = "#C2A5CF",
  "Enh E7.5 union"    = "#E7D4E8",
  "H3K4me3 E7.5 Ect"  = "#C51B7D",
  "H3K4me3 E7.5 End"  = "#DE77AE",
  "H3K4me3 E7.5 Mes"  = "#F1B6DA",
  "H3K4me3 E7.5 common" = "#FDE0EF",
  "ESC p300"           = "#8C510A",
  "ESC DHS"            = "#BF812D"
)

# ---------------------------------------------------------------------------
# Ecker mouse brain
# ---------------------------------------------------------------------------

## 3 cell classes -- colorblind-safe triad
ecker_cell_class_pal <- c(
  "Exc"  = "#0072B2",
  "Inh"  = "#009E73",
  "NonN" = "#D55E00"
)

## 25 major types -- Paired extended with Set3
ecker_major_type_pal <- c(
  "ANP"       = "#A6CEE3",
  "ASC"       = "#1F78B4",
  "CGE-Lamp5" = "#B2DF8A",
  "CGE-Vip"   = "#33A02C",
  "CLA"       = "#FB9A99",
  "CT-L6"     = "#E31A1C",
  "EC"        = "#FDBF6F",
  "IT-L23"    = "#FF7F00",
  "IT-L4"     = "#CAB2D6",
  "IT-L5"     = "#6A3D9A",
  "IT-L6"     = "#FFFF99",
  "L6b"       = "#B15928",
  "MGC"       = "#8DD3C7",
  "MGE-Pvalb" = "#BEBADA",
  "MGE-Sst"   = "#FB8072",
  "MSN-D2"    = "#80B1D3",
  "NP-L6"     = "#FDB462",
  "ODC"       = "#B3DE69",
  "OPC"       = "#FCCDE5",
  "PAL-Inh"   = "#D9D9D9",
  "PC"        = "#BC80BD",
  "PT-L5"     = "#CCEBC5",
  "Unc5c"     = "#1B9E77",
  "VLMC"      = "#D95F02",
  "VLMC-Pia"  = "#7570B3"
)

## Ecker genomic annotations: grouped by biological category
## - Gene/promoter: blues
## - Repeats: greens
## - Activating marks (H3K4me3, H3K4me1, H3K27ac): oranges/warm
## - Repressive marks (H3K9me3, H3K27me3): reds/cool
ecker_annotation_pal <- c(
  "Genes"    = "#2166AC",
  "Promoters" = "#67A9CF",
  "LINEs"    = "#1B7837",
  "SINEs"    = "#7FBC41",
  "H3K4me3"  = "#E08214",
  "H3K9me3"  = "#B2182B",
  "H3K27me3" = "#D6604D",
  "H3K4me1"  = "#F4A582",
  "H3K27ac"  = "#FDDBC7"
)

# ---------------------------------------------------------------------------
# CRC
# ---------------------------------------------------------------------------

## 6 tumour locations -- qualitative, ordered from normal to distant metastasis
crc_location_pal <- c(
  "NC" = "#1B9E77",
  "PT" = "#D95F02",
  "LN" = "#7570B3",
  "ML" = "#E7298A",
  "MP" = "#66A61E",
  "MO" = "#E6AB02"
)

crc_patient_pal <- c(
  "CRC01" = "#A6CEE3",
  "CRC02" = "#1F78B4",
  "CRC04" = "#B2DF8A",
  "CRC10" = "#33A02C",
  "CRC11" = "#FB9A99",
  "CRC13" = "#E31A1C",
  "CRC15" = "#FDBF6F"
)

crc_patient_shapes <- c(
  "CRC01" = 0L,
  "CRC02" = 1L,
  "CRC04" = 2L,
  "CRC10" = 5L,
  "CRC11" = 6L,
  "CRC13" = 3L,
  "CRC15" = 4L
)

## CRC genomic annotations: grouped by biological category
## - Chromatin states (HMM): greys/muted for numbered states
## - Gene/CpGi/promoter: blues
## - Histone marks: oranges/reds
## - Structural (lamin, PMDs, HMDs): browns/greens
## - Repeats: greens
## - Copy number: purples
crc_annotation_pal <- c(
  "0 Enh."        = "#DADAEB",
  "1 Transc. "    = "#BCBDDC",
  "2 Enh."        = "#9E9AC8",
  "3 Quies."      = "#807DBA",
  "4 Transc."     = "#6A51A3",
  "5 Reg. perm."  = "#54278F",
  "6 Low conf."   = "#D9D9D9",
  "7 Reg. perm."  = "#BDBDBD",
  "8 Quies."      = "#969696",
  "9 Cons. het."  = "#737373",
  "10 Quies."     = "#525252",
  "11 Prom."      = "#252525",
  "12 Prom."      = "#3F007D",
  "13 Cons. het"  = "#4A1486",
  "CpGi"         = "#4292C6",
  "H3K4me3"      = "#E08214",
  "Genes"        = "#2166AC",
  "H3K27me3"     = "#D6604D",
  "H3K9me3"      = "#B2182B",
  "Lamin B1"     = "#8C510A",
  "PMDs"         = "#BF812D",
  "LINEs"        = "#1B7837",
  "SINEs"        = "#7FBC41",
  "HMDs"         = "#DFC27D",
  "SCNA loss"    = "#762A83",
  "SCNA neutral" = "#9970AB",
  "SCNA gain"    = "#C2A5CF"
)

# ---------------------------------------------------------------------------
# Shared: driver categorization (adjH vs adjS)
# ---------------------------------------------------------------------------

driver_pal <- c(
  "across-cell driven" = "#0072B2",
  "within-cell driven" = "#D55E00",
  "both"               = "#009E73",
  "neither"            = "#999999"
)

driver_shapes <- c(
  "across-cell driven" = 17L,
  "within-cell driven" = 15L,
  "both"               = 16L,
  "neither"            = 4L
)

## Argelaguet-specific mark type shapes (enhancer/promoter scatter)
mark_type_shapes <- c(
  "enhancer" = 16L,
  "promoter" = 17L
)

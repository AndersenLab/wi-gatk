# This script will download, index, and prep reference genomes.

# Prep N2 reference genome
RELEASE=WS274
PROJECT=PRJNA13758

mkdir -p "${PROJECT}.${RELEASE}" && cd "${PROJECT}.${RELEASE}"

# Annotation file
wget ftp://ftp.ensembl.org/pub/current_gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.99.gff3.gz

# Download and prepare FASTA Genome
wget "ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/${PROJECT}/sequence/genomic/c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz"

bwa index c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz
gunzip -kfc c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz > c_elegans.${PROJECT}.${RELEASE}.genomic.fa

# This script will download, index, and prep reference genomes.


# RELEASE=WS274
# PROJECT=PRJNA13758

RELEASE=WS274
PROJECT=PRJEB28388 # PD1074

cd "$(git rev-parse --show-toplevel)" || raise "This is not a git repo"
mkdir -p "genome/${PROJECT}_${RELEASE}"
cd "genome/${PROJECT}_${RELEASE}" || raise "Unable access dir"

# Download and prepare FASTA Genome
wget "ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/${PROJECT}/sequence/genomic/c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz"

bwa index c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz
gunzip -kfc c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz > c_elegans.${PROJECT}.${RELEASE}.genomic.fa

# Download Annotations
wget "ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/${PROJECT}/gff/c_elegans.${PROJECT}.${RELEASE}.annotations.gff3.gz"
tabix "c_elegans.${PROJECT}.${RELEASE}.annotations.gff3.gz"

# Unzip and rezip using BGZF
gunzip "c_elegans.${PROJECT}.${RELEASE}.annotations.gff3.gz"
bgzip "c_elegans.${PROJECT}.${RELEASE}.annotations.gff3"

# Index
bedtools sort -i "c_elegans.${PROJECT}.${RELEASE}.annotations.gff3.gz" | bgzip > "c_elegans.${PROJECT}.${RELEASE}.annotations.sorted.gff3.gz"
tabix -p gff "c_elegans.${PROJECT}.${RELEASE}.annotations.sorted.gff3.gz"
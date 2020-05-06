#!/usr/bin/bash
# This script will download, index, and prepare reference genomes for testing purposes
set -e
PARENT_DIR="$(git rev-parse --show-toplevel)"

# Prep N2 reference genome
RELEASE="WS276"
PROJECT="PRJNA13758"

# Create release directory
mkdir -p "${PARENT_DIR}/genomes/${PROJECT}_${RELEASE}"
cd "${PARENT_DIR}/genomes/${PROJECT}_${RELEASE}" || return 1

# Annotation file ~ Note this is different from what is used on QUEST!
wget ftp://ftp.ensembl.org/pub/release-100/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.100.gff3.gz

# Download and prepare FASTA Genome
wget "ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/${PROJECT}/sequence/genomic/c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz"

bwa index c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz
gunzip -kfc c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz > c_elegans.${PROJECT}.${RELEASE}.genomic.fa

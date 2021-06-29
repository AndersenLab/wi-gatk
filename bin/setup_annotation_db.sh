#!/usr/bin/env bash
# Specify wormbase build
# Fetch genome path
genome_path="`brew info snpEff | grep '/data' | cut -f 7 -d ' '`"
# Create directory
mkdir -p ${genome_path}/$1
# Update config file
echo "${1}.genome : C. elegans" >> ${genome_path}/../snpEff.config
echo "${1}.MtDNA.codonTable : Invertebrate_Mitochondrial" >> ${genome_path}/../snpEff.config
# Download genome

wget -O ${genome_path}/${1}/sequences.fa.gz ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.${1}.genomic.fa.gz
# Extract sequence
zcat ${genome_path}/${1}/sequences.fa.gz > ${genome_path}/${1}/sequences.fa
# Download and extract protein fasta file
wget -O \${genome_path}/${1}/protein.fa.gz ftp://ftp.wormbase.org/pub/wormbase/releases/${1}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.${1}.protein.fa.gz 
zcat ${genome_path}/${1}/protein.fa.gz > ${genome_path}/${1}/protein.fa
# Download gtf
wget -O ${genome_path}/${1}/genes.gtf.gz ftp://ftp.wormbase.org/pub/wormbase/releases/${1}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.${1}.canonical_geneset.gtf.gz
# Build genome
snpEff build -gtf22 -v ${1}

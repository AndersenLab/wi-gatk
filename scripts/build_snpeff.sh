# Use this script to build the necessary SnpEff database!
# MUST USE A BUILD >= WS254
# Specify wormbase build
build="WS258"
# Fetch genome path
genome_path="`brew info snpEff | grep '/data' | cut -f 7 -d ' '`"
# Create directory
mkdir -p ${genome_path}/${build}
# Update config file
echo "${build}.genome : C. elegans" >> `dirname $genome_path`/snpEff.config
# Download genome
wget -O ${genome_path}/${build}/sequences.fa.gz ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.WS253.genomic.fa.gz
# Extract sequence
zcat ${genome_path}/${build}/sequences.fa.gz > ${genome_path}/${build}/sequences.fa
# Download gtf
wget -O ${genome_path}/${build}/genes.gtf.gz ftp://ftp.wormbase.org/pub/wormbase/releases/${build}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.${build}.canonical_geneset.gtf.gz
# Build genome
snpEff build -gtf22 -v ${build}
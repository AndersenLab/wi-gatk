#!/bin/bash
# Fetch Genome
# Author: Daniel E. Cook
# Use this script to fetch a reference genome
set -e

# Modify this variable to download the appropriate reference genome
REFERENCE="c_elegans/PRJNA13758/"

# Fetch test data
PARENT_DIR=$(git rev-parse --show-toplevel)
mkdir -p "${PARENT_DIR}/genomes"

# Setup human reference genome
error() {
    >&2 echo -e "\n\t$(tput setaf 1)${1}$(tput sgr0)\n"
    exit 1
}

msg() {
    >&2 echo -e "$(tput setaf 3)${1}$(tput sgr0)"
}

# Copy reference genome from quest
quest_host=$(grep -o 'Host quest' ~/.ssh/config | wc -l)
if [[ ! ${quest_host} -eq 1 ]]; then
    error 'You do not have a quest host configured in your profile.'
fi;

msg "Downloading reference genome: ${REFERENCE}"
# Set rsync defaults using an alias here.
mkdir -p "${PARENT_DIR}/genomes/${REFERENCE}"
rsync -rauL --progress --copy-links "quest:/projects/b1059/data/genomes/${REFERENCE}/*" "${PARENT_DIR}/genomes/${REFERENCE}"


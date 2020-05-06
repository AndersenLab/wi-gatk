![Build Docker (env/gatk4.Dockerfile)](https://github.com/AndersenLab/wi-gatk/workflows/Build%20Docker%20(env/gatk4.Dockerfile)/badge.svg)

# wi-gatk

The new GATK-based pipeline for wild isolate C. elegans strains


# Running Locally

1. Fetch and prepare the reference genome by running `scripts/setup_genome.sh`, adjusting the `$PROJECT` and `$BUILD` as appropriate. This script will create a `genomes` folder, download a reference file, generate a BWA index, and create an unzipped version.

# Software Used

* [Somalier](https://github.com/brentp/somalier)
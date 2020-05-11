![Build Docker (env/gatk4.Dockerfile)](https://github.com/AndersenLab/wi-gatk/workflows/Build%20Docker%20(env/gatk4.Dockerfile)/badge.svg)

# wi-gatk

The new GATK-based pipeline for wild isolate C. elegans strains


# Running Locally

1. Fetch and prepare the reference genome by running `scripts/setup_genome.sh`, adjusting the `$PROJECT` and `$BUILD` as appropriate. This script will create a `genomes` folder, download a reference file, generate a BWA index, and create an unzipped version. Software requirements:

* `picard`
* `samtools`
* `bwa`
* `wget`

2. Install [nextflow](http://www.nextflow.io)
3. Install docker and pull the container ([`andersenlab/gatk4`](https://www.dockerhub.com/andersenlab/gatk4))
4. Run `scripts/download_genome.sh` to fetch a reference genome.

The repo includes test data and is configured to use a design file.

```bash
nextflow run main.nf -resume -profile local --debug
```

__Tip__

* Use [`ctop`](https://github.com/bcicen/ctop) to monitor docker containers locally.
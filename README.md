![Build Docker (env/gatk4.Dockerfile)](https://github.com/AndersenLab/wi-gatk/workflows/Build%20Docker%20(env/gatk4.Dockerfile)/badge.svg)

# wi-gatk

The new GATK-based pipeline for wild isolate Caenorhabditis strains

# Pipeline Overview

```

     _______ _______ _______ __  __      _______ _______ 
    |     __|   _   |_     _|  |/  |    |    |  |    ___|
    |    |  |       | |   | |     <     |       |    ___|
    |_______|___|___| |___| |__|\__|    |__|____|___|    
                                              

To run the pipeline:

nextflow main.nf --help
nextflow main.nf --debug
nextflow main.nf --sample_sheet=/path/sample_sheet.txt --species c_elegans --bam_location=/path/to/bams --gvcf_location=/path/to/gvcfs

    parameters                 description                           Set/Default
    ==========                 ===========                           ========================
    --debug                    Use --debug to indicate debug mode    false
    --species                  Species to call variants from         null
    --sample_sheet             Sample sheet                          null
    --bam_location             Directory of BAM files                {dataDir}/{species}/WI/alignments
    --gvcf_location            Directory of gVCF files               {dataDir}/{species}/WI/gVCFs
    --mito_name                Contig not to polarize hetero sites   MtDNA
    --username                                                       {user}

    Reference Genome
    --------------- 
    --reference                The fa.gz reference file to use       {dataDir}/{species}/genomes/{project}/{ws_build}/{species}.{project}.{ws_build}.genome.fa.gz

    Variant Filters         
    ---------------           
    --min_depth                Minimum variant depth                 3
    --qual                     Variant QUAL score                    30.0
    --strand_odds_ratio        SOR_strand_odds_ratio                 5.0
    --quality_by_depth         QD_quality_by_depth                   20.0
    --fisherstrand             FS_fisher_strand                      100.0
    --high_missing             Max % missing genotypes               0.95
    --high_heterozygosity      Max % max heterozygosity              0.10

```

## Software Requirements

* The latest update requires Nextflow version 24.10.0+. On Rockfish, you can access this version by loading the `nf24_env` conda environment prior to running the pipeline command:

```
ml anaconda
conda activate /data/eande106/software/conda_envs/nf24_env
```

## Typical use for debugging:
```
nextflow run -latest andersenlab/wi-gatk --debug
```

## Typical use for new bam files:
```
nextflow run -latest andersenlab/wi-gatk --sample_sheet=/path/sample_sheet_GATK.txt --species c_elegans --bam_location=/path/to/bams --gvcf_location=/path/to/gvcfs
```

# Parameters

## -profile ( optional )

__default__ = rockfish

There are three configuration profiles for this pipeline.

* `rockfish` - Used for running on Rockfish
* `quest`    - Used for running on Quest
* `local`    - Used for local development

## --sample_sheet

The sample sheet is a list of strains to jointly call variants from, one strain per line

|          |
|:---------|
| AB1      |
| AB4      |
| BRC20067 |

## --bam_location (optional if model worm species specified)

Path to directory holding all the alignment files for strains in the analysis. Defaults to `/vast/eande106/data/{species}/WI/alignments/`

## --gvcf_location (optional if model worm species specified)

Path to directory holding all the gVCF files for strains. If one or more strain doesn't have a gVCF file, for example a new strain, one is created. Defaults to `/vast/eande106/data/{species}/WI/alignments/`

## --species

__default__ = null

Options: c_elegans, c_briggsae, or c_tropicalis (for default reference, bam, and gvcf paths) or any other with reference, bam directory, and gvcf directory specified

## --project (optional, only used with model worm species)

__default__ = PRJNA13758

WormBase project ID for selected species. Choose from some examples [here](https://github.com/AndersenLab/genomes-nf/blob/master/bin/project_species.tsv)

## --ws_build (optional, only used with model worm species)

__default__ = WS283

WormBase version to use for reference genome.

## --reference (optional if model worm species specified)

A fasta reference indexed with BWA. On Rockfish, the reference is available here:

```
/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa.gz
```

>[!Note]
>If running on Rockfish, instead of changing the `reference` parameter, opt to change the `species`, `project`, and `ws_build` for other reference build (and then the reference will change automatically) 

## --mito_name (optional)

__default__ = MtDNA

Name of contig to skip het polarization. Might need to change for other species besides c_elegans if the mitochondria contig is named differently

## -output-dir (optional)

__default__ = WI-{today's date} where the date is formatted as `YYYYMMDD` 

A directory in which to output results

>[!Note]
>This option is a nextflow parameter and so only uses a single dash to specify it

# Output

The final output directory looks like this:

```

├── variation
│   ├── *.hard-filter.vcf.gz
│   ├── *.hard-filter.vcf.tbi
│   ├── *.hard-filter.stats.txt
│   ├── *.hard-filter.filter_stats.txt
│   ├── *.soft-filter.vcf.gz
│   ├── *.soft-filter.vcf.tbi
│   ├── *.soft-filter.stats.txt
│   └── *.soft-filter.filter_stats.txt
└── report
    ├── multiqc.html
	└── multiqc_data
        └── multiqc_*.json
   
```

# Relevant Docker Images
* `quay.io-biocontainers-samtools-1.21--h50ea8bc_0` ([link](https://quay.io/biocontainers/samtools)): Docker image maintained by biocontainers for samtools
* `quay.io-biocontainers-bcftools-1.16--hfe4b78e_1` ([link](https://quay.io/biocontainers/bcftools)): Docker image maintained by biocontainers for bcftools
* `quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0` ([link](https://quay.io/biocontainers/gatk4)): Docker image maintained by biocontainers for gatk
* `quay.io-biocontainers-multiqc-1.8--py_1` ([link](https://quay.io/biocontainers/multiqc)): Docker image maintained by biocontainers for multiqc
* `andersenlab-hetpolarization-1.10` ([link](https://hub.docker.com/r/andersenlab/hetpolarization)): Docker image is created within this pipeline using GitHub actions. Whenever a change is made to `env/hetpolarization.Dockerfile` or `.github/workflows/build_docker.yml` GitHub actions will create a new docker image and push if successful
* `andersenlab-r_packages-v0.7` ([link](https://hub.docker.com/r/andersenlab/r_packages)): Docker image is created manually, code can be found in the [dockerfile](https://github.com/AndersenLab/dockerfile/tree/master/r_packages) repo.

Make sure that you have followed the [Nextflow configuration](https://andersenlab.org/dry-guide/latest/rockfish/rf-nextflow/) described in the dry-guide prior to running the workflow.
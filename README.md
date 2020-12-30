![Build Docker (env/gatk4.Dockerfile)](https://github.com/AndersenLab/wi-gatk/workflows/Build%20Docker%20(env/gatk4.Dockerfile)/badge.svg)

# wi-gatk

The new GATK-based pipeline for wild isolate C. elegans strains


### Typical use for debugging:
```
nextflow main.nf --debug
```

### Typical use for new bam files:
```
nextflow main.nf --sample_sheet=/path/sample_sheet_GATK.tsv --bam_location=/projects/b1059/workflows/alignment-nf/
```

### Parameters
    parameters                 description                           Set/Default
    ==========                 ===========                           ========================
    --debug                    Use --debug to indicate debug mode    (Optional)
    --output                   Release Directory                     alignment-YYYYMMDD
    --sample_sheet             Sample sheet                          (Required)
    --bam_location             Directory of bam files                (Optional but a full path will help Nextflow find files inside of Docker container)
    --reference                Reference Genome                      (Rquired. Default is in quest.config)
    --mito_name                Contig not to polarize hetero sites   MtDNA                                                     
    Reference Genome (these 4 parameters are building parts of params.reference)
    --------------- 
    --reference_base           Location of ref genomes               (Optional. Default is in quest.config)
    --species/project/build    These 4 params form --reference       (Optional. Default is in quest.config)

    Variant Filters         
    ---------------           
    --min_depth                Minimum variant depth                 5
    --qual                     Variant QUAL score                    30.0
    --strand_odds_ratio        SOR_strand_odds_ratio                 5.0 
    --quality_by_depth         QD_quality_by_depth                   20.0 
    --fisherstrand             FS_fisher_strand                      100.0
    --high_missing             Max % missing genotypes               0.95
    --high_heterozygosity      Max % max heterozygosity              0.10



### Notes

1. This pipeline uses a Docker container ([`andersenlab/gatk4`](https://www.dockerhub.com/andersenlab/gatk4)). So on Quest need to do `module add singularity` before running. 

2. The dockerfile can be found at `wi-gatk/env/gatk4.Dockerfile` in this repo. An automatic Gtihub action is set up so that the Docker image gets re-built with each push to Github. In other words, editing the dockerfile and push the change will automatically update the container. 

3. The docker image is usually stored in a Nextflow cache folder (typically specified in config). Each time the Docker image is updated, the existing image in the cache folder needs to be manually deleted, otherwise Nextflow will NOT automatically pull the latest Docker image. Pulling Docker images is best run on Quest login node.

4. GATK requires the genome fasta to be bgzip-ped. It's taken care of in genomes-nf


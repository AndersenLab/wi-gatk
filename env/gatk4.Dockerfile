# docker build --file gatk4.Dockerfile -t swantonlab/manta .
FROM continuumio/miniconda3
RUN apt-get update && apt-get install -y procps && \
    apt-get clean
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge
RUN conda create -n gatk4 \
                        bioconda::gatk4=4.1.4.0 \
                        bioconda::bcftools=1.10 \
                        bioconda::samtools=1.10 \
                        snpeff=4.3.1t \
                        vcflib=1.0.0_rc3 \
                        multiqc=1.8 \
                        parallel=20200322 \
    && conda clean -a
ENV PATH /opt/conda/envs/gatk4/bin:$PATH
# Use libhts.so from conda
ENV LD_LIBRARY_PATH /opt/conda/envs/gatk4/lib
RUN conda env export --name gatk4 > gatk4.yml

# Add het_polarization (pre-built for linux)
ADD het_polarization /usr/local/bin

# Install vcfanno
RUN wget -O /opt/conda/envs/gatk4/bin/vcfanno https://github.com/brentp/vcfanno/releases/download/v0.3.2/vcfanno_linux64 && \
    chmod +x /opt/conda/envs/gatk4/bin/vcfanno

LABEL Name="gatk4" Author="Daniel Cook"

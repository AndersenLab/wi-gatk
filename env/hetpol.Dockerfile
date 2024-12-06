FROM continuumio/miniconda3
RUN apt-get update && apt-get install -y procps && \
    apt-get clean
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge
RUN conda create -n hetpol -y \
                        bioconda::bcftools=1.10 \
                        vcflib=1.0.0_rc3 \
    && conda clean -a
ENV PATH=/opt/conda/envs/hetpol/bin:$PATH
# Use libhts.so from conda
ENV LD_LIBRARY_PATH=/opt/conda/envs/hetpol/lib
RUN conda env export --name hetpol > hetpol.yml

# Add het_polarization (pre-built for linux)
ADD het_polarization /usr/local/bin

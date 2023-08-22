FROM continuumio/miniconda3
LABEL authors="Michael J. Bale" \
      description="Docker image containing all software requirements for the michaeljbale/jlabflow pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/JLabFlow-2.0.0b/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name jlabflow-2.0.0b > jlabflow-2.0.0b.yml

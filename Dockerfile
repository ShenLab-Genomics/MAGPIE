# Before using this docker image, you should apply for all database listed inAdditional requirements
# and download files using `bash download.sh`.

FROM ubuntu
RUN apt-get update && \
    apt-get install -y wget unzip

# Install MATLAB Runtime.
RUN mkdir /mcr-install && \
    cd /mcr-install && \
    wget https://ssd.mathworks.com/supportfiles/downloads/R2022b/Release/7/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2022b_Update_7_glnxa64.zip && \
    unzip MATLAB_Runtime_R2022b_Update_7_glnxa64.zip && \
    ./install -mode silent -agreeToLicense yes

RUN rm -rf /mcr-install
ENV LD_LIBRARY_PATH /usr/local/MATLAB/MATLAB_Runtime/R2022b/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/R2022b/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/R2022b/sys/os/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/R2022b/extern/bin/glnxa64

# Install Miniconda.
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/miniconda/bin:${PATH}

# Copy files.
WORKDIR /app
COPY . /app

# Create conda environments
RUN conda env create -n magpie -f magpie.yml && \
    conda env create -n spliceai -f spliceai.yml 
ENV PATH /opt/conda/envs/magpie/bin:$PATH
ENV PATH /opt/conda/envs/spliceai/bin:$PATH
RUN echo "source activate magpie" > ~/.bashrc

RUN /bin/bash -c "source activate magpie && bash magpie.sh --mode pred --test_file data/datasets/test.csv --visualization"
# RUN /bin/bash -c "source activate magpie && bash magpie.sh --mode train --input_file data/datasets/denovo.csv"
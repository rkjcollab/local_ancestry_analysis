FROM ubuntu:16.04

RUN apt-get update -y && apt-get install -y \
    wget \
    unzip \
    python3 \
    g++
    
# shapeit
RUN wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz  && \
    tar -zxvf shapeit.v2.r904.glibcv2.12.linux.tar.gz  && \
    mv shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit /usr/local/bin  && \
    rm shapeit.v2.r904.glibcv2.12.linux.tar.gz  && \
    rm -r shapeit.v2.904.2.6.32-696.18.7.el6.x86_64

# Genetic map required by shapeit
COPY shapeit_input.tar.gz .
RUN gzip -d shapeit_input.tar.gz && \
    tar -xvf shapeit_input.tar && \
    mv shapeit_input /home/ && \
    rm -r shapeit_input.tar

# RFMix
RUN wget https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip  && \
    unzip RFMix_v1.5.4.zip  && \
    cd RFMix_v1.5.4/PopPhased/  && \
    g++ -Wall -O3 -ftree-vectorize -fopenmp main.cpp getdata.cpp randomforest.cpp crfviterbi.cpp windowtosnp.cpp -o RFMix_PopPhased  && \
    mv RFMix_PopPhased /usr/local/bin  && \
    cd ../..  && \
    rm RFMix_v1.5.4.zip  && \
    rm -r RFMix_v1.5.4 

# Install scripts
COPY run_shapeit.sh /usr/local/bin
COPY run_rfmix.sh /usr/local/bin
RUN mkdir /home/rfmix_file_creation_scripts
COPY *.py /home/rfmix_file_creation_scripts/


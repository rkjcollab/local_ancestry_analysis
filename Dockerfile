FROM python:3.8-slim-buster

RUN apt-get update -y && apt-get install -y \
    wget \
    unzip \
    r-base \
    tabix \
    bcftools \
    vcftools \
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

# CrossMap
RUN pip3 install CrossMap

# admixture
RUN wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz  && \
    tar -zxvf admixture_linux-1.3.0.tar.gz  && \
    mv ./dist/admixture_linux-1.3.0/admixture /usr/local/bin  && \
    rm admixture_linux-1.3.0.tar.gz  && \
    rm -r ./dist

# plink
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip  && \
    unzip plink_linux_x86_64_20200616.zip -d plink_linux_x86_64_20200616/  && \
    mv plink_linux_x86_64_20200616/plink /usr/local/bin  && \
    rm plink_linux_x86_64_20200616.zip  && \
    rm -r plink_linux_x86_64_20200616

# Install scripts
COPY run_*.sh /usr/local/bin/
RUN chmod a+x /usr/local/bin/run_*.sh
RUN mkdir /home/rfmix_file_creation_scripts
COPY ./rfmix_file_creation_scripts/* /home/rfmix_file_creation_scripts/
RUN mkdir /home/shapeit_formatting_scripts
COPY ./shapeit_formatting_scripts/* /home/shapeit_formatting_scripts/


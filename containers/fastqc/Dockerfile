FROM openjdk:18.0.1.1-oraclelinux7

# File Author / Maintainer
MAINTAINER Zarnaz Hadi <bs21zh@leeds.ac.uk>

#Install extra tools
RUN yum -y install wget unzip perl

# Download FASTQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip 
RUN rm fastqc_v0.11.9.zip 
RUN chmod 755 /FastQC/fastqc

# Create a stub directory for /nobackup on ARC
RUN mkdir /nobackup 

# Set working directory
WORKDIR	/FastQC/


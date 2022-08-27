FROM openjdk:19-jdk-oraclelinux7

# File Author / Maintainer
MAINTAINER Zarnaz Hadi <bs21zh@leeds.ac.uk>

# Install extra tools
RUN yum -y install wget unzip

# Download picard using wget
RUN wget https://github.com/broadinstitute/picard/archive/refs/tags/2.27.4.zip && \
    unzip 2.27.4.zip && \
    rm 2.27.4.zip && \
    wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar

RUN chmod +x ./picard.jar

# Create a stub directory for /nobackup on ARC3
RUN mkdir /nobackup

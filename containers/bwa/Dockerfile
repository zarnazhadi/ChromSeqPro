FROM ubuntu

# File Author / Maintainer
MAINTAINER Zarnaz Hadi <bs21zh@leeds.ac.uk>

# Update essential libraries and install bwa
RUN apt-get update && apt-get install -y bwa

# Create a stub directory for /nobackup on ARC3
RUN mkdir /nobackup 

# Entrypoint
ENTRYPOINT ["bwa"]

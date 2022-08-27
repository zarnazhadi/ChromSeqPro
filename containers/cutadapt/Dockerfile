FROM ubuntu

# File Author / Maintainer
MAINTAINER Zarnaz Hadi <bs21zh@leeds.ac.uk>

# Update repository and install essential libraries
RUN apt-get update && apt-get install --yes build-essential

# Install pip and cutadapt required libraries
RUN apt-get install --yes python-pip libpython2.7-dev 
Run apt-get install --yes pip

# Install cutadapt
RUN pip install cutadapt

# Create a stub directory for /nobackup on ARC3
RUN mkdir /nobackup 

# Entrypoint
ENTRYPOINT ["cutadapt"]

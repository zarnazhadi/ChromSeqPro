FROM python:3.7

# File Author / Maintainer
MAINTAINER Zarnaz Hadi <bs21zh@leeds.ac.uk>

# Install python and MACS2
RUN pip install --upgrade pip
RUN pip install numpy 
RUN pip install MACS2

ENTRYPOINT ["macs2"]
CMD ["-h"]

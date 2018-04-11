FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
    apt-get install -y build-essential wget unzip python3 \
    python3-dev python3-pip bats awscli git \
    libcurl4-openssl-dev make gcc zlib1g-dev curl

# Install some helper code 
# TODO pin to release
RUN pip3 install -e git://github.com/FredHutch/aws-batch-helpers.git#egg=aws-batch-helpers

# Run tests and then remove the folder
ADD tests /usr/local/tests
RUN bats /usr/local/tests/ && rm -r /usr/local/tests/

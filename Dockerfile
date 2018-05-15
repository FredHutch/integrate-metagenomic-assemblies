FROM ubuntu:16.04
LABEL MAINTAINER=sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
    apt-get install -y build-essential wget unzip python3 \
    python3-dev python3-pip bats awscli git \
    libcurl4-openssl-dev make gcc zlib1g-dev curl \
    cmake g++ libfile-slurp-perl jq

# Install some helper code (pinned to a commit)
RUN pip3 install -e git://github.com/FredHutch/aws-batch-helpers.git@b06bcaf56479c24532ddcf87f4c2f3722e260144#egg=aws-batch-helpers

# Install BioPython
RUN pip3 install biopython==1.70 tables bucket_command_wrapper==0.3.0

# Install DIAMOND v0.9.10
RUN cd /usr/local/bin && \
    wget -q https://github.com/bbuchfink/diamond/releases/download/v0.9.10/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    rm diamond-linux64.tar.gz

# Install MMseqs2
RUN cd /usr/local/bin && \
    git clone https://github.com/soedinglab/MMseqs2.git && \
    cd MMseqs2 && \
    git checkout 2-23394 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. && \
    make && \
    make install
ENV PATH="/usr/local/bin/MMseqs2/build/bin:${PATH}"

# Make the /share directory and set as the default working directory
RUN mkdir /share
WORKDIR /share

# Add the entire folder
ADD . /usr/local/ima
# Link the main script to the PATH
RUN ln -s /usr/local/ima/integrate_assemblies.py /usr/local/bin/

# Run tests and then remove the folder
RUN bats /usr/local/ima/tests/ && rm -r /usr/local/ima/tests/

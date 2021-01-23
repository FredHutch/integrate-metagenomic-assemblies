FROM ubuntu:20.04
LABEL MAINTAINER=sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" \
    apt-get install -y build-essential wget unzip python3 \
    python3-dev python3-pip bats awscli git \
    libcurl4-openssl-dev make gcc zlib1g-dev curl \
    cmake g++ libfile-slurp-perl jq

# Install BioPython
RUN pip3 install biopython==1.70 tables bucket_command_wrapper==0.3.0

# Install DIAMOND v2.0.6
RUN cd /usr/local/bin && \
    wget -q https://github.com/bbuchfink/diamond/releases/download/v2.0.6/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    rm diamond-linux64.tar.gz

# Install MMseqs2
RUN cd /usr/local/bin && \
    git clone https://github.com/soedinglab/MMseqs2.git && \
    cd MMseqs2 && \
    git checkout 12-113e3 && \
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
RUN ln -s /usr/local/ima/integrate_assemblies.py /usr/local/bin/ && \
    ln -s /usr/local/ima/cluster_proteins.py /usr/local/bin/

# Run tests and then remove the folder
RUN bats /usr/local/ima/tests/ && rm -r /usr/local/ima/tests/

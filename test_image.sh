#!/bin/bash

docker run --rm \
    -v $PWD:/share \
    -v ~/.aws/credentials:/root/.aws/credentials \
    ima \
    python3 /share/integrate_assemblies.py \
    --gff-folder /share/tests \
    --prot-folder /share/tests \
    --output-name TEST \
    --output-folder /share/tests


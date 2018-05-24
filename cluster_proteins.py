#!/usr/bin/env python3

import os
import io
import sys
import gzip
import json
import uuid
import boto3
import shutil
import logging
import argparse
import subprocess
import pandas as pd
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
from batch_helpers.helpers import run_cmds, exit_and_clean_up


def run_mmseqs2(
    input_fp,
    temp_folder,
    min_identity,
    min_coverage=0.5
):
    # Output is a list of dicts, with:
    # "protein_id"
    # "sequence"
    # "members"

    # Make sure the paramters are allowed
    assert isinstance(min_identity, float)
    assert min_identity <= 1.0
    assert min_identity >= 0.0
    assert isinstance(min_coverage, float)
    assert min_coverage <= 1.0
    assert min_coverage >= 0.0

    logging.info("Clustering all proteins with MMSeqs2")
    # Write the database to the temporary folder
    db_name = os.path.join(temp_folder, "db")
    cluster_name = os.path.join(temp_folder, 'mmseqs.cluster')
    tsv_name = os.path.join(temp_folder, 'mmseqs.tsv')
    rep_name = os.path.join(temp_folder, 'mmseqs.rep')
    rep_fasta_name = os.path.join(temp_folder, 'mmseqs.rep.fasta')

    # Make the MMSeqs2 database
    run_cmds(["mmseqs", "createdb", input_fp, db_name])
    # Cluster the protein sequences
    run_cmds([
        "mmseqs", "cluster", db_name, cluster_name, temp_folder,
        "--min-seq-id", str(min_identity),
        "--max-seqs", str(100000), # Don't limit the size of clusters
        "-c", str(min_coverage)])
    # Make TSV output for clustering
    run_cmds(["mmseqs", "createtsv", db_name, db_name, cluster_name, tsv_name])
    # Get the representative sequences
    run_cmds(["mmseqs", "result2repseq", db_name, cluster_name, rep_name])
    run_cmds(["mmseqs", "result2flat", db_name, db_name, rep_name, rep_fasta_name, "--use-fasta-header"])

    # Read in the representative sequences
    with open(rep_fasta_name, "rt") as f:
        rep_seqs = {
            header: seq
            for header, seq in SimpleFastaParser(f)
        }

    # Read in the cluster membership
    cluster_members = defaultdict(list)
    with open(tsv_name, "rt") as f:
        for line in f:
            if len(line) <= 1:
                continue
            # Tab delimited lines
            cluster, member, _ = line.rstrip("\n").split("\t", 2)
            # Make sure the cluster has a representative sequence
            assert cluster in rep_seqs
            cluster_members[cluster].append(member)

    # Format the output
    output = []
    # Loop over the representative sequences
    for cluster_name, seq in rep_seqs.items():
        # Make sure the cluster has members from the TSV
        assert cluster_name in cluster_members
        output.append({
            "cluster": cluster_name,
            "sequence": seq,
            "members": cluster_members[cluster_name]
        })
    logging.info("Read in {:,} clusters".format(len(output)))

    return output


def return_results(
    clustered_proteins,
    temp_folder,
    output_json
):
    """Upload all of the final result files to the output folder."""

    # Write the JSON to a file
    temp_json_path = os.path.join(temp_folder, "output.json")
    logging.info("Writing output to " + temp_json_path)
    with open(temp_json_path, "wt") as fo:
        json.dump(clustered_proteins, fo)
    run_cmds(["gzip", temp_json_path])
    temp_json_path = temp_json_path + ".gz"

    # Check to see if this is an S3 folder to upload to
    if output_json.startswith("s3://"):
        logging.info("Uploading to S3 path: " + output_json)
        run_cmds(["aws", "s3", "cp", temp_json_path, output_json])
    else:
        logging.info("Copying to local path: " + output_json)
        run_cmds(["mv", temp_json_path, output_json])
    logging.info("Copied output to {}".format(output_json))


def get_file(file_path, dest_folder):
    if dest_folder.endswith("/") is False:
        dest_folder = dest_folder + "/"

    assert os.path.exists(dest_folder)

    if file_path.startswith("s3://"):
        run_cmds(["aws", "s3", "cp", file_path, dest_folder])
    else:
        assert os.path.exists(file_path)
        run_cmds(["cp", file_path, dest_folder])

    new_path = os.path.join(dest_folder, file_path.split("/")[-1])
    return new_path


def cluster_proteins(
    input_fasta,
    output_json,
    min_identity,
    temp_folder="/share"
):
    # Make a temporary folder for all files to be placed in
    temp_folder = os.path.join(temp_folder, str(uuid.uuid4())[:8])
    assert os.path.exists(temp_folder) is False
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = os.path.join(temp_folder, "temp.log")
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [cluster_proteins.py] %(message)s'
    )
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to file
    fileHandler = logging.FileHandler(log_fp)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Also write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    assert output_json.endswith(".json.gz")

    # Fetch the input FASTA
    try:
        input_fp = get_file(input_fasta, temp_folder)
    except:
        exit_and_clean_up(temp_folder)

    if input_fp.endswith(".gz"):
        logging.info("Decompressing " + input_fp)
        try:
            run_cmds(["gunzip", input_fp])
        except:
            exit_and_clean_up(temp_folder)
        input_fp = input_fp.replace(".gz", "")

    # Deduplicate protein sequences by clustering
    try:
        # Output is a list of dicts, with:
        # "protein_id"
        # "sequences"
        # "members"
        clustered_proteins = run_mmseqs2(
            input_fp,
            temp_folder,
            min_identity
        )
    except:
        exit_and_clean_up(temp_folder)

    # Upload all of the results
    try:
        return_results(
            clustered_proteins,
            temp_folder,
            output_json
        )
    except:
        exit_and_clean_up(temp_folder)

    logging.info("All done, cleaning up")
    shutil.rmtree(temp_folder)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Cluster a set of proteins.
    """)

    parser.add_argument("--input",
                        type=str,
                        required=True,
                        help="""Input file in FASTA format""")
    parser.add_argument("--output-json",
                        type=str,
                        required=True,
                        help="""Output file in compressed JSON format""")
    parser.add_argument("--min-identity",
                        type=float,
                        required=True,
                        help="""Minimum sequence identity for clustering""")
    parser.add_argument("--temp-folder",
                        type=str,
                        default="/share",
                        help="""Folder for temporary files""")

    # No arguments were passed in
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()

        cluster_proteins(
            args.input,
            args.output_json,
            args.min_identity,
            temp_folder=args.temp_folder
        )

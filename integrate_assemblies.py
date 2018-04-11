#!/usr/bin/env python3

import os
import io
import sys
import gzip
import uuid
import boto3
import logging
import argparse
import subprocess
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
from batch_helpers.helpers import run_cmds, exit_and_clean_up


def split_s3_bucket_prefix(path):
    # Get the bucket name
    bucket = path[5:].split("/", 1)[0]
    # Get the prefix path
    prefix = "/".join(path[5:].split("/")[1:])

    return bucket, prefix


def list_folder(folder):
    """List the contents of a folder, either local or on S3."""
    # Make sure the folder ends with a trailing slash
    if folder.endswith("/") is False:
        folder = folder + "/"

    # Check to see if it is an S3 URL
    if folder.startswith("s3://"):
        bucket, prefix = split_s3_bucket_prefix(folder)
        # Set up the connection
        conn = boto3.client('s3')
        # Iterate over the contents of the folder
        contents = conn.list_objects_v2(
            Bucket=bucket,
            Prefix=prefix
        )
        for key in contents['Contents']:
            # Yield the filepath within the folder
            yield key['Key'].replace(prefix, "")
        while "NextContinuationToken" in contents:
            contents = conn.list_objects_v2(
                Bucket=bucket,
                Prefix=prefix,
                ContinuationToken=contents["NextContinuationToken"]
            )
            for key in contents['Contents']:
                # Yield the filepath within the folder
                yield key['Key'].replace(prefix, "")

    else:
        # Treat as as local path
        assert os.path.exists(folder), "Folder does not exist ({})".format(folder)
        for fp in os.listdir(folder):
            yield fp


def open_handle(path, mode="rt"):
    """Open a file handle for a path, silently handling gzip and S3."""
    if path.startswith("s3://"):
        s3 = boto3.client('s3')
        bucket, key = split_s3_bucket_prefix(path)
        retr = s3.get_object(Bucket=bucket, Key=key)
        if path.endswith(".gz"):
            bytestream = io.BytesIO(retr['Body'].read())
            return io.StringIO(gzip.GzipFile(None, 'rb', fileobj=bytestream).read().decode('utf-8'))
        else:
            return io.StringIO(retr['Body'].read())

    else:
        if path.endswith(".gz"):
            return gzip.open(path, mode)
        else:
            return open(path, mode)


def deduplicate_proteins(
    prot_folder,
    temp_folder,
    output_name,
    prot_suffix="fastp",
    min_seq_id=0.9,
    min_coverage=0.8
):
    # Output is a list of dicts, with:
    # "protein_id"
    # "sequence"
    # "members"

    # Make sure the paramters are allowed
    assert isinstance(min_seq_id, float)
    assert min_seq_id <= 1.0
    assert min_seq_id >= 0.0
    assert isinstance(min_coverage, float)
    assert min_coverage <= 1.0
    assert min_coverage >= 0.0

    # Keep track of all of the sample names
    sample_names = set([])
    # Keep track of which protein is from which sample
    protein_membership = {}

    # Make a single file with all of the proteins
    all_prot_fp = os.path.join(temp_folder, output_name + ".fastp")
    assert os.path.exists(all_prot_fp) is False
    with open(all_prot_fp, "wt") as fo:
        # Rename the proteins to include the sample name
        for fp in list_folder(prot_folder):
            if fp.endswith(("." + prot_suffix, "." + prot_suffix + ".gz")):
                # Get the sample name from the file path
                if fp.endswith(".gz"):
                    sample_name = fp.replace("." + prot_suffix + ".gz", "")
                else:
                    sample_name = fp.replace("." + prot_suffix, "")
                # Make sure this is not a duplicated name
                assert sample_name not in sample_names, "Duplicate name: " + sample_name
                sample_names.add(sample_name)
                # Read the file
                logging.info("Reading from " + fp)
                # Count the number of proteins from each sample
                n_prots = 0
                # Parse as a FASTA
                for header, seq in SimpleFastaParser(
                    # Open the file, regardless of gzip or s3
                    open_handle(os.path.join(prot_folder, fp))
                ):
                    n_prots += 1
                    # Strip the white space from the header
                    protein_name = header.split(" ")[0]
                    # Add the sample name to the protein name
                    protein_name = sample_name + "_" + protein_name
                    # Make sure that this name is unique
                    assert protein_name not in protein_membership
                    # Record which sample this protein is from
                    protein_membership[protein_name] = sample_name

                    fo.write(">{}\n{}\n".format(protein_name, seq))
                logging.info("Read in {:,} sequences from {}".format(n_prots, fp))
    logging.info("Read in a total of {:,} protein sequences".format(len(protein_membership)))

    logging.info("Clustering all proteins with MMSeqs2")
    # Write the database to the temporary folder
    db_name = os.path.join(temp_folder, output_name)
    cluster_name = os.path.join(temp_folder, output_name + '.cluster')
    tsv_name = os.path.join(temp_folder, output_name + '.tsv')
    rep_name = os.path.join(temp_folder, output_name + '.rep')
    rep_fasta_name = os.path.join(temp_folder, output_name + '.rep.fasta')

    # Make the MMSeqs2 database
    run_cmds(["mmseqs", "createdb", all_prot_fp, db_name])
    # Cluster the protein sequences
    run_cmds([
        "mmseqs", "cluster", db_name, cluster_name, temp_folder,
        "--min-seq-id", str(min_seq_id),
        "--max-seqs", str(len(sample_names) * 10), # Don't limit the size of clusters
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
            "sequences": seq,
            "members": cluster_members[cluster_name]
        })
    logging.info("Read in {:,} clusters".format(len(output)))

    return output


def return_results(
    temp_folder,
    output_name,
    output_folder
):
    """Make two files (FASTA + JSON) and write to the output folder."""
    pass


def summarize_proteins(
    gff_folder,
    dedup_proteins
):
    """Output is a list describing all proteins, as in the Readme."""
    return {}


def index_proteins(
    output_name,
    temp_folder
):
    """Make a DIAMOND database for the proteins."""
    pass


def write_results(
    dedup_proteins,
    prot_summary,
    temp_folder,
    output_name,
):
    """Write all of the results to a set of files."""
    pass


def integrate_assemblies(
    gff_folder,
    prot_folder,
    output_name,
    output_folder,
    gff_suffix="gff",
    prot_suffix="fastp",
    temp_folder="/share"
):
    # Make a temporary folder for all files to be placed in
    temp_folder = os.path.join(temp_folder, str(uuid.uuid4())[:8])
    assert os.path.exists(temp_folder) is False
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = os.path.join(temp_folder, output_name + ".log")
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [integrate_assemblies.py] %(message)s'
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

    # Deduplicate protein sequences by clustering
    try:
        # Output is a list of dicts, with:
        # "protein_id"
        # "sequences"
        # "members"
        dedup_proteins = deduplicate_proteins(
            prot_folder,
            temp_folder,
            output_name,
            prot_suffix=prot_suffix
        )
    except:
        exit_and_clean_up(temp_folder)

    # Make a summary JSON with the description of each deduplicated protein
    try:
        # Output is a list describing all proteins, as in the Readme
        prot_summary = summarize_proteins(
            gff_folder,
            dedup_proteins
        )
    except:
        exit_and_clean_up(temp_folder)

    # Write out the results
    try:
        write_results(
            dedup_proteins,
            prot_summary,
            temp_folder,
            output_name,
        )
    except:
        exit_and_clean_up(temp_folder)

    # Make a DIAMOND database for the final proteins
    try:
        index_proteins(
            output_name,
            temp_folder
        )
    except:
        exit_and_clean_up(temp_folder)

    # Upload all of the results
    try:
        return_results(
            temp_folder,
            output_name,
            output_folder
        )
    except:
        exit_and_clean_up(temp_folder)

    exit_and_clean_up(temp_folder)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Integrate the information from a set of assemblies.
    """)

    parser.add_argument("--gff-folder",
                        type=str,
                        required=True,
                        help="""Folder containing GFF files""")
    parser.add_argument("--prot-folder",
                        type=str,
                        required=True,
                        help="""Folder containing protein sequences (FASTA format)""")
    parser.add_argument("--output-name",
                        type=str,
                        required=True,
                        help="""Prefix for output files""")
    parser.add_argument("--output-folder",
                        type=str,
                        required=True,
                        help="""Folder for output files""")
    parser.add_argument("--gff-suffix",
                        type=str,
                        default="gff",
                        help="""Suffix used for GFF files""")
    parser.add_argument("--prot-suffix",
                        type=str,
                        default="fastp",
                        help="""Suffix used for protein sequences (FASTA format)""")
    parser.add_argument("--temp-folder",
                        type=str,
                        default="/share",
                        help="""Folder for temporary files""")

    # No arguments were passed in
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()

        integrate_assemblies(
            args.gff_folder,
            args.prot_folder,
            args.output_name,
            args.output_folder,
            gff_suffix=args.gff_suffix,
            prot_suffix=args.prot_suffix,
            temp_folder=args.temp_folder
        )

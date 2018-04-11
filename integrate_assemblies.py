#!/usr/bin/env python3

import os
import sys
import uuid
import boto3
import logging
import argparse
from batch_helpers.helpers import run_cmds, exit_and_clean_up


def list_folder(folder):
    """List the contents of a folder, either local or on S3."""
    # Make sure the folder ends with a trailing slash
    if folder.endswith("/") is False:
        folder = folder + "/"

    # Check to see if it is an S3 URL
    if folder.startswith("s3://"):
        # Get the bucket name
        bucket = folder[5:].split("/", 1)[0]
        # Get the prefix path
        prefix = "/".join(folder[5:].split("/")[1:])
        # Set up the connection
        conn = boto3.client('s3')
        # Iterate over the contents of the folder
        for key in conn.list_objects(
            Bucket=bucket,
            Prefix=prefix
        )['Contents']:
            # Yield the filepath within the folder
            yield key['Key'].replace(prefix, "")
    else:
        # Treat as as local path
        assert os.path.exists(folder), "Folder does not exist ({})".format(folder)
        for fp in os.listdir(folder):
            yield fp


def deduplicate_proteins(
    prot_folder,
    temp_folder
):
    # Output is a list of dicts, with:
    # "protein_id"
    # "sequences"
    # "members"

    # Make a single file with all of the proteins
    # Rename the proteins to include the sample name
    for fp in list_folder(prot_folder):
        print(fp)
    return {}


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
            temp_folder
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

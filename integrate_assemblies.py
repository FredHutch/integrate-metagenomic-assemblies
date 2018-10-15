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


def read_gff(fp):
    """Read in a GFF file as a Pandas DataFrame."""

    # Names of the columns
    columns = [
        "seqname", "source", "feature",
        "start", "end", "score", "strand",
        "frame", "attribute"
    ]
    n_columns = len(columns)

    # Explicitly coerce the data in each column
    column_types = [str, str, str, int, int, str, str, int, str]

    # Read in the whole file as a list
    gff = []

    for line in open_handle(fp):
        # Skip the commented lines
        if line[0] == '#':
            continue
        
        # Split on tabs
        line = line.rstrip("\n").split("\t")
        
        # Subset to the lines with the right number of columns
        if len(line) != n_columns:
            continue

        keep_line = True
        for v, t in zip(line, column_types):
            try:
                v = t(v)
            except:
                keep_line = False

        if keep_line:
            gff.append(line)

    # Make a DataFrame
    gff = pd.DataFrame(gff, columns=columns)

    return gff


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
    protein_fps,
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
        for fp in protein_fps:
            if fp.endswith(("." + prot_suffix, "." + prot_suffix + ".gz")):
                # Get the sample name from the file path
                if fp.endswith(".gz"):
                    sample_name = fp.split("/")[-1].replace("." + prot_suffix + ".gz", "")
                else:
                    sample_name = fp.split("/")[-1].replace("." + prot_suffix, "")
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
                    open_handle(fp)
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
            "sequence": seq,
            "members": cluster_members[cluster_name]
        })
    logging.info("Read in {:,} clusters".format(len(output)))

    return output


def return_results(
    all_output_files,
    output_folder
):
    """Upload all of the final result files to the output folder."""
    # Make sure the output folder ends with "/"
    if output_folder.endswith("/") is False:
        output_folder = output_folder + "/"

    # Make sure the local files exist    
    for fp in all_output_files:
        assert os.path.exists(fp), "File does not exist locally: " + fp

    # Check to see if this is an S3 folder to upload to
    if output_folder.startswith("s3://"):
        logging.info("Uploading to S3 folder: " + output_folder)
        for fp in all_output_files:
            run_cmds(["aws", "s3", "cp", fp, output_folder])
    else:
        logging.info("Copying to local folder: " + output_folder)
        assert os.path.exists(output_folder)
        for fp in all_output_files:
            run_cmds(["mv", fp, output_folder])
    logging.info("Copied {} files to {}".format(len(all_output_files), output_folder))


def summarize_proteins(
    gff_fps,
    dedup_proteins,
    gff_suffix="gff"
):
    """Output is a list describing all proteins, as in the Readme."""
    # Read in all of the protein annotations from the GFF files
    # For each, record all information about each individual protein,
    # as well as UPSTREAM and DOWNSTREAM genes

    # Key the output by cluster, but ultimately save it as a flat list of the dictionary values
    output = {}

    # Iterate over the files in the gff_folder
    sample_names = set([])
    for fp in gff_fps:
        # Check to see if they are GFF files
        if fp.endswith(("." + gff_suffix, "." + gff_suffix + ".gz")):
            # Get the sample name from the file path
            if fp.endswith(".gz"):
                sample_name = fp.split("/")[-1].replace("." + gff_suffix + ".gz", "")
            else:
                sample_name = fp.split("/")[-1].replace("." + gff_suffix, "")
            # Make sure this is not a duplicated name
            assert sample_name not in sample_names, "Duplicate name: " + sample_name
            sample_names.add(sample_name)
            # Read the file
            logging.info("Reading from " + fp)
            # Read in as a pandas DataFrame
            prots = read_gff(fp)
            logging.info("Read in {:,} total annotations".format(prots.shape[0]))
            # Subset to only the CDS features
            prots = prots.loc[prots["feature"] == "CDS"]
            # Reset the index
            prots.reset_index(drop=True, inplace=True)

            logging.info("Read in {:,} CDS anotations".format(prots.shape[0]))

            # Add columns with the protein attributes
            prots = pd.concat([
                prots,
                pd.DataFrame([
                    dict([
                        field.split("=", 1)
                        for field in attr.split(";")
                    ])
                    for attr in prots["attribute"].values
                ])
            ], axis=1)

            # Delete the attributes column
            del prots["attribute"]

            # Make sure that the ID is an attribute
            assert "ID" in prots.columns, "ID not found in GFF attributes"

            # Add the sample name to the ID
            prots["ID"] = sample_name + "_" + prots["ID"]

            # Map to the protein clusters
            cluster_dict = {
                member: cluster["cluster"]
                for cluster in dedup_proteins
                for member in cluster["members"]
            }
            # Add a column with the protein cluster to the DataFrame
            prots["cluster"] = prots["ID"].apply(cluster_dict.get)

            # Iterate over every CDS feature
            for ix, r in prots.iterrows():
                # Set up the format for every record
                if r["cluster"] not in output:
                    output[r["cluster"]] = {
                        "protein_id": r["cluster"],
                        "members": [],
                        "annotation": []
                    }
                # Add information for this member
                member_annot = r.dropna().to_dict()
                # Remove extraneous data
                for k in ["cluster", "score", "source", "inference", "feature", "frame"]:
                    del member_annot[k]
                output[r["cluster"]]["members"].append(
                    member_annot
                )

                # Add to the cluster-level annotation
                if "product" in member_annot:
                    if member_annot["product"] not in output[r["cluster"]]["annotation"]:
                        output[
                            r["cluster"]
                        ]["annotation"].append(
                            member_annot["product"].rstrip("\n")
                        )

    # Return a list with entries for each unique protein cluster
    return list(output.values())


def write_hdf5_summary(
    dat, 
    fp_out, 
    temp_folder, 
    chunksize=10000, 
    gene_positions_headers=[
        'ID',
        'seqname',
        'start',
        'end',
        'strand'
        'Name',
        'gene',
        'locus_tag',
        'eC_number',
        'product'
    ]
):
    """Write out a summary of the protein clusters in HDF5 format"""
    cluster_members = []
    gene_positions = []
    orf_clusters = {}

    for ix, cluster in enumerate(dat):

        # Iterate over each member of the cluster
        for member in cluster["members"]:
            # Which ORFs are in which ORF clusters
            cluster_members.append({
                "cluster": cluster["protein_id"],
                "member": member["ID"]
            })

            # What are the coordinates for each ORF
            gene_positions.append(member)

            # Which ORFs are in which clusters
            orf_clusters[member["ID"]] = cluster["protein_id"]

        if ix % chunksize == 0 and ix > 0:
            logging.info("Reformatted data for {:,} protein clusters for HDF5 format".format(ix))
    logging.info("Reformatted data for {:,} protein clusters for HDF5 format".format(ix + 1))

    cluster_members = pd.DataFrame(cluster_members)

    store = pd.HDFStore(fp_out)
    logging.info("Writing 'cluster_members' to HDF5")
    cluster_members.to_hdf(
        store,
        'cluster_members',
        format="table",
        data_columns=True
    )
    
    logging.info("Writing 'gene_positions' to HDF5")
    # Make a DataFrame
    gene_positions = pd.DataFrame(gene_positions)

    # Enforce particular headers for the `gene_positions` table
    gene_positions = gene_positions.reindex(columns=gene_positions_headers)

    # Add the annotation of the particular cluster each gene is in        
    gene_positions["cluster"] = gene_positions["ID"].apply(orf_clusters.get)

    for k in ["seqname", "cluster"]:
        msg = "{} not found in headers: {}".format(k, ", ".join(gene_positions.columns.values))
        assert k in gene_positions.columns.values, msg
        null_ix = gene_positions[k].isnull()
        if null_ix.any():
            logging.info("{:,} / {:,} of records being dropped for missing a cluster".format(
                null_ix.sum() / null_ix.shape[0]
            ))
            gene_positions = gene_positions.loc[~null_ix]
    
    try:
        gene_positions.to_hdf(
            store,
            'gene_positions',
            format="table",
            data_columns=["seqname", "cluster"]
        )
    except:
        logging.info("Problem writing gene positions")
        exit_and_clean_up(temp_folder)

    store.close()
    logging.info("Done writing to HDF5")


def write_results(
    dedup_proteins,
    prot_summary,
    temp_folder,
    output_name,
    gene_positions_headers=[
        'ID',
        'seqname',
        'start',
        'end',
        'strand'
        'Name',
        'gene',
        'locus_tag',
        'eC_number',
        'product'
    ]
):
    """Write all of the results to a set of files."""
    # Write out the centroids in FASTA format
    centroid_fasta = os.path.join(temp_folder, output_name + ".fastp")
    logging.info("Writing out " + centroid_fasta)
    with open(centroid_fasta, "wt") as fo:
        logging.info("Writing out centroids")
        for cluster in dedup_proteins:
            fo.write(">{}\n{}\n".format(
                cluster["cluster"],
                cluster["sequence"]
            ))
    # Compress the file
    run_cmds(["gzip", centroid_fasta])
    centroid_fasta = centroid_fasta + ".gz"

    # Make a DIAMOND database of the centroids
    logging.info("Making a DIAMOND database for the centroids")
    run_cmds([
        "diamond", "makedb", 
        "--db", os.path.join(temp_folder, output_name),
        "--in", centroid_fasta
    ])
    dmnd_fp = os.path.join(temp_folder, output_name + ".dmnd")
    assert os.path.exists(dmnd_fp)

    # Write out the protein structure information in JSON format
    summary_json = os.path.join(temp_folder, output_name + ".json")
    logging.info("Writing out " + summary_json)
    with open(summary_json, "wt") as fo:
        json.dump(prot_summary, fo)
    # Compress the file
    run_cmds(["gzip", summary_json])
    summary_json = summary_json + ".gz"

    # Write out the protein structure in HDF5 format
    summary_hdf5 = os.path.join(temp_folder, output_name + ".hdf5")
    logging.info("Writing out " + summary_hdf5)
    write_hdf5_summary(
        prot_summary, 
        summary_hdf5, 
        temp_folder, 
        gene_positions_headers=gene_positions_headers
    )

    all_output_files = [
        centroid_fasta,
        summary_json,
        summary_hdf5,
        dmnd_fp
    ]

    return all_output_files


def overlapping_protein_gff_files(
    prot_folder,
    gff_folder,
    prot_suffix="fastp",
    gff_suffix="gff"
):
    # Key by sample name, file type, and path
    all_fps = defaultdict(lambda: defaultdict(dict))

    for folder, suffix, file_type in [
        (prot_folder, prot_suffix, "protein_fasta"),
        (gff_folder, gff_suffix, "gff")
    ]:
        for fp in list_folder(folder):
            if fp.endswith(("." + suffix, "." + suffix + ".gz")):
                # Get the sample name from the file path
                if fp.endswith(".gz"):
                    sample_name = fp.replace("." + suffix + ".gz", "")
                else:
                    sample_name = fp.replace("." + suffix, "")
                # Make sure this is not a duplicated name
                assert file_type not in all_fps[sample_name], "Duplicate name: " + sample_name
                all_fps[sample_name][file_type] = os.path.join(folder, fp)

    # Subset to those samples with both protein and GFF files available
    all_fps = {
        sample_name: sample_fps
        for sample_name, sample_fps in all_fps.items()
        if len(sample_fps) == 2
    }
    assert len(all_fps) > 0, "Did not find any data in the indicated folders"
    logging.info("Found {:,} samples with both protein FASTA and GFF information".format(len(all_fps)))
    prot_fps = [
        sample_fps["protein_fasta"]
        for sample_fps in all_fps.values()
    ]
    gff_fps = [
        sample_fps["gff"]
        for sample_fps in all_fps.values()
    ]
    return prot_fps, gff_fps


def integrate_assemblies(
    gff_folder,
    prot_folder,
    output_name,
    output_folder,
    gff_suffix="gff",
    prot_suffix="fastp",
    temp_folder="/share",
    gene_positions_headers=[
        'ID',
        'seqname',
        'start',
        'end',
        'strand'
        'Name',
        'gene',
        'locus_tag',
        'eC_number',
        'product'
]
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

    # Get the set of samples that have both proteins and GFF files
    protein_fps, gff_fps = overlapping_protein_gff_files(
        prot_folder,
        gff_folder,
        gff_suffix=gff_suffix,
        prot_suffix=prot_suffix
    )

    # Deduplicate protein sequences by clustering
    try:
        # Output is a list of dicts, with:
        # "protein_id"
        # "sequences"
        # "members"
        dedup_proteins = deduplicate_proteins(
            protein_fps,
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
            gff_fps,
            dedup_proteins,
            gff_suffix=gff_suffix
        )
    except:
        exit_and_clean_up(temp_folder)

    # Write out the results
    try:
        all_output_files = write_results(
            dedup_proteins,
            prot_summary,
            temp_folder,
            output_name,
            gene_positions_headers=gene_positions_headers
        )
    except:
        exit_and_clean_up(temp_folder)

    # Add the log file to the output files
    all_output_files.append(log_fp)

    # Upload all of the results
    try:
        return_results(
            all_output_files,
            output_folder
        )
    except:
        exit_and_clean_up(temp_folder)

    logging.info("All done, cleaning up")
    shutil.rmtree(temp_folder)

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

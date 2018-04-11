# Integrate Metagenomic Assemblies

Combine a set of metagenomic assemblies into a common set of references

**Inputs**

To use this software, generate a set of assemblies with predicted protein-coding
genes in a single FASTA file (e.g. *.fastp.gz) and the annotations in GFF format.

All of those files should be found in a single folder, with the name of the file
matching the name of the sample that the assembly originated from. 


**Outputs**

After running this code, two files should be generated:

  1. A FASTA file with the deduplicated protein-coding gene sequences
  2. A JSON file describing the identity and physical relationship between
  those protein-coding sequences

The format of the JSON will be as follows:

``` json
[
    {
        "protein_id": "<name of deduplicated reference>",
        "members": [
            "<ids of assembled proteins within the group>"
        ],
        "annotation": "<annotation of protein>",
        "neighbors": {
            "<upstream or downstream>": {
                "<protein id of upstream neighbor>": "<number of assemblies with connection>"
            }
        }
    },
    {
        "etc"
    }
]
```

**Invocation**

```
integrate_assemblies.py \
    --gff-folder "<folder containing GFF files>" \
    --prot-folder "<folder containing protein FASTA files>" \
    --output-name "<base name for output files>" \
    --output-folder "<folder for output files>"
```

`integrate_assemblies.py -h` for more options.

**Installation**

To use this code, we recommend using the Docker image, because it contains all the
needed dependencies and is validated with a set of tests. If that is not a satisfying
option, see the Dockerfile to see how to set up the appropriate environment.

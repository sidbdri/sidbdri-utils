Download and process genome sequences and annotations from Ensembl
==========

This script downloads and processes genome sequences and annotations for a particular species from an Ensembl release, e.g.

    download_ensembl mouse 90 <your_email>

The script:
1) Downloads top-level sequences for the species' genome in FASTA format (primary assembly sequences for human and mouse),
2) Downloads gene annotations in GTF format,
3) Creates a genome index for the STAR mapping tool,
4) Builds Salmon and Kallisto indexes for the transcripts defined in the GTF file,
5) Downloads gene information (gene id, name, chromosome, description, Entrez ID) into a
file 'genes.tsv'. Only genes on primary assembly sequences are retained (i.e.
those on patch sequences are removed),
6) Downloads information on mouse and human orthologous genes to files
'<species>_orthologs.tsv'.

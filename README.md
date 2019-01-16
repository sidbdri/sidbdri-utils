Download and process genome sequences and annotations from Ensembl
==========

This script downloads and processes genome sequences and annotations for a particular species from an Ensembl release, e.g.

    download_ensembl mouse 90 <your_email>

The script:
1) Downloads top-level sequences for the species' genome in FASTA format
(primary assembly sequences for human and mouse),
2) Downloads gene annotations in GTF format,
3) Creates a genome index for the STAR mapping tool,
4) Builds Salmon and Kallisto indexes for the transcripts defined in the GTF file,
5) Downloads gene information (gene id, name, chromosome, description, Entrez
ID, gene biotype) into a file 'genes.tsv'. Only genes on primary assembly
sequences are retained (i.e.  those on patch sequences are removed),
6) Downloads transcript information (transcript id, transcript biotype, gene
id, chromosome) into a file 'transcripts.tsv'. Only transcripts for genes on
primary assembly sequences are retained,
7) Downloads information on mouse and human orthologous genes to files
'<species>_orthologs.tsv'.


Now the genome folder will look like this: (example from mouse ensembl 94)
```
|-- BOWTIE2_indices
|   |-- primary_assembly_2.3.4
|   `-- primary_assembly_2.3.4.3
|-- KALLISTO_indices
|   |-- primary_assembly -> KALLISTO_indices/primary_assembly_kallisto0.45.0/kallisto_index
|   |-- primary_assembly_0.43.1
|   `-- primary_assembly_0.45.0
|-- Mus_musculus.GRCm38.94.gtf
|-- README
|-- SALMON_indices
|   |-- primary_assembly -> SALMON_indices/primary_assembly_salmon0.12.0
|   |-- primary_assembly_0.12.0
|   `-- primary_assembly_0.8.2
|-- STAR_indices
|   |-- primary_assembly -> STAR_indices/primary_assembly_STAR2.6.1d
|   |-- primary_assembly_2.5.3a
|   `-- primary_assembly_2.6.1d
|-- download_ensembl.sh
|-- mouse_primary_assembly.fa
|-- primary_assembly
|   |-- 1.fa
|   |-- 10.fa
|   |-- .........
`-- transcripts_ref
    |-- transcripts
    |-- transcripts.chrlist
    |-- transcripts.grp
    |-- transcripts.idx.fa
    |-- transcripts.n2g.idx.fa
    |-- transcripts.seq
    |-- transcripts.ti
    `-- transcripts.transcripts.fa

```

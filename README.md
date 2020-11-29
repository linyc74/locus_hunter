# Locus Hunter
**Mining Genomic Loci of Prokaryotes**

Given a query protein fasta file, in a set of genbank files as the database,
this package finds genomic loci that contains any of the query proteins.
It outputs a genbank file containing all loci, and a figure in both PDF and PNG format.

### Usage
    git clone https://github.com/linyc74/locus_hunter.git
    python locus_hunter -q PATH/TO/QUERY.FAA -g PATH/TO/GENBANK_DIR

Genbank filenames in the folder `GENBANK_DIR` should not contain any blank space.

### Dependency

Download and install [Anaconda](https://www.anaconda.com/products/individual) on either Mac or Linux.
If you are a Windows user, [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) works as well.

Once Anaconda is set up, install the following packages in the terminal:

    conda install -c bioconda dna_features_viewer blast CD-HIT

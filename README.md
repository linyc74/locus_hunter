# Locus Hunter
**Mining Genomic Loci of Prokaryotes**

Given a query protein fasta file, in a set of genbank files (the database),
this package finds genomic loci that contains any of the query proteins.
It outputs a genbank file containing all loci, and a figure in both PDF and PNG format.

### Usage

```bash
git clone https://github.com/linyc74/locus_hunter.git
python locus_hunter -q PATH/TO/QUERY.FAA -g PATH/TO/GENBANK_DIR
```

Genbank filenames in the folder `GENBANK_DIR` should not contain any blank space, e.g. `E coli.gbk` is not allowed.

### Dependency

Download and install [Anaconda](https://www.anaconda.com/products/individual) on either Mac or Linux.
For Windows users, [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) works as well.

Once Anaconda is set up, install the following packages in the terminal:

```bash
pip install numpy pandas scipy ngslite dna_features_viewer
conda install -c bioconda blast cd-hit
```

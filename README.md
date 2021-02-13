# VPsero
rapid serotyping of  Vibrio parahaemolyticus from whole genome sequencing data using serogroup-specific genes

## Introduction

`VPsero` is a software for serotype prediction of Vibrio parahaemolyticus based on Next-Generation 
Sequencing technology. By inputting the strain genome assembly file or prokka genome annotation result, 
it can predict strain’s O/K serotype and determine whether it is a new serotype combination.

## Installation

In order to download `VPsero`, you should clone this repository via the commands

```
git clone https://github.com/shengzheBian/VPsero.git
cd VPsero
chmod 777 ./scripts_of/blastall
chmod 777 ./scripts_of/formatdb
```
In order to install the `Python dependencies`, 
you will need the Anaconda Python distribution and package manager. After installing Anaconda, run the following commands to create an
enviroment with VPsero’s dependencies:

```
conda env create -f enviroment.yaml
source activate VPsero
```

In order to perform a complete analysis process from strain genome assembly file,
you will need to install `prokka` referring to https://github.com/tseemann/prokka. If you have the prokka result, this installation step is not necessary.
## Getting Started
#### In order to get started with the software quickly, you can run the `example data` provided by `VPsero`.

* Perform a complete analysis process from strain genome assembly file:
```
python program.py -i example_data/genome_assembly_seq -o my_out_put_i  -n 5
```
* If you have prokka results, you can set -p parameter to skip genome annotation step：
```
python program.py -p example_data/prokka_result -o my_out_put_p -n 5
```

#### Command line options
```
-h  This help
-i  a director that contains all genome assemble fasta
-p  a director that contains all prokka results
-o  a director that generate analyze result
-n  set the thread number when genome annotate
```

## Output Files:
The key output file is `{output director name}/serotype_predict/04.predict_result/all_strain_predict_result.xlsx`.
The meaning of each column is as following:

| Column Name | Description |
| --------- | ----------- |
| .gff | This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV. |
| .gbk | This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence. |
| .fna | Nucleotide FASTA file of the input contig sequences. |
| .faa | Protein FASTA file of the translated CDS sequences. |
| .ffn | Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA) |

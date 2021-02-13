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
| strain_name | The name of the input strain genome assembly file. |
| O_coaD_contig - O_hldD_direct | The information about O-serogroup gene cluster border genes. |
| K_hldD_contig - K_glpX_direct | The information about K-serogroup gene cluster border genes. |
| O_Spec_Gene | The specific genes found in O-serogroup gene cluster. If the suffix is _a or _b, it means that the O-serogroup needs to be identified by multiple genes |
| K_Spec_Gene | The specific genes found in K-serogroup gene cluster. If the suffix is _a or _b, it means that the K-serogroup needs to be identified by multiple genes. |
| Predict_O_sero | The predicted O-serogroup."One" means that the O-serogoup gene cluster didn't been extracted; "Ont" means that the it may be other known O-serogroup not included in VPsero or OUT. |
| Predict_K_sero | The predicted K-serogroup. "Kne" and "Knt" are similar as "One" and "Ont". |
| New_serotype | "New" means that VPsero predicted new serotype combination not in GB 4789.7-2013; "Exist" means that VPsero predicted existing serotype; "NULL" means that VPsero predicted serotype containing One/Kne or Ont/Knt. |

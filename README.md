# VPsero
#### rapid serotyping of  Vibrio parahaemolyticus from whole genome sequencing data using serogroup-specific genes

## Introduction

`VPsero` is a software for serotype prediction of *Vibrio parahaemolyticus* from genomic sequences generated especially from high throughput Sequencing. 
By inputting the strain genome assembly file or prokka genome annotation result, 
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
-i  a directory that contains all genome assemble fasta
-p  a directory that contains all prokka results
-o  a directory that generate analyze result
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
| O_Spec_Gene | The specific genes found in O-serogroup gene cluster. If the suffix is `_a` or `_b`, it means that the O-serogroup needs to be identified by multiple genes. |
| K_Spec_Gene | The specific genes found in K-serogroup gene cluster. If the suffix is `_a` or `_b`, it means that the K-serogroup needs to be identified by multiple genes. |
| Predict_O_sero | The predicted O-serogroup `"One"` means that the O-serogoup gene cluster didn't been extracted; Most of `"Ont"` might be the serogroups uncovered by VPsero or the sub-popluation of certain serogroup or novel serogroup populations. The prefix `"p"` means that the prediction robustness of this O-serogroup is limited by strain number. |
| Predict_K_sero | The predicted K-serogroup. `"Kne"`, `"Knt"` and `"p"` are similar as `"One"`, `"Ont"` and `"pOx"`. |
| New_serotype | `"New"` means that VPsero predicted new serotype which is combined by known O and K serogroups and not list in China National Food Safety Standard GB 4789.7-2013(see below table S1); `"Exist"` means that VPsero predicted existing serotype; `"NULL"` means that VPsero predicted serotype containing `One/Kne` or `Ont/Knt`. |

  
Supplemental table 1 (China National Food Safety Standard GB 4789.7-2013)
| O serogroup | K serogroup|
|-------------|------------|
|1|1, 5, 20, 25, 26, 32, 38, 41, 56, 58, 60, 64, 69|
|2|3, 28|
|3|4, 5, 6, 7, 25, 29, 30, 31, 33, 37, 43, 45, 48, 54, 56, 57, 58, 59, 72, 75|
|4|4, 8 ,9,10,11,12,13,34,42,49,53,55,63,67,68,73|
|5|15,17,30,47,60,61,68|
|6|18,46|
|7|19|
|8|20,21,22,39,41,70,74|
|9|23,44|
|10|24,71|
|11|19,36,40,46,50,51,61|
|12|19,52,61,66|
|13|65|


## Details of sensitivity and specifity of O and K serogroup
The `sensitivity` and `specifity` information of O and K serogroup is helpful to evaluating the prediction results from `VPsero`


Supplemental table 2
|No.|O serogroup|Sensitivity|Specifity|Report serogroup|
|---|----------|-----------|----------|----------------|
|1|O4|0.910|0.970|O4|
|2|O1|0.930|0.972|O4|
|3|O3|0.890|0.982|O3|
|4|O5|0.950|0.993|O5|
|5|O2|0.940|0.993|O2|
|6|O10|0.740|0.997|O10|
|7|O8|0.790|0.996|O8|
|8|O11|0.950|0.999|O11|
|9|O6|1.000|0.999|O6|
|10|O9|1.000|1.000|O9|
|11|O12|1.000|0.994|pO12|
|12|O7|1.000|0.993|pO7|
  
Supplemental table 3
|No.|K serogroup|Sensitivity|Specificity|Report serogroup|
|--|--|--|--|---|
|1|K6|0.980 |0.983 |K6|
|2|K8|0.970 |0.998 |K8|
|3|K56|0.970 |0.997 |K56|
|4|K36|0.970 |1.000 |K36|
|5|K68|1.000 |1.000 |K68|
|6|K9|1.000 |0.998 |K9|
|7|K25|1.000 |0.992 |K25|
|8|K12|1.000 |0.982 |K12|
|9|K17|1.000 |0.987 |K17|
|10|K3|1.000 |1.000 |K3|
|11|K28|1.000 |1.000 |K28|
|12|K29|1.000 |1.000 |K29|
|13|K18|1.000 |0.998 |K18|
|14|K42|1.000 |0.998 |K42|
|15|K60|1.000 |1.000 |K60|
|16|K11|1.000 |1.000 |K11|
|17|K44|1.000 |1.000 |K44|
|18|K63|1.000 |0.997 |K63|
|19|K5|1.000 |0.998 |K5|
|20|K34|1.000 |0.998 |K34|
|21|K41|0.900 |0.998 |K41|
|22|K13|0.850 |1.000 |K13|
|23|K20|0.710 |1.000 |K20|
|24|K32|1.000 |0.997 |pK32|
|25|K33|1.000 |1.000 |pK33|
|26|K4|1.000 |1.000 |pK4|
|27|K1|1.000 |1.000 |pK1|
|28|K30|1.000 |0.998 |pK30|
|29|K19|0.670 |1.000 |pK19|
|30|K58|0.670 |1.000 |pK58|
|31|K69|0.670 |1.000 |pK69|
|32|K15|1.000 |0.997 |pK15|
|33|K70|1.000 |0.998 |pK70|
|34|K21|1.000 |1.000 |pK21|
|35|K31|1.000 |1.000 |pK31|
|36|K38|1.000 |0.998 |pK38|
|37|K39|1.000 |1.000 |pK39|
|38|K47|1.000 |1.000 |pK47|
|39|K48|1.000 |1.000 |pK48|
|40|K49|1.000 |1.000 |pK49|
|41|K71|1.000 |1.000 |pK71|
|42|K55|1.000 |0.994 |pK55|
|43|K23|1.000 |0.992 |pK23|
|44|K37|-|-|-|
|45|K10|-|-|-|
|46|K53|-|-|-|
## Citiation

## Author

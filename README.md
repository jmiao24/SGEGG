# SGEGG (Simple_GxE_GxG)


## Introduction

SGEGG(Simple_GxE_GxG) is an in-house software to run linear regression for genome-wide GxE or GxG analysis. 

## Updates
- Jan 6, 2022: Initial release. Release the codes for Genome-wide GxG and GxE analysis.

## Prerequisites

The software is developed using R and tested in Linux environments. The statistical computing software R (>=3.5.1) and the following R packages are required:

* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) (>=1.11.8)
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html) (>=1.6.6)
* [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) (>=3.5.1)
* [BEDMatrix](https://cran.r-project.org/web/packages/BEDMatrix/index.html) (>=2.0.3)

## Download SGEGG

You can download SGEGG by:

```
$ git clone https://github.com/jmiao24/SGEGG
$ cd ./SGEGG
```


## Input Data Format
### Phenotype file

The input phenotype file need to be an n x 3 table, where n is the sample size. The columns in order are the FID, IID, and the phenotype. Here is an example of the phenotype file `pheno.txt`:
```
FID	IID	bmi
1	1	31.0920070855422
2	2	34.9722947258896
3	3	26.651358296193
4	4	11.3265202767703
5	5	28.1784063613473
...
...
```


### Covariate file

The input Covariate file needs to be an n x (m+2) table, where m is the number of covariates to include. The columns in order are the FID, IID, and the covariates. Here is an example of the covariate file `covar.txt`:

```
FID	IID	sex	age	pc1
1	1	0	40	0.0496566810517587
2	2	1	41	-0.30689382604118
3	3	1	42	0.0662231532694345
4	4	0	43	0.653681021333332
5	5	1	44	0.50856163868585
...
...
```

### Environment/SNP file for interaction analysis

The input Covariate file needs to be an n x 3 table, where n is the sample size. The columns in order are the FID, IID, and the Environment/SNP for interaction test **(Note: the third colnams of the table must be `I`)**. Here is an example of the phenotype file `interaction.txt`:

```
FID	IID	I
1	1	1
2	2	0
3	3	1
4	4	1
5	5	1
...
...
```


### Genotype file

The input genotype file needs to be in the plink bed/bim/fam format. The path only includes the prefix, not the suffix. For example, the path to the input genotype file is `geno` where the genotype files are `geno.bed, geno.bim, geno.fam`.



## Genome-wide GxE or GxG analysis


#### Example:
The following script perform genome-wide GxE interaction analysis from the `1`-`1000` SNPs in plink genotype file `geno` for phenotype `pheno.txt`.  This analysis adjusted the covaraites in `covar.txt` and `5` cores are used for parallel computing. 
```bash
$ Rscript SGEGG.R \
  --pheno pheno.txt \
  --geno geno \
  --covar covar.txt \
  --interaction interaction.txt \
  --output output_1-100.txt \
  --num_cores 5 \
  --type GxE \
  --start 1 \
  --end 1000
```
where the inputs are

| Flag | Description |
|-----|------------------------------------------------------------------------|
| pheno     | The path to phenotype file|
| geno         | The path of the genotype file following the [PLINK format](https://www.cog-genomics.org/plink/1.9). |
| covar        | The path to the covariate file |
| interaction       | The path to the variable of interaction file |  
| output     | The path to the output summary statistics file |
| num_cores        | Number of cores for parallel computing |
| type        | GxE or GxG analysis (This will only change the columns names for output)  |
| start          | (Optional) Index of SNP that starts computing |
| end       | (Optional) Index of SNP that ends computing |

Note that 

* The script is designed to run on chromosome segments to facilitate parallel computation on the cluster. If `--start` or `--end` is not specified, the script will perform the analysis on all SNPs in the plink `test` file.



#### Explanation of Output

The final result of GxE analysis has the following fields:

| Column | Description |
|-----|-------------|
| CHR | Chromosome |
| SNP | rs ID |
| BP | Base pair position |                                                 
| A1 | Allele 1 (effect allele) |
| A2 | Allele 2 (non-effect allele) |
| FREQ | allele frequency of allele 1 |
| BETA_G | The estimated effect size for genetic main effects |
| SE_G | The estimated standard error of BETA |
| P_G | The P-value for testing variance effects |
| BETA_E | The estimated effect size |
| SE_E | The estimated standard error of BETA |
| P_E | The P-value for testing variance effects |
| BETA_I | The estimated effect size |
| SE_I | The estimated standard error of BETA |
| P_I | The P-value for testing variance effects |
| N | Sample size |

The final result of GxG analysis has the following fields:

| Column | Description |
|-----|-------------|
| CHR | Chromosome |
| SNP | rs ID |
| BP | Base pair position |                                                 
| A1 | Allele 1 (effect allele) |
| A2 | Allele 2 (non-effect allele) |
| FREQ | allele frequency of allele 1 |
 | BETA_G1 | The estimated effect size for SNPs in the geno file |
| SE_G1 | The estimated standard error of BETA_G1 for SNPs in the geno file |
| P_G1 | The P-value for SNPs in the geno file |
| BETA_G2 | The estimated effect size for variable of interaction main effects |
| SE_G2 | The estimated standard error of BETA_G2 for variable of interaction main effects  |
| P_G2 | The P-value for variable of interaction main effects  |
| BETA_I | The estimated effect size of interaction |
| SE_I | The estimated standard error of interaction |
| P_I | The P-value for testing interaction effects |
| N | Sample size |

## Citation

If you use SGEGG, please cite

Miao, J., Lin, Y., Wu, Y., Zheng, B., Schmitz, L. L., Fletcher, J. M., & Lu, Q. (2021). [A quantile integral linear model to quantify genetic effects on phenotypic variability.](https://www.biorxiv.org/content/10.1101/2021.04.14.439847v1) bioRxiv, 2021.2004.2014.439847



## License

All rights reserved for [Lu-Laboratory](http://qlu-lab.org/)

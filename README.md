[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)


## Overview of mediation analysis Workflow

![](graphs/mediation_work_flow.PNG)

## Installation

- __Running environment__: 
    - The workflow was constructed based on the __Linux system__ running the [R v4.1.1](https://cran.r-project.org/).

- __Required software packages__: 
    - [Rstudio](https://www.rstudio.com/products/rstudio/download/)
    - `install.packages(c("data.table", "glmnet", "MASS", "rrBLUP", "parallel", "doParallel", "CompQuadForm"))`
        
```R
install.packages(c("data.table", "glmnet", "MASS", "rrBLUP", "parallel", "doParallel", "CompQuadForm"))
install.packages(devtools)
devtools::install_github("jyanglab/GMA")
```

## Input Data

#### Required input data:

- __y matrix__ (`input/y_matrix.txt`): 
  - A `n x 1` matrix of phenotype values; `n` is the number of individuals.
- __Z matrix__ (`input/Z_matrix.txt`): 
  - The `n x m` genotype matrix of  data; `m` is the number of bi-allelic SNP markers, coded as `-1, 0, 1`.
- __X matrix__ (`input/X_matrix.txt`): 
  - The `n x k` intermediate Omics matrix; `k` is the number of Omics traits. 
  In the example, gene expression data (i.e., RNA-seq read counts of `k` genes) was used as the intermediate traits .


#### Optional input data:
- __X0 matrix__ (`input/X0_matrix.txt`): 
  - A matrix of confounding effects. In the example, three principal components calculated from the Z matrix were used to control population structure.
  


## Major steps

#### Step 1: running the 1_mediation_demo.R to conduct mediation analysis
- Note that you have to adjust the path at the begining in R script; and the path, ntasks, mem, time in the shell script.

```
sbatch workflow/1_mediation_demo.sh
```

#### Step 2: Visualize the results

- You can plot the results yourself using the below R code.
- Note that you have to adjust the path at the begining in R script;

```
2_circos_visual.R
```


## Expected results

The outputs of the example data:  

- `output/res_fixed_bic_trait_V1.csv, res_fixed_eq_trait_V1.csv, res_mixed_linear_trait_V1.csv, res_mixed_shrink_trait_V1.csv` : summaries of the results for different methods of mediation: MedFix_BIC, MedFix_0.5, MedMixed_Linear, MedMixed_Shrink (pmed.pure: proportion of the variance mediated, v.tot: variance of total effect, v.med: variance of indirect effect, v.dir: variance of direct effect, n.med: number of significant mediators after adjustment to control the false discovery proportion (FDP) in mediator selection, pval.cut: the threshold for deciding significance, n.direct: number of direct snps.)
- `output/mediators_fixed_bic_trait_V1.csv, mediators_fixed_eq_trait_V1.csv, mediators_mix_linear_trait_V1.csv, mediators_mix_shrink_trait_V1.csv` : non-adjusted mediators of different methods of mediation: MedFix_BIC, MedFix_0.5, MedMixed_Linear, MedMixed_Shrink (id: mediator gene id, e2m: p-value of effect from exposure to mediator, m2y: p-value of effect from mediator to outcome, e2m2y: maximum value between e2m and m2y, padj: adjusted p-value, coef: product of effect from exposure to mediator and effect from mediator to outcome? )
- `output/dsnps_fixed_bic_trait_V1.csv, dsnps_fixed_eq_trait_V1.csv` : direct SNPs identified (only MedFix methods will report direct SNPs; snp: direct SNP, pval: p-value of effect from exposure to outcome, coef: effect of from exposure to outcome)
- `output/isnps_fixed_bic_trait_V1.csv, isnps_fixed_eq_trait_V1.csv` : indirect SNPs identified (only MedFix methods will report direct SNPs; medi: mediator gene, snps_for_medi: indirect SNPs for the corresponding mediator, coef: effect from exposure to mediator)


Visualization of outputs of MedFix_BIC and GWAS:

![](graphs/circos.PNG)

- In the circos plots, the outermost circular track represents the ten maize chromosomes; the next inner track shows the GWAS results, with two circular blue dashed lines indicating -log(p-value) of 5 and 10 and the red lines denoting the position of direct SNPs;
the next inner track shows the relative positions of identified mediator genes with different genes represented by different colors; the lines in the innermost circle connects mediators with their corresponding indirect SNPs.

## License
It is a free and open source software, licensed under []() (choose a license from the suggested list:  [GPLv3](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/gpl-3.0.txt), [MIT](https://github.com/github/choosealicense.com/blob/gh-pages/LICENSE.md), or [CC BY 4.0](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/cc-by-4.0.txt)).

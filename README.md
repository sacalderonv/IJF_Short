# Replication repository for "Bayesian estimation of a multivariate TAR model when the noise process distribution belongs to the class of Gaussian variance mixtures"

**Vanegas, LH., Calderón, SA., & Rondón, LM.**

Contact: [sacalderonv@unal.edu.co](mailto\:sacalderonv@unal.edu.co)

Preprint available at [arXiv](https://www.arxiv.org/pdf/2503.04593). Conditionally accepted at the [International Journal of Forecasting](https://forecasters.org/ijf).

---

## Overview

This repository contains the replication materials for the paper:

> "Bayesian estimation of a multivariate TAR model when the noise process distribution belongs to the class of Gaussian variance mixtures"

The structure includes scripts, simulation results, and data to replicate the key tables and figures of the manuscript.

### Repository contents

- `R Files for M1 model structure 2 reg/`: R scripts and `.rds` files (1000 replications) for Table 10 and Student-t columns of Tables 1, 2, 5 and 9.
- `R files for M2 structure model 3 reg/`: R scripts and `.rds` files (1000 replications for each distribution) for Tables 3,4 and 6. Initially, it is set to show results for Gaussian distribution.
- `data/`: Contains the empirical application dataset on river flow in Colombia.

---

## Reproducibility instructions

### Step 1: Clone the repository and set working directory
We may choose the folder for M1 structure
```r
git clone https://github.com/sacalderonv/IJF.git
setwd("path/to/IJF/R Files for M1 model structure 2 reg")
```
or for M2 structure

```r
git clone https://github.com/sacalderonv/IJF.git
setwd("path/to/IJF/R files for M2 structure model 3 reg")
```

### Step 2: Install required packages with specific versions

```r
install.packages("devtools")
devtools::install_version("GIGrvg", version = "0.8", dependencies=TRUE)
devtools::install_version("Formula", version = "1.2.5",dependencies=TRUE)
devtools::install_version("Rfast", version = "2.1.5.1", dependencies=TRUE)
devtools::install_version("tsDyn", version = "11.0.5.2" ,dependencies=TRUE)
devtools::install_version("foreach", version = "1.5.2",dependencies=TRUE)
devtools::install_version("doParallel", version = "1.0.17",dependencies=TRUE)
devtools::install_version("ltsa", version = "1.4.6.1",dependencies=TRUE)
devtools::install_version("mtarm",version="0.1.2",dependencies=TRUE)
```
**Remark:** The R version used to run the simulations was 4.3.3, however you can use the recent R version 4.5.1. Additionally, you can also use the most recent version of the package `mtarm`(0.1.6) with the same version of the other packages.

### Step 3: Run full simulations (optional)

You may choose to re-run the simulations from scratch for M1 structure:

```r
source("summarymtar_simulation.R")
source("IJFSimulChequeoDistribution2regbaseFinalParalelizar.R")
```
and for M2 struture

```r
source("summarymtar_simulation.R")
source("SimulayEstimaMtar_Replicas3Reg.R")
```

Or use the saved `.rds` replication results to extract tables directly.

### Step 4: Generate tables

For structure M1
```r
source("Resumen_Replicas_Student.R")
```

For structure M2
```r
source("Resumen_Replicas_3Reg_Final.R")
```

Tables will be printed in the R console. You may redirect output to files if desired.


**Remark 1:** Note that the reference examples for simulating in the context of M1 structure is for student-t distribution error, and for structure 2 is for slash distribution error.

**Remark 2:** For visuslizing results, by default it is considered student-t distribution as true distribution for errors in M1 structure. On the other hand, for M2 structure, we consider Gaussian distribution for errors.

---

## Notes on reproducibility

- No fixed seed were specified in the scripts for simulations and data analysis; results may vary slightly with each execution.
- To ensure reproducibility, you can modify the scripts to include `set.seed(1234)` before any random generation.
- Simulations take approximately 1.37 minutes per replication for M1 and 57 seconds per replication for M2 on an iMac 2019 (3.0GHz i5, 8GB RAM).
- 75 minutes was the approximate time that took the script to generate table 11.

---

## Data description

The data correspond to an empirical application involving two daily river flows and precipitation data in Colombia. Originally analyzed in Calderón & Nieto (2016) [doi](https://doi.org/10.1080/03610926.2014.990758), the dataset is available in `data/riverflows.rda`.

The data set for financial time series in Online Appendix can be obtained from `mtarm` package using


```r
library(mtarm)
data(returns)
```

This data set was used in article Romero and Calderón(2019) [doi](https://doi.org/10.1080/03610926.2019.1669807) and was taken from Banco de la República de Colombia web page and investing platform. 


---

## Wrapper scripts for automation

You may create and use a wrapper script like `run_all_M1.R` with the following content:

```r
setwd("path/to/IJF/R Files for M1 model structure 2 reg")
library(devtools)
library(mtarm); library(GIGrvg); library(Formula); library(Rfast); library(tsDyn)

# Optional: Set reproducibility seed
set.seed(1234)

source("summarymtar_simulation.R")
source("IJFSimulChequeoDistribution2regbaseFinalParalelizar.R")
source("Resumen_Replicas_Student.R")
```

Adjust the file names and parameters as needed for other distributions or model structures.
If you require any other details, refers to  pdf file "[DetailsforScript.pdf](DetailsforScript.pdf)".

---

## License

This work is distributed for academic reproducibility. Please cite the original paper if you use or adapt any part of the code or data.

---

## Changelog

- **2025-07-16**: Improved README for reproducibility, package versions added, clearer script instructions, and wrapper script included.


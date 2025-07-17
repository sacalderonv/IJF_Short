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

- `R Files for M1 model structure 2 reg/`: R scripts and `.rds` files (1000 replications) for Table 10 and Student-t columns of Tables 1, 2, and 5.
- `R files for M2 structure model 3 reg/`: R scripts and `.rds` files (1000 replications for each distribution) for Tables 3 and 6.
- `data/`: Contains the empirical application dataset on river flow in Colombia.

---

## Reproducibility instructions

### Step 1: Clone the repository and set working directory

```r
git clone https://github.com/sacalderonv/IJF.git
setwd("path/to/IJF")
```

### Step 2: Install required packages with specific versions

```r
install.packages("devtools")
devtools::install_version("mtarm", version = "0.1.2", repos = "http://cran.us.r-project.org")
devtools::install_version("GIGrvg", version = "0.7", repos = "http://cran.us.r-project.org")
devtools::install_version("Formula", version = "1.2.4", repos = "http://cran.us.r-project.org")
devtools::install_version("Rfast", version = "2.1.0", repos = "http://cran.us.r-project.org")
devtools::install_version("tsDyn", version = "11.0.4.1", repos = "http://cran.us.r-project.org")
devtools::install_version("foreach", version = "1.5.1", repos = "http://cran.us.r-project.org")
devtools::install_version("doParallel", version = "1.0.16", repos = "http://cran.us.r-project.org")
```

### Step 3: Run full simulations (optional)

You may choose to re-run the simulations from scratch for M1 structure:

```r
source("R Files for M1 model structure 2 reg/summarymtar_simulation.R")
source("R Files for M1 model structure 2 reg/IJFSimulChequeoDistribution2regbaseFinalParalelizar.R")
```
and for M2 struture

```r
source("/R files for M2 structure model 3 reg/summarymtar_simulation.R")
source("/R files for M2 structure model 3 reg/SimulayEstimaMtar_Replicas3Reg.R")
```

Or use the saved `.rds` replication results to extract tables directly.

### Step 4: Generate tables

```r
source("R Files for M1 model structure 2 reg/Resumen_Replicas_Student.R")
source("R files for M2 structure model 3 reg/Resumen_Replicas_3Reg_Final.R")
```

Tables will be printed in the R console. You may redirect output to files if desired.

---

## Notes on reproducibility

- No fixed seed is specified in the simulation scripts; results may vary slightly with each execution.
- To ensure reproducibility, you can modify the scripts to include `set.seed(1234)` before any random generation.
- Simulations take approximately 1.37 minutes per replication for M1 and 57 seconds per replication for M2 on an iMac 2019 (3.0GHz i5, 8GB RAM).

---

## Data description

The data correspond to an empirical application involving two daily river flows and precipitation data in Colombia. Originally analyzed in Calderón & Nieto (2016) [DOI](https://doi.org/10.1080/03610926.2014.990758), the dataset is available in `data/riverflows.rda`.

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

---

## License

This work is distributed for academic reproducibility. Please cite the original paper if you use or adapt any part of the code or data.

---

## Changelog

- **2025-07-16**: Improved README for reproducibility, package versions added, clearer script instructions, and wrapper script included.


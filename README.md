# pairwiseCorrelation

This R package provides functions for calculating pairwise correlations efficiently. It's built with an underlying C implementation to handle large datasets efficiently.

## Features

- Efficient computation of correlation matrices
- Support for datasets with missing values
- Helper functions to extract correlations

## Installation

### Prerequisites

Make sure you have the following prerequisites installed:

- R (version 3.4 or higher)
- Rtools (for Windows users)
- Development tools such as `gcc` for compiling C code (for Mac/Linux)
- `devtools` package for R

#### Linux
- Install libomp (e.g. `yum install libomp`, `dnf install libomp`, etc)


#### macOS
- Install llvm and libomp (`brew install llvm libomp`)
- You may also need to modify ~/.R/Makevars to include the following to override Apple's C compiler:

```
CC=gcc-14
CXX=c++-14
```


### Installing from GitHub

You can install the development version of this package directly from GitHub using the `devtools` package.

1. **Install the `devtools` package** (if you haven't already):

    ```r
    install.packages("devtools")
    ```

2. **Install the package from GitHub**:

    ```r
    devtools::install_github("schaffman5/pairwiseCorrelation")
    ```

## Examples

```r
library(pairwiseCorrelation)

# compute correlations
corrs<-cormat(rbind(c(2.5,4.2,1.6,6.2,3.1), c(5.2,3.5,8.2,7.4,5.3), c(0.3, 0.4, 1.4, 2.7, 3.1) ))

corrs

            1           2
2 -0.07970577            
3  0.35220623  0.42153557

# get correlations for index 1
get.all.dists(1, corrs)

[1]  1.00000000 -0.07970577  0.35220623
```

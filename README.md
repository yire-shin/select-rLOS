

# modselrLOS

**Model Selection for r-Largest Order Statistics (rLOS) Models**

The `modselrLOS` package provides two automated algorithms to select a good r-largest order statistics (rLOS) model with an appropriate value of *r*.  
The algorithms differ in the order of selecting the best model and determining the optimal *r*.

---

## Overview

Two algorithms are implemented:

1. **Algorithm 1: Survival**
   - Selects the **best model first**, then determines an appropriate *r* for the fixed model.

2. **Algorithm 2: r_median**
   - Determines an **appropriate *r*** first, then selects the best model for the fixed *r*.

---

## Available Models

The pool of candidate models for rLOS includes:

- `rK4D`  
- `rK3D`  
- `rGEV`  
- `rGLO`  
- `rLD`  
- `rGGD`  
- `rGD`

---

## Installation

You can install the modselrLOS package by manually downloading the source file (.tar.gz) from GitHub.

Download the file below from the GitHub repository:
modselrLOS_0.0.0.9000.tar.gz

In R or RStudio, run the following command:

```r
install.packages("C:/path_to_file/modselrLOS_0.0.0.9000.tar.gz",
                 repos = NULL,
                 type = "source")
```

※ Replace "C:/path_to_file/" with the actual directory where the file was saved.

## Getting Started

Before using the main functions, load the required packages:

```r
library(lmomco)
library(eva)
library(Rsolnp)
library(modselrLOS)
```

### Main Program 1: sel.rmod.surv

The Survival algorithm selects the best model first, then determines an appropriate r for that model.
This approach focuses on model-based stability, meaning it chooses the most suitable distribution family first and then optimizes the number of order statistics r.

```r
data(pohang)
sel.rmod.surv(pohang)
```

### Main Program 2: sel.rmod.rmed

The r_median algorithm determines an appropriate r first, then selects the best model for that fixed r.
This approach emphasizes data-driven selection, where the stability and information content of r are determined before fitting multiple models.

```r
data(pohang)
sel.rmod.rmed(pohang)
```


## Additional Information

please refer to the supplementary document README_modsel.pdf



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

To install the latest **development version** from GitHub, use:

```r
# Install the {remotes} package if not already installed
install.packages("remotes")

# Install modselrLOS from GitHub
remotes::install_github("yire-shin/select-rLOS")

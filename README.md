# msrlos

**Model Selection for r-Largest Order Statistics (rLOS) Models**

`msrlos` provides model-selection procedures for r-largest order statistics models.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yire-shin/select-rLOS")
library(msrlos)
```

## Input data

`xdat` should be a numeric matrix in which rows represent independent blocks and columns contain order statistics in decreasing order:

```text
x(1) >= x(2) >= ... >= x(R)
```

The package includes the `pohang` dataset, containing the 10 largest daily precipitation observations by year for Pohang, Korea, from 1949 to 2022.

```r
data(pohang)
dim(pohang)
head(pohang)
```

## Usage

### Example: Pohang rainfall data

```r
library(msrlos)
data(pohang)

fit <- ms.rlos(pohang)
fit
```

### General usage

```r
fit <- ms.rlos(
  xdat,
  sig.ed = 0.1,
  dmin = 0.03,
  numh = 8
)
```

- `xdat`: matrix of r-largest order statistics
- `sig.ed`: significance level for ED-based selection
- `dmin`: minimum spacing for fixed-`h` rK3D candidates
- `numh`: number of initial candidate `h` values

## Candidate models

- `rGLO`: generalized logistic
- `rGGD`: generalized Gumbel
- `rGEV`: generalized extreme value
- `rK4D`: four-parameter kappa
- `rLD`: logistic
- `rGD`: Gumbel
- `rK3D`: three-parameter kappa with fixed `h`

## Results

```r
fit$surv
fit$rmed
fit$alt.surv
fit$alt.rmed
```

- `surv`: Survival algorithm
- `rmed`: r-Median algorithm
- `alt.surv`: alternating algorithm initialized by Survival
- `alt.rmed`: alternating algorithm initialized by r-Median

## References

- Bader, B., Yan, J., and Zhang, X. (2017). Automated selection of r for the r largest order statistics approach with adjustment for sequential testing. *Statistics and Computing*, 27, 1435--1451.
- Hosking, J. R. M. (1994). The four-parameter kappa distribution. *IBM Journal of Research and Development*, 38, 251--258.
- Shin, Y. and Park, J. S. (2023). Modeling climate extremes using the four-parameter kappa distribution for r-largest order statistics. *Weather and Climate Extremes*, 39, 100533.
- Shin, Y. and Park, J. S. (2024). Generalized logistic model for r largest order statistics, with hydrological application. *Stochastic Environmental Research and Risk Assessment*, 38, 1567--1581.
- Shin, Y. and Park, J. S. (2025). Generalized Gumbel model for r-largest order statistics, with an application to peak streamflow. *Scientific Reports*, 15, 7614.
- Tawn, J. A. (1988). An extreme-value theory model for dependent observations. *Journal of Hydrology*, 101, 227--250.

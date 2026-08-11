# msrlos

**Automatic Model Selection in the r Largest Order Statistics Approach**

`msrlos` provides automatic model-selection procedures for jointly selecting a probability model and the number of order statistics, `r`, in the r largest order statistics (rLOS) approach.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yire-shin/select-rLOS")
library(msrlos)
```

## Input data

`xdat` should be a numeric matrix or data frame of the r largest order statistics. Each row represents an independent block, typically one year, and the columns contain the order statistics in decreasing order. Thus, the first column contains the block maxima, the second column contains the second-largest observations, and so on.

The package includes the `pohang` dataset, containing the 10 largest daily precipitation observations for each year in Pohang, Korea, from 1949 to 2022.

```r
data(pohang)

dim(pohang)
head(pohang)
```

## Usage

### Example: Pohang rainfall data

`ms.rlos()` runs the complete model-selection procedure.

```r
data(pohang)

fit <- ms.rlos(
  xdat = pohang,
  sig.ed = 0.1,
  dmin = 0.03,
  numh = 8
)

fit
```

| Argument | Description |
|---|---|
| `xdat` | Matrix or data frame of the r largest order statistics |
| `sig.ed` | Significance level for the entropy difference (ED) test |
| `dmin` | Minimum distance used when constructing candidate rK3D models |
| `numh` | Number of candidate `h` values generated for the rK3D models |


## Candidate models

`msrlos` considers the following probability models for the r largest order statistics:

| Model | Description | Parameters |
|---|---|---:|
| `rGLO` | generalized logistic distribution | 3 |
| `rGGD` | generalized Gumbel distribution | 3 |
| `rGEV` | generalized extreme value distribution | 3 |
| `rK4D` | four-parameter kappa distribution | 4 |
| `rLD` | logistic distribution | 2 |
| `rGD` | Gumbel distribution | 2 |
| `rK3D` | three-parameter kappa distribution with fixed `h` | 3 |

The candidate `h` values for the rK3D models are generated automatically from the fitted rK4D model.

## Results

`ms.rlos()` returns the results from four model-selection procedures:

```r
fit$surv
fit$rmed
fit$alt.surv
fit$alt.rmed
```

- `surv`: Survival algorithm, which selects the probability model first and then determines `r`
- `rmed`: r-Median algorithm, which first determines a common value of `r` and then selects the probability model
- `alt.surv`: alternating procedure initialized with the Survival solution
- `alt.rmed`: alternating procedure initialized with the r-Median solution

For example, the selected model, selected value of `r`, and parameter estimates from the Survival algorithm can be obtained by

```r
fit$surv$best.model
fit$surv$rstar
fit$surv$theta
```

The exact returned components vary slightly among the four procedures.

## References

- Bader, B., Yan, J., and Zhang, X. (2017). Automated selection of r for the r largest order statistics approach with adjustment for sequential testing. *Statistics and Computing*, 27, 1435–1451.
- Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values. Springer, London.
- G'Sell, M. G., Wager, S., Chouldechova, A., and Tibshirani, R. (2016). Sequential selection procedures and false discovery rate control. *Journal of the Royal Statistical Society: Series B*, 78, 423–444.
- Shin, Y., and Park, J. S. (2023). Modeling climate extremes using the four-parameter kappa distribution for r-largest order statistics. *Weather and Climate Extremes*, 39, 100533.
- Shin, Y., and Park, J. S. (2024). Generalized logistic model for r largest order statistics, with hydrological application. *Stochastic Environmental Research and Risk Assessment*, 38, 1567–1581.
- Shin, Y., and Park, J. S. (2025). Generalized Gumbel model for r-largest order statistics, with an application to peak streamflow. *Scientific Reports*, 15, 7614.
- Shin, Y., Park, J., and Park, J.-S. (2026). *Automatic model selection in the r largest order statistics approach*. Manuscript.
- Tawn, J. A. (1988). An extreme-value theory model for dependent observations. *Journal of Hydrology*, 101, 227–250.

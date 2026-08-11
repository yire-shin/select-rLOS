# Methods and code map

This file summarizes the statistical procedures implemented in `msrlos` and the main functions associated with them.

## Workflow

```r
fit <- ms.rlos(xdat, sig.ed = 0.1, dmin = 0.03, numh = 8)
```

`ms.rlos()` performs the common ED/AIC calculations and returns four results:

```text
ms.rlos()
  |
  +-- redtest.all()
  |
  +-- sel.surv()   -> surv
  +-- sel.rmed()   -> rmed
  +-- Alteria.new(start = "surv") -> alt.surv
  +-- Alteria.new(start = "rmed") -> alt.rmed
```

Main files:

- `R/ms-rlos.R`: top-level interface
- `R/redtest-all.R`: ED/AIC calculations shared by the procedures
- `R/selection-survival.R`: Survival procedure and fixed-`h` candidate construction
- `R/selection-rmedian.R`: r-Median procedure
- `R/alternating.R`: alternating procedures and cycle resolution

## Candidate models

| Internal name | Model | Parameters |
|---|---|---:|
| `rglo` | generalized logistic | 3 |
| `rggd` | generalized Gumbel | 3 |
| `rgev` | generalized extreme value | 3 |
| `rk4d` | four-parameter kappa | 4 |
| `rld` | logistic | 2 |
| `rgd` | Gumbel | 2 |
| `rk3d` | three-parameter kappa with fixed `h` | 3 |

The fitting routines are in `R/fit-*.R`.

## Fixed-h rK3D candidates

`h.select()` fits the rK4D, obtains an ED-based order, and constructs candidate `h` values around the fitted second shape parameter. Candidates too close to the fitted value, `h = -1`, or `h = 0`, and values violating the rK4D constraints are removed. `dmin` controls the minimum spacing.

Main functions:

- `h.select()`
- `rk4d.lik.park.fisher()`

## Selection of r

The entropy-difference (ED) procedure compares consecutive orders and is applied across candidate values of `r`. ForwardStop and StrongStop quantities are also calculated for sequential testing.

Main functions:

- `multi.rEdtest.park()`
- model-specific `*Ed.park()` functions
- `one_optr()` and related helpers

Files: `R/entropy-difference-tests.R`, `R/optimal-r.R`.

## AIC comparison

Models are compared by AIC only at a common value of `r`.

Main functions:

- `com.aic()`
- `com.aic.fast()`
- `AIC1.model.rfix()`

Files: `R/aic.R`, `R/aic-fixed-r.R`.

## Model-selection procedures

### Survival

`sel.surv()` starts with all candidate models and eliminates the model with the largest AIC as `r` increases. The selected model is then paired with its ED-selected order.

### r-Median

`sel.rmed()` obtains an ED-selected `r` for each candidate model, takes their median as a common order, and selects the model with the smallest AIC at that order.

### Alternating procedures

`Alteria.new()` alternates between model selection at the current `r` and selection of `r` for the current model. It is initialized from either Survival or r-Median. `solve.cycle()` handles cycling solutions.

## Parameter convention

The common four-element parameter representation is

```text
(mu, sigma, k, h)
```

The implementation generally uses the Hosking-style sign convention for the GEV/Kappa shape parameter. Some internal likelihood routines use the opposite Coles-style sign, with conversions performed where required.

## References

- Bader, B., Yan, J., and Zhang, X. (2017). Automated selection of r for the r largest order statistics approach with adjustment for sequential testing. *Statistics and Computing*, 27, 1435--1451.
- G'Sell, M. G., Wager, S., Chouldechova, A., and Tibshirani, R. (2016). Sequential selection procedures and false discovery rate control. *Journal of the Royal Statistical Society: Series B*, 78, 423--444.
- Hosking, J. R. M. (1994). The four-parameter kappa distribution. *IBM Journal of Research and Development*, 38, 251--258.
- Shin, Y. and Park, J. S. (2023). Modeling climate extremes using the four-parameter kappa distribution for r-largest order statistics. *Weather and Climate Extremes*, 39, 100533.
- Shin, Y. and Park, J. S. (2024). Generalized logistic model for r largest order statistics, with hydrological application. *Stochastic Environmental Research and Risk Assessment*, 38, 1567--1581.
- Shin, Y. and Park, J. S. (2025). Generalized Gumbel model for r-largest order statistics, with an application to peak streamflow. *Scientific Reports*, 15, 7614.
- Tawn, J. A. (1988). An extreme-value theory model for dependent observations. *Journal of Hydrology*, 101, 227--250.

Two algorithms to select a good r-largest order statistics (rLOS) model with an appropriate r are implemented. The rst one is the algorithm Survival. The second is the algorithm r_median.
The rst algorithm selects the best model rst, then determines an appropriate r for the fixed model. The second algorithm determines an appropriate r, then selects the best model for the
fixed r. The pool of models for rLOS are rK4D, rK3D, rGEV, rGLO, rLD, rGGD, and rGD.

1. Installation :
Install the latest development version (on GitHub) via {remotes}:

remotes::install_github("yire-shin/select-rLOS")

Getting started :

Main program-1: sel.rmod.surv
The algorithm Survival selects the best model first, then determines an appropriate r for the fixed model.

Example :
library(lmomco)
library(eva)
library(Rsolnp)
library(modselrLOS)
data(pohang)
sel.rmod.surv(pohang)

Main program-2: sel.rmod.rmed

The algorithm r_median determines an appropriate r, then selects the best model for the fixed r.

Example :
library(lmomco)
library(eva)
library(Rsolnp)
library(modselrLOS)
data(pohang)
sel.rmod.rmed(pohang)

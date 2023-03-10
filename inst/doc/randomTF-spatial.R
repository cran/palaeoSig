## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(palaeoSig)
library(rioja)
library(sf)
library(gstat)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(ggplot2)
set.seed(42) # for reproducibility 

## ---- results='hide', warning=FALSE-------------------------------------------
# load data
data(Atlantic)
meta <- c("Core", "Latitude", "Longitude", "summ50")
Atlantic <- as.data.frame(Atlantic) # prevents rowname warnings

# pseudocore as no fossil foram data in palaeoSig
fosn <- Atlantic |>
  filter(between(summ50, 5, 10)) |>
  slice_sample(n = 20)

# remaining samples as training set
Atlantic <- Atlantic |>
  anti_join(fosn, by = "Core") |> 
  slice_sample(n = 300) # random subset to speed analysis up

Atlantic_meta <- Atlantic |>
  select(one_of(meta)) # to keep rdist.earth happy
Atlantic <- Atlantic |> # species
  select(-one_of(meta))

fos <- fosn |>
  select(-one_of(meta))

## ---- message=FALSE, results = 'hide', fig.cap = "Figure 2: Semi-variogram fitted to detrended residuals of a weighted averaging model."----
Atlantic_meta <- st_as_sf(
  x = Atlantic_meta,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

## ---- fig.cap = "Figure 1: Semi-variogram model (Matérn class) fitted to the Atlantic summer sea temperature at 50 m depth."----
# Estimate the variogram model for the environmental variable of interest
ve <- variogram(summ50 ~ 1, data = Atlantic_meta)
vem <- fit.variogram(
  object = ve,
  model = vgm(40, "Mat", 5000, .1, kappa = 1.8)
)
plot(ve, vem)
vem

## -----------------------------------------------------------------------------
# Simulating environmental variables
sim <- krige(sim ~ 1,
  locations = Atlantic_meta,
  dummy = TRUE,
  nsim = 100,
  beta = mean(Atlantic_meta$"summ50"),
  model = vem,
  newdata = Atlantic_meta
)

# convert sf back to a regular data.frame
sim <- sim |> st_drop_geometry()

## -----------------------------------------------------------------------------
rtf_auto <- randomTF(
  spp = Atlantic,
  env = Atlantic_meta$summ50,
  fos = fos,
  autosim = sim,
  fun = MAT,
  col = "MAT.wm"
)

plot(rtf_auto)

## -----------------------------------------------------------------------------
rtf_ind <- randomTF(
  spp = Atlantic,
  env = Atlantic_meta$summ50,
  fos = fos,
  fun = MAT,
  col = "MAT.wm"
)

plot(rtf_ind)


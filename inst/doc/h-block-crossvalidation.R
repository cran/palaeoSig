## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE------------------------------------------------------
library(palaeoSig)
library(rioja)
library(sp)
library(gstat)
library(tidyverse)

# suppress warnings from MAT about duplicate samples (presumably 100% N. pachyderma)
MAT <- function(...){
  suppressWarnings(rioja::MAT(...))
}

## ---- results='hide', warning=FALSE--------------------------------------
#load data
data(Atlantic)
meta <- c("Core", "Latitude", "Longitude", "summ50")

#N Atlantic
N_Atlantic <- Atlantic %>% 
  filter(Latitude > 3)
N_Atlantic_meta <- N_Atlantic %>% 
  select(one_of(meta)) %>% 
  as.data.frame() # to keep rdist.earth happy
N_Atlantic <- N_Atlantic %>% 
  select(-one_of(meta)) 

#S Atlantic
S_Atlantic <- Atlantic %>% 
  filter(Latitude < -3)
S_Atlantic_meta <- S_Atlantic %>% 
  select(one_of(meta))
S_Atlantic <- S_Atlantic %>% 
  select(-one_of(meta))


#calculating distances among the sampled points in the North Atlantic foraminifera data set 
geodist <- fields::rdist.earth(select(N_Atlantic_meta, Longitude, Latitude), miles = FALSE)

#values of h for which h-block cross-validation is calculated
threshs <- c(0.01, 100, 200, 400, 600, 700, 800, 900, 1000, 1200, 1400, 1600)

#h-block cross-validation of the NA foraminifera dataset for different values of h
res_h <- map_df(threshs, function(h){
  mod <- MAT(N_Atlantic, N_Atlantic_meta$summ50, k = 5, lean = FALSE)
  mod <- crossval(mod, cv.method = "h-block", h.dist = geodist, h.cutoff = h)
  
  tibble(h = h, 
         RMSE = performance(mod)$crossval["N05", "RMSE"],
         R2 = performance(mod)$crossval["N05", "R2"]
  )
})

## ---- warning = FALSE, results = 'hide', tidy = TRUE---------------------
#Leave-one-out cross-validated RMSEP using MAT with k = 5
round(res_h[1, 'RMSE'], 2)
#Predicting the Sout Atlantic test set
mod.NA <- MAT(N_Atlantic, N_Atlantic_meta$summ50, k = 5)
pred.SA <- predict(mod.NA, newdata = S_Atlantic)$fit
#Determining RMSEP of the SA test set
rmse.mat <- sqrt(mean((pred.SA[, 1] - S_Atlantic_meta$summ50)^2))
#RMSEP of the SA test set using MAT with k = 5
round(rmse.mat, 2)

## ---- echo=FALSE, results = 'hide', fig.cap = "Figure 1: Root mean square error of prediction (RMSEP) as a function of removal distance h. Dashed horizontal line indicates RMSEP found on a spatially independent test set."----
est_h <- approx(y = res_h$h, x = res_h$RMSE, xout = rmse.mat)$y
seg_dat <- tibble(x = c(-Inf, est_h), 
              xend = c(est_h, est_h),
              y = c(rmse.mat, rmse.mat), 
              yend = c(rmse.mat, -Inf)
)

ggplot(res_h, aes(x = h, y = RMSE)) + 
  geom_point() +
  geom_line() +
  geom_segment(aes(x = x, y= y, xend = xend, yend = yend), data = seg_dat, colour = "red") +
  labs(x = "h [km]", y = "RMSEP [°C]") 

## ---- message=FALSE, results = 'hide', fig.cap = "Figure 2: Semi-variogram fitted to detrended residuals of a weighted averaging model."----
#WA model
modwa <- crossval(WA(sqrt(N_Atlantic), N_Atlantic_meta$summ50, mono = TRUE))
#residuals of the WA model
wa.resid <- residuals(modwa, cv = TRUE)
#detrend to remove edge effects (loess with span = 0.1) 
detrended_resid <- resid(loess(wa.resid[, 1] ~ N_Atlantic_meta$summ50, span = 0.1))
#copy meta data and add coordinate system
N_Atlantic_meta_c <- N_Atlantic_meta
coordinates(N_Atlantic_meta_c) <- ~ Longitude + Latitude
proj4string(N_Atlantic_meta_c) <- CRS("+proj=longlat +datum=WGS84")
#variogram of the detrended residuals of the WA model 
v <- variogram(detrended_resid ~ 1, data = N_Atlantic_meta_c)
#Fitting a spherical variogram (partial sill, range and nugget are approximately estimated from the empirical variogram)
vm <- fit.variogram(v, vgm(psill = 2, "Sph", range = 1500, nugget =  0.5))
plot(v, vm)

## ---- fig.cap = "Figure 3: Semi-variogram model (Matérn class) fitted to the North Atlantic summer sea temperature at 50 m depth."----
#Estimate the variogram model for the environmental variable of interest
ve <- variogram(summ50 ~ 1, data = N_Atlantic_meta_c)
vem <- fit.variogram(ve, vgm(40, "Mat", 5000, .1, kappa = 1.8))
plot(ve, vem)

#Simulating environmental variables
sim <- krige(sim ~ 1, locations = N_Atlantic_meta_c, dummy = TRUE, nsim = 100, beta = mean(N_Atlantic_meta[,"summ50"]), model = vem, newdata = N_Atlantic_meta_c)

#convert spatialpointsdataframe back to a regular data.frame
sim <- as.data.frame(sim) %>% 
  select(-Longitude, -Latitude)

## ---- message = FALSE, warning= FALSE,results='hide', fig.cap="Figure 4: Histogram of squared correlation coefficients between simulated variables and the environmental variable of interest."----

#Function for h-block cross-validating several simulations at a time
mat.h1 <- function (y, x, noanalogues, geodist, thresh)
{
    if (!inherits(y, "dist")) {
        if (is.data.frame(y) || !(ncol(y) == nrow(y) & sum(diag(y)) ==
            0)) {
            y <- dist(sqrt(y))^2#squared chord distance
        }
    }
    y <- as.matrix(y)
    diag(y) <- Inf
    if (inherits(geodist, "dist"))
        geodist = as.matrix(geodist)
    sapply(seq_len(nrow(y)), function(n) {browser()
        exneigh <- geodist[n, ] >= thresh
        x2 <- x[exneigh, ]
        y2 <- y[n, ][exneigh]
        analogues <- which(rank(y2, ties.method = "random") <= noanalogues)
        colMeans(x2[analogues, ])
    })
}


#h-block cross-validation of the simulated variables
simhr <- sapply(threshs, function(h) {
  hn <- mat.h1(N_Atlantic, sim, noanalogues = 5, geodist = geodist, thresh = h)
    diag(cor(t(hn), sim)^2)
})
 
#Estimating squared correlation between environmental varible of interest and simulated variables  
sim.obs.r2 <- sapply(sim, cor, N_Atlantic_meta$summ50)^2
#Calculating sum of squares between the two squared correlations
so.squares <- apply(simhr, 2, function(x){
  sum((x - sim.obs.r2) ^ 2)
  })

## ----message = FALSE, warning= FALSE,results='hide', fig.cap = "Figure 5: Scatterplot of squared correlation coefficients between simulated variables and the environmental variable of interest and transfer function r^2^."----

simhr %>% 
  as.data.frame() %>% 
  set_names(threshs) %>% 
  mutate(sim.obs.r2 = sim.obs.r2) %>% 
  gather(key = h, value  = value, -sim.obs.r2) %>% 
  mutate(h = factor(h, levels = threshs)) %>% 
  ggplot(aes(x = sim.obs.r2, y = value)) +
  geom_point() + 
  geom_abline() + 
  facet_wrap(~h) + 
  labs( x = "Simulated-observed environmental r²", y = "Transfer function r²")

## ----fig.cap =  "Figure 6: Relationship between the sum of squares between the two r^2^ as function of distance *h*."----
tibble(threshs, so.squares) %>% 
  ggplot(aes(x = threshs, y = so.squares)) +
  geom_point() + 
  labs(x = "h km", y = "Sum of squares")


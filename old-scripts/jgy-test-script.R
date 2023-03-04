# Install all the packages--------------------------------------------------------

install.packages("tidyverse")
install.packages("sf")
install.packages("fields")
install.packages("furrr")
install.packages(
    "INLA",
    repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), 
    dep=TRUE
)

# Load packages and test reading a shapefile--------------------------------------

library(INLA)
library(tidyverse)
library(sf)
library(fields)
library(furrr)

nc <- st_read(system.file("shape/nc.shp", package="sf"))
print(paste0("I read a shapefile with ", nrow(nc), " rows"))

# Test INLA script from website---------------------------------------------------

n = 100; a = 1; b = 1; tau = 100
z = rnorm(n)
eta = a + b*z

scale = exp(rnorm(n))
prec = scale*tau
y = rnorm(n, mean = eta, sd = 1/sqrt(prec))


data = list(y=y, z=z)
formula = y ~ 1+z
result = inla(formula, family = "gaussian", data = data)

summary(result)
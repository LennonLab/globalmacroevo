## functions.R
## Ford Fishman
## Common functions used across scripts
require("png")
require("grid")
require("tidyr")
require("ggplot2")
require("viridis")
require("scales")
require("here")

## Setup environment
# rm(list=ls()) # removes all objects in the given environment
wd <- here()
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)

## Transformation for plot
magnify_trans <- function(intercept, reducer) {
  
  trans <- function(x, i = intercept, r = reducer) {
    sapply(x, function(x) {
      if (x < i) x
      else x / r + i
    })
  }
  
  inv <- function(x, i = intercept, r = reducer) {
    sapply(x, function(x) {
      if(!is.na(x)) {
        if (x < i) x
        else (x - i) * r
      }
    })
  }
  
  trans_new(name = 'custom',
            transform = trans,
            inverse = inv
  )
}

# reverse log axis tranformation
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}


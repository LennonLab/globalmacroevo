## Ford Fishman
## December 16, 2019
## Common functions used across scripts
require("png")
require("grid")
require("tidyr")
require("ggplot2")
require("viridis")
require("scales")

## Setup environment
rm(list=ls()) # removes all objects in the given environment
wd <- "~/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)

## import Rohde & Muller extinction events data
# Originally from MicroSpeciation-Logistic3

# Say x is diversity and y is extinction intensity for an era lasting 3 years
# Current dataset has data like the following for an era lasting 3 My:
# Time  Diversity   ExtinctionIntensity
# t     x           y 
# t+1   x           y
# t+2   x           y

# I'm altering it to look like this to spread out extinction over
# the entire length of the era (avg_ext):

# Time  ExtinctionIntensity
# t     y/3 
# t+1   y/3
# t+2   y/3

# Alternatively (peak_ext):

# Time  ExtinctionIntensity
# t     0 
# t+1   0
# t+2   y
importRohdeMuller <- function(path = paste(data_dir, "Extinction_Rohde_Muller_raw.csv", sep = "")){
  ext <- read.csv(path)
  original <- ext$ExtinctionIntensity # Untouched extinction values, for the purposes of making graphs
  peak_ext <- c() # Only consider extinction events at end of era
  avg_ext <- c()  # Average out exitinction percentage over entire time period
  era_lengths <- c() # vector of era lengths
  era_means <- c() #vector of mean era extinction intensities
  
  for (i in 1:nrow(ext)){
    
    if (i == 1){
      # First entry in the dataframe starts a paleontologic era
      era <- c(ext$ExtinctionIntensity[i])
      
      # For all entries after the first: 
    } else {
      # I assume that when either ExtinctionIntensity or Diversity changes,
      # then a new era has begun, otherwise, an era continues
      if (ext$Diversity[i] == ext$Diversity[i-1] & ext$ExtinctionIntensity[i] == ext$ExtinctionIntensity[i-1]) {
        era <- c(era, ext$ExtinctionIntensity[i])
        if (i == nrow(ext)) { # if the last entry is its own era
          
          # integer(x) makes a vector of zeros of length x
          # the mean of extinction intensity is evenly distributed over an era's duration
          avg_ext <- c(avg_ext, integer(length(era)) + mean(era)/length(era))
          a <- integer(length(era))
          a[1] <- mean(era)
          peak_ext <- c(peak_ext, a)
          era_lengths <- c(era_lengths, length(era))
          era_means <- c(era_means, a[1])
        }
        
      } else { # if a new era is beginning
        # Add previous era into vector
        avg_ext <- c(avg_ext, integer(length(era)) + mean(era)/length(era))
        a <- integer(length(era))
        a[1] <- mean(era)
        peak_ext <- c(peak_ext, a)
        era_lengths <- c(era_lengths, length(era))
        era_means <- c(era_means, a[1])
        era <- c(ext$ExtinctionIntensity[i]) # begin recording new era
        
        if (i == nrow(ext)) { # if the last entry is part of another era
          avg_ext <- c(avg_ext, integer(length(era)) + mean(era)/length(era))
          a <- integer(length(era))
          a[1] <- mean(era)
          peak_ext <- c(peak_ext, a)
          era_lengths <- c(era_lengths, length(era))
          era_means <- c(era_means, a[1])
        }
      }
    }
  }
  ext$time <- as.character(4000 - ext$Mya)
  names(peak_ext) <- ext$time
  peak_ext <- peak_ext[order(names(peak_ext))]/100  
  names(avg_ext) <- ext$time
  avg_ext <- avg_ext[order(names(avg_ext))]/100
  names(original) <- ext$time
  original <- original[order(names(original))]/100
  return(list(original=original,peak_ext=peak_ext,avg_ext=avg_ext,era_lengths=era_lengths,era_means=era_means))
}

## Used for binning into colors for heatmap
rounding_scheme <- function(S){
  logS <-  log10(S)
  # roundS <- ifelse(logS > 30, 30, logS) # Anything >10^30 is one color
  roundS <- ifelse(logS > 23, 23, logS) # Anything <10^30 and >10^20 is another color
  roundS <- ifelse(roundS < 4, 4, roundS) # can't be less than 0
  roundS <- floor(roundS) # round down order of magnitude
  return(data.frame(logS, roundS))
}

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



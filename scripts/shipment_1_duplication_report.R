#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(lubridate)
library(geosphere)
library(googlesheets)
library(geonames)
library(rebus)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))
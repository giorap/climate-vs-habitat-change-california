######################################################
##### -- climate-vs-habitat-change-california -- #####
######################################################
##################### ANALYSIS #######################
##### -- Set things up -- #####
### Load packages
my_packages <- c("dplyr", "magrittr", "knitr", ### data-handling
                 "ggplot2", "fields", "RColorBrewer", "shiny", "grid", ### graphics
                 "raster", "rgdal", "sp", "rgeos", "maptools", ### geospatial
                 "foreach" ### parallel
)
lapply(my_packages, require, character.only = TRUE)
### Install and load library multispeciesPP
# NOTE: RUN THE LINE BELOW ONLY THE FIRST TIME YOU RUN THIS SCRIPT
# install.packages("multispeciesPP_1.0.tar.gz", repos = NULL, type="source")
library(multispeciesPP)
### Load custom functions
source("src/climate-vs-habitat-change-california-functions.R")

#### -- Load presence-absence data -- ####
### NOTE: Number of sites is less than total number available because all sites not overlapping with historical vegetation survey coverage have been excluded
### Birds
t1_pa_bird <- readRDS("data/t1_pa_bird.rds")
t2_pa_bird <- readRDS("data/t2_pa_bird.rds")
### Mammals
t1_pa_mamm <- readRDS("data/t1_pa_mamm.rds")
t2_pa_mamm <- readRDS("data/t2_pa_mamm.rds")
#### -- Load presence-only data -- ####
### Birds
t1_po_bird <- readRDS("data/t1_po_bird.rds")
t2_po_bird <- readRDS("data/t2_po_bird.rds")
### Mammals
t1_po_mamm <- readRDS("data/t1_po_mamm.rds")
t2_po_mamm <- readRDS("data/t2_po_mamm.rds")

#### -- Load background data -- ####
### NOTE: The spatial grain of the environmental datasets is 30-arcsec (~ 800-m) and the spatial extent is limited by the extent of the historic vegetatation survey (~ 2/3 of California)
### Historic (T1)
t1_bg <- readRDS("data/t1_bg.rds")
### Modern (T2)
t2_bg <- readRDS("data/t2_bg.rds")

#### -- Load habitat associations -- ####
habitat_associations_bird <- readRDS("data/habitat_associations_bird.rds")
habitat_associations_mamm <- readRDS("data/habitat_associations_mamm.rds")

#### -- Load final species sets -- ####
final_species_set_bird <- readRDS("data/final_species_set_bird.rds")
final_species_set_mamm <- readRDS("data/final_species_set_mamm.rds")

#### -- Build models -- ####
source("src/climate-vs-habitat-change-california-models.R")

#### -- Extract model output -- ####
mPP_out <- multispeciesPP_output()
mPP_out <- readRDS("mPP_out.rds")
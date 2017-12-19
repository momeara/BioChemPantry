# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(plyr)
library(dplyr)
library(readr)
library(BioChemPantry)
library(SEAR)


staging_directory <- BioChemPantry::get_staging_directory("chembl21")
load(paste0(staging_directory, "/data/full_chembl_data.Rdata"))




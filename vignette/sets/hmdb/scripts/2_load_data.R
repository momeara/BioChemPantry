# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(assertthat)
library(stringr)
library(readr)
library(BioChemPantry)

schema <- "hmdb_180122"
pantry <- BioChemPantry::get_pantry(schema)
staging_directory <- BioChemPantry::get_staging_directory(schema)

m_data <- readr::read_tsv(paste0(staging_directory, "/data/metabolites.tsv"))
c_data <- readr::read_tsv(paste0(staging_directory, "/data/metabolite_concentrations.tsv"))

m_tbl <- pantry %>% dplyr::copy_to(
	m_data,
	name=c("metabolites"),
	temporary=F,
	indexes=list("accession"))

c_tbl <- pantry %>% dplyr::copy_to(
	c_data,
	name=c("concentrations"),
	temporary=F,
	indexes=list("accession"))


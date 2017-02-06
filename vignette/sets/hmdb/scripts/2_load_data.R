# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(assertthat)
library(stringr)
library(readr)
library(BioChemPantry)

pantry <- get_pantry()
pantry %>% create_schema('hmdb_170128')
pantry %>% set_schema('hmdb_170128')


z_data <- readr::read_tsv(paste0(staging_directory, "/dump/hmdbendo_substances_170201.tsv")) %>%
	dplyr::select(
		zinc_id,
		supplier_code,
		preferred_name = substance.preferred_name,
		zinc_smiles = substance.smiles)

m_data <- readr::read_tsv(paste0(staging_directory, "/data/metabolites.tsv"))

m_tbl <- pantry %>% dplyr::copy_to(
	m_data,
	name=c("metabolites"),
	temporary=F,
	indices=list("accession"))

z_tbl <- pantry %>% dplyr::copy_to(
	z_data,
	name=c("zinc_map"),
	temporary=F,
	indices=list("zinc_id"))


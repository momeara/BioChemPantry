# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(assertthat)
library(stringr)
library(readr)
library(BioChemPantry)

pantry <- BioChemPantry::get_pantry()
#pantry %>% BioChemPantry::create_schema('hmdb_170128')
pantry %>% BioChemPantry::set_schema('hmdb_170128')
staging_directory <- BioChemPantry::get_staging_directory('hmdb')

z_data <- readr::read_tsv(paste0(staging_directory, "/dump/hmdbendo_substances_170201.tsv")) %>%
	dplyr::select(
		zinc_id,
		supplier_code,
		preferred_name = substance.preferred_name,
		zinc_smiles = substance.smiles)

m_data <- readr::read_tsv(paste0(staging_directory, "/data/metabolites.tsv"))

kingdom_is_chemical_entities <- m_data$taxonomy_kingdom == "Chemical entities"
m_data <- m_data %>%
	dplyr::mutate(
		taxonomy_kingdom = ifelse(kingdom_is_chemical_entities, taxonomy_super_class, taxonomy_kingdom),
		taxonomy_super_class = ifelse(kingdom_is_chemical_entities, taxonomy_class, taxonomy_super_class),
		taxonomy_class = ifelse(kingdom_is_chemical_entities, taxonomy_sub_class, taxonomy_class),
		taxonomy_sub_class = ifelse(kingdom_is_chemical_entities, NA, taxonomy_sub_class))

m_tbl <- pantry %>% BioChemPantry::copy_to.src_postgres(
	m_data,
	name=c("metabolites"),
	temporary=F,
	fast=T,
	indices=list("accession"))




z_tbl <- pantry %>% BioChemPantry::copy_to.src_postgres(
	z_data,
	name=c("zinc_map"),
	temporary=F,
	fast=T,
	indices=list("zinc_id"))


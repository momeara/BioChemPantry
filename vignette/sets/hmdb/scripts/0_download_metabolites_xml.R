# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(curl)
library(BioChemPantry)
library(Zr)
library(plyr)
library(dplyr)

schema <- paste0("hmdb_180122")
staging_directory <- BioChemPantry::get_staging_directory(schema)

dir.create(paste0(staging_directory, "/dump"), recursive=TRUE)
dump_fname <- paste0(staging_directory, "/dump/hmdb_metabolites.zip")
curl::curl_download("http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip", dump_fname)

unzip(dump_fname, exdir=paste0(staging_directory, "/dump"))



hmdbendo_substances <- Zr::catalog_items(
	"hmdbendo",
	output_fields=c(
		"zinc_id",
		"supplier_code",
		"substance.preferred_name",
		"substance.smiles"),
	result_batch_size=10000,
	verbose=T) %>%
	dplyr::select(
		accession=supplier_code,
		zinc_id,
		zinc_name = substance.preferred_name,
		zinc_smiles = substance.smiles)

hmdbendo_substances %>% readr::write_tsv(
		paste0(staging_directory, "/data/hmdbendo_substances.tsv"))

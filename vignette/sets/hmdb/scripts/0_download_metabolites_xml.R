# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(curl)
library(BioChemPantry)
library(Zr)

staging_directory <- BioChemPantry::get_staging_directory("hmdb")

dir.create(paste0(staging_directory, "/dump"), recursive=TRUE)
dump_fname <- paste0(staging_directory, "/data/hmdb_metabolites.zip")
curl::curl_download("http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip", dump_fname)

unzip(dump_fname, exdir=paste0(staging_directory, "/dump")))
unlink(dump_fname)


metabolites <- Zr::catalog_items(
	"hmdbendo",
	output_fields=c(
		"zinc_id",
		"supplier_code",
		"substance.preferred_name",
		"substance.smiles"),
	result_batch_size=1000,
	verbose=T)

metabolites %>% readr::write_tsv("data/hmdbendo_substances_170201.tsv")

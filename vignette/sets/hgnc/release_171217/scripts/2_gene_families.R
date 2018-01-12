# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(jsonlite)
library(BioChemPantry)
library(assertthat)

staging_directory <- BioChemPantry::get_staging_directory("hgnc/release_171217")

system(paste0("
cd ", staging_directory, "
wget https://www.genenames.org/cgi-bin/genefamilies/download-all/json
mv json gene_families.json"))

json_fname <- paste0(staging_directory, "/gene_families.json")
assertthat::assert_that(file.exists(json_fname), msg="The HGNC data was downloaded correctly.")
gene_families <- jsonlite::fromJSON(txt=json_fname)


pantry %>% dplyr::copy_to(
	gene_families,
	"gene_families",
	temporary=F,
	indices=list(
		"symbol",
		"gene_family_id"))


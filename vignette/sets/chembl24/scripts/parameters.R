# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(BioChemPantry)
schema <- "chembl24"
sea_schema <- "sea_chembl24"

chembl_version <- "24"

staging_directory <- BioChemPantry::get_staging_directory(schema)
dump_dir <- paste0(staging_directory, "/dump/postgres")



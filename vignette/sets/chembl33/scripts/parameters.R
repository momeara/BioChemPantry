library(BioChemPantry)
schema <- "chembl33"
sea_schema <- "sea_chembl33"

chembl_version <- "33"

staging_directory <- BioChemPantry::get_staging_directory(schema)
dump_dir <- paste0(staging_directory, "/dump/postgres")


user_agent <- "httr mattjomeara@gmail.com"
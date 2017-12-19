# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(BioChemPantry)

schema <- "csfmetabolome"


pantry <- get_pantry()
pantry %>% create_schema(schema)

pantry_login <- get_pantry_config()$login
staging_directory <- get_staging_directory(schema)

xml_url <- "http://www.csfmetabolome.ca/system/downloads/current/csf_metabolites.zip"
dump_dir <- paste0(staging_directory, "/dump/")
zip_fname <- paste0(dump_dir, "csf_metabolites.zip")
xml_fname <- paste0(dump_dir, "csf_metabolites.xml")


cat("
Downloading CSF Metabolome XML data


   url: ", xml_url, "
   dump_url: ", dump_dir, "

and loading it into the postgres database

  host: ", pantry_login$host, "
  port: ", pantry_login$port, "
  username: ", pantry_login$user, "
  password: **********
  schema: ", schema, "

", sep="")


# download and unzip data
paste0("
mkdir -p ", dump_dir, "
cd ", dump_dir, "
wget ", xml_url, "
unzip ", zip_fname, "
") %T>% cat %>% system




# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(assertr)
library(BioChemPantry)

# this downloads chembl and uploads into the pantry database in the given schema.

schema <- "chembl22"
db_url <- "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_22_1_postgresql.tar.gz"

pantry_login <- get_pantry_config()$login
staging_directory <- get_staging_directory(schema)

dump_dir <- paste0(staging_directory, "/dump/postgres")
tar_fname <- paste0(dump_dir, "/", db_url %>% stringr::str_extract("[^/]+$"))
dump_fname <- paste0(
	tar_fname %>% stringr::str_replace(".tar.gz$", ""), "/",
	tar_fname %>% stringr::str_extract("[^/]+$") %>% stringr::str_replace("_postgresql.tar.gz", ".pgdump.sql"))

cat("
Downloading ChEMBL postgres database

  url: ", db_url, "
  dump_dir: ", dump_dir, "

and loading it into the postgres database

  host: ", pantry_login$host, "
  port: ", pantry_login$port, "
  username: ", pantry_login$user, "
  password: **********
  schema: ", schema, "

", sep="")


# download and unzip SQL dump
paste0("
mkdir -p ", dump_dir, "
cd ", dump_dir, "
wget ", db_url, "
tar xzvf ", tar_fname, "
") %T>% cat %>% system

# set schema

paste("
sed -i'' -e 's/SET search_path = public, pg_catalog;/SET search_path = ", schema, ", pg_catalog;/g' ", dump_fname, "
", sep="") %T>% cat %>% system

# load into database
# you may have to type in password
paste(
		"psql ",
		"--host=", pantry_login$host, " ",
		"--port=", pantry_login$port, " ",
		"--username=", pantry_login$user, " ",
		"--password=", pantry_login$password, " ",
		"< ", dump_fname, "
", sep="") %T>% cat %>% system


# delete dump
unlink(dump_dir)

# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(assertr)
library(BioChemPantry)

# this downloads chembl and uploads into the pantry database in the given schema.

schema <- "chembl25"
db_url <- "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_25_postgresql.tar.gz"

pantry <- get_pantry()
pantry %>% create_schema(schema)

# for loading restoring the pg_dump
pantry_login <- get_pantry_config()$login
staging_directory <- get_staging_directory(schema)

dump_dir <- paste0(staging_directory, "/dump/postgres")
#dump_dir <- "/mnt/nfs/work/momeara/pantry_sets/chembl25/dump/postgres"

tar_fname <- paste0(dump_dir, "/", db_url %>% stringr::str_extract("[^/]+$"))
#tar_fname <- "/mnt/nfs/work/momeara/pantry_sets/chembl25/dump/postgres/chembl_25_postgresql.tar.gz"

dump_fname <- paste0(dump_dir, "/chembl_25/chembl_25_postgresql/chembl_25_postgresql.dmp")


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
# and don't create the plpgsql extension because it already exists and creating it requires super-user privileges
# This assumes that the schema name in the sql is the same as the one on phi server
paste0(
    "sed ",
        "-i'' ",
            “-e ‘s/SELECT pg_catalog.set_config('search_path', '', false);/SELECT pg_catalog.set_config('search_path',” schema, “ , false);/g’ “,
        "-e 's/CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;//g' ",
        "-e "s/COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';//g" ",
        "-e 's/Owner:  /Owner: ", pantry_login$user, "/g' ",

dump_fname) %T>% cat %>% system


#Old Matt's code, no longer work with ChEMBL25 sql format
paste0(
	"sed ",
		"-i'' ",
	  "-e 's/SET search_path = public, pg_catalog;/SET search_path = ", schema, ", pg_catalog;/g' ",
		"-e 's/CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;//g' ",
		"-e \"s/COMMENT ON EXTENSION plpgsql IS 'PL\/pgSQL procedural language';//g\" ",
		"-e 's/Owner: user/Owner: ", pantry_login$user, "/g' ",
		"-e 's/OWNER TO \"user\"/OWNER TO \"", pantry_login$user, "\"/g' ",
		dump_fname) %T>% cat %>% system

# load into database
# you may have to type in password
paste(
		"psql ",
		"--host=", pantry_login$host, " ",
		"--port=", pantry_login$port, " ",
		"--username=", pantry_login$user, " ",
		"--password ",
		"-v ON_ERROR_STOP=1",
		"< ", dump_fname, "
", sep="") %T>% cat %>% system


# delete dump
unlink(dump_dir)

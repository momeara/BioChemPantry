library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(assertr)
library(BioChemPantry)

source("scripts/parameters.R")
# this downloads chembl and uploads into the pantry database in the given schema

dump_dir <- paste0(staging_directory, "/dump/postgres")
if(!dir.exists(dump_dir)) {
  cat("Creating folder to load chembl data from '", dump_dir, "'\n", sep = "")
  dir.create(dump_dir, recursive = TRUE)
}


db_url <- paste0(
  "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/",
  "chembl_", chembl_version, "_postgresql.tar.gz")

pantry <- BioChemPantry::get_pantry()
pantry |> BioChemPantry::create_schema(schema)

# for loading restoring the pg_dump
pantry_login <- BioChemPantry::get_pantry_config()$login
staging_directory <- BioChemPantry::get_staging_directory(schema)

tar_fname <- paste0(dump_dir, "/", db_url |> stringr::str_extract("[^/]+$"))
dump_fname <- paste0(
  dump_dir, "/chembl_", chembl_version, "/chembl_", chembl_version,
  "_postgresql/chembl_", chembl_version, "_postgresql.dmp")


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
cmd <- paste0("
cd ", dump_dir, "
wget ", db_url, "
tar xzvf ", tar_fname, "
")
cat(cmd)
system(cmd)

# set schema
# and don't create the plpgsql extension because it already exists and creating
# it requires super-user privileges
cmd <- paste0(
	"sed ",
		"-i'' ",
#	  "-e 's/SET search_path = public, pg_catalog;/SET search_path = ", schema, ", pg_catalog;/g' ",
#    "-e 's/SELECT pg_catalog.set_config(\'search_path\', '', false);/SELECT pg_catalog.set_config(\'search_path\', '", schema, "', false);/g' ",	
#		"-e 's/CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;//g' ",
#		"-e \"s/COMMENT ON EXTENSION plpgsql IS 'PL\\/pgSQL procedural language';//g\" ",
		"-e 's/Owner: user/Owner: ", pantry_login$user, "/g' ",
		"-e 's/OWNER TO \"user\"/OWNER TO \"", pantry_login$user, "\"/g' ",
		dump_fname)
cat(cmd)
system(cmd)

# load into database
# you may have to type in password
cmd <- paste(
		"psql ",
		"--host=", pantry_login$host, " ",
		"--port=", pantry_login$port, " ",
		"--username=", pantry_login$user, " ",
		"--password ",
#		"-v ON_ERROR_STOP=1",
		"< ", dump_fname, "
", sep="")
cat(cmd)
system(cmd)

# delete dump
unlink(dump_dir)

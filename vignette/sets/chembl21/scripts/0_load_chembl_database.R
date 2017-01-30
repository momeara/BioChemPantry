# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(assertr)
library(BioChemPantry)

# this downloads chembl and uploads into the pantry database in the given schema.

schema <- "chembl21"

pantry_login <- get_pantry_config()$login
staging_directory <- get_staging_directory("chembl21/dump/chembl_21_postgresql")


# download and unzip SQL dump
system(paste("
cd ", staging_directory, "
wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_21_postgresql.tar.gz
tar xzvf chembl21_postgresql.tar.gz
", sep=""))

# set schema to chembl21
system(paste("
cd ", staging_directory, "
sed -i'' -e 's/SET search_path = public, pg_catalog;/SET search_path = ", schema, ", pg_catalog;/g' chembl_21.pgdump.sql
", sep=""))

# load into database
# you may have to type in password
system(paste(
		"psql ",
		"--host=", pantry_login$host, " ",
		"--port=", pantry_login$port, " ",
		"--username=", pantry_login$username, " ",
		"--password ", pantry_login$password, " ",
		"< ", staging_directory, "/chembl_21.pgdump.sql
", sep=""))


# delete dump
unlink(staging_directory)

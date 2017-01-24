# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(assertr)
source("~/work/sea/scripts/data_repo.R")

data_repo <- get_data_repo()
data_repo %>% create_schema("chembl21")


# download and unzip SQL dump
system("
mkdir dump/
cd dump
wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_21_postgresql.tar.gz
tar xzvf chembl21_postgresql.tar.gz
")

# set schema to chembl21
system("
cd dump/chembl_21_postgresql/
sed -i'' -e 's/SET search_path = public, pg_catalog;/SET search_path = chembl21, pg_catalog;/g' chembl_21.pgdump.sql
")

# load into database
# you may have to type in password
system("
psql --host=phi.cluster.ucsf.bkslab.org --port=PORT --username=momeara --password momeara < dump/chembl_21_postgresql/chembl_21.pgdump.sql
")

# delete dump
system("
rm -rf dump
")

# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(plyr)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(magrittr)
library(BioChemPantry)


pantry <- BioChemPantry::get_pantry()
pantry %>% BioChemPantry::create_schema('uniprot_2017_2')
pantry %>% BioChemPantry::set_schema('uniprot_2017_2')
staging_directory <- BioChemPantry::get_staging_directory("uniprot_2017_2")


system(paste0("
mkdir -p ", staging_directory, "/dump
mkdir -p ", staging_directory, "/data
cd ", staging_directory, "/dump

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Reference_Proteomes_2017_12.tar.gz
tar -xzvf Reference_Proteomes_2017_12.tar.gz
rm -rf Reference_Proteomes_2017_12.tar.gz
"))




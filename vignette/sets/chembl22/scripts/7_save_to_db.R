# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(BioChemPantry)

staging_directory <- get_staging_directory("chembl21")
target_info <- read_tsv(paste0(staging_directory, "/data/target_info.tsv"))
load(paste0(staging_directory, "/data/chembl_compounds.Rdata"))
load(paste0(staging_directory, "/data/chembl_activities.Rdata"))
scores <- read_csv(paste0(staging_directory, "/data/chembl21.scores.csv"))

pantry <- get_pantry()
pantry %>% create_schema('sea_chembl21')
pantry %>% set_schema('sea_chembl21')

pantry %>%
	copy_to(
		target_info,
		"targets",
		temporary=F,
		indices=list("uniprot_accn", "uniprot_entry", "gene_name", "entrez_id"))

z <- pantry %>%
	copy_to(
		chembl_compounds,
		"compounds",
		temporary=F,
		indices=list("zinc_id", "chembl_id"),
		fast=T)

z <- pantry %>%
	copy_to(
		chembl_activities,
		"activities",
		temporary=F,
		indices=list(
			c("uniprot_entry", "chembl_id")),
		fast=T)

pantry %>%
	copy_to(
		scores,
		"scores",
		temporary=F,
		indices=list(
			"target1",
			"target2",
			c("target1", "target2")),
		fast=T)


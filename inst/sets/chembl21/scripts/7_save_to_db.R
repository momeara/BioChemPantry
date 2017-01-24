# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
source("~/work/sea/scripts/data_repo.R")


target_info <- read_tsv("data/target_info.tsv")
load("data/chembl_compounds.Rdata")
load("data/chembl_activities.Rdata")
scores <- read_csv("data/chembl21.scores.csv")

data_repo <- get_data_repo()
data_repo %>% create_schema('sea_chembl21')
data_repo %>% set_schema('sea_chembl21')

data_repo %>%
	copy_to(
		target_info,
		"targets",
		temporary=F,
		indices=list("uniprot_accn", "uniprot_entry", "gene_name", "entrez_id"))

z <- data_repo %>%
	copy_to(
		chembl_compounds,
		"compounds",
		temporary=F,
		indices=list("zinc_id", "chembl_id"),
		fast=T)

z <- data_repo %>%
	copy_to(
		chembl_activities,
		"activities",
		temporary=F,
		indices=list(
			c("uniprot_entry", "chembl_id")),
		fast=T)

data_repo %>%
	copy_to(
		scores,
		"scores",
		temporary=F,
		indices=list(
			"target1",
			"target2",
			c("target1", "target2")),
		fast=T)


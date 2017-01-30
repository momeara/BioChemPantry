# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(BioChemPantry)
source("~/work/sets/uniprot_idmapping/scripts/uniprot_id_map.R")



pantry <- get_pantry("pfam30")
staging_directory <- get_staging_directory("pfam30/dump")

system(paste0("
cd ", staging_directory, "
#
# Don't load at this time
#wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.full.uniprot.gz
#gunzip Pfam-A.full.uniprot.gz
#
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.clans.tsv.gz
gunzip Pfam-A.clans.tsv.gz
#
wget wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.regions.tsv.gz
gunzip Pfam-A.regions.tsv.gz
"))

pfam_regions <- readr::read_tsv(paste0(staging_directory, "/Pfam-A.regions.tsv")) %>%
	dplyr::select(
		uniprot_accn = pfamseq_acc,
		pfamA_acc,
		seq_start,
		seq_end)

pfam_clans <- read_tsv(
	paste0(staging_directory, "/Pfam-A.clans.tsv"),
	col_names=c("pfamA_acc", "clan_acc", "pfamA_alt_id", "pfamA_id", "pfamA_description"))



pantry %>% copy_to(
	pfam_regions,
	"regions",
	temporary=F,
	indices=list(
		"uniprot_accn"),
	fast=T)

pantry %>% copy_to(
	pfam_clans,
	"clans",
	temporary=F,
	indices=list(
			"pfamA_acc",
			"pfamA_id"),
	fast=T)

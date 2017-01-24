# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(plyr)
library(dplyr)
library(readr)
library(stringr)

source("~/work/sets/uniprot_idmapping/scripts/uniprot_id_map.R")
source("~/work/sea/scripts/data_repo.R")


data_repo <- get_data_repo()
data_repo %>% create_schema("pfam30")
data_repo %>% set_schema("pfam30")



system("
mkdir dump
cd dump
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
")

pfam_regions <- readr::read_tsv("dump/Pfam-A.regions.tsv") %>%
	dplyr::select(
		uniprot_accn = pfamseq_acc,
		pfamA_acc,
		seq_start,
		seq_end)

pfam_clans <- read_tsv(
	"dump/Pfam-A.clans.tsv",
	col_names=c("pfamA_acc", "clan_acc", "pfamA_alt_id", "pfamA_id", "pfamA_description"))



data_repo %>% copy_to(
	pfam_regions,
	"regions",
	temporary=F,
	indices=list(
		"uniprot_accn"),
	fast=T)

data_repo %>% copy_to(
	pfam_clans,
	"clans",
	temporary=F,
	indices=list(
			"pfamA_acc",
			"pfamA_id"),
	fast=T)

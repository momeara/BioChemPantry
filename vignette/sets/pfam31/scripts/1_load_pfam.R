# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(plyr)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(magrittr)
library(ape)
library(httr)
library(data.tree)
library(BioChemPantry)

pfam_pantry <- BioChemPantry::get_pantry()
#pfam_pantry %>% BioChemPantry::create_schema('pfam31')
pfam_pantry %>% BioChemPantry::set_schema('pfam31')
staging_directory <- BioChemPantry::get_staging_directory("pfam31")

system(paste0("
mkdir -p ", staging_director, "/dump
mkdir -p ", staging_directory, "/data
cd ", staging_directory, "/dump
#
# Don't load at this time
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.full.uniprot.gz
gunzip Pfam-A.full.uniprot.gz
#
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.clans.tsv.gz
gunzip Pfam-A.clans.tsv.gz
#
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.regions.tsv.gz
gunzip Pfam-A.regions.tsv.gz
#
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/trees.tgz
tar -xzf trees.tgz
mv nfs/production/xfam/pfam/sequpd/31.0/Release/31.0/trees trees
rm -rf nfs
rm -rf trees.tgz
"))

uniprot_target_columns <- paste(
			"id",
			"entry name",
			"organism",
			"organism-id")

uniprot_targets <- httr::GET(
	url="http://www.uniprot.org/uniprot/",
	httr::user_agent("httr mattjomeara@gmail.com"),
	query=list(
		query="database:(type:pfam)",
		format='tab',
		compress='yes',
		columns=uniprot_target_columns)) %>%
	httr::content() %>%
	rawConnection() %>%
	gzcon() %>%
	readr::read_tsv()


pfam_regions <- readr::read_tsv(
	paste0(staging_directory, "/dump/Pfam-A.regions.tsv"),
	col_types=cols(
	  pfamseq_acc = col_character(),
	  seq_version = col_double(),
	  crc64 = col_character(),
	  md5 = col_character(),
	  pfamA_acc = col_character(),
	  seq_start = col_double(),
	  seq_end = col_double())) %>%
	dplyr::select(
		uniprot_accn = pfamseq_acc,
		pfamA_acc,
		seq_start,
		seq_end)
save("pfam_regions", file=paste0(staging_directory, "/data/pfam_regions.Rdata"))
#load(paste0(staging_directory, "/data/pfam_regions.Rdata"))

pfam_clans <- read_tsv(
	paste0(staging_directory, "/dump/Pfam-A.clans.tsv"),
	col_names=c("pfamA_acc", "clan_acc", "pfamA_alt_id", "pfamA_id", "pfamA_description"),
	col_types=cols(
	  pfamA_acc = col_character(),
	  clan_acc = col_character(),
	  pfamA_alt_id = col_character(),
	  pfamA_id = col_character(),
	  pfamA_description = col_character()))
save("pfam_clans", file=paste0(staging_directory, "/data/pfam_clans.Rdata"))
#load(paste0(staging_directory, "/data/pfam_clans.Rdata"))

pfam_trees <- list.files(
	path=paste0(staging_directory, "/dump/trees"),
	pattern="*.tree",
	all.files=TRUE,
	full.name=TRUE) %>%
	plyr::ldply(function(tree_fname){
		# cat("reading tree ", phylo_fname, " ...\n", sep="")
		pfamA_acc <- tree_fname %>%
			stringr::str_match("(PF[0-9]+)[.]tree$") %>%
			magrittr::extract(,2)
		newick <- readr::read_file(tree_fname)
		tree <- ape::read.tree(tree_fname)
		tibble::data_frame(pfamA_acc=pfamA_acc, newick=newick, tree=list(tree))
	})
save("pfam_trees", file=paste0(staging_directory, "/data/pfam_trees.Rdata"))
#load(paste0(staging_directory, "/data/pfam_trees.Rdata"))

pfam_pantry %>% copy_to(
	pfam_regions,
	"regions",
	temporary=FALSE,
	indices=list(
		"uniprot_accn",
		"pfamA_acc",
		"pfamA_id"),
	analyze=TRUE)

pfam_pantry %>% copy_to(
	pfam_clans,
	"clans",
	temporary=FALSE,
	indices=list(
			"pfamA_acc",
			"pfamA_id"),
	analyze=TRUE)

pfam_pantry %>% copy_to(
	pfam_trees %>%
		dplyr::select(pfamA_acc, newick),
	"trees",
	temporary=FALSE,
	indices=list("pfamA_acc"),
	analyze=TRUE)

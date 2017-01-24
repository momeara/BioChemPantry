# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(seqinr)
library(stringr)
library(Bethany)
source("~/work/sets/uniprot_idmapping/scripts/uniprot_id_map.R")
source("~/work/sea/scripts/data_repo.R")

data_repo <- get_data_repo("hgnc_151123")

genes <- data_repo %>% tbl("genes") %>% collect
hgnc_human_fname <- "data/uniprot_151123/hgnc_human.fasta"


scores <- blastp(hgnc_human_fname, hgnc_human_fname, "hgnc_human")

scores <- scores %>%
	dplyr::rename(
		uniprot_accn1 = ref_target,
		uniprot_accn1 = query_target) %>%
	dplyr::arrange(EValue) %>%
	dplyr::distinct(target1, target2)

scores <- scores %>%
	rename(
		uniprot_accn1 = target1,
		uniprot_accn2 = target2) %>%
	left_join(
		genes %>%
			select(
				symbol1 = symbol,
				uniprot_accn1 = uniprot_accn,
				uniprot_entry1 = uniprot_entry),
		by="uniprot_accn1") %>%
	left_join(
		genes %>%
			select(
				symbol2 = symbol,
				uniprot_accn2 = uniprot_accn,
				uniprot_entry2 = uniprot_entry),
		by="uniprot_accn2")

data_repo %>% copy_to(
	df=scores,
	name="blastp_target_vs_target",
	temporary=F,
	indexes=list(
		"symbol1",
		"uniprot_accn1",
		"uniprot_entry1",
		"symbol2",
		"uniprot_accn2",
		"uniprot_entry2",
		c("symbol1", "symbol2"),
		c("uniprot_accn1", "uniprot_accn2"),
		c("uniprot_entry1", "uniprot_entry2")),
	fast=T)

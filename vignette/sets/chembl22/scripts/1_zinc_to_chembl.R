# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(assertr)
library(Zr)

source("~/work/sea/scripts/zinc_tools.R")

staging_directory <- get_staging_directory("chembl21")

zinc_to_chembl <- catalog_items(
	catalog_short_name="chembl21",
	output_fields=c(
		"zinc_id",
		"supplier_code",
		"substance.preferred_name",
		"substance.smiles",
		"substance.purchasable",
		"substance.purchasability",
#		"substance.gene_names",
		"substance.rb",
		"substance.reactive",
		"substance.features"),
	count='all',
	result_batch_size=1,
	page=914931,
	verbose=T)

zinc_to_chembl <- zinc_to_chembl_raw %>%
	dplyr::mutate(
#		# zinc sub.stance.gene_names looks like "[u'GENEID1', u'GENEID2', ...]"
#		# split this into an R list
#		gene_names = substance.gene_names %>%
#			stringr::str_replace("^[\\[]", "") %>%
#			stringr::str_replace("[\\]]$", "") %>%
#			stringr::str_replace("^u'", "") %>%
#			stringr::str_replace("'$", "") %>%
#			strsplit("', u'", fixed=T)) %>%
	dplyr::mutate(
		smiles = substance.smiles,
		preferred_name = substance.preferred_name,
		n_genes = gene_names %>% vapply(length, 1L),
		rotatable_bonds = as.numeric(substance.rb),
		reactivity = as.numeric(substance.reactive)) %>%
	dplyr::mutate(
#		gene_names = gene_names %>% llply(function(x) paste(x, collapse=";")) %>% unlist,
		aggregator = !is.na(substance.features) & str_detect(substance.features, 'aggregator'),
		purchasable_code = as.numeric(substance.purchasable),
		purchasable_level = substance.purchasability,
		drug_code =
			ifelse(str_detect(substance.features, 'fda'), 210,
			ifelse(str_detect(substance.features, 'world'), 211,
			ifelse(str_detect(substance.features, 'investigationa'), 212,
			ifelse(str_detect(substance.features, 'in-man'), 213,
			ifelse(str_detect(substance.features, 'in-vivo'), 213,
			ifelse(str_detect(substance.features, 'in-cells'), 215,
			ifelse(str_detect(substance.features, 'in_vitro'), 216,
			NA))))))),
		drug_level =
			ifelse(str_detect(substance.features, 'fda'), 'fda',
			ifelse(str_detect(substance.features, 'world'), 'world',
			ifelse(str_detect(substance.features, 'investigationa'), 'investigational',
			ifelse(str_detect(substance.features, 'in-man'), 'in-man',
			ifelse(str_detect(substance.features, 'in-vivo'), 'in-vivo',
			ifelse(str_detect(substance.features, 'in-cells'), 'in-cells',
			ifelse(str_detect(substance.features, 'in_vitro'), 'in-vitro',
			NA))))))),
		biological_code =
			ifelse(str_detect(substance.features, 'endogenous'), 201,
			ifelse(str_detect(substance.features, 'metabolite'), 202,
			ifelse(str_detect(substance.features, 'biogenic'), 203,
			NA))),
		biological_level =
			ifelse(str_detect(substance.features, 'endogenous'), 'endogenous',
			ifelse(str_detect(substance.features, 'metabolite'), 'metabolite',
			ifelse(str_detect(substance.features, 'biogenic'), 'biogenenic',
			NA)))) %>%
	dplyr::select(
		chembl_id = supplier_code,
		zinc_id,
		preferred_name = preferred_name,
		smiles,
		rotatable_bonds,
		reactivity,
#		gene_names,
#		n_genes,
		aggregator,
		purchasable_code,
		purchasable_level,
		drug_code,
		drug_level,
		biological_code,
		biological_level)


zinc_to_chembl <- zinc_to_chembl %>%
	dplyr::rename(chembl_id = supplier_code)

zinc_to_chembl %>%
	assertr::verify(zinc_id %>% str_detect("^ZINC[0-9]{12}$")) %>%
	assertr::verify(chembl_id %>% str_detect("^CHEMBL[0-9]{1,7}$"))

zinc_to_chembl %>%
	readr::write_tsv(paste0(staging_directory, "/data/zinc_to_chembl.tsv"))


aggregators <- catalog_items(
	catalog_short_name="aggregators",
	output_fields=c("zinc_id", "smiles"),
	count='all',
	result_batch_size=1000,
	verbose=T)

aggregators %>%
	readr::write_tsv(paste0(staging_directory, "/data/aggregators_160703.tsv"))


## chembl compounds similar ot aggregators
chembl_to_aggregators <- query_compounds_near_reference(
	aggregators,
	compounds,
	cutoff=.7,
	ref_compound="chembl_id",
	ref_smiles="chembl_smiles",
	query_compound="zinc_id",
	query_smiles="smiles",
	verbose=T)

chembl_to_aggregators <- chembl_to_aggregators %>%
	dplyr::select(
		chembl_id = compound,
		MaxTC_aggregator = MaxTC)

chembl_to_aggregators %>%
	readr::write_tsv(paste0(staging_directory, "/data/chembl_to_aggregators_160704.tsv"))




#
#zinc_to_chembl <- zinc_to_chembl_raw %>%
#	dplyr::mutate(
#		# zinc sub.stance.gene_names looks like "[u'GENEID1', u'GENEID2', ...]"
#		# split this into an R list
#		gene_names = substance.gene_names %>%
#			stringr::str_replace("^[\\[]", "") %>%
#			stringr::str_replace("[\\]]$", "") %>%
#			stringr::str_replace("^u'", "") %>%
#			stringr::str_replace("'$", "") %>%
#			strsplit("', u'", fixed=T)) %>%
#	dplyr::mutate(
#		smiles = substance.smiles,
#		preferred_name = substance.preferred_name,
#		n_genes = gene_names %>% vapply(length, 1L),
#		rotatable_bonds = as.numeric(substance.rb),
#		reactivity = as.numeric(substance.reactive)) %>%
#	dplyr::mutate(
#		gene_names = gene_names %>% llply(function(x) paste(x, collapse=";")) %>% unlist,
#		aggregator = !is.na(substance.features) & str_detect(substance.features, 'aggregator'),
#		purchasable_code = as.numeric(substance.purchasable),
#		purchasable_level = substance.purchasability,
#		drug_code =
#			ifelse(str_detect(substance.features, 'fda'), 210,
#			ifelse(str_detect(substance.features, 'world'), 211,
#			ifelse(str_detect(substance.features, 'investigationa'), 212,
#			ifelse(str_detect(substance.features, 'in-man'), 213,
#			ifelse(str_detect(substance.features, 'in-vivo'), 213,
#			ifelse(str_detect(substance.features, 'in-cells'), 215,
#			ifelse(str_detect(substance.features, 'in_vitro'), 216,
#			NA))))))),
#		drug_level =
#			ifelse(str_detect(substance.features, 'fda'), 'fda',
#			ifelse(str_detect(substance.features, 'world'), 'world',
#			ifelse(str_detect(substance.features, 'investigationa'), 'investigational',
#			ifelse(str_detect(substance.features, 'in-man'), 'in-man',
#			ifelse(str_detect(substance.features, 'in-vivo'), 'in-vivo',
#			ifelse(str_detect(substance.features, 'in-cells'), 'in-cells',
#			ifelse(str_detect(substance.features, 'in_vitro'), 'in-vitro',
#			NA))))))),
#		biological_code =
#			ifelse(str_detect(substance.features, 'endogenous'), 201,
#			ifelse(str_detect(substance.features, 'metabolite'), 202,
#			ifelse(str_detect(substance.features, 'biogenic'), 203,
#			NA))),
#		biological_level =
#			ifelse(str_detect(substance.features, 'endogenous'), 'endogenous',
#			ifelse(str_detect(substance.features, 'metabolite'), 'metabolite',
#			ifelse(str_detect(substance.features, 'biogenic'), 'biogenenic',
#			NA)))) %>%
#	dplyr::select(
#		chembl_id = supplier_code,
#		zinc_id,
#		preferred_name = preferred_name,
#		smiles,
#		rotatable_bonds,
#		reactivity,
#		gene_names,
#		n_genes,
#		aggregator,
#		purchasable_code,
#		purchasable_level,
#		drug_code,
#		drug_level,
#		biological_code,
#		biological_level)
#


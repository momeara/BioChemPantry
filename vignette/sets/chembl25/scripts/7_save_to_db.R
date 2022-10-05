# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(BioChemPantry)

staging_directory <- get_staging_directory("chembl25")
load(paste0(staging_directory, "/data/full_chembl_data.Rdata"))
target_info <- read_tsv(paste0(staging_directory, "/data/target_info.tsv"))
scores <- readr:read_csv(paste0(staging_directory, "/data/chembl25_rdkit_ecfp4.scores.csv"))
load(paste0(staging_directory, "/data/chembl_compound_images.Rdata"))

load(paste0(staging_directory, "/data/chembl_compounds.Rdata"))
load(paste0(staging_directory, "/data/chembl_activities.Rdata"))




pantry <- get_pantry('sea_chembl25')
#pantry %>% create_schema('sea_chembl25')
#pantry %>% set_schema('sea_chembl25')

##########################################33
pantry %>%
	copy_to(
		target_info,
		"targets",
		temporary=F,
		indices=list("uniprot_accn", "uniprot_entry", "gene_name", "entrez_id"))


##########################################
compounds <- full_chembl_data %>%
	dplyr::distinct(zinc_id, .keep_all=T) %>%
	dplyr::select(
		zinc_id,
		chembl_id,
		smiles,
		chembl_smiles,
		preferred_name,
		molecular_weight,
		ALogP,
		rotatable_bonds,
		reactivity,
		PSA,
		HBA,
		HBD,
		most_basic_pKa,
		most_acidic_pKa,
		molecular_species,
		n_heavy_atoms,
		aggregator,
		aggregator_zinc_id,
		tc_to_aggregator,
		n_genes,
		purchasable_code,
		purchasable_level,
		drug_code,
		drug_level,
		biological_code,
		biological_level)


compound_n_test_targets <- full_chembl_data %>%
	dplyr::distinct(uniprot_entry, zinc_id, .keep_all) %>%
	dplyr::count(zinc_id) %>%
	dplyr::transmute(zinc_id, n_tested_targets=n)


compound_n_potent_targets <- full_chembl_data %>%
	dplyr::arrange(standard_value) %>%
	dplyr::distinct(uniprot_entry, zinc_id, .keep_all=T) %>%
	dplyr::group_by(zinc_id) %>%
	dplyr::summarize(n_potent_targets=sum(standard_value < 10000))

compound_n_active_targets <- full_chembl_data %>%
	dplyr::arrange(activity) %>%
	dplyr::distinct(uniprot_entry, zinc_id, .keep_all=T) %>%
	dplyr::group_by(zinc_id) %>%
	dplyr::summarize(n_active_targets=sum(active))

compounds <- compounds %>%
#	dplyr::left_join(compound_n_test_targets, by="zinc_id") %>%
#	dplyr::left_join(compound_n_potent_targets, by="zinc_id") %>%
	dplyr::left_join(compound_n_active_targets, by="zinc_id")


# for some reason I'm getting errors when I try to insert the whole table at once
# it does seem to work to insert the first 200k rows and the rest of it
z <- pantry %>%
	copy_to(
		compounds %>% head(200000),
		"compounds",
		temporary=F,
		indices=list("zinc_id", "chembl_id", "smiles"),
		fast=T)
dplyr::db_begin(pantry$con)
	dplyr::db_insert_into(pantry$con, "compounds", compounds %>% slice(200001:474591))
	dplyr::db_analyze(pantry$con, "compounds")
dplyr::db_commit(pantry$con)

#####################################################
full_chembl_data <- full_chembl_data %>%
	dplyr::mutate(gene_names = vapply(gene_names, function(g) paste0(g, collapse=";"), "A"))

z <- pantry %>%
	copy_to(
		full_chembl_data %>% head(100000),
		"activities",
		temporary=F,
		indices=list(
			c("uniprot_entry", "chembl_id")),
		fast=T)
dplyr::db_begin(pantry$con)
	dplyr::db_insert_into(pantry$con, "activities", full_chembl_data %>% slice(100001:828327))
	dplyr::db_analyze(pantry$con, "activities")
dplyr::db_commit(pantry$con)

###################################################333
pantry %>%
	copy_to(
		scores %>% head(300000),
		"scores_rdkit_ecfp4",
		temporary=F,
		indices=list(
			"target1",
			"target2",
			c("target1", "target2")),
		fast=T)
dplyr::db_begin(pantry$con)
	dplyr::db_insert_into(pantry$con, "scores_rdkit_ecfp4", scores %>% slice(300001:653579))
	dplyr::db_analyze(pantry$con, "scores_rdkit_ecfp4")
dplyr::db_commit(pantry$con)




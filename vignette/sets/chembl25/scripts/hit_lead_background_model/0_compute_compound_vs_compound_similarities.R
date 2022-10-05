# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(DBI)
library(RSQLite)
library(dbplyr)


library(BioChemPantry)
library(SEAR)



staging_directory <- BioChemPantry::get_staging_directory("chembl23")


load(paste0(staging_directory, "/data/full_chembl_data.Rdata"))

chembl_smi <- full_chembl_data %>%
	dplyr::filter(active) %>%
	dplyr::distinct(zinc_id, .keep_all=T) %>%
	dplyr::mutate(compound_id = row_number()) %>%
	dplyr::select(zinc_id, smiles, compound_id)



#tcs <- SEAR::tc_matrix(
#	ref_fp=data.frame(
#		smiles=c(
#			"c1cc(c(cc1CC[NH3+])O)O",
#			"c1cc2c(cc1O)c(c[nH]2)CC[NH3+]",
#			"C[NH2+]C[C@@H](c1ccc(c(c1)O)O)O"),
#		compound=c(
#			"dopamine",
#			"serotonin",
#			"adrenaline")),
#	query_fp=data.frame(
#		smiles=c(
#			"COc1cc(ccc1O)CC(=O)[O-]",
#			"c1cc(c(cc1C[C@@H](C(=O)[O-])[NH3+])O)O"),
#		compound=c(
#			"homovanillic acid",
#			"levodopa")))


# TODO use indices rather than ZINC ids to save on disk io
tcs <- SEAR::tc_matrix(
	ref_fp = chembl_smi,
	query_fp = chembl_smi,
	cutoff = .28,
	ref_compound = "zinc_id",
	ref_smiles = "smiles",
	query_compound = "zinc_id",
	query_smiles = "smiles",
	verbose=TRUE)

tcs <- tcs %>%
	dplyr::rename(
		ref_zinc_id=ref_cid,
		query_zinc_id=query_cid) %>%
	dplyr::left_join(
		chembl_smi %>%
			dplyr::select(ref_zinc_id=zinc_id, ref_cid=compound_id),
		by=c("ref_zinc_id")) %>%
	dplyr::left_join(
		chembl_smi %>%
			dplyr::select(query_zinc_id=zinc_id, query_cid=compound_id),
		by=c("query_zinc_id")) %>%
	dplyr::arrange(ref_cid) %>%
	dplyr::select(ref_cid, query_cid, tc)

save("tcs", file=paste0(staging_directory, "/data/tcs_rdkit_ecfp4.Rdata"))


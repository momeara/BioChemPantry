# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(assertr)
library(Zr)
library(BioChemPantry)
library(SEAR)

staging_directory <- get_staging_directory("chembl27")

# don't just get all the info for the catalog items because it times
# out instead get the basic info and then look up more detailed data
# in batches
zinc_to_chembl <- Zr::catalog_items(
	catalog_short_name="chembl27",
	count='all',
	verbose=T) %>%
	dplyr::select(
		zinc_id,
		chembl_id = supplier_code)

zinc_to_chembl <- readr::read_delim(paste0(staging_directory, "/dump/zinc_to_chembl.txt"), delim=' ')

zinc_to_chembl_info <- Zr::substance_info(
	zinc_ids=zinc_to_chembl$zinc_id,
	output_fields=c(
		"zinc_id",
		"preferred_name",
		"smiles",
		"purchasable",
		"purchasability",
		"gene_names",
		"rb",
		"reactive",
		"features"),
	batch_size=10000,
	verbose=T)
zinc_to_chembl_info <- readr::read_tsv(paste0(staging_directory, "/dump/chembl27_to_zinc.tsv"))

processed_zinc_to_chembl_info <- zinc_to_chembl_info %>%
	Zr::process_substance_info() %>%
	dplyr::select(
		zinc_id,
		preferred_name,
		smiles,
		rotatable_bonds=rb,
		reactivity=reactive,
		gene_names,
		n_genes,
		aggregator,
		purchasable_code=purchasable,
		purchasable_level=purchasability,
		drug_code,
		drug_level,
		biological_code,
		biological_level)

zinc_to_chembl_info %>%
	readr::write_tsv(
		paste0(staging_directory, "/dump/chembl27_to_zinc_info.tsv"))

missing_info <- zinc_to_chembl %>%
	dplyr::anti_join(zinc_to_chembl_info)
missing_info %>% glimpse

zinc_to_chembl <- zinc_to_chembl %>%
	dplyr::inner_join(zinc_to_chembl_info, by="zinc_id")

zinc_to_chembl %>%
	assertr::verify(zinc_id %>% str_detect("^ZINC[0-9]{12}$")) %>%
	assertr::verify(chembl_id %>% str_detect("^CHEMBL[0-9]{1,7}$"))

BioChemPantry::summarize_map(
	zinc_to_chembl,
	"zinc_id",
	"smiles")
#X<-[zinc_id]:
#  |X|: 2062965
#  count*size: 2062965*1
#Y<-[smiles]:
#  |Y|: 2062881
#  count*size: 2062797*1, 84*2
#[X U Y]:
#  |X U Y|: 2062965
#  count*size: 1924578*1, 565*2, 132093*4, 43*6, 3645*9, 28*12, ... 2*324, 1*361, 1*462, 1*528, 1*552, 1*675
#[X @ Y]:
#  |X ~ Y|: 2062797
#  |X:X < Y|, |Y:X < Y|: 0, 0
#  |X:X > Y|, |Y:X > Y|: 168, 84


save(zinc_to_chembl, file=paste0(staging_directory, "/data/zinc_to_chembl.Rdata"))

aggregators <- Zr::catalog_items(
	catalog_short_name="aggregators",
	output_fields=c("zinc_id", "smiles"),
	count='all',
	result_batch_size=1000,
	verbose=T) %>%
	dplyr::distinct("zinc_id", "smiles")

aggregators %>%
	readr::write_tsv(paste0(staging_directory, "/data/aggregators_170818.tsv"))


## chembl compounds similar ot aggregators
chembl_to_aggregators <- tc_matrix(
	ref_fp=aggregators,
	query_fp=zinc_to_chembl %>%
		dplyr::distinct(zinc_id, .keep_all=TRUE),
	cutoff=.7,
	ref_compound="zinc_id",
	ref_smiles="smiles",
	query_compound="zinc_id",
	query_smiles="smiles",
	fp_format='sea_native',
	fp_type='rdkit_ecfp',
	verbose=T) %>%
	rename(aggregator_zinc_id = ref, chembl_zinc_id = query)

chembl_to_aggregators %>%
	readr::write_tsv(paste0(staging_directory, "/data/chembl_to_aggregators_170813.tsv"))




library(dplyr)
library(readr)
library(stringr)
library(assertr)
library(Zr)
library(BioChemPantry)
library(SEAR)

source("scripts/parameters.R")
staging_directory <- BioChemPantry::get_staging_directory(schema)



# don't just get all the info for the catalog items because it times
# out instead get the basic info and then look up more detailed data
# in batches
zinc_to_chembl <- Zr::catalog_items(
	catalog_short_name = "chembl30",
	page = 1,
	result_batch_size = 3000,
	output_fields = c("zinc_id", "supplier_code"),
	verbose = TRUE) |>
	dplyr::select(
		zinc_id,
		chembl_id = supplier_code)

## For me this took lots or re-starts so, load from tempfiles:
zinc_to_chembl <- data.frame(
  page_fname = list.files(
    path = staging_directory,
    pattern = "*_zinc_query_[0-9]+.csv",
    full.names = TRUE)) |>
  dplyr::rowwise() |>
  dplyr::do({
    page_fname <- .$page_fname[1]
    readr::read_csv(
      file = page_fname,
      show_col_types = FALSE)
  }) |>
  dplyr::ungroup() |>
  dplyr::select(
    zinc_id,
    chembl_id = supplier_code)

zinc_to_chembl |> 
  readr::write_tsv(
    file = paste0(staging_directory, "/dump/zinc_to_chembl.tsv"))


n_batches <- 5000
zinc_to_chembl_info <- zinc_to_chembl |>
  dplyr::mutate(batch_index = dplyr::row_number() %% n_batches) |>
  dplyr::group_by(batch_index) |>
  dplyr::do({
    batch <- .
    cat(
      "getting batch '", .$batch_index[1], ", ",
      "nrow:", nrow(batch), "\n", sep = "")
    tictoc::tic()
    info <- Zr::substance_info(
	    zinc_ids = batch$zinc_id,
	    output_fields = c(
		    "zinc_id",
		    "preferred_name",
		    #"smiles",
		    "purchasable",
		    "purchasability",
		    #"gene_names",
		    "rb",
		    "reactive"),
		    #"features"),
	    verbose = TRUE)
    tictoc::toc()
    info |> readr::write_tsv(
      paste0(
        staging_directory,
        "/zinc_to_chembl_info_", batch, "_of_", batch_size, ".tsv"))
    info
  }) |>
  dplyr::ungroup()

zinc_to_chembl_info <- zinc_to_chembl_info |>
	Zr::process_substance_info() |>
	dplyr::select(
		zinc_id,
		preferred_name,
		smiles,
		rotatable_bonds,
		reactivity,
		gene_names,
		n_genes,
		aggregator,
		purchasable_code,
		purchasable_level,
		drug_code,
		drug_level,
		biological_code,
		biological_level)

zinc_to_chembl_info |>
	readr::write_tsv(
		paste0(staging_directory, "/dump/chembl33_to_zinc_info.tsv"))

missing_info <- zinc_to_chembl |>
	dplyr::anti_join(zinc_to_chembl_info)
missing_info |> glimpse()

zinc_to_chembl <- zinc_to_chembl |>
	dplyr::inner_join(zinc_to_chembl_info, by = "zinc_id")

zinc_to_chembl |>
	assertr::verify(
	  zinc_id |>
	    stringr::str_detect("^ZINC[0-9]{12}$")) |>
	assertr::verify(
	  chembl_id |>
	    stringr::str_detect("^CHEMBL[0-9]{1,7}$"))

BioChemPantry::summarize_map(
	zinc_to_chembl,
	"zinc_id",
	"smiles")

save(
  zinc_to_chembl,
  file=paste0(staging_directory, "/data/zinc_to_chembl.Rdata"))

aggregators <- Zr::catalog_items(
	catalog_short_name = "aggregators",
	output_fields = c("zinc_id", "smiles"),
	count='all',
	result_batch_size = 1000,
	verbose = TRUE) |>
	dplyr::distinct("zinc_id", "smiles")

aggregators |>
	readr::write_tsv(
	  paste0(
	    staging_directory, "/data/aggregators_",
	    BioChemPantry::date_code(), ".tsv"))


## chembl compounds similar to aggregators
chembl_to_aggregators <- SEAR::tc_matrix(
	ref_fp = aggregators,
	query_fp = zinc_to_chembl |>
		dplyr::distinct(zinc_id, .keep_all=TRUE),
	cutoff = .7,
	ref_compound = "zinc_id",
	ref_smiles = "smiles",
	query_compound = "zinc_id",
	query_smiles = "smiles",
	fp_format = 'sea_native',
	fp_type = 'rdkit_ecfp',
	verbose = TRUE) |>
	dplyr::rename(aggregator_zinc_id = ref, chembl_zinc_id = query)

chembl_to_aggregators |>
	readr::write_tsv(
	  paste0(
	    staging_directory, "/data/chembl_to_aggregators_",
	    BioChemPantry::date_code(), ".tsv"))


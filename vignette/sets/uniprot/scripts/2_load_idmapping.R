# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(plyr)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(magrittr)
library(BioChemPantry)

schema <- "uniprot_2017_2"
pantry <- BioChemPantry::get_pantry(schema)
staging_directory <- BioChemPantry::get_staging_directory(schema)

id_types <- tibble::tribble(
	~id_type, ~idmapping_type,
	"uniprot_entry", "UniProtKB-ID",
	"gene_id", "GeneID",
	"hgnc", "HGNC",
	"tax_id", "NCBI_TaxID") %>%
	dplyr::mutate(
		table_name=paste0("to_", id_type))

# rebuild tables
pantry %>% dplyr::db_begin()
tryCatch({
	pantry %>%
		BioChemPantry::get_tables() %>%
		dplyr::filter(tablename %>% stringr::str_detect("^to_")) %>%
		plyr::a_ply(1, function(df){
			pantry %>% BioChemPantry::drop_table(df$tablename[1])
		})

	id_types %>%
		plyr::a_ply(1, function(df){
				dplyr::db_create_table(
					con=pantry,
					table=df$table_name[1],
					types=c(uniprot_accn="TEXT", id_value="TEXT"),
					temporary=FALSE)
	})
	pantry %>% dplyr::db_commit()
}, error = function(err) {
	cat("Transaction failed: ", err$message, "\n", sep="")
	pantry %>% dplyr::db_rollback()
})



kingdom <- "Eukaryota"
ancestor <- 7711 # chordate
taxa <- BioChemPantry::uniprot_taxa(ancestor=ancestor)

idmapping_files <- list.files(
	path=paste0(staging_directory, "/dump/", kingdom),
	pattern=".*idmapping.gz$",
	full.names=TRUE) %>%
	tibble::data_frame(fname=.) %>%
	dplyr::mutate(taxon = fname %>%
		stringr::str_match("UP[0-9]+_([0-9]+)[.]idmapping.gz") %>%
		magrittr::extract(,2) %>%
		as.numeric()) %>%
	dplyr::left_join(taxa, by="taxon")


pantry %>% dplyr:::db_begin()
tryCatch({
	idmapping_files %>%
		plyr::a_ply(1,
			function(df) {
				cat("Reading ids for '", df$scientific_name[1], "' taxon='", df$taxon[1], "' ...\n", sep="")
				idmapping <- readr::read_tsv(
					file=df$fname[1],
					col_names=c(
						"uniprot_accn",
						"idmapping_type",
						"id_value"),
					col_types = readr::cols(
						uniprot_accn = readr::col_character(),
						idmapping_type = readr::col_character(),
						id_value = readr::col_character())) %>%
					dplyr::inner_join(id_types, by="idmapping_type")

				idmapping %>%
					plyr::d_ply(c("id_type"), function(df){
						id_type <- df$id_type[1]
						table_name <- df$table_name[1]
						cat(
							"\tCopying ", nrow(df), " ids for '",
							id_type, "' to pantry table '", schema, ".", table_name, "' ...\n", sep="")

						pantry %>%
							DBI::dbWriteTable(
								name=table_name,
								value=df %>% dplyr::select(uniprot_accn, id_value),
								append=TRUE,
								overwrite=FALSE)
					})
			})

	pantry %>% dplyr::db_commit()
}, error = function(err){
	cat("Transaction failed: ", err$message, "\n", sep="")
	pantry %>% dplyr::db_rollback()
})


id_types %>%
	plyr::a_ply(1, function(df){
#			dplyr::db_create_indexes(
#				con=pantry,
#				table=df$table_name[1],
#				indexes=list("uniprot_accn", "id_value"),
#				unique = FALSE)
			dplyr::db_analyze(
				con=pantry,
				table=df$table_name[1])
	})

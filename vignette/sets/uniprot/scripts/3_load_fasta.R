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


# build tables
pantry %>% dplyr::db_begin()
tryCatch({
	pantry %>% BioChemPantry::drop_table("protein_sequences")
	pantry %>%
		dplyr::db_create_table(
			table="protein_sequences",
			types=c(uniprot_entry="TEXT", sequence="TEXT"),
			temporary=FALSE)
	pantry %>% dplyr::db_commit()
}, error = function(err) {
	cat("failed to create protein_sequences_table: ", err$message, "\n", sep="")
	pantry %>% dplyr::db_rollback()
})


kingdom <- "Eukaryota"
fasta_files <- list.files(
	path=paste0(staging_directory, "/dump/", kingdom),
	pattern="[0-9][.]fasta.gz$",
	full.names=TRUE) %>%
	tibble::data_frame(fname=.) %>%
	dplyr::mutate(taxon = fname %>%
		stringr::str_match("UP[0-9]+_([0-9]+)[.]fasta.gz") %>%
		magrittr::extract(,2) %>%
		as.numeric())

pantry %>% dplyr:::db_begin()
tryCatch({
	fasta_files %>%
		plyr::a_ply(1,
			function(df) {
				cat("Reading fasta for taxon='", df$taxon[1], "' ...\n", sep="")

				sequences <- seqinr::read.fasta(df$fname[1], seqtype="AA")
				sequences <-
					tibble::data_frame(
						uniprot_entry = sequences %>%
							seqinr::getName() %>%
							stringr::str_match("([^|]+)$") %>%
							magrittr::extract(,2),
						sequence = sequences %>%
							seqinr::getSequence() %>%
							lapply(function(seq) paste0(seq, collapse="")) %>%
							unlist())

				cat(
					"\tCopying ", nrow(sequences), " sequences to pantry table ",
					"'", schema, ".protein_sequences' ...\n", sep="")

				pantry %>%
					DBI::dbWriteTable(
						name="protein_sequences",
						value=sequences,
						append=TRUE,
						overwrite=FALSE)
			})

	pantry %>% dplyr::db_commit()
}, error = function(err){
	cat("Transaction failed: ", err$message, "\n", sep="")
	pantry %>% dplyr::db_rollback()
})


pantry %>% dplyr::db_create_indexes(
	table="protein_sequences",
	indexes=list("uniprot_entry"),
	unique=TRUE)

pantry %>% dplyr::db_analyze(
	table="protein_sequences")

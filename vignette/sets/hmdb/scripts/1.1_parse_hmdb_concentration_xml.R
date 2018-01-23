# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(xml2)
library(BioChemPantry)

#####################################
# Analysis parameters

staging_directory <- BioChemPantry::get_staging_directory("hmdb")

# input paths
hmdb_metabolites_xml <- paste0(staging_directory, "/dump/hmdb_metabolites.xml") %>%
	xml2::read_xml()


# concentrations
concentration_data <- hmdb_metabolites_xml %>%
	xml2::xml_children() %>%
	plyr::ldply(function(metabolite, .parallel=TRUE){
			metabolite %>%
				xml_child_by_name("normal_concentrations") %>%
				xml2::xml_children() %>%
				plyr::ldply(function(concentration){
					get_value <- function(name){
						node <- concentration %>% xml_child_by_name(name)
						if(is.na(node)) return(NA)
						contents <- node %>% xml2::xml_contents()
						if(length(contents) == 0) return(NA)
						contents %>% as.character
					}

					tryCatch({
						reference <- concentration %>%
							xml_child_by_name("references") %>%
							xml_child_by_name("reference")
					}, error=function(e){
						reference <<- NA
					})

					if(reference %>% is.na){
							reference_text = NA
							reference_pubmed_id = NA
					} else {
						tryCatch({
							reference_text = reference %>%
								xml_child_by_name("reference_text") %>%
								xml2::xml_contents()
							reference_text <- ifelse(
								length(reference_text) == 0, NA, as.character(reference_text))
						}, error=function(e){
							reference_text <<- NA
						})

						tryCatch({
							reference_pubmed_id = reference %>%
								xml_child_by_name("pubmed_id") %>%
								xml2::xml_contents()
							reference_pubmed_id <- ifelse(
								length(reference_pubmed_id) == 0, NA, as.character(reference_pubmed_id))
						}, error=function(e){
							reference_pubmed_id <<- NA
						})
					}

					data.frame(
						accession = metabolite %>% xml_child_by_name("accession") %>% xml2::xml_contents() %>% as.character,
						biofluid = get_value("biofluid"),
						concentration_value = get_value("concentration_value"),
						concentration_units = get_value("concentration_units"),
						subject_age = get_value("subject_age"),
						subject_sex = get_value("subject_sex"),
						subject_condition = get_value("subject_condition"),
						reference_text = reference_text,
						reference_pubmed_id = reference_pubmed_id)
			})
	})

concentration_data <- concentration_data %>%
	dplyr::mutate(
		concentration_value_sd = concentration_value %>%
			stringr::str_extract("[0-9]+\\.[0-9]+$") %>%
			as.numeric(),
		concentration_value = concentration_value %>%
			stringr::str_extract("^[0-9]+\\.[0-9]*") %>%
			as.numeric(),
		subject_age = subject_age %>%
			stringr::str_replace("&gt;", ">") %>%
			stringr::str_replace("&lt;", "<"),
		subject_condition = subject_condition %>%
			stringr::str_replace("[ ]+$", ""),
		subject_condition = ifelse(
			subject_condition == "normal",
			"Normal", subject_condition)) %>%
		dplyr::filter(!(is.na(concentration_value))

concentration_data %>% readr::write_tsv(
   paste0(staging_directory, "/data/metabolite_concentrations.tsv"))

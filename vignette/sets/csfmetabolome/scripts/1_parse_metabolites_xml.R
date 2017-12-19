# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(xml2)
library(xlsx)
library(BioChemPantry)

#####################################
# Analysis parameters

staging_directory <- BioChemPantry::get_staging_directory("csfmetabolome")

# input paths
metabolites_xml <- paste0(staging_directory, "/dump/csf_metabolites.xml") %>%
	xml2::read_xml()


xml_child_by_name <- function(xml, name){
	child_index <- xml %>%
		xml2::xml_children() %>%
		xml2::xml_name() %>%
		magrittr::equals(name) %>%
		which
	if(length(child_index) == 0) return(NA)
	xml %>%
		xml2::xml_children() %>%
		magrittr::extract2(child_index)
}

m_data <- metabolites_xml %>%
	xml2::xml_children() %>%
	plyr::ldply(function(metabolite){
		get_value <- function(name){
			node <- metabolite %>% xml_child_by_name(name)
			if(is.na(node)) return(NA)
			contents <- node %>% xml2::xml_contents()
			if(length(contents) == 0) return(NA)
			contents %>% as.character
		}
		get_taxon <- function(level){
			taxonomy <- metabolite %>% xml_child_by_name("taxonomy")
			if(is.na(taxonomy)) return(NA)
			node <- taxonomy %>% xml_child_by_name(level)
			if(is.na(node)) return(NA)
			contents <- node %>% xml2::xml_contents()
			if(length(contents) == 0) return(NA)
			contents %>% as.character
		}
		ontology_status <- function(){
			ontology_node <- metabolite %>% xml_child_by_name("ontology")
			if(is.na(ontology_node)) return(NA)
			status <- ontology_node %>% xml_child_by_name("status")
			if(is.na(status)) return(NA)
			contents <- status %>% xml2::xml_contents()
			if(length(contents) == 0) return(NA)
			contents %>% as.character
		}
		from_origin <- function(source){
			ontology_node <- metabolite %>% xml_child_by_name("ontology")
			if(is.na(ontology_node)) return(NA)
			origins_nodeset <- ontology_node %>% xml_child_by_name("origins")
			if(is.na(origins_nodeset) || length(origins_nodeset) == 0) return(NA)
			origins <- origins_nodeset %>% xml2::as_list() %>% purrr::flatten
			source %in% origins
		}

		tibble::data_frame(
			accession = get_value("accession"),
			name = get_value("name"),
			smiles = get_value("smiles"),
			description = get_value("description"),
			taxonomy_description = get_taxon("description"),
			taxonomy_kingdom = get_taxon("kingdom"),
			taxonomy_super_class = get_taxon("super_class"),
			taxonomy_class = get_taxon("class"),
			taxonomy_sub_class = get_taxon("sub_class"),
			taxonomy_molecular_framework = get_taxon("molecular_framework"),
			ontology_status = ontology_status(),
			origin_food = from_origin("Food"),
			origin_endogenous = from_origin("Endogenous"),
			origin_microbial = from_origin("Microbial"),
			origin_drug = from_origin("Drug"),
			origin_drug_or_steroid_metabolite = from_origin("Drug or Steroid Metabolite"),
			origin_drug_metabolite = from_origin("Drug Metabolite"),
			origin_plant = from_origin("Plant"),
			origin_toxin_pollutant = from_origin("Toxin/Pollutant"),
			origin_cosmetic = from_origin("Cosmetic"))
	})

# concentrations
concentration_data <- metabolites_xml %>%
	xml2::xml_children() %>%
	plyr::ldply(function(metabolite){
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
			"Normal", subject_condition))

concentration_data <- concentration_data %>%
	dplyr::filter(concentration_value != "1")

dir.create(paste0(staging_directory, "/data"))
m_data %>% readr::write_tsv(paste0(staging_directory, "/data/metabolites.tsv"))


prevelent_csf_metabolites <- m_data %>%
	dplyr::filter(
		taxonomy_super_class != "Inorganic compounds") %>%
	dplyr::inner_join(
		concentration_data %>%
			dplyr::filter(
				concentration_value > 1,
				concentration_units == "uM",
				biofluid == "Cerebrospinal Fluid (CSF)") %>%
			dplyr::arrange(desc(concentration_units)) %>%
			dplyr::group_by(accession) %>%
			dplyr::slice(1) %>%
			dplyr::ungroup()) %>%
	dplyr::arrange(
		taxonomy_kingdom,
		taxonomy_super_class,
		taxonomy_class,
		taxonomy_sub_class,
		taxonomy_molecular_framework) %>%
	dplyr::select(
		accession,
		name,
		concentration_value,
		taxonomy_super_class,
		taxonomy_class,
		taxonomy_sub_class) %>%
	data.frame

prevelent_csf_metabolites %>%
	readxl::write_xlsx(paste0(staging_directory, "data/prevelent_csf_metabolites.xlsx"))


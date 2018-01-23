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

schema <- "hmdb_180122"
staging_directory <- BioChemPantry::get_staging_directory(schema)

dir.create(paste0(staging_directory, "/data"))

z_data <- readr::read_tsv(paste0(staging_directory, "/data/hmdbendo_substances.tsv"))

# input paths
hmdb_metabolites_xml <- paste0(staging_directory, "/dump/hmdb_metabolites.xml") %>%
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

# this takes several hours to run
m_data <- hmdb_metabolites_xml %>%
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


# HMDB classification, some compounds have an extra chemical entities layer at the top that provides no information
kingdom_is_chemical_entities <- m_data$taxonomy_kingdom == "Chemical entities"
m_data <- m_data %>%
	dplyr::mutate(
		taxonomy_kingdom = ifelse(kingdom_is_chemical_entities, taxonomy_super_class, taxonomy_kingdom),
		taxonomy_super_class = ifelse(kingdom_is_chemical_entities, taxonomy_class, taxonomy_super_class),
		taxonomy_class = ifelse(kingdom_is_chemical_entities, taxonomy_sub_class, taxonomy_class),
		taxonomy_sub_class = ifelse(kingdom_is_chemical_entities, NA, taxonomy_sub_class))


m_data <- m_data %>%
	dplyr::left_join(
		z_data,
		by=c("accession"))

m_data %>% readr::write_tsv(paste0(staging_directory, "/data/metabolites.tsv"))

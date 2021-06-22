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

schema <- "drugbank"
staging_directory <- BioChemPantry::get_staging_directory(schema)

dir.create(paste0(staging_directory, "/data"))

z_data <- readr::read_tsv(paste0(staging_directory, "/data/dball_substances.tsv"))

# input paths
drugbank_xml <- paste0(staging_directory, "/data/db_database.xml") %>%
  xml2::read_xml()


xml_child_by_name <- function(xml, name){
  child_index <- xml %>%
    xml2::xml_children() %>%
    xml2::xml_name() %>%
    magrittr::equals(name) %>%
    which
  if(length(child_index) == 0) return(NA)
  if(length(child_index) > 1) {
    child_index <- child_index[1]
  }
  xml %>%
    xml2::xml_children() %>%
    magrittr::extract2(child_index)
}

get_value <- function(name, root){
  node <- root %>% xml_child_by_name(name)
  if(is.na(node)) return(NA)
  contents <- node %>% xml2::xml_contents()
  if(length(contents) == 0) return(NA)
  contents %>% as.character
}

# for single data. Takes about 1 hour
d_data <- drugbank_xml %>%
  xml2::xml_children() %>%
  plyr::ldply(function(drug){
    get_value <- function(name){
      node <- drug %>% xml_child_by_name(name)
      if(is.na(node)) return(NA)
      contents <- node %>% xml2::xml_contents()
      if(length(contents) == 0) return(NA)
      contents %>% as.character
    }
    get_group <- function(name){
      group <- drug %>% xml_child_by_name(name)
      if(is.na(group)) return(NA)
      group_child <- xml_child_by_name(drug, name)
      name <- substr(name, 1, nchar(name) - 1)
      node <- group_child %>% xml_child_by_name(name)
      if(is.na(node)) return(NA)
      siblings <- node %>% xml_siblings()
      contents <- node %>% xml2::xml_text()
      if(length(contents) == 0) return(NA)
      if(length(siblings) > 0) {
        siblings <- siblings %>% xml_text()
        contents <- append(contents, siblings)
        contents <- paste(contents, collapse = ", ")
      }
      else {contents}
    }
    tibble::data_frame(
      name = get_value("name"), drugbank_id = get_value("drugbank-id"),
      description = get_value("description"),
      cas_number = get_value("cas-number"), state = get_value("state"),
      groups = get_group("groups"), indication = get_value("indication"),
      pharmacodynamics = get_value("pharmacodynamics"),
      mechanism_action = get_value("mechanism-of-action"),
      toxicity = get_value("toxicity"), metabolism = get_value("metabolism"),
      absorption = get_value("absorption"), half_life = get_value("half-life"),
      protein_binding = get_value("protein-binding"),
      route_elim = get_value("route-of-elimination"),
      vol_of_dist = get_value("volume-of-distribution"),
      food_int = get_group("food-interactions")
    )
  })

#for drug interactions. Takes about 8 hours
drug_interactions <- drugbank_xml %>%
  xml2::xml_children() %>%
  plyr::ldply(function(drug){
    get_interact <- function(){
      drug_ints <- drug %>% xml_child_by_name("drug-interactions")
      if(is.na(drug_ints)) return(NA)
      # df for drug-interact
      interactions <- drug_ints %>% xml2::xml_children() %>%
      plyr::ldply(function(drug_int){
        tibble::data_frame(
          child_id = get_value("drugbank-id", drug_int),
          child_name = get_value("name", drug_int),
          child_description = get_value("description", drug_int))
      })
    }
    tibble::data_frame(parent_id = get_value("drugbank-id", drug),
                       parent_name = get_value("name", drug),
                       interaction <- get_interact())
  })

#for experimental properties. takes <1 hour
experimental_properties <- drugbank_xml %>%
  xml2::xml_children() %>%
  plyr::ldply(function(drug){
    get_props <- function(){
      exp_props <- drug %>% xml_child_by_name("experimental-properties")
      if(is.na(exp_props)) return(NA)
      # df for each exp prop
      properties <- exp_props %>% xml2::xml_children() %>%
      plyr::ldply(function(exp_prop){
        tibble::data_frame(
          child_type = get_value("kind", exp_prop),
          child_value = get_value("value", exp_prop))
        })
    }
    tibble::data_frame(parent_id = get_value("drugbank-id", drug),
                       parent_name = get_value("name", drug),
                       exp_props <- get_props())
  })

#drug pathway. takes <1 hour 
pathways <- drugbank_xml %>%
  xml2::xml_children() %>%
  plyr::ldply(function(drug){
    get_pathways <- function(){
      pathways <- drug %>% xml_child_by_name("pathways")
      if(is.na(pathways)) return(NA)
      # df for pathways
      drug_pathways <- pathways %>% xml2::xml_children() %>%
        plyr::ldply(function(path){
          tibble::data_frame(
            smpdb_id = get_value("smpdb-id", path),
            pathway = get_value("name", path),
            category = get_value("category", path))
        })
    }
    # retrieves ID's of enzymes in pathway
    get_enzymes <- function(){
      pathways <- drug %>% xml_child_by_name("pathways")
      if(is.na(pathways)) return(NA)
      pathway <- pathways %>% xml_child_by_name("pathway")
      if(is.na(pathway)) return(NA)
      enzymes <- pathway %>% xml_child_by_name("enzymes")
      if(is.na(enzymes)) return(NA)
      enzyme <- enzymes %>% xml_child_by_name("uniprot-id")
      if(is.na(enzyme)) return(NA)
      siblings <- enzyme %>% xml_siblings()
      contents <- enzyme %>% xml2::xml_text()
      if(length(contents) == 0) return(NA)
      if(length(siblings) > 0) {
        siblings <- siblings %>% xml_text()
        contents <- append(contents, siblings)
        contents <- paste(contents, collapse = ", ")
      }
      else {contents}
    }
    tibble::data_frame(parent_id = get_value("drugbank-id", drug),
                       parent_name = get_value("name", drug),
                       pathway <- get_pathways(), enzymes = get_enzymes())
  })

#drug reactions. takes <1 hour
reactions <- drugbank_xml %>%
  xml2::xml_children() %>%
  plyr::ldply(function(drug){
    # retrieves ID's of enzymes per reaction
    get_enzymes <- function(node){
      enzyme <- node %>% xml_child_by_name("enzyme")
      if(is.na(enzyme)) return(NA)
      uniprot_id <- enzyme %>% xml_child_by_name("uniprot-id")
      if(is.na(uniprot_id)) return(NA)
      siblings <- enzyme %>% xml_siblings()
      contents <- uniprot_id %>% xml2::xml_text()
      if(length(contents) == 0) return(NA)
      if(length(siblings) > 0) {
        siblings <- as_list(siblings)
        siblings_ids <- lapply(siblings, '[[', "uniprot-id")
        siblings_ids <- lapply(siblings_ids, '[[', 1)
        contents <- append(contents, siblings_ids)
        contents <- paste(contents, collapse = ", ")
      }
      else {contents}
    }
    get_reactions <- function(){
      reactions <- drug %>% xml_child_by_name("reactions")
      if(is.na(reactions)) return(NA)
      # df for reactions
      drug_reactions <- reactions %>% xml2::xml_children() %>%
        plyr::ldply(function(react){
          left <- react %>% xml_child_by_name("left-element")
          right <- react %>% xml_child_by_name("right-element")
          enzymes <- react %>% xml_child_by_name("enzymes")
          tibble::data_frame(
            sequence = get_value("sequence", react), 
            left_element = get_value("name", left), 
            right_element = get_value("name", right), 
            enzymes = get_enzymes(enzymes))
          })
    }
    tibble::data_frame(parent_id = get_value("drugbank-id", drug),
                       parent_name = get_value("name", drug),
                       reacts <- get_reactions())
  })

#drug targets. takes <1 hour
targets <- drugbank_xml %>%
  xml2::xml_children() %>%
  plyr::ldply(function(drug){
    get_targets <- function(){
      targets <- drug %>% xml_child_by_name("targets")
      if(is.na(targets)) return(NA)
      # df for targets
      drug_targets <- targets %>% xml2::xml_children() %>%
        plyr::ldply(function(targ){
          tibble::data_frame(
            target_id = get_value("drugbank-id", targ),
            target_name = get_value("name", targ),
            organism = get_value("organism", targ), 
            known_action = get_value("known-action", targ))
        })
    }
    tibble::data_frame(parent_id = get_value("drugbank-id", drug),
                       parent_name = get_value("name", drug),
                       target <- get_targets())
  })

#change column in z_data for merge
names(z_data)[names(z_data) == "accession"] <- "parent_id"

#merge each df with z_data
d_data_final <- d_data %>%
  dplyr::left_join(
    z_data,
    by=c("parent_id"))

drug_int_final <- drug_interactions %>%
  dplyr::left_join(
    z_data,
    by=c("parent_id"))

exp_prop_final <- experimental_properties %>%
  dplyr::left_join(
    z_data,
    by=c("parent_id"))

pathways_final <- pathways %>%
  dplyr::left_join(
    z_data,
    by=c("parent_id"))

reactions_final <- reactions %>%
  dplyr::left_join(
    z_data,
    by=c("parent_id"))

targets_final <- targets %>%
  dplyr::left_join(
    z_data,
    by=c("parent_id"))

#write each df to a tsv
d_data_final %>% readr::write_tsv(paste0(staging_directory, "/data/db_short_data.tsv"))
drug_int_final %>% readr::write_tsv(paste0(staging_directory, "/data/db_drug_int.tsv"))
exp_prop_final %>% readr::write_tsv(paste0(staging_directory, "/data/db_exp_prop.tsv"))
pathways_final %>% readr::write_tsv(paste0(staging_directory, "/data/db_pathways.tsv"))
reactions_final %>% readr::write_tsv(paste0(staging_directory, "/data/db_reactions.tsv"))
targets_final %>% readr::write_tsv(paste0(staging_directory, "/data/db_targets.tsv"))

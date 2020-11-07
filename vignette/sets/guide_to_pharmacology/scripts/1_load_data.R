# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(assertr)
library(BioChemPantry)
library(Zr)
library(curl)


httr::set_config( httr::config( ssl_verifypeer = 0L ) )
schema <- "guide_to_pharmacology_2020_03"
pantry <- get_pantry()
pantry %>% BioChemPantry::create_schema(schema)
pantry %>% BioChemPantry::set_schema(schema)
staging_directory <- BioChemPantry::get_staging_directory(schema)
dir.create(paste0(staging_directory, "/dump"), recursive=TRUE)
dir.create(paste0(staging_directory, "/data"), recursive=TRUE)


#########################
HGNC_mapping_fname <- paste0(staging_directory, "/dump/GtP_to_HGNC_mapping.csv")
httr::GET(url="http://www.guidetopharmacology.org/DATA/GtP_to_HGNC_mapping.csv",
  httr::write_disk(HGNC_mapping_fname, overwrite=TRUE))

uniprot_mapping_fname <- paste0(staging_directory, "/dump/GtP_to_UniProt_mapping.csv")
httr::GET(
	url="http://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv",
	httr::write_disk(uniprot_mapping_fname, overwrite=TRUE))


targets_fname <- paste0(staging_directory, "/dump/targets_and_families.csv")
httr::GET(
	url="http://www.guidetopharmacology.org/DATA/targets_and_families.csv",
	httr::write_disk(targets_fname, overwrite=TRUE))

ligands_fname <- paste0(staging_directory, "/dump/ligands.csv")
httr::GET(
	url="http://www.guidetopharmacology.org/DATA/ligands.csv",
	httr::write_disk(ligands_fname, overwrite=TRUE))

peptides_fname <- paste0(staging_directory, "/dump/peptides.csv")
httr::GET(
	url="http://www.guidetopharmacology.org/DATA/peptides.csv",
  httr::write_disk(peptides_fname, overwrite=TRUE))

interactions_fname <- paste0(staging_directory, "/dump/interactions.csv")
httr::GET(
	url="http://www.guidetopharmacology.org/DATA/interactions.csv",
	httr::write_disk(interactions_fname, overwrite=TRUE))


zinc_ids <- Zr::catalog_items(
	catalog_short_name="iuphar",
	output_fields=c(
		"zinc_id",
		"supplier_code",
		"substance.preferred_name",
		"substance.smiles"),
	result_batch_size=100,
	verbose=T) %>%
	dplyr::select(
		zinc_id,
		ligand_id=supplier_code,
		preferred_name = substance.preferred_name,
		zinc_smiles = substance.smiles)
zinc_ids %>% readr::write_tsv(
	paste0(staging_directory, "/data/zinc_ids.tsv"))


##################################

uniprot_mapping <- readr::read_csv(
	file=uniprot_mapping_fname,
	col_types= readr::cols(
  	uniprot_id = readr::col_character(),
  	species = readr::col_character(),
  	iuphar_name = readr::col_character(),
  	iuphar_id = readr::col_double(),
  	gtp_url = readr::col_character())) %>%
	dplyr::rename(
		uniprot_accn = uniprot_id)

uniprot_mapping <- pantry %>%
	dplyr::copy_to(
		uniprot_mapping,
		"uniprot_mapping",
		temporary=TRUE) %>%
	dplyr::left_join(
		pantry %>%
			BioChemPantry::schema_tbl("uniprot_2017_2.to_uniprot_entry") %>%
			dplyr::select(uniprot_accn, uniprot_entry=id_value),
		by=c("uniprot_accn")) %>%
	dplyr::collect(n=Inf)

#pantry %>% BioChemPantry::drop_table("uniprot_mapping")

############################################33

targets_1 <- readr::read_csv(
	file=targets_fname,
	col_types=readr::cols_only(
  	Type = col_character(),
  	`Target id` = col_double(),
  	`Target name` = col_character(),
  	`Family id` = col_double(),
  	`Family name` = col_character(),
  	`Human SwissProt` = col_character(),
  	`Human Entrez Gene` = col_character(),
  	`Rat SwissProt` = col_character(),
  	`Rat Entrez Gene` = col_character(),
  	`Mouse SwissProt` = col_character(),
  	`Mouse Entrez Gene` = col_character()))

targets <- rbind(
	targets_1 %>%
		dplyr::mutate(organism='human') %>%
		dplyr:::select(
			target_id=`Target id`,
			target_name=`Target name`,
			type=Type,
			family_id=`Family id`,
			family_name=`Family name`,
			organism,
			gene_symbol=`Human Entrez Gene`,
			uniprot_accn=`Human SwissProt`),
	targets_1 %>%
		dplyr::mutate(organism='rat') %>%
		dplyr:::select(
			target_id=`Target id`,
			target_name=`Target name`,
			type=Type,
			family_id=`Family id`,
			family_name=`Family name`,
			organism,
			gene_symbol=`Rat Entrez Gene`,
			uniprot_accn=`Rat SwissProt`),
	targets_1 %>%
		mutate(organism='mouse') %>%
		dplyr:::select(
			target_id=`Target id`,
			target_name=`Target name`,
			type=Type,
			family_id=`Family id`,
			family_name=`Family name`,
			organism,
			gene_symbol=`Mouse Entrez Gene`,
			uniprot_accn=`Mouse SwissProt`)) %>%
	dplyr::left_join(
		uniprot_mapping,
		by=c("uniprot_accn"))

targets <- targets %>%
	dplyr::mutate(
		target_name = str_replace_all(target_name, "<[^>]+>", ""),
		target_name = str_replace_all(target_name, "&alpha;", "a"),
		target_name = str_replace_all(target_name, "&beta;", "b"),
		target_name = str_replace_all(target_name, "&delta;", "d"),
		target_name = str_replace_all(target_name, "&Delta;", "D"),
		target_name = str_replace_all(target_name, "&epsilon;", "e"),
		target_name = str_replace_all(target_name, "&gamma;", "g"),
		target_name = str_replace_all(target_name, "&kappa;", "k"),
		target_name = str_replace_all(target_name, "&lambda;", "l"),
		target_name = str_replace_all(target_name, "&eta;", "n"),
		target_name = str_replace_all(target_name, "&pi;", "p"),
		target_name = str_replace_all(target_name, "&rho;", "r"),
		target_name = str_replace_all(target_name, "&sigma;", "s"),
		target_name = str_replace_all(target_name, "&theta;", "t"),
		target_name = str_replace_all(target_name, "&mu;", "u"),
		target_name = str_replace_all(target_name, "&zeta;", "z"),
		target_name = str_replace_all(target_name, "&uuml;", "u"),
		family_name = str_replace_all(family_name, "<[^>]+>", ""),
		family_name = str_replace_all(family_name, "&alpha;", "a"),
		family_name = str_replace_all(family_name, "&beta;", "b"),
		family_name = str_replace_all(family_name, "&delta;", "d"),
		family_name = str_replace_all(family_name, "&Delta;", "D"),
		family_name = str_replace_all(family_name, "&epsilon;", "e"),
		family_name = str_replace_all(family_name, "&gamma;", "g"),
		family_name = str_replace_all(family_name, "&kappa;", "k"),
		family_name = str_replace_all(family_name, "&lambda;", "l"),
		family_name = str_replace_all(family_name, "&eta;", "n"),
		family_name = str_replace_all(family_name, "&pi;", "p"),
		family_name = str_replace_all(family_name, "&rho;", "r"),
		family_name = str_replace_all(family_name, "&sigma;", "s"),
		family_name = str_replace_all(family_name, "&theta;", "t"),
		family_name = str_replace_all(family_name, "&mu;", "u"),
		family_name = str_replace_all(family_name, "&zeta;", "z"),
		family_name = str_replace_all(family_name, "&uuml;", "u"))

targets %>%
	dplyr::filter(
		target_name %>% str_detect("[^a-zA-Z0-9/\\- .(),;+':]")) %>%
		dplyr::select(target_name, uniprot_entry) %>%
		data.frame %>%
		head(10)

targets %>%
	dplyr::filter(
		family_name %>% str_detect("[^a-zA-Z0-9/\\- .(),;+':]")) %>%
		dplyr::select(family_name, uniprot_entry) %>%
		data.frame %>%
		head(10)

############################################################3

ligands_1 <- readr::read_csv(
	file=ligands_fname,
	col_types = readr::cols_only(
	  `Ligand id` = readr::col_double(),
	  Name = readr::col_character(),
	  Approved = readr::col_character(),
	  Withdrawn = readr::col_character(),
	  Labelled = readr::col_character(),
	  Radioactive = readr::col_character(),
	  Type = readr::col_character(),
	  `PubChem SID` = readr::col_double(),
	  `PubChem CID` = readr::col_double(),
	  SMILES = col_character()))

ligands <- ligands_1 %>%
	dplyr::mutate(
		approved = ifelse(Approved == "yes" & Withdrawn != "yes", TRUE, NA),
		labelled = ifelse(Labelled == "yes", TRUE, NA),
		radioactive = ifelse(Radioactive == "yes", TRUE, NA)) %>%
	dplyr:::select(
		ligand_id=`Ligand id`,
		ligand_type=Type,
		name=Name,
		smiles=SMILES,
		approved,
		labelled,
		radioactive,
		pubchem_sid = `PubChem SID`,
		pubchem_cid = `PubChem CID`) %>%
	dplyr::left_join(zinc_ids, by=c("ligand_id"))

ligands <- ligands %>%
	dplyr::mutate(
		name = str_replace_all(name, " +$", ""),
		name = str_replace_all(name, "<[^>]+>", ""),
		name = str_replace_all(name, "&alpha;", "a"),
		name = str_replace_all(name, "&beta;", "b"),
		name = str_replace_all(name, "&gamma;", "g"),
		name = str_replace_all(name, "&delta;", "d"),
		name = str_replace_all(name, "&Delta;", "d"),
		name = str_replace_all(name, "&epsilon;", "e"),
		name = str_replace_all(name, "&gamma;", "g"),
		name = str_replace_all(name, "&kappa;", "k"),
		name = str_replace_all(name, "&lambda;", "l"),
		name = str_replace_all(name, "&eta;", "n"),
		name = str_replace_all(name, "&pi;", "p"),
		name = str_replace_all(name, "&rho;", "r"),
		name = str_replace_all(name, "&sigma;", "s"),
		name = str_replace_all(name, "&theta;", "t"),
		name = str_replace_all(name, "&mu;", "u"),
		name = str_replace_all(name, "&zeta;", "z"),
		name = str_replace_all(name, "&uuml;", "u"),
		name = str_replace_all(name, "&plusmn;", "+/-"))

ligands %>%
	dplyr::filter(
		name %>% str_detect("[^a-zA-Z0-9/\\- .(),;+':\\[\\]]\\{\\}")) %>%
		dplyr::select(name, ligand_id) %>%
		data.frame %>%
		head(10)



#########################################3

interactions_1 <- readr::read_csv(
	file=interactions_fname,
	col_types=readr::cols_only(
			target_id = readr::col_integer(),
			ligand_id = readr::col_integer(),
			type = readr::col_character(),
			action = readr::col_character(),
			action_comment = readr::col_character(),
			endogenous = readr::col_character(),
			primary_target = readr::col_character(),
			original_affinity_units = readr::col_character(),
			original_affinity_median_nm = readr::col_double(),
			original_affinity_relation = readr::col_character(),
			assay_description = readr::col_character(),
			receptor_site = readr::col_character(),
			ligand_context = readr::col_character(),
			pubmed_id = readr::col_character()))

interactions <- interactions_1 %>%
	dplyr::mutate(
		endogenous = ifelse(endogenous == "t", TRUE, FALSE),
		primary_target = ifelse(primary_target == "t", TRUE, FALSE)) %>%
	dplyr::rename(
		interaction_type = type,
		interaction_action = action) %>%
	dplyr::left_join(
		ligands %>% dplyr::select(ligand_id, name, smiles, zinc_id),
		by="ligand_id") %>%
	dplyr::left_join(
		targets %>% dplyr::select(target_id, target_name, uniprot_entry),
		by="target_id")



#####

pantry %>% dplyr::copy_to(
	targets,
	"targets",
	temporary=FALSE,
	indexes = list(
		"target_id",
		"uniprot_accn"))

pantry %>% dplyr::copy_to(
	ligands,
	"ligands",
	temporary=FALSE,
	indexes = list(
		"ligand_id",
		"zinc_id",
		"smiles"))

pantry %>% dplyr::copy_to(
	interactions,
	"interactions",
	temporary=FALSE,
	indexes = list(
		"target_id", "uniprot_entry",
		"ligand_id",
		c("target_id", "ligand_id"),
		c("uniprot_entry", "zinc_id")))



metabolic_interactions <- interactions %>%
#	dplyr::filter(endogenous) %>%
	dplyr::select(
		target_id,
		ligand_id,
		zinc_id,
		affinity_median) %>%
	dplyr::left_join(
		ligands_1 %>%
			dplyr::select(
				ligand_id =`Ligand id`,
				name=`Name`,
				synonyms=`Synonyms`),
		by="ligand_id") %>%
	left_join(
		targets %>%
			dplyr::select(
				target_id, target, uniprot_accn, target_type=type, family_name),
		by=c("target_id")) %>%
	dplyr::arrange(desc(affinity_median)) %>%
	distinct(target_id, ligand_id)

signaling_metabolites <- metabolic_interactions %>%
	dplyr::distinct(ligand_id, target_type) %>%
	dplyr::count(ligand_id, sort=T) %>%
	dplyr::filter(n>1)

signaling_metabolite_interactions <- signaling_metabolites %>%
	left_join(
		metabolic_interactions,
		by="ligand_id")

signaling_metabolite_interactions %>%
	write.xlsx(paste0(output_path, "/signaling_metabolite_interactions.xlsx"))

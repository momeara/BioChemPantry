# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(assertr)
library(plyr)
library(dplyr)
library(readr)
library(assertthat)
library(BioChemPantry)

source("scripts/4__activity_types.R")
staging_directory <- get_staging_directory("chembl21")
chembl_db <- get_pantry("chembl21")


target_info <- read_tsv(paste0(staging_directory, "/data/target_info.tsv"))
zinc_to_chembl <- read_tsv(paste0(staging_directory, "/data/zinc_to_chembl.tsv"))
chembl_to_aggregators <- read_tsv(
	paste0(staging_directory, "/data/chembl_to_aggregators_160704.tsv")



c("assays", "activities", "target_dictionary", "target_components",
	"component_sequences", "chembl_id_lookup", "compound_structures",
	"source") %in% src_tbls(chembl_db) %>% all %>%
	assert_that

assays_tbl <- chembl_db %>% tbl("assays")
activities_tbl <- chembl_db %>% tbl("activities")
docs_tbl <- chembl_db %>% tbl("docs")
chembl_id_lookup_tbl <- chembl_db %>% tbl("chembl_id_lookup")
component_sequences_tbl <- chembl_db %>% tbl("component_sequences")
compound_properties_tbl <- chembl_db %>% tbl("compound_properties")
compound_structures_tbl <- chembl_db %>% tbl("compound_structures")
source_tbl <- chembl_db %>% tbl("source")
target_dictionary_tbl <- chembl_db %>% tbl("target_dictionary")
target_components_tbl <- chembl_db %>% tbl("target_components")




### explore data criteria
target_dictionary_tbl %>% count(target_type, sort=T) %>% collect(n=Inf) %>% cdata.frame


# any duplicate smiles?
compound_structures_tbl %>%
	dplyr::count(canonical_smiles) %>% dplyr::ungroup() %>%
	dplyr::count(n, sort=T) %>%	dplyr::ungroup() %>%
	data.frame()

#assays_relationship_types
assays_tbl %>%
	dplyr::count(relationship_type) %>%
	dplyr::left_join( tbl(chembl_db, "relationship_type"), by="relationship_type") %>%
	dplyr::arrange(desc(n)) %>%
	dplyr::collect(n=Inf) %>%
	data.frame()

#assays_confidence_scores
bassays_tbl %>%
	dplyr::count(confidence_score) %>%
	dplyr::arrange(confidence_score) %>%
	dplyr::left_join(tbl(chembl_db, "confidence_score_lookup"), by=c("confidence_score")) %>%
	dplyr::collect(n=Inf) %>%
	data.frame()

#assays_types
assays_tbl %>%
	dplyr::count(assay_type) %>%
	dplyr::left_join(tbl(chembl_db, "assay_type"), by="assay_type") %>%
	dplyr::arrange(desc(n)) %>%
	dplyr::collect(n=Inf) %>%
	data.frame()



activities_tbl %>%
	dplyr::count(potential_duplicate) %>%
	dplyr::collect(n=Inf) %>%
	data.frame()

activities_tbl %>%
	dplyr::count(data_validity_comment) %>%
	dplyr::collect(n=Inf) %>%
	data.frame()

x <- activities_tbl %>%
	dplyr::count(activity_comment, sort=T) %>%
	dplyr::collect(n=Inf) %>%
	data.frame()

activities_tbl %>%
	dplyr::filter(!(standard_type %in% activity_types)) %>%
	dplyr::count(standard_type, sort=T) %>%
	dplyr::filter(n > 250)

activities_tbl %>%
	dplyr::filter(!(standard_units %in% activity_units)) %>%
	dplyr::count(standard_units, sort=T) %>%
	dplyr::filter(n > 1)
#
############## prepare primary data ############

targets <- target_dictionary_tbl %>%
	dplyr::filter(
		target_type == "PROTEIN COMPLEX" |
		target_type == "SINGLE PROTEIN") %>%
	dplyr::left_join(
		target_components_tbl %>%
			dplyr::select(
				tid,
				component_id),
		by="tid") %>%
	dplyr::left_join(
		component_sequences_tbl %>%
			dplyr::select(
				component_id,
				uniprot_accn = accession,
				description),
		by="component_id") %>%
	dplyr::select(
		tid,
		target_type,
		uniprot_accn,
		target_description = description) %>%
	dplyr::collect(n=Inf) %>%
	dplyr::inner_join(
		target_info %>%
			dplyr::select(
				uniprot_accn,
				uniprot_entry,
				gene_name = `Gene names  (primary )`),
		by=c("uniprot_accn"))

assays <- assays_tbl %>%
	dplyr::filter(
		# D -> direct protein target assigned
		# H -> homologous protein target assigned
		# M -> Molecular target other than protein assigned
		relationship_type %in% c('D', 'H', 'M'),
		# Confidence score, indicating how accurately the assigned target(s)
		# represents the actually assay target. 1=low confidence; 9=high confidence
		confidence_score >= 7,
		# F -> functional
		# B -> binding
		assay_type %in% c('F', 'B')) %>%
	dplyr::select(
		assay_id,
		assay_chembl_id=chembl_id,
		tid,
		src_id, src_assay_id,
		assay_description = description,
		assay_test_type,
		assay_category,
		assay_organism,
		assay_relationship_type=relationship_type,
		assay_confidence_score=confidence_score,
		assay_curated_by=curated_by) %>%
	dplyr::collect(n=Inf)


activities <- activities_tbl %>%
	dplyr::filter(
		potential_duplicate == 0,
		(is.na(data_validity_comment) | (data_validity_comment == "Manually validated")),
		(is.na(activity_comment) | !(activity_comment %in% c("inactive", "Not Active", "inconclusive", "Inconclusive"))),
		!is.na(standard_units), standard_units %in% activity_units,
		!is.na(standard_type), standard_type %in% activity_types,
		!is.na(standard_relation), standard_relation %in% c('<', '=')) %>%
	dplyr::select(
		assay_id,
		doc_id,
		molregno,
		standard_type,
		standard_units,
		standard_value,
		data_validity_comment,
		activity_comment) %>%
	dplyr::left_join(
		docs_tbl %>%
			dplyr::mutate(pubmed_id = as.integer(pubmed_id)) %>%
			dplyr::select(doc_id, journal, year, doi, pubmed_id),
		by="doc_id") %>%
	dplyr::collect(n=Inf) %>%
	dplyr::mutate(
		activity = normalize_activity_value(standard_type, standard_units, standard_value))

compounds <- compound_structures_tbl %>%
	dplyr::filter(
		nchar(canonical_smiles) < 1024) %>%
	dplyr::left_join(
		chembl_id_lookup_tbl %>%
			dplyr::filter(entity_type == "COMPOUND"),
		by=c("entity_id"="molregno")) %>%
	dplyr::left_join(compound_properties_tbl, by="molregno") %>%
	dplyr::select(
		molregno,
		chembl_id,
		chembl_smiles=canonical_smiles,
		molecular_weight = mw_freebase,
		ALogP=alogp,
		HBA=hba,
		HBD=hbd,
		PSA=psa,
		RTB=rtb,
		most_acidic_pKa=acd_most_apka,
		most_basic_pKa=acd_most_bpka,
		molecular_species,
		n_heavy_atoms=heavy_atoms) %>%
	dplyr::collect(n=Inf) %>%
	dplyr::filter(
		!(chembl_smiles %>% stringr::str_detect("[\\[]11CH3[\\]]")))



# note there are 389 smiles that map to more than one chembl id e.g.
# Nc1ccc2ccccc2c1N=Nc3ccccc3  CHEMBL86177
# Nc1ccc2ccccc2c1N=Nc3ccccc3  CHEMBL1734418

# chembl_id doesn't map 1-1 to zinc_id, for example
#   chembl distinguishes salt forms, but zinc does not
#   zinc distinguishes steochemistry, but in some cases chembl does not
# We'd like to use zinc ids as that interfaces better with our other tools
#   for each chembl_id select a single zinc_id
compounds <- compounds %>%
	dplyr::left_join(chembl_to_aggregators, by="chembl_id") %>%
	dplyr::inner_join(
		zinc_to_chembl %>%
			dplyr::arrange(zinc_id) %>%
			group_by(chembl_id) %>%
				dplyr::summarize(
					zinc_id = zinc_id[1],
					alt_zinc_ids = paste(zinc_id, collapse=";")) %>%
			dplyr::ungroup() %>%
			dplyr::distinct(chembl_id, .keep_all=T),
		by=c("chembl_id")) %>%
	dplyr::distinct(chembl_smiles, .keep_all=T)


full_chembl_data <- targets %>%
	dplyr::inner_join(assays, by=c("tid")) %>%
	dplyr::inner_join(activities, by=c("assay_id")) %>%
	dplyr::inner_join(compounds, by=c("molregno"))

# remove duplicate activies
# preferentially taking activities by the following criteria
# and then remove low activity measurements
# and measurements from those that are thougth to be highly promiscuous
full_chembl_data <- full_chembl_data %>%
	dplyr::arrange(
		target_type == "SINGLE_PROTEIN",
		assay_relationship_type == 'D',
		assay_curated_by == "Expert",
		assay_confidence_score,
		is.na(activity_comment),
		activity,
		desc(year),
		zinc_id) %>%
	dplyr::distinct(uniprot_accn, chembl_smiles, .keep_all=T) %>%
	dplyr::mutate(
		active = (activity <= 10000) &
			(activity > 0) &
			(is.na(MaxTC_aggregator) | MaxTC_aggregator <= .70 | ALogP <= 3.5 | activity <= 1000))






##### Save data ####

full_chembl_data %>% write_tsv(paste0(staging_directory, "/data/full_chembl_data.tsv"))
save("full_chembl_data", file=paste0(staging_directory, "/data/full_chembl_data.Rdata"))





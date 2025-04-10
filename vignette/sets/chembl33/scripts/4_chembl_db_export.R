
library(assertr)
library(plyr)
library(dplyr)
library(readr)
library(assertthat)
library(BioChemPantry)

source("scripts/4__activity_types.R")
source("scripts/parameters.R")
staging_directory <- BioChemPantry::get_staging_directory(schema)
pantry <- BioChemPantry::get_pantry(schema)

target_info <- readr::read_tsv(
  paste0(staging_directory, "/data/target_info.tsv"),
  show_col_types = FALSE)

load(paste0(staging_directory, "/data/zinc_to_chembl.Rdata"))
chembl_to_aggregators <- readr::read_tsv(
	paste0(staging_directory, "/data/chembl_to_aggregators_170813.tsv"))

c("assays",
  "activities",
  "target_dictionary",
  "target_components",
	"component_sequences",
  "chembl_id_lookup",
  "compound_structures",
	"source") %in% DBI::dbListTables(pantry) |>
  all() |>
	assertthat::assert_that()

assays_tbl <- pantry |> dplyr::tbl("assays")
activities_tbl <- pantry |> dplyr::tbl("activities")
docs_tbl <- pantry |> dplyr::tbl("docs")
chembl_id_lookup_tbl <- pantry |> dplyr::tbl("chembl_id_lookup")
component_sequences_tbl <- pantry |> dplyr::tbl("component_sequences")
compound_properties_tbl <- pantry |> dplyr::tbl("compound_properties")
compound_structures_tbl <- pantry |> dplyr::tbl("compound_structures")
source_tbl <- pantry |> dplyr::tbl("source")
target_dictionary_tbl <- pantry |> dplyr::tbl("target_dictionary")
target_components_tbl <- pantry |> dplyr::tbl("target_components")


### explore data criteria
target_dictionary_tbl |>
  dplyr::count(target_type, sort = TRUE) |>
  dplyr::collect(n = Inf) |> 
  data.frame()


# any duplicate smiles?
compound_structures_tbl |>
	dplyr::count(canonical_smiles) |>
  dplyr::ungroup() |>
	dplyr::count(n, sort = TRUE) |>
  dplyr::ungroup() |>
	data.frame()

#assays_relationship_types
assays_tbl |>
	dplyr::count(relationship_type) |>
	dplyr::left_join(
	  dplyr::tbl(pantry, "relationship_type"),
	  by = "relationship_type") |>
	dplyr::arrange(desc(n)) |>
	dplyr::collect(n = Inf) |>
	data.frame()

#assays_confidence_scores
assays_tbl |>
	dplyr::count(confidence_score) |>
	dplyr::left_join(
	  dplyr::tbl(pantry, "confidence_score_lookup"),
	  by = c("confidence_score")) |>
	dplyr::collect(n = Inf) |>
  dplyr::arrange(confidence_score) |>
	data.frame()

#assays_types
assays_tbl |>
	dplyr::count(assay_type) |>
	dplyr::left_join(
	  dplyr::tbl(pantry, "assay_type"),
	  by = "assay_type") |>
	dplyr::arrange(dplyr::desc(n)) |>
	dplyr::collect(n = Inf) |>
	data.frame()

activities_tbl |>
	dplyr::count(potential_duplicate) |>
	dplyr::collect(n = Inf) |>
	data.frame()

activities_tbl |>
	dplyr::count(data_validity_comment) |>
	dplyr::collect(n = Inf) |>
	data.frame()

activities_tbl |>
	dplyr::count(activity_comment, sort = T) |>
	dplyr::collect(n = Inf) |>
	data.frame()

activities_tbl |>
	dplyr::filter(!(standard_type %in% activity_types)) |>
	dplyr::count(standard_type, sort = TRUE) |>
	dplyr::filter(n > 250)

activities_tbl |>
	dplyr::filter(!(standard_units %in% activity_units)) |>
	dplyr::count(standard_units, sort = TRUE) |>
	dplyr::filter(n > 1)
#
############## prepare primary data ############

targets <- target_dictionary_tbl |>
	dplyr::filter(
		target_type == "PROTEIN COMPLEX" |
		target_type == "SINGLE PROTEIN") |>
	dplyr::left_join(
		target_components_tbl |>
			dplyr::select(
				tid,
				component_id),
		by="tid") |>
	dplyr::left_join(
		component_sequences_tbl |>
			dplyr::select(
				component_id,
				uniprot_accn = accession,
				description),
		by="component_id") |>
	dplyr::select(
		tid,
		target_type,
		uniprot_accn,
		target_description = description) |>
	dplyr::collect(n = Inf) |>
	dplyr::inner_join(
		target_info |>
			dplyr::select(
				uniprot_accn,
				uniprot_entry,
				gene_name),
		by = c("uniprot_accn"))

assays <- assays_tbl |>
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
		assay_type %in% c('F', 'B')) |>
	dplyr::select(
		assay_id,
		assay_chembl_id=chembl_id,
		tid,
		src_id, src_assay_id,
		assay_description = description,
		assay_test_type,
		assay_category,
		assay_organism,
		assay_relationship_type = relationship_type,
		assay_confidence_score = confidence_score,
		assay_curated_by = curated_by) |>
	dplyr::collect(n = Inf)


activities <- activities_tbl |>
	dplyr::filter(
		potential_duplicate == 0,
		(is.na(data_validity_comment) | (data_validity_comment == "Manually validated")),
		(is.na(activity_comment) | !(activity_comment %in% c("inactive", "Not Active", "inconclusive", "Inconclusive"))),
		!is.na(standard_units), standard_units %in% activity_units,
		!is.na(standard_type), standard_type %in% activity_types,
		!is.na(standard_relation), standard_relation %in% c('<', '=')) |>
	dplyr::select(
		assay_id,
		doc_id,
		molregno,
		standard_type,
		standard_units,
		standard_value,
		data_validity_comment,
		activity_comment) |>
	dplyr::left_join(
		docs_tbl |>
			dplyr::mutate(pubmed_id = as.integer(pubmed_id)) |>
			dplyr::select(doc_id, journal, year, doi, pubmed_id),
		by = "doc_id") |>
	dplyr::collect(n = Inf) |>
	dplyr::mutate(
		activity = normalize_activity_value(standard_type, standard_units, standard_value))

chembl_compounds <- compound_structures_tbl |>
	dplyr::filter(
		nchar(canonical_smiles) < 1024) |>
	dplyr::left_join(
		chembl_id_lookup_tbl |>
			dplyr::filter(entity_type == "COMPOUND"),
		by = c("molregno"="entity_id")) |>
	dplyr::left_join(compound_properties_tbl, by = "molregno") |>
	dplyr::select(
		molregno,
		chembl_id,
		chembl_smiles = canonical_smiles,
		molecular_weight = mw_freebase,
		ALogP = alogp,
		HBA = hba,
		HBD = hbd,
		PSA = psa,
		RTB = rtb,
		#most_acidic_pKa = acd_most_apka,
		#most_basic_pKa = acd_most_bpka,
		molecular_species,
		n_heavy_atoms = heavy_atoms) |>
	dplyr::collect(n = Inf) |>
	dplyr::filter(
		!(chembl_smiles |> stringr::str_detect("[\\[]11CH3[\\]]")))



# note there are 389 smiles that map to more than one chembl id e.g.
# Nc1ccc2ccccc2c1N=Nc3ccccc3  CHEMBL86177
# Nc1ccc2ccccc2c1N=Nc3ccccc3  CHEMBL1734418

# chembl_id doesn't map 1-1 to zinc_id, for example
#   chembl distinguishes salt forms, but zinc does not
#   zinc distinguishes steochemistry, but in some cases chembl does not
# We'd like to use zinc ids as that interfaces better with our other tools
#   for each chembl_id select a single zinc_id
compounds <- chembl_compounds |>
	dplyr::left_join(zinc_to_chembl, by = "chembl_id") |>
	dplyr::left_join(
		chembl_to_aggregators |>
			rename(tc_to_aggregator = Tc) |>
			dplyr::group_by(chembl_zinc_id) |>
				dplyr::arrange(desc(tc_to_aggregator)) |>
				dplyr::slice(1) |>
			dplyr::ungroup(),
			by=c("zinc_id"="chembl_zinc_id"))

# there are now 5 potential ids: molregno, chembl_id, chembl_smiles, zinc_id, smiles
# verify these are redundant:
compounds |> BioChemPantry::summarize_map("zinc_id", "smiles")
#   168 examples with same smiles but different zinc_id
compounds |> BioChemPantry::summarize_map("molregno", "chembl_id")
#   molregno and chembl_id are 1-1
compounds |> BioChemPantry::summarize_map("chembl_id", "chembl_smiles")
#   464 examples with same chembl_smiles but different chembl_ids

compounds <- compounds |> filter(!is.na(zinc_id))

#	dplyr::inner_join(
#		zinc_to_chembl |>
#			dplyr::arrange(zinc_id) |>
#			dplyr::group_by(chembl_id) |>
#				dplyr::summarize(
#					zinc_id = zinc_id[1],
#					alt_zinc_ids = paste(zinc_id, collapse=";")) |>
#			dplyr::ungroup() |>
#			dplyr::distinct(chembl_id, .keep_all=T),
#		by=c("chembl_id")) |>
#	dplyr::distinct(chembl_smiles, .keep_all=T)


full_chembl_data <- targets |>
	dplyr::inner_join(
	  assays,
	  by = c("tid"),
	  relationship = "many-to-many") |>
	dplyr::inner_join(
	  activities,
	  by = c("assay_id"),
	  relationship = "many-to-many") |>
	dplyr::inner_join(
	  compounds,
	  by = c("molregno"),
	  relationship = "many-to-many")

# remove duplicate activities
# preferentially taking activities by the following criteria
# and then remove low activity measurements
# and measurements from those that are thougth to be highly promiscuous
full_chembl_data <- full_chembl_data |>
	dplyr::arrange(
		target_type == "SINGLE_PROTEIN",
		assay_relationship_type == 'D',
		assay_curated_by == "Expert",
		assay_confidence_score,
		is.na(activity_comment),
		activity,
		desc(year),
		zinc_id) |>
	dplyr::distinct(uniprot_accn, chembl_smiles, .keep_all = T) |>
	dplyr::mutate(
		active = (activity <= 10000) &
			(activity > 0))
# TODO Get more stringent filtering for potential aggregators working again...
# &
#			(is.na(tc_to_aggregator) |
#			   tc_to_aggregator > .70 |
#			   ALogP <= 3.5 |
#			   activity <= 1000))

# now uniquely map to zinc_id
full_chembl_data <- full_chembl_data |>
	dplyr::arrange(
		active,
		target_type == "SINGLE_PROTEIN",
		assay_relationship_type == 'D',
		assay_curated_by == "Expert",
		assay_confidence_score,
		is.na(activity_comment),
		activity,
		desc(year),
		zinc_id) |>
		#desc(n_genes),
		#desc(purchasable_code),
		#desc(drug_code),
		#desc(biological_code)) |>
	dplyr::distinct(uniprot_accn, chembl_smiles, .keep_all = TRUE) |>
	dplyr::distinct(uniprot_accn, zinc_id, .keep_all = TRUE)



##### Save data ####

#full_chembl_data |> dplyr::write_tsv(paste0(staging_directory, "/data/full_chembl_data.tsv"))
save(
  "full_chembl_data",
  file = paste0(staging_directory, "/data/full_chembl_data.Rdata"))





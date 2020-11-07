# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(BioChemPantry)
library(SEAR)


staging_directory <- BioChemPantry::get_staging_directory("chembl23")
load(paste0(staging_directory, "/data/full_chembl_data.Rdata"))


chembl_sets <- full_chembl_data %>%
	dplyr::filter(active) %>%
	dplyr::transmute(
		target=uniprot_entry,
		name=gene_name,
		compound=zinc_id,
		affinity=5,
		description=target_description) %>%
	dplyr::mutate(
		name=ifelse(name %>% is.na, "", name),
		description=ifelse(description %>% is.na, "", description))
save("chembl_sets", file=paste0(staging_directory, "/data/chembl_sets.csv"))

# this takes ~6 days to run
chembl_compound_images <- full_chembl_data %>%
	dplyr::filter(active) %>%
	dplyr::distinct(zinc_id, .keep_all=T) %>%
	dplyr::select(smiles, zinc_id) %>%
	SEAR::generate_images(
		smiles="smiles",
		compound="zinc_id",
		batch_size=256,
		verbose=TRUE)
save("chembl_compound_images", file=paste0(staging_directory, "/data/chembl_compound_images.Rdata"))
load(file=paste0(staging_directory, "/data/chembl_compound_images.Rdata"))


# smiles file with data
chembl_smi <- full_chembl_data %>%
	dplyr::filter(active) %>%
	dplyr::distinct(zinc_id, .keep_all=T) %>%
	dplyr::left_join(
		chembl_compound_images %>%
			dplyr::select(
				zinc_id,
				image),
		by=c("zinc_id")) %>%
	plyr::adply(1, function(row){
		data <- row %>% with(list(
			preferred_name=preferred_name,
			physical_properties = list(
				molecular_weight=molecular_weight,
				ALogP=ALogP,
				rotatable_bonds=rotatable_bonds,
				reactivity=reactivity,
				polar_surface_area=PSA,
				n_hydrogen_bond_acceptors=HBA,
				n_hydrogen_bond_donors=HBD,
				pKa_range = c(most_basic_pKa, most_acidic_pKa),
				molecular_species=molecular_species,
				n_heavy_atoms=n_heavy_atoms),
			closest_aggregator = list(
				zinc_id=aggregator_zinc_id,
				tc=tc_to_aggregator),
			gene_names=gene_names,
			n_genes=n_genes,
			purchasable_code=purchasable_code,
			purchasable_level=purchasable_level,
			drug_code=drug_code,
			drug_level=drug_level,
			biological_code=biological_code,
			biological_level=biological_level,
			image=image) %>%
			jsonlite::toJSON(auto_unbox=TRUE) %>%
			as.character())
		row %>% with(
			data.frame(
				zinc_id=zinc_id,
				smiles=smiles,
				data=data))
	}) %>%
	dplyr::select(zinc_id, smiles, data)

save("chembl_smi", file=paste0(staging_directory, "/data/chembl.smi"))
load(file=paste0(staging_directory, "/data/chembl.smi"))
###############

make_library <- function(
	library_tag,
	fingerprint_type,
	fpcore_config,
	description,
	verbose=TRUE){
	library_fname <- paste0(staging_directory, "/data/chembl23_", library_tag, ".sea")
	pack.library(
		molecules=chembl_smi,
		targets=chembl_sets,
		library_fname=library_fname,
		fingerprint_type=fingerprint_type,
		name=description,
		config_files=list(fpcore=fpcore_config),
		write.smi_args=list(compound="zinc_id", smiles="smiles", data="data", verbose=TRUE),
		verbose=T)

	### Fit the library
	SEAR::plot_fits.library(library_fname, verbose=T)
	# check plot ...
	SEAR::set_fit.library(library_fname, threshold=.28, verbose=T)
}

make_library(
	library_tag="rdkit_ap",
	fingerprint_type="rdkit_ap",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with RDKit Atom Pair (AP) fingerprints")

make_library(
	library_tag="rdkit_avalon",
	fingerprint_type="rdkit_avalon",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with Avalon Toolkit fingerprints")

make_library(
	library_tag="rdkit_ecfc",
	fingerprint_type="rdkit_ecfc",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint (Path-based, Counts) fingerprints")

make_library(
	library_tag="rdkit_fcfc",
	fingerprint_type="rdkit_fcfc",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint (Feature-based, Counts) fingerprints")

make_library(
	library_tag="rdkit_fcfp4",
	fingerprint_type="rdkit_fcfp",
	fpcore_config="[rdkit_fcfp]\ncircle_radius=2\n",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint with diameter 4 (Feature-based, Hashed) fingerprints")

make_library(
	library_tag="rdkit_fcfp6",
	fingerprint_type="rdkit_fcfp",
	fpcore_config="[rdkit_fcfp]\ncircle_radius=3\n",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint with diameter 6 (Feature-based, Hashed) fingerprints")

make_library(
	library_tag="rdkit_ecfp4",
	fingerprint_type="rdkit_ecfp",
	fpcore_config="[rdkit_ecfp]\ncircle_radius=2\n",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint with diameter 4 (Hashed) fingerprints")

make_library(
	library_tag="rdkit_ecfp6",
	fingerprint_type="rdkit_ecfp",
	fpcore_config="[rdkit_ecfp]\ncircle_radius=3\n",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint with diameter 6 (Hashed) fingerprints")

make_library(
	library_tag="rdkit_hashap",
	fingerprint_type="rdkit_hashap",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with RDKit Atom Pair (AP, Hashed) fingerprints")

make_library(
	library_tag="rdkit_hashtt",
	fingerprint_type="rdkit_hashtt",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with RDKit Topological Torsion (Hashed) fingerprints")

make_library(
	library_tag="rdkit_maccs",
	fingerprint_type="rdkit_maccs",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with Molecular ACCess System keys")

make_library(
	library_tag="rdkit_path",
	fingerprint_type="rdkit_path",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with RDKit Fingerprints (Path-based, Hashed)")

make_library(
	library_tag="rdkit_tt",
	fingerprint_type="rdkit_tt",
	fpcore_config="",
	description="ChEMBL23 by Uniprot Entry to Zinc ID with RDKit Topological Torsion Fingerprints")











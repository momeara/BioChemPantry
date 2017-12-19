# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(BioChemPantry)
library(SEAR)


staging_directory <- BioChemPantry::get_staging_directory("chembl21")
load(paste0(staging_directory, "/data/full_chembl_data.Rdata"))


chembl_sets <- full_chembl_data %>%
	dplyr::filter(active) %>%
	dplyr::transmute(
		target=uniprot_entry,
		name=gene_name,
		compound=zinc_id,
		affinity=5,
		description=description) %>%
	dplyr::mutate(
		name=ifelse(name %>% is.na, "", name),
		description=ifelse(description %>% is.na, "", description))

chembl_smi <- full_chembl_data %>%
	dplyr::filter(active) %>%
	dplyr::distinct(zinc_id, .keep_all=T) %>%
	plyr::adply(1, function(row){
		row %>% with(
			data.frame(
				zinc_id=zinc_id,
				smiles=smiles,
				data=list(
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
					biological_level=biological_level) %>%
						jsonlite::toJSON(auto_unbox=TRUE) %>%
						as.character))})


###############

make_library <- function(fingerprint_type, fpcore_config, name, verbose=TRUE){
	library_fname <- paste0(staging_directory, "/data/chembl21_", fingerprint_type, ".sea")
	SEAR::pack.library(
		molecules=chembl_smi,
		targets=chembl_sets,
		library_fname=library_fname,
		fingerprint_type=fingerprint_type,
		name=name,
		config_files=list(fpcore=fpcore_config),
		verbose=T)

	### Fit the library
	SEAR::plot_fits.library(library_fname, verbose=T)
	# check plot ...
	SEAR::set_fit.library(library_fname, threshold=.28, verbose=T)
}

make_library(
	fingerprint_type="rdkit_ap",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with RDKit Atom Pair (AP) fingerprints")

make_library(
	fingerprint_type="rdkit_avalon",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with Avalon Toolkit fingerprints")

make_library(
	fingerprint_type="rdkit_ecfc",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint (Path-based, Counts) fingerprints")

make_library(
	fingerprint_type="rdkit_fcfc",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint (Feature-based, Counts) fingerprints")

make_library(
	fingerprint_type="rdkit_fcfp4",
	fpcore_config="[rdkit_ecfp]\ncircle_radius=2\n",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint with diameter 4 (Feature-based, Hashed) fingerprints")

make_library(
	fingerprint_type="rdkit_fcfp6",
	fpcore_config="[rdkit_ecfp]\ncircle_radius=3\n",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with Extended Connectivity Fingerprint with diameter 6 (Feature-based, Hashed) fingerprints")

make_library(
	fingerprint_type="rdkit_hashap",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with RDKit Atom Pair (AP, Hashed) fingerprints")

make_library(
	fingerprint_type="rdkit_hashtt",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with RDKit Topological Torsion (Hashed) fingerprints")

make_library(
	fingerprint_type="rdkit_maccs",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with Molecular ACCess System keys")

make_library(
	fingerprint_type="rdkit_path",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with RDKit Fingerprints (Path-based, Hashed)")

make_library(
	fingerprint_type="rdkit_tt",
	fpcore_config="",
	name="ChEMBL21 by Uniprot Entry to Zinc ID with RDKit Topological Torsion Fingerprints")











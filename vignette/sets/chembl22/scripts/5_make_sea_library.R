# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)

source("~/work/sea/scripts/sea.R")

staging_directory <- get_staging_directory("chembl21")


load(paste0(staging_directory, "/data/full_chembl_data.Rdata"))

library_fname <-paste0(staging_directory, "/data/chembl21.sea")
fingerprint_type <- "rdkit_ecfp"
name <- "ChEMBL21_by_uniprot_entry_to_zinc_id"


chembl_sets <- full_chembl_data %>%
	filter(active) %>%
	dplyr::transmute(
		target=uniprot_entry,
		name=gene_name,
		compound=zinc_id,
		affinity=5,
		description=target_description) %>%
	mutate(
		name=ifelse(name %>% is.na, "", name),
		description=ifelse(description %>% is.na, "", description))

chembl_smi <- full_chembl_data %>%
	filter(active) %>%
	distinct(zinc_id, .keep_all=T) %>%
	dplyr::select(
		compound=zinc_id,
		smiles=chembl_smiles)

pack.library(
	molecules=chembl_smi,
	targets=chembl_sets,
	library_fname=library_fname,
	fingerprint_type=fingerprint_type,
	name=name,
	verbose=T)


#### Check what targets/compounds actually got read into the library
lib <- unpack.library(library_fname, verbose=T)

missing_sets <- chembl_sets %>% anti_join(lib$targets, by=c("target", "compound"))


test_sets <- chembl_sets %>% filter(compound == "ZINC000169746649")
write.set(test_sets, paste0("/scratch/momeara/test_set.set")
test_smi <- chembl_smi %>% filter(compound == "ZINC000169746649")
write.smi(test_smi, "/scratch/momeara/test.csv")
pack.library(
	molecules=test_smi,
	targets=test_sets,
	library_fname="/scratch/momeara/test_library.sea",
	fingerprint_type=fingerprint_type,
	name=name,
	verbose=T)



chembl_targets <- target_info %>%
	semi_join(lib$targets, by=c("uniprot_entry" = "target"))

chembl_targets %>% write_tsv(paste0(staging_directory, "/data/chembl_targets.tsv"))
save("chembl_targets", file=paste0(staging_directory, "/data/chembl_targets.Rdata"))

chembl_compounds <- full_chembl_data %>%
	distinct(zinc_id, .keep_all=T) %>%
	semi_join(lib$molecules, by=c("zinc_id" = "compound")) %>%
	dplyr::select(
		zinc_id,
		alt_zinc_ids,
		chembl_id,
		chembl_smiles,
		MaxTC_aggregator,
		molecular_weight,
		ALogP,
		HBA,
		HBD,
		PSA,
		RTB,
		most_acidic_pKa,
		most_basic_pKa,
		molecular_species,
		n_heavy_atoms)

chembl_compounds %>% write_tsv(paste0(staging_directory, "/data/chembl_compounds.tsv"))
save("chembl_compounds", file=paste0(staging_directory, "/data/chembl_compounds.Rdata"))


chembl_activities <- full_chembl_data %>%
	distinct(uniprot_accn, chembl_id, .keep_all=T) %>%
	semi_join(lib$targets, by=c("uniprot_entry" = "target", "zinc_id" = "compound")) %>%
	dplyr::select(
		assay_id,
		tid,
		uniprot_accn,
		uniprot_entry,
		gene_name,
		zinc_id,
		chembl_id,
		chembl_smiles,
		activity,
		standard_type,
		assay_description,
		assay_test_type,
		assay_category,
		assay_organism,
		assay_relationship_type,
		assay_confidence_score,
		doc_id,
		journal,
		year,
		doi,
		pubmed_id)

chembl_activities %>% write_tsv(paste0(staging_directory, "/data/chembl_activities.tsv"))
save("chembl_activities", file=paste0(staging_directory, "/data/chembl_activities.Rdata"))



### Fit the library
plot_fits.library(library_fname, verbose=T)
# check plot ...
set_fit.library(library_fname, threshold=.28, verbose=T)





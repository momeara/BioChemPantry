# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(data.table)
library(stringr)
library(sqldf)

source("~/work/sets/uniprot_idmapping/scripts/uniprot_id_map.R")
source("~/work/sea/scripts/sea.R")

#########################
guide_to_pharmacology_targets_table_fname <- "data/targets_and_families_140708.csv"
guide_to_pharmacology_interactions_fname <- "data/interactions_140708.csv"
guide_to_pharmacology_ligands_fname <- "data/ligands_140708.csv"

targets <- read.table(guide_to_pharmacology_targets_table_fname, header=T, sep=",", quote="\"") %>%
	select(
			target_id=Target.id,
			family=Type,
			gene_symbol=Human.Entrez.Gene,
			uniprot_entry=Human.SwissProt)

interactions <- read.table(guide_to_pharmacology_interactions_fname, header=T, sep=",", quote="\"") %>%
	mutate(target = str_replace_all(target, "<[^>]+>", "")) %>%
	mutate(target = str_replace_all(target, "&alpha;", "a")) %>%
	mutate(target = str_replace_all(target, "&beta;", "b")) %>%
	mutate(target = str_replace_all(target, "&gamma;", "g")) %>%
	mutate(target = str_replace_all(target, "&deta;", "d")) %>%
	mutate(ligand = str_replace_all(ligand, " +$", "")) %>%
	mutate(ligand = str_replace_all(ligand, "<[^>]+>", "")) %>%
	mutate(ligand = str_replace_all(ligand, "&alpha;", "a")) %>%
	mutate(ligand = str_replace_all(ligand, "&beta;", "b")) %>%
	mutate(ligand = str_replace_all(ligand, "&gamma;", "g")) %>%
	mutate(ligand = str_replace_all(ligand, "&delta;", "d")) %>%
	mutate(ligand = str_replace_all(ligand, "&Delta;", "d")) %>%
	mutate(ligand = str_replace_all(ligand, "&plusmn;", "+/-")) %>%
	mutate(ligand = str_replace_all(ligand, " +$", ""))

ligands <- read.table(guide_to_pharmacology_ligands_fname, header=T, sep=",", quote="\"")

ligands2 <- ligands %>%
	select(
		ligand_id=Ligand.id,
		ligand_type=Type,
		smiles=SMILES)


bridging_metabolites <- interactions %>%
	left_join(targets, by="target_id") %>%
	left_join(ligands2, by="ligand_id") %>%
	filter(endogenous=="t") %>%
	filter(ligand_type != "Inorganic", ligand_type != "Peptide") %>%
	filter(
		ligand != "d9-tetrahydrocannabinol",
		ligand != "PGD3",
		ligand != "abnormal cannabidiol") %>%
	group_by(target_id, ligand_id) %>% filter(row_number(target_id) == 1) %>% ungroup() %>%
	select(ligand_id, ligand, smiles, target_id, target, uniprot_entry, target_family=family) %>%
	arrange(ligand, target_family, target) %>%
	as.data.frame()

write.table(bridging_metabolites, "~/work/phenologs/sea_ChEMBL18_HMDB_131030/product/known_metabolite_activities_IUPHAR_140902.xls",
						sep="\t", row.names=F, quote=F)

# example: make a sea set from interaction data
GPR55 <- interactions %>%
	filter(target=="GPR55") %>%
	left_join(ligands2, by="ligand_id") %>%
	left_join(targets, by="target_id") %>%
	group_by(target_id, ligand_id) %>%
	filter(row_number(ligand_id) == 1) %>%
	ungroup() %>%
	select(target=target_id, name=target, family, compound=ligand_id, ligand, endogenous, uniprot_entry, ligand_type, smiles) %>%	
	as.data.frame()
write.smi(GPR55, "~/work/phenologs/sea_GPR55_CB2_140708/data/GPR55_IUPHAR.smi")
write.set(GPR55, "~/work/phenologs/sea_GPR55_CB2_140708/data/GPR55_IUPHAR.set")



cross_gpcr <- bridging_metabolites %>%
	filter(family == "gpcr") %>%
	group_by(target_id, ligand_id) %>%
	filter(row_number(target) == 1) %>%
	ungroup() %>%
	group_by(ligand_id) %>%
	filter(n() != 1) %>%
	left_join(interactions2)


cross_gpcr %>%
	as.data.frame() %>%
	arrange(ligand) %>%
	head(200)

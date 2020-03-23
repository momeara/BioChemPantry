
library(plyr)
library(tidyverse)
library(BioChemPantry)

# new filtered set
load("intermediate_data/nCov_host_ppi_stringent_200318.Rdata")
load("intermediate_data/nCov_host_ppi_filtered_200318.Rdata")

####################
# Load IUPHAR Data #
####################

pantry <- BioChemPantry::get_pantry("guide_to_pharmacology_2020_03")
iuphar_targets <- pantry %>%
    dplyr::tbl("targets") %>%
    dplyr::collect(n=Inf)

iuphar_ligands <- pantry %>%
    dplyr::tbl("ligands") %>%
    dplyr::collect(n=Inf)

iuphar_interactions <- pantry %>%
    dplyr::tbl("interactions") %>%
    dplyr::collect(n=Inf)

################
# Merge IUPHAR #
################

proteomic_iuphar_ligands_filtered_200318 <- nCov_host_ppi_filtered_200318 %>%
    dplyr::select(Bait, Preys, PreyGene, uniprot_entry, protein_name) %>%
    dplyr::left_join(
        iuphar_targets %>%
            dplyr::select(
                target_id,
                uniprot_entry),
        by=c("uniprot_entry"="uniprot_entry")) %>%
    dplyr::filter(!is.na(target_id)) %>%
    dplyr::inner_join(
        iuphar_interactions %>%
           dplyr::select(
               target_id, ligand_id,
               compound_name=name,
               zinc_id,
               smiles,
               primary_target,
               assay_description,
               interaction_type,
               assay_units = original_affinity_units,
               affinity_nM = original_affinity_median_nm,
               pubmed_id),
        by=c("target_id")) %>%
    dplyr::left_join(
        iuphar_ligands %>%
            dplyr::select(
                ligand_id,
                ligand_type,
                approved,
                radioactive),
        by=c("ligand_id")) %>%
    dplyr::distinct(uniprot_entry, compound_name, .keep_all=TRUE) %>%
    dplyr::mutate(
        is_stringent_prey = PreyGene %in% nCov_host_ppi_stringent_200318$PreyGene)

save(proteomic_iuphar_ligands_filtered_200318, file="intermediate_data/proteomic_iuphar_ligands_filtered_200318.Rdata")

proteomic_iuphar_ligands_final <- proteomic_iuphar_ligands_filtered_200318 %>%
    dplyr::mutate(Bait = dplyr::case_when(
        Bait == "SARS-CoV2 protein14" ~ "SARS-CoV-2 Orf9c",
        Bait == "SARS-CoV2 orf8"      ~ "SARS-CoV-2 Orf8",  
        Bait == "SARS-CoV2 nsp6"      ~ "SARS-CoV-2 Nsp6",   
        Bait == "SARS-CoV2 E"         ~ "SARS-CoV-2 E",   
        Bait == "SARS-CoV2 nsp13"     ~ "SARS-CoV-2 Nsp13",
        Bait == "SARS-CoV2 nsp7"      ~ "SARS-CoV-2 Nsp7",  
        Bait == "SARS-CoV2 M"         ~ "SARS-CoV-2 M",   
        Bait == "SARS-CoV2 N"         ~ "SARS-CoV-2 N",      
        Bait == "SARS-CoV2 orf9b"     ~ "SARS-CoV-2 Nrf9b",    
        Bait == "SARS-CoV2 nsp9"      ~ "SARS-CoV-2 Nsp9",  
        Bait == "SARS-CoV2 nsp10"     ~ "SARS-CoV-2 Nsp10",   
        Bait == "SARS-CoV2 nsp2"      ~ "SARS-CoV-2 Nsp2",  
        Bait == "SARS-CoV2 nsp14"     ~ "SARS-CoV-2 Nsp14",   
        Bait == "SARS-CoV2 nsp5"      ~ "SARS-CoV-2 Nsp5",  
        Bait == "SARS-CoV2 orf3a"     ~ "SARS-CoV-2 Orf3a",  
        Bait == "SARS-CoV2 nsp4"      ~ "SARS-CoV-2 Nsp4",  
        Bait == "SARS-CoV2 nsp12"     ~ "SARS-CoV-2 Nsp12",  
        Bait == "SARS-CoV2 nsp8"      ~ "SARS-CoV-2 Nsp8",  
        Bait == "SARS-CoV2 nsp1"      ~ "SARS-CoV-2 Nsp1",   
        Bait == "SARS-CoV2 orf10"     ~ "SARS-CoV-2 Orf10")) %>%    
    dplyr::rename(
        bait=Bait,
        prey_uniprot_accn = Preys,
        prey_uniprot_entry = uniprot_entry,
        prey_gene_name = PreyGene,
        prey_protein_name = protein_name)

proteomi_chembl_ligands_final %>%
    readr::write_tsv("product/proteomic_iuphar_ligands_final.tsv")

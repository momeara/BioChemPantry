
library(plyr)
library(tidyverse)
library(BioChemPantry)
library(xlsx)

load("intermediate_data/nCov_host_ppi_stringent_200318.Rdata")
load("intermediate_data/nCov_host_ppi_filtered_200318.Rdata")

pantry <- BioChemPantry::get_pantry("sea_chembl25")
chembl_targets <- pantry %>%
    dplyr::tbl("targets") %>%
    dplyr::collect(n=Inf)

chembl_compounds <- pantry %>%
    dplyr::tbl("compounds") %>%
    dplyr::collect(n=Inf)

chembl_activities <- pantry %>%
    dplyr::tbl("activities") %>%
    dplyr::filter(!is.na(gene_name)) %>%
    dplyr::filter(active) %>%
    dplyr::collect(n=Inf)

chembl_targets_summary <- chembl_activities %>%
    dplyr::filter(!is.na(gene_name), active) %>%
    dplyr::distinct(uniprot_accn, uniprot_entry, gene_name, zinc_id) %>%
    dplyr::count(uniprot_accn, uniprot_entry, gene_name) %>%
    chembl_targets %>% dplyr::rename(n_active_compound = n)        
chembl_targets_summary %>%
    readr::write_tsv("product/chembl25_targets.tsv")

proteomic_chembl_ligands_200318 <- nCov_host_ppi_filtered_200318 %>%
    dplyr::select(Bait, Preys, PreyGene,  prey_uniprot_entry = uniprot_entry, protein_name) %>%
    dplyr::left_join(
        chembl_activities %>%
            dplyr::select(
                uniprot_accn,
                gene_name,
                target_type,
                assay_description,
                standard_type,
                standard_units,
                standard_value,
                doi,
                pubmed_id,
                chembl_id,
                zinc_id,
                preferred_name,
                smiles,
                n_genes,
                purchasable_level,
                drug_code,
                biological_level,
                aggregator,
                tc_to_aggregator,
                active),
        by=c("PreyGene"="gene_name")) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
        is_stringent_prey = PreyGene %in% nCov_host_ppi_stringent_200318$PreyGene)

proteomic_chembl_ligands_200318 %>%
    readr::write_tsv("intermediate_data/proteomic_chembl_ligands_200318.tsv")

proteomic_chembl_ligands_200318_top <- proteomic_chembl_ligands_200318 %>%
    dplyr::arrange(PreyGene, standard_value) %>%
    plyr::ddply("PreyGene", function(df){
        cat("PreyGene: ", df$PreyGene[1], " nrow: ", nrow(df), "\n")
        z <- df %>% head(min(nrow(df), 100))
        cat("   now nrow: ", nrow(z), "\n")
        z
    })
proteomic_chembl_ligands_200318_top %>%
    readr::write_tsv("intermediate_data/proteomic_chembl_ligands_200318_top.tsv")



proteomi_chembl_ligands_final <- proteomi_chembl_ligands_200318 %>%
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
        prey_gene_name=PreyGene,
        prey_uniprot_accn=uniprot_accn,
        prey_protein_name=protein_name) %>%
    dplyr::select(-Preys)

proteomi_chembl_ligands_final %>%
    readr::write_tsv("product/proteomic_chembl_ligands_final.tsv")

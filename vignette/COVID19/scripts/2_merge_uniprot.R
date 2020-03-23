
library(plyr)
library(tidyverse)
library(BioChemPantry)

source("scripts/uniprot.R")



load("intermediate_data/nCov_host_ppi_filtered_200317.Rdata")
user_agent_arg <- httr::user_agent("httr <fill in your email here>")

pantry <- BioChemPantry::get_pantry()

nCov_host_ppi_filtered_200317 <- pantry %>%
    dplyr::copy_to(
        nCov_host_ppi_filtered_200317,
        "nCov_host_ppi_filtered_200317",
        temporary=TRUE, overwrite=TRUE) %>%
    dplyr::left_join(
        pantry %>%
            BioChemPantry::schema_tbl("hgnc_171217.genes") %>%
                dplyr::select(uniprot_accn, uniprot_entry, PreyGene=symbol, protein_name=name, gene_family),
        by=c("PreyGene")) %>%
    dplyr::collect(n=Inf)

# there a few that didn't map, so do them manually
nCov_host_ppi_filtered_200317z <- nCov_host_ppi_filtered_200317 %>%
    dplyr::mutate(
        PreyGene = ifelse(Preys == "Q5VT66", "MTARC1", PreyGene),
        uniprot_accn = ifelse(
            PreyGene %in% c("MTARC1", "POGLUT3", "POGLUT2", "DIPK1B"),
            Preys, uniprot_accn),
        uniprot_entry = ifelse(
            PreyGene %in% c("MTARC1", "POGLUT3", "POGLUT2", "DIPK1B"),
            c("MARC1_HUMAN", "PLGT3_HUMAN", "PLGT2_HUMAN", "DIK1B_HUMAN"),
            uniprot_accn),
        protein_name = ifelse(
            PreyGene %in% c("MTARC1", "POGLUT3", "POGLUT2", "DIPK1B"),
            c("Mitochondrial amidoxime-reducing component 1",
              "Protein O-glucosyltransferase 3",
              "Protein O-glucosyltransferase 2",
              "Divergent protein kinase domain 1B"),
            protein_name))        
          


proteomic_uniprot_targets <- nCov_host_ppi_filtered_200315 %>%
    magrittr::extract2("uniprot_entry") %>%
    uniprot_entry_web_lookup(
        columns=c(
            "id",
            "entry name",
            "comment(SUBUNIT)",
            "families",
            "interactor",
            "comment(DEVELOPMENTAL STAGE)",
            "comment(INDUCTION)",
            "comment(TISSUE SPECIFICITY)",
            "comment(FUNCTION)",
            "comment(PATHWAY)",
            "go(biological process)",
            "go(molecular function)",
            "go(cellular component)",
            "comment(DISEASE)",
            "comment(DISRUPTION PHENOTYPE)",
            "comment(SUBCELLULAR LOCATION)",
            "feature(INTRAMEMBRANE)",
            "comment(DOMAIN)",
            "3d"),
        user_agent_arg=user_agent_arg,
        verbose=TRUE)



proteomic_uniprot_targets <- nCov_host_ppi_filtered_200315 %>%
    dplyr::left_join(
        proteomic_uniprot_targets %>%
            dplyr::mutate(uniprot_accn=Entry) %>%
            dplyr::select(-uniprot_entry, -`Entry name`, -`Entry`),
        by=c("Preys"="uniprot_accn"))

proteomic_uniprot_targets %>%
    readr::write_tsv("product/proteomic_uniprot_targets_200315.tsv")




library(plyr)
library(tidyverse)
library(BioChemPantry)
library(readxl)


nCov_host_ppi_stringent_200318 <- readxl::read_excel(
    "raw_data/2020-03-18_Krogan_SARSCoV2_27baits.xlsx")

nCov_host_ppi_filtered_200318 <- readxl::read_excel(
    "raw_data/2020-03-18_Krogan_SARSCoV2_27baits_LowThreshold.xlsx")

pantry <- BioChemPantry::get_pantry("hgnc_171217")
nCov_host_ppi_stringent_200318 <- pantry %>%
    dplyr::copy_to(
        nCov_host_ppi_stringent_200318,
        "nCov_host_ppi_filtered_200318",
        temporary=TRUE, overwrite=TRUE) %>%
    dplyr::left_join(
        pantry %>%
            dplyr::tbl("genes") %>%
                dplyr::select(
                    uniprot_entry,
                    PreyGene=symbol,
                    protein_name=name,
                    gene_family),
        by=c("PreyGene")) %>%
    dplyr::collect(n=Inf)
save(nCov_host_ppi_stringent_200318, file="intermediate_data/nCov_host_ppi_stringent_200318.Rdata")

pantry <- BioChemPantry::get_pantry("hgnc_171217")
nCov_host_ppi_filtered_200318 <- pantry %>%
    dplyr::copy_to(
        nCov_host_ppi_filtered_200318,
        "nCov_host_ppi_filtered_200318",
        temporary=TRUE, overwrite=TRUE) %>%
    dplyr::left_join(
        pantry %>%
            dplyr::tbl("genes") %>%
                dplyr::select(
                    uniprot_entry,
                    PreyGene=symbol,
                    protein_name=name,
                    gene_family),
        by=c("PreyGene")) %>%
    dplyr::collect(n=Inf)
save(nCov_host_ppi_filtered_200318, file="intermediate_data/nCov_host_ppi_filtered_200318.Rdata")




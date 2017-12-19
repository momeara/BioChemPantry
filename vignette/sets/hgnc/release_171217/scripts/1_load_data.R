# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(jsonlite)
library(readr)
library(httr)
library(BioChemPantry)
library(assertthat)
source("~/work/sea/scripts/id_utils.R")

source("~/work/sets/uniprot_idmapping/scripts/uniprot_id_map_web.R")

# collect genes with gene products from HGNC
# Each gene symbol corresponds to at most one uniprot accn
# but some uniprot entries correspond to multiple gene ids most of these are pretty weird or uncharacterized gene products, so filter all these out

staging_directory <- BioChemPantry::get_staging_directory("hgnc/release_171217")


system(paste0("
cd ", staging_directory, "
wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/locus_types/gene_with_protein_product.json"))

json_fname <- paste0(staging_directory, "/gene_with_protein_product.json")
assertthat::assert_that(file.exists(json_fname), msg="The HGNC data was downloaded correctly.")
genes <- jsonlite::fromJSON(txt=json_fname)$response$docs

genes %>% BioChemPantry::summarize_map("hgnc_id", "symbol")
# hgnc_id, gene_symbol are 1-1


# work around for bug https://github.com/hadley/tidyr/issues/141
my_unnest <- function(data, col_name){
	data %>% adply(1, function(row){
		if(is.null(row[,col_name][[1]])){
			row[,col_name] <- NA
		} else {
			unnested <- row[,col_name][[1]]
			row <- row[rep(1,length(unnested)),]
			row[,col_name] <- unnested
		}
		row
	})
}

genes2 <- genes %>%
	my_unnest("uniprot_ids") %>%
	rename(uniprot_accn = uniprot_ids)


genes %>% summarize_map("hgnc_id", "symbol")
# hgnc_id, gene_symbol are 1-1

genes %>% summarize_map("hgnc_id", "uniprot_accn")
#X<-[hgnc_id]:
#  |X|: 18989
#  count*size: 18989*1
#Y<-[uniprot_accn]:
#  |Y|: 18834 (30 NA)
#  count*size: 18765*1, 53*2, 6*3, 3*4, 3*5, 1*7, 1*10, 1*12, 1*14
#[X U Y]:
#  |X U Y|: 18959 (30 NA)
#  count*size: 18959*1
#[X @ Y]:
#  |X ~ Y|: 18765
#  |X:X < Y|, |Y:X < Y|: 0, 0
#  |X:X > Y|, |Y:X > Y|: 194, 69

#just get rid of them,
genes <- genes %>%
	dplyr::anti_join(
		genes %>%
			count(uniprot_accn) %>%
			filter(n>1),
		by="uniprot_accn")


genes <- genes %>%
	mutate(iuphar = iuphar %>% str_replace("objectId:","") %>% as.numeric) %>%
	dplyr::select(
		symbol,
		gene_family,
		gene_family_id,
		name,
		uniprot_accn,
		entrez_id,
		ucsc_id,
		ensembl_gene_id,
		refseq_accession,
		omim_id,
		iuphar)


uniprot_data <- httr::GET(
	url="http://www.uniprot.org/uniprot/",
	httr::user_agent("httr http::/github.com/momeara/BioChemPantry"),
	query=list(
		query="database:(type:hgnc)",
		format='tab',
		compress='yes')) %>%
	httr::content() %>%
	rawConnection() %>%
	gzcon() %>%
	readr::read_tsv()

uniprot_data <- uniprot_data %>% rename(
	uniprot_accn = `Entry`,
	uniprot_entry = `Entry name`,
	protein_names = `Protein names`,
	status = `Status`,
	uniprot_gene_names = `Gene names`,
	organism = `Organism`,
	length = `Length`)



genes2 <- genes %>%
	dplyr::left_join(
		uniprot_data,
		by="uniprot_accn")


pantry <- get_pantry()
pantry %>% create_schema("hgnc_171217")
pantry %>% set_schema("hgnc_171217")

pantry %>% dplyr::copy_to(
	genes2,
	"genes",
	temporary=F,
	indices=list(
		"symbol",
		"uniprot_accn",
		"uniprot_entry"))

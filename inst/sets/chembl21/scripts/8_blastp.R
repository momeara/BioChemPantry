# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(Bethany)
library(seqinr)
source("~/work/sea/scripts/data_repo.R")

data_repo <- get_data_repo("sea_chembl21")

chembl_target_sequences <- data_repo %>%
	dplyr::tbl("targets") %>%
	dplyr::select(uniprot_entry, sequence) %>%
	dplyr::collect(n=Inf)

chembl_sequences <- chembl_target_sequences %>%
	plyr::alply(1, function(df){
		seqinr::as.SeqFastaAA(df$sequence[1], name=df$uniprot_entry[1])
	})
names(chembl_sequences) = chembl_target_sequences$uniprot_entry

# compute target_vs_target blast scores
scores <- Bethany::blastp(
	ref=chembl_sequences,
	query=chembl_sequences,
	run_id="blastp_chembl21_vs_chembl21",
	verbose=T) %>%
	rename(
		target1 = ref_target,
		target2 = query_target) %>%
	arrange(EValue) %>%
	distinct(target1, target2, .keep_all=T)

data_repo %>% copy_to(
	df=scores,
	name="blastp_target_vs_target",
	temporary=F,
	indexes=list(
		"target1",
		"target2",
		c("target1", "target2")),
	fast=T)


# compute target_vs_human blast scores
scores_human <- blastp(
	ref=paste0(getwd(), "/../gencode24/dump/gencode.v24.pc_translations.fa"),
	query=chembl_sequences,
	run_id="blastp_chembl21_vs_gencode_v24_human",
	verbose=T) %>%
	rename(
		chembl_target_uniprot_entry = ref_target,
		gencode_ids = query_target)
gencode_ids <- str_split_fixed(scores_human$gencode_ids, "[|]", 8)

# work around for "incorrect number of dimensions" bug?
b <- gencode_ids

scores_human <- scores_human %>%
	dplyr::mutate(
		gencode_ensembl_protein_id =  b[,1],
		gencode_ensembl_transcript_id = b[,2],
		gencode_ensembl_gene_id = b[,3],
		gencode_hgnc_gene_symbol = b[,7],
		gencode_entrez_gene_id = as.integer(b[,8])) %>%
	dplyr::select(
		chembl_target_uniprot_entry,
		gencode_ensembl_protein_id,
		gencode_ensembl_transcript_id,
		gencode_ensembl_gene_id,
		gencode_hgnc_gene_symbol,
		gencode_entrez_gene_id,
		bit_score,
		EValue)

scores_human <- scores_human %>%
	arrange(EValue) %>%
	distinct(chembl_target_uniprot_entry, gencode_hgnc_gene_symbol, .keep_all=T)



data_repo %>% copy_to(
	df=scores_human,
	name="blastp_target_vs_human",
	temporary=F,
	indexes=list(
		"chembl_target_uniprot_entry",
		"gencode_hgnc_gene_symbol",
		c("chembl_target_uniprot_entry", "gencode_hgnc_gene_symbol")),
	fast=T)

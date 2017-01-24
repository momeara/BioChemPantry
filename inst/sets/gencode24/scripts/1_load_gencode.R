# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(seqinr)

data_repo <- get_data_repo()
data_repo %>% create_schema("gencode24")
data_repo %>% set_schema("gencode24")

system("
mkdir dump
cd dump
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.pc_translations.fa.gz
gunzip gencode.v24.pc_translations.fa.gz
")

protein_sequences <- seqinr::read.fasta(
	"dump/gencode.v24.pc_translations.fa",
	seqtype="AA",
	as.string=T) %>%
	plyr::adply(1, function(df){
		name <- df %>% attr("name") %>% stringr::str_split_fixed("[|]", 8)
		data_frame(
		ensembl_protein_id =  name[1],
		ensembl_transcript_id = name[2],
		ensembl_gene_id = name[3],
		hgnc_gene_symbol = name[7],
		entrez_id = as.integer(name[8]),
		sequence = as.character(df))
	}) %>%
	dplyr::select(-X1)


data_repo %>% copy_to(
	df=protein_sequences,
	name="protein_sequences",
	temporary=F,
	indexes=list(
		"ensembl_protein_id",
		"entrez_id"))



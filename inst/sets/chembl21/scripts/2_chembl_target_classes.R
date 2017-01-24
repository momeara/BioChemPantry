# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(reshape2)
library(readr)
source("~/work/sea/scripts/data_repo.R")

# chembl stores protein classes as tree structured ontology this
# converts it into a flat format where the classes with higher numbers
# are more specific or repeated.

#  chembl_target_classes %>% str
#  Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	7063 obs. of  9 variables:
#   $ uniprot_id      : chr  "O09028" "P02708" "P04637" "P04757" ...
#   $ protein_class_id: int  1173 422 12 422 11 11 323 11 424 559 ...
#   $ short_name      : chr  "GABAA" "CHRN alpha" "Transcription Factor" "CHRN alpha" ...
#   $ class_1         : chr  "Ion channel" "Ion channel" "Transcription Factor" "Ion channel" ...
#   $ class_2         : chr  "LGIC" "LGIC" "Transcription Factor" "LGIC" ...
#   $ class_3         : chr  "GABAA" "ACH" "Transcription Factor" "ACH" ...
#   $ class_4         : chr  "GABAA" "CHRN alpha" "Transcription Factor" "CHRN alpha" ...
#   $ class_5         : chr  "GABAA" "CHRN alpha" "Transcription Factor" "CHRN alpha" ...
#   $ class_6         : chr  "GABAA" "CHRN alpha" "Transcription Factor" "CHRN alpha" ...


chembl_db <- get_data_repo("chembl21")

protein_class_ids <- chembl_db %>%
	tbl("protein_classification") %>%
	collect %>%
	mutate(parent_id = ifelse(parent_id == 0, protein_class_id, parent_id)) %>%
	select(
		protein_class_id,
		parent_id,
		short_name)

class_id_lookup <- 1:nrow(protein_class_ids)
names(class_id_lookup) <- protein_class_ids$protein_class_id

parents <-
	lapply(1:nrow(protein_class_ids), function(i){
		class_ids = protein_class_ids$protein_class_id[i]
		curr_id = class_ids[1]
		while(T){
			parent_index <- class_id_lookup[as.character(curr_id)]
			parent_id <- protein_class_ids$parent_id[parent_index]
			if(is.na(parent_id) | parent_id == curr_id){
				break
			} else {
				class_ids <- c(parent_id, class_ids)
				curr_id <- parent_id
			}
		}
		class_ids
	})

max_depth <- parents %>% sapply(length) %>% max

# fill in hierarchy with terminal class and turn into a data.frame
parents <- parents %>%
	lapply(function(p) {
		c(
			p[length(p)],
			p,
			rep(
				p[length(p)], max_depth-length(p)))
	}) %>%
	unlist %>%
	matrix(ncol=max_depth+1, byrow=T) %>%
	as.data.frame

colnames(parents) <- c("protein_class_id", paste("class_", 1:max_depth, sep=""))

protein_class_ids <- protein_class_ids %>%
	select(-parent_id) %>%
	left_join(parents, by="protein_class_id")

# convert parent_ids to names
protein_class_ids <- protein_class_ids %>%
	melt(id.vars=c("protein_class_id", "short_name")) %>%
		rename(parent_class_level = variable, parent_class_id = value) %>%
	left_join(
		protein_class_ids %>% select(protein_class_id, parent_name = short_name),
		by=c("parent_class_id" = "protein_class_id")) %>%
	dcast(protein_class_id + short_name ~ parent_class_level, value.var="parent_name")

chembl_target_classes <-
	chembl_db %>% tbl("component_sequences") %>%
	select(
		component_id,
		uniprot_accn = accession) %>%
	left_join(
		chembl_db %>% tbl("component_class") %>%
		select(
			component_id,
			protein_class_id),
		by="component_id") %>%
	select(-component_id) %>%
	collect %>%
	left_join(protein_class_ids, by="protein_class_id")

chembl_target_classes %>%
	write_tsv(
		"data/chembl_target_classes.tsv")

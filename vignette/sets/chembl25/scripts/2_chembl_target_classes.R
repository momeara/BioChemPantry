# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(reshape2)
library(readr)
library(BioChemPantry)

# chembl stores protein classes as tree structured ontology this
# converts it into a flat format where the classes with higher numbers
# are more specific or repeated.

#  chembl_target_classes %>% str
#  Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	7063 obs. of  9 variables:
#   $ uniprot_id      : chr  "O09028" "P02708" "P04637" "P04757" ...
#   $ protein_class_id: int  1173 422 12 422 11 11 325 11 424 559 ...
#   $ short_name      : chr  "GABAA" "CHRN alpha" "Transcription Factor" "CHRN alpha" ...
#   $ class_1         : chr  "Ion channel" "Ion channel" "Transcription Factor" "Ion channel" ...
#   $ class_2         : chr  "LGIC" "LGIC" "Transcription Factor" "LGIC" ...
#   $ class_3         : chr  "GABAA" "ACH" "Transcription Factor" "ACH" ...
#   $ class_4         : chr  "GABAA" "CHRN alpha" "Transcription Factor" "CHRN alpha" ...
#   $ class_5         : chr  "GABAA" "CHRN alpha" "Transcription Factor" "CHRN alpha" ...
#   $ class_6         : chr  "GABAA" "CHRN alpha" "Transcription Factor" "CHRN alpha" ...


pantry <- get_pantry("chembl25")

staging_directory <- get_staging_directory("chembl25")

protein_class_ids <- pantry %>%
	dplyr::tbl("protein_classification") %>%
	dplyr::collect() %>%
	dplyr::mutate(parent_id = ifelse(parent_id == 0, protein_class_id, parent_id)) %>%
	dplyr::select(
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
	dplyr::select(-parent_id) %>%
	dplyr::left_join(parents, by="protein_class_id")

# convert parent_ids to names
protein_class_ids <- protein_class_ids %>%
	reshape2::melt(id.vars=c("protein_class_id", "short_name")) %>%
		rename(parent_class_level = variable, parent_class_id = value) %>%
	dplyr::left_join(
		protein_class_ids %>% select(protein_class_id, parent_name = short_name),
		by=c("parent_class_id" = "protein_class_id")) %>%
	reshape2::dcast(protein_class_id + short_name ~ parent_class_level, value.var="parent_name")

chembl_target_classes <-
	pantry %>% dplyr::tbl("component_sequences") %>%
	dplyr::select(
		component_id,
		uniprot_accn = accession) %>%
	dplyr::left_join(
		pantry %>% tbl("component_class") %>%
		select(
			component_id,
			protein_class_id),
		by="component_id") %>%
	dplyr::select(-component_id) %>%
	dplyr::collect() %>%
	dplyr::left_join(protein_class_ids, by="protein_class_id")

dir.create(paste0(staging_directory, "/data"), showWarnings=FALSE)
chembl_target_classes %>%
	readr::write_tsv(
		paste0(staging_directory, "/data/chembl_target_classes.tsv"))

# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(fdrtool)
source("~/work/sea/scripts/sea.R")

target_info <- read_tsv("data/target_info.tsv")

sea_library_fname <- paste0(getwd(), "/data/chembl21.sea")
run_base <- paste0(getwd(), "/data/runs/chembl21_vs_chembl21")
sea_wrapper_fname <- paste0(getwd(), "/scripts/sea-wrapper.sh")

scores_fname <- paste0(getwd(), "/data/chembl21.scores.csv")


####
sea_library <- unpack.library(
	library_fname, verbose=T)

sea_library <- sea_library$targets %>%
	left_join(
		sea_library$molecules,
		by=c("compound"))

#### prepare inputs

outputs_path <- paste0(run_base, "/outputs")
unlink(outputs_path, recursive=T, force=T)
dir.create(outputs_path)

logs_path <-paste0(run_base, "/logs")
unlink(logs_path, recursive=T, force=T)
dir.create(logs_path)

inputs_path <- paste0(run_base, "/inputs")
unlink(inputs_path, recursive=T, force=T)
dir.create(inputs_path, recursive=T)
sea_library %>%
	d_ply(c("target"), function(df){
		target <- df$target[1]
		target_fname <- paste0(inputs_path, "/", target, ".csv")
		df %>%
			select(compound, smiles) %>%
			write.table(
				target_fname,
				quote=F,
				sep=",",
				row.names=F,
				col.names=T)
	})

target_range <- paste0("1-", sea_library %>% distinct(target) %>% nrow)

cmd <- paste0("
cd ", logs_path, "
qsub -t ", target_range, " ", sea_wrapper_fname, " ", run_base, " ", sea_library_fname, "\n")
cat(cmd)

system(paste0(cmd, "\n"))



# wait for it to finish...

scores <- list.files(paste0(run_base, "/outputs")) %>%
	plyr::ldply(function(target){
		df <- readr::read_csv(paste0(run_base, "/outputs/", target, "/", target, ".csv.out.csv")) %>%
			dplyr::transmute(
				target1 = `Target ID`,
				target2 = `Query ID`,
				affinity = `Affinity Threshold (nM)`,
				Pvalue = `P-Value`,
				MaxTC = `Max Tc`,
				Zscore = `Z-Score`)
	})

fdr <- expand.grid(
	target1 = target_info$uniprot_entry,
	target2 = target_info$uniprot_entry) %>%
	dplyr::left_join(
		scores %>% dplyr::select(target1, target2, Pvalue),
		by=c("target1", "target2")) %>%
	dplyr::mutate(
		Pvalue = ifelse(is.na(Pvalue), runif(n()), Pvalue))

a <- fdrtool(fdr$Pvalue, statistic="pvalue")
fdr <- fdr %>%
	dplyr::mutate(
		Qvalue=a$qval,
		Evalue=Pvalue * ((nrow(target_info) * nrow(target_info))/2 - nrow(target_info)) ) %>%
	dplyr::select(-Pvalue)

# this is slow because there are lot of fdr values
scores <- scores %>%
	left_join(fdr, by=c("target1", "target2"))

scores <- scores %>%
	left_join(
		target_info %>%
			dplyr::select(
				target1 = uniprot_entry,
				entrez_id1 = entrez_id,
				gene_name1 = gene_name,
				description1 = description,
				target1_class_1 = class_1,
				target1_class_2 = class_2,
				target1_class_3 = class_3,
				target1_class_4 = class_4,
				target1_class_5 = class_5,
				target1_class_6 = class_6),
		by=c("target1")) %>%
	left_join(
		target_info %>%
			dplyr::select(
				target2 = uniprot_entry,
				entrez_id2 = entrez_id,
				gene_name2 = gene_name,
				description2 = description,
				target2_class_1 = class_1,
				target2_class_2 = class_2,
				target2_class_3 = class_3,
				target2_class_4 = class_4,
				target2_class_5 = class_5,
				target2_class_6 = class_6),
		by=c("target2")) %>%
	dplyr::select(
		target1,
		target2,
		MaxTC,
		Qvalue,
		entrez_id1,
		gene_name1,
		description1,
		target1_class_1,
		target1_class_2,
		target1_class_3,
		target1_class_4,
		target1_class_5,
		target1_class_6,
		entrez_id2,
		gene_name2,
		description2,
		target2_class_1,
		target2_class_2,
		target2_class_3,
		target2_class_4,
		target2_class_5,
		target2_class_6,
		Zscore,
		Pvalue,
		Evalue)


scores %>% write_csv(scores_fname)

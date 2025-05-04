# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(fdrtool)
library(BioChemPantry)
library(SEAR)

staging_directory <- get_staging_directory("chembl27")

runs_dir <- paste0(staging_directory, "/data/runs")
if(!dir.exists(runs_dir)){
	dir.create(runs_dir)
}

target_info <- readr::read_tsv(
	paste0(staging_directory, "/data/target_info.tsv"),
	col_types = readr::cols(
		.default = col_character(),
		component_id = col_integer(),
		tax_id = col_integer(),
		iuphar_ids = col_integer(),
		entrez_id = col_integer()))

sea_wrapper_fname <- paste0(getwd(), "/scripts/sea-wrapper.sh")

#' Submit a cluster job to compute all target vs. target SEA scores
submit_compute_sea_network <- function(
	library_id,
	fp_type,
	config_files=NULL,
	target_range="all",
	verbose=TRUE,
	dry_run=FALSE
){

	sea_library_fname <- paste0(staging_directory, "/data/", library_id, ".sea")
	sea_library <- unpack.library(
		sea_library_fname, verbose=T)
	sea_library <- sea_library$targets %>%
		left_join(
			sea_library$molecules,
			by=c("compound"))

	run_base <- paste0(
		runs_dir, "/",
		library_id, "_vs_", library_id)
	unlink(run_base, recursive=T, force=T)
	dir.create(run_base)

	if(verbose){
		cat(
			"Computing SEA target vs. target for all targets in library:\n",
			"\tlibrary_id: ", library_id, "\n",
			"\tfp_type: ", fp_type, "\n",
			"\tconfig_files: ", paste(names(config_files), config_files, collapse=", ", sep=": "), "\n",
			"\trun_base: ", run_base, "\n",
			"\tsea_library_fname: ", sea_library_fname, "\n",
			"\tdry_run: ", ifelse(dry_run, "yes", "no"), "\n",
			sep="")
	}

	n_targets <- sea_library %>% distinct(target) %>% nrow
	if(target_range == 'all'){
		target_range <- paste0("1-", n_targets)
	}

	if(verbose){
		cat(
			"\tn_targets: ", n_targets, "\n",
			"\ttarget_range: ", target_range, "\n",
			sep="")
	}

	#### prepare inputs
	if(verbose){
		cat("Preparing run directories and data...\n")
	}

	if(!is.null(config_files)){
		do.call(SEAR::write.config_files, c(dir=run_base, config_files))
	}

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

	cmd <- paste0(
		sea_wrapper_fname, " ",
		run_base, " ",
		sea_library_fname, " ",
		fp_type)

	script <- paste0("
	cd ", logs_path, "
	qsub -t ", target_range, " ", cmd, "\n")

  if(verbose){
		cat("script:\n", script, sep="")
	}

	if(!dry_run){
		system(paste0(script, "\n"))
		if(verbose){
			cat("Run submitted, to check if it is done do 'qstat'\n\n", sep="")
		}
	} else {
		if(verbose){
			cat("To submit run copy and paste in shell, then to check if it is done do 'qstat'\n\n", sep="")
		}
	}
	run_base
}

collect_compute_sea_network <- function(
	staging_directory,
	library_id,
	run_base,
	target_info
){
	scores <- list.files(paste0(run_base, "/outputs")) %>%
		plyr::ldply(function(target){
			df <- readr::read_csv(
				paste0(run_base, "/outputs/", target, "/", target, ".csv.out.csv"),
				col_types=cols(
					`Query ID` = col_character(),
					`Target ID` = col_character(),
					`Affinity Threshold (nM)` = col_integer(),
					`P-Value` = col_double(),
					`Cut Sum` = col_double(),
					`Max Tc` = col_double(),
					`Z-Score` = col_double(),
					Name = col_character(),
					Description = col_character())) %>%
				dplyr::transmute(
					target1 = `Target ID`,
					target2 = `Query ID`,
					affinity = `Affinity Threshold (nM)`,
					Pvalue = `P-Value`,
					CutSum = `Cut Sum`,
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

	a <- fdrtool::fdrtool(fdr$Pvalue, statistic="pvalue")
	fdr <- fdr %>%
		dplyr::mutate(
			Qvalue=a$qval,
			Evalue=Pvalue * ((nrow(target_info) * nrow(target_info))/2 - nrow(target_info)) ) %>%
		dplyr::select(-Pvalue)

	# this is slow because there are lot of fdr values
	scores <- scores %>%
		left_join(fdr, by=c("target1", "target2"))

	scores <- scores %>%
		dplyr::left_join(
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
		dplyr::left_join(
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
			CutSum,
			Zscore,
			Pvalue,
			Evalue)

	scores_fname <- paste0(staging_directory, "/data/", library_id, ".scores.csv")
	scores %>% readr::write_csv(scores_fname)
	s <- scores
}



#######
run_base <- submit_compute_sea_network(
	library_id="chembl27_rdkit_ecfp4",
	target_range='all',
	fp_type='rdkit_ecfp',
	config_files=list(fpcore="[rdkit_ecfp]\ncircle_radius=2\n"),
	verbose=TRUE,
	dry_run=TRUE)

scores <- collect_compute_sea_network(
	staging_directory=staging_directory,
	library_id="chembl27_rdkit_ecfp4",
	run_base=run_base,
	target_info=target_info)

##############
run_base <- submit_compute_sea_network(
	library_id="chembl27_ECFC6",
	target_range='all',
	fp_type='rdkit_ecfc',
	config_files=list(fpcore="[rdkit_ecfc]\ncircle_radius=3\n"),
	verbose=TRUE,
	dry_run=TRUE)

collect_compute_sea_network(
	staging_directory=staging_directory,
	library_id="chembl27_ECFC6",
	run_base=run_base,
	target_info=target_info)




#######
run_base <- submit_compute_sea_network(
	library_id="chembl27_APDP",
	target_range='all',
	fp_type='rdkit_apdp',
	config_files=NULL,
	verbose=TRUE,
	dry_run=TRUE)

collect_compute_sea_network(
	staging_directory=staging_directory,
	library_id="chembl27_ECFC6",
	run_base=run_base,
	target_info=target_info)



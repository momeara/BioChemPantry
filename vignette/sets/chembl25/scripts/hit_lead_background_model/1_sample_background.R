# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(DBI)
library(RSQLite)
library(dbplyr)
library(ggplot2)

library(BioChemPantry)
library(SEAR)

staging_directory <- BioChemPantry::get_staging_directory("chembl23")

# this is used to match the number of targets
# and make comparison plots
load(paste0(staging_directory, "/data/full_chembl_data.Rdata"))

generate_model <- function(
	n_targets,
	n_hits_mu,
	n_hits_size_v,
	n_hits_prob_mu,
	n_hits_prob_v,
	n_leads_mu,
	n_leads_size_v,
	n_leads_prob_mu,
	n_leads_prob_v){

	n_hits_size_mu <- n_hits_mu / n_hits_prob_mu
	n_hits_size_shape <- n_hits_size_mu * n_hits_size_mu / n_hits_size_v
	n_hits_size_scale <- n_hits_size_v / n_hits_size_mu

	n_leads_size_mu <- n_leads_mu / n_leads_prob_mu
	n_leads_size_shape <- n_leads_size_mu * n_leads_size_mu / n_leads_size_v
	n_leads_size_scale <- n_leads_size_v / n_leads_size_mu
	model <- data.frame(tid=1:n_targets) %>%
		dplyr::mutate(
			n_hits_size = rgamma(
				n=n_targets,
				shape=n_hits_size_shape,
				scale=n_hits_size_scale) %>%
				floor,
			n_hits_prob = rbeta(
				n=n_targets,
				shape1=n_hits_prob_mu*n_hits_prob_v,
				shape2=(1-n_hits_prob_mu)*n_hits_prob_v),
			n_hits = rbinom(
				n=n_targets,
				size=n_hits_size,
				prob=n_hits_prob)) %>%
		dplyr::filter(n_hits>0) %>%
		dplyr::rowwise() %>%
		dplyr::do({
			data.frame(
				tid = rep(.$tid, times=.$n_hits),
				n_hits_size = rep(.$n_hits_size, times=.$n_hits),
				n_hits_prob = rep(.$n_hits_prob, times=.$n_hits),
				n_hits = rep(.$n_hits, times=.$n_hits),
				hit_cid = sample.int(n=n_compounds, size=.$n_hits, replace=FALSE),
				hit_id = 1:.$n_hits,
				n_leads_size = rgamma(
					n=.$n_hits,
					shape=n_leads_size_shape,
					scale=n_leads_size_scale) %>%
					floor,
				n_leads_prob = rbeta(
					n=.$n_hits,
					shape1=n_leads_prob_mu*n_leads_prob_v,
					shape2=(1-n_leads_prob_mu)*n_leads_prob_v)) %>%
				dplyr::mutate(
					n_leads = rbinom(n=.$n_hits, size=n_leads_size, prob=n_leads_prob))
		}) %>%
		dplyr::ungroup()
}


plot_model <- function(
	full_chembl_data,
	model
){

	data <- rbind(
		chembl_target_sizes <- full_chembl_data %>%
			dplyr::mutate(tid=uniprot_entry) %>%
			dplyr::group_by(tid) %>%
			dplyr::summarize(n_actives=n()) %>%
			dplyr::mutate(source="ChEMBL23"),
		model_target_sizes <- model %>%
			dplyr::group_by(tid) %>%
			dplyr::summarize(n_actives=sum(n_leads+1)) %>%
			dplyr::mutate(source="Model"))

	p <- ggplot(data=data) + theme_bw() +
		geom_density(
			mapping=aes(x=n_actives, color=source),
			size=2) +
		scale_x_log10("Set Size", breaks=c(1,3,10,30,100,300,1000,3000)) +
		scale_colour_discrete("Source")

	fname_base <- paste0(staging_directory, "/product/hit_lead_background_model/set_sizes")
	ggsave(paste0(fname_base, ".pdf"))
	ggsave(paste0(fname_base, ".png"))

	data %>% group_by(source) %>% summarize(mean(n_actives), sd(n_actives))
}

sample_compounds <- function(
	chembl_smi,
	tcs,
	threshold,
	model
){
	chembl_smi <- chembl_smi %>% dplyr::select(zinc_id, smiles, compound_id)
	n_compounds <- chembl_smi %>% nrow
	tcs <- tcs %>%
		dplyr::select(ref_cid, query_cid, tc) %>%
		dplyr::filter(tc >= threshold)
	model <- model %>% dplyr::select(tid, hit_cid, n_leads)

	Rcpp::sourceCpp("hit_to_lead.cpp", rebuild=TRUE)
	model <- hit_to_lead(n_compounds, tcs, model)
	# model now has columns [tid, hit_cid, lead_cid]
}

save_model_as_library <- function(
	chembl_smi,
	model,
	library_fname,
	...
){

	model_sets <- rbind(
		#hits
		model %>%
			dplyr::distinct(tid, hit_cid) %>%
			dplyr::left_join(
				chembl_smi %>%
					dplyr::select(
						compound=zinc_id,
						hit_cid=compound_id),
				by=c("hit_cid")) %>%
			dplyr::transmute(
				target=tid,
				name=tid,
				affinity=5,
				description="",
				compound=compound),
		#leads
		model %>%
			dplyr::left_join(
				chembl_smi %>%
					dplyr::select(
						compound=zinc_id,
						lead_cid=compound_id),
				by=c("lead_cid")) %>%
			dplyr::transmute(
				target=tid,
				name=tid,
				affinity=5,
				description="",
				compound=compound))

	model_smi <- model_sets %>%
		dplyr::distinct(compound) %>%
		dplyr::left_join(
			chembl_smi %>%
				dplyr::select(
					compound=zinc_id,
					smiles),
			by=c("compound"))

	SEAR::pack.library(
		molecules=model_smi,
		targets=model_sets,
		library_fname=library_fname,
		...)
}

###########################################################################################33

run_base <- SEAR::submit_compute_sea_network(
	library_fname=paste0(staging_directory, "/data/chembl23_rdkit_ecfp4.sea"),
	library_id="chembl23_rdkit_ecfp4_mmd",
	fp_type="rdkit_ecfp4",
	runs_dir=paste0(staging_directory, "/data/runs"),
	config_files=list(
		seacore="[run]\npvalue_cutoff=infinity\nscoretype=mmd\n",
		fpcore="[rdkit_ecfp]\ncircle_radius=2\n"),
	verbose=TRUE,
	dry_run=TRUE)

scores <- SEAR::collect_compute_sea_network(
	run_base=run_base,
	verbose=TRUE)

scores %>%
	readr::write_csv(
		paste0(staging_directory, "/data/chembl23_rdkit_ecfp4.scores.all.csv"))




# hit lead model 2
model <- generate_model(
	n_targets = full_chembl_data %>% distinct(uniprot_entry) %>% nrow,
	n_hits_mu = 8,
	n_hits_size_v = 3000,
	n_hits_prob_mu = .1,
	n_hits_prob_v = 3,
	n_leads_mu = 15,
	n_leads_size_v = 5000,
	n_leads_prob_mu = .5,
	n_leads_prob_v = 3)

plot_model(
	full_chembl_data,
	model)

model <- sample_compounds(
	chembl_smi,
	tcs,
	threshold=.35,
	model)

save_model_as_library(
	chembl_smi,
	model,
	library_fname=paste0(staging_directory, "/data/hit_lead_background_2.sea"),
	fingerprint_type="rdkit_ecfp",
	name="hit_lead_background_2",
	config_files=NULL,
	verbose=TRUE)

# use the fit from the normal library to fit the background
SEAR::set_fit.library(
	library_fname=paste0(staging_directory, "/data/hit_lead_background_2.sea"),
	threshold=.28,
	diagnostics_fname=paste0(staging_directory, "/data/chembl23_rdkit_ecfp4.sea"),
	verbose=TRUE)

run_base <- SEAR::submit_compute_sea_network(
	library_fname=paste0(staging_directory, "/data/hit_lead_background_2.sea"),
	library_id="hit_lead_background_2",
	fp_type="rdkit_ecfp4",
	runs_dir=paste0(staging_directory, "/data/runs"),
	config_files=list(
		seacore="[run]\npvalue_cutoff='infinity'\n",
		fpcore="[rdkit_ecfp]\ncircle_radius=2\n"),
	verbose=TRUE,
	dry_run=TRUE)

scores <- SEAR::collect_compute_sea_network(
	run_base=run_base,
	verbose=TRUE)

scores %>%
	readr::write_csv(
		paste0(staging_directory, "/data/hit_lead_background_2.scores.all.csv"))


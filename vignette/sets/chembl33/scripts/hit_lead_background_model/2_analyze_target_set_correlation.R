# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(BioChemPantry)

staging_directory <- BioChemPantry::get_staging_directory("chembl23")

pantry <- BioChemPantry::get_pantry("sea_chembl23")

date_code <- function(d=NA){
  # reference http://www.r-cookbook.com/node/17
  if(is.na(d)) d <- Sys.Date()
  pattern <- '20([[:digit:]]{2})-([[:digit:]]{2})-([[:digit:]]{2})'
  paste(
    sub(pattern, '\\1', d), sub(pattern, '\\2', d), sub(pattern, '\\3', d),
    sep="")
}

collect_per_target_stats <- function(library_fname, tc_cutoff, source){
	stats_fname <- paste0(tempfile(), "_", source, "_stats.csv")
	cmd <- paste0("
LD_PRELOAD=/nfs/home/momeara/work/tools/anaconda2/lib/libmkl_core.so
python compute_target_set_correlation.py ", library_fname, " ", tc_cutoff, " library ", stats_fname)
	cat("cmd: ", cmd, "\n")
	results <- system(cmd)
	stats <- readr::read_csv(stats_fname) %>%
		dplyr::mutate(source = source)

}

#' stats should be the result of `rbind`ing one or more calls to collect_stats together
plot_per_target_stats <- function(
	scores,
	stats,
	output_fname_base
){
	point_size <- .5
	point_alpha <- .8
	scale_x_set_size <- list(
		scale_x_continuous("Set Size", breaks=c(1, 3, 10, 30, 100, 300, 1000, 3000)))
	plot_height <- 4
	plot_width <- 6

	save_plot <- function(plot, plot_tag){
		sources <- stats %>%
			dplyr::distinct(source) %>%
			magrittr::extract2("source") %>%
			paste0(collapse="_")

		dir.create(file.path(output_fname_base, sources), showWarnings = FALSE)
		fname_base <- paste0(
			output_fname_base, "/",
			sources, "/",
			plot_tag, "",
			date_code())

		ggplot2::ggsave(
			filename=paste0(fname_base, ".pdf"),
			plot=plot,
			height=plot_height,
			width=plot_width)

		ggplot2::ggsave(
			filename=paste0(fname_base, ".png"),
			plot=plot,
			height=plot_height,
			width=plot_width)
	}

#	scores %>%
#		dplyr::left_join(
#			stats %>%
#				dplyr::mutate(...),
#			by=...) %>%
#		dplyr::left_join(
#			stats %>%
#				dplyr::mutate(...),
#			by=...) %>%
#		dplyr::muate(
#			mean_CutSetSum = (ref_CutSetSum + query_CutSetSum)/2)
#
#	p <- ggplot() +
#		theme_bw() +
#		ggtitle("Z-score by mean(CutSetSum)") +
#		geom_point(
#			data=scores,
#			mapping=aes(x=mean_CutSetSum, y=Zscore, color=source)
#		scale_y_continuous("Zscore") +
#		scale_x_continuous("Mean(CutSetSum)")
#	save_plot(p, "Zscore_by_meanCutSetSum")

	# cut set sum
	p <- ggplot() +
		theme_bw() +
		ggtitle("Set self similarity: CutSetSum") +
		geom_point(
			data=stats,
			mapping=aes(x=ref_set_size, y=cut_set_sum, color=source),
			size=.5,
			alpha=.8) +
		coord_trans(x="log10", y="log10") +
		scale_y_continuous("CutSetSum", breaks=c(1,10,100,1000,10000, 100000)) +
		scale_x_set_size
	save_plot(p, "CutSetSum")


	# CutSetSum - set_size
	p <- ggplot() +
		theme_bw() +
		ggtitle("Set self similarity: CutSetSum - SetSize") +
		geom_point(
			data=stats,
			mapping=aes(
				x=ref_set_size,
				y=cut_set_sum-ref_set_size,
				color=source),
			size=.5,
			alpha=.8) +
		coord_trans(x="log10", y="log10") +
		scale_y_continuous("CutSetSum - SetSize") +
		scale_x_set_size
	save_plot(p, "CutSetSum_mSetSize")

	# gini coefficient
	p <- ggplot() +
		theme_bw() +
		ggtitle("Set self similarity: Gini coefficient") +
		geom_point(
			data=stats,
			mapping=aes(x=ref_set_size, y=gini_coeffient, color=source),
			size=.5,
			alpha=.8) +
		coord_trans(x="log10", y=scales::exp_trans(10)) +
		scale_y_continuous("Gini coefficient", ) +
		scale_x_set_size
	save_plot(p, "gini_coefficient")

	# largest eigenvalue
	p <- ggplot() + theme_bw() +
		ggtitle("Set self similarity: largest eigenvalue") +
		geom_point(
			data=stats,
			mapping=aes(x=ref_set_size, y=largest_eigenvalue, color=source),
			size=.5,
			alpha=.8) +
		coord_trans(x="log10",y="log10") +
		scale_y_continuous("largest eigenvalue", breaks=c(1,3,10,30,100,300)) +
		scale_x_set_size
	save_plot(p, "largest_eigenvalue")

	# smallest largest eigenvalue
	p <- ggplot() +
		theme_bw() +
		ggtitle("Set self similarity: smallest eigenvalue") +
		geom_point(
			data=stats,
			mapping=aes(x=ref_set_size, y=smallest_eigenvalue, color=source),
			size=.5,
			alpha=.8) +
		coord_trans(x="log10") +
		scale_y_continuous("smallest eigenvalue") +
		scale_x_set_size
	save_plot(p, "smallest_eigen_value")


#	# MP distribution
#	p <- ggplot() +
#		theme_bw() +
#		ggtitle("Max MP Distribution prob)") +
#		geom_point(
#			data=stats,
#			mapping=aes(x=ref_set_size, y=ref_mp_max_prob, color=source),
#			size=.5,
#			alpha=.8) +
#		scale_y_continuous("MP Distribution prob)") +
#		scale_x_set_size
#	save_plot(p, "mp_max_prob")
}


collect_target_target_scores <- function(
	scores_fname,
	library_fname,
	stats,
	source
) {

	scores <- readr::read_csv(scores_fname)

	lib <- SEAR::unpack.library(library_fname, verbose=T)
	fit <- lib$fit
	scores <- scores %>%
		dplyr::left_join(
				stats %>%
					dplyr::select(
						target1=ref_set_id,
						ref_set_size,
						ref_cut_set_sum=cut_set_sum),
				by="target1") %>%
		dplyr::left_join(
				stats %>%
					dplyr::select(
						target2=query_set_id,
						query_set_size,
						query_cut_set_sum=cut_set_sum),
				by="target2") %>%
		dplyr::mutate(
			cut_set_sum = CutSum,
			cut_set_sum_fit_mu = fit$mu[1] * ((ref_set_size*query_set_size) ** fit$mu[2]) + fit$mu[3],
			cut_set_sum_fit_stddev = fit$sigma[1] * ((ref_set_size*query_set_size) ** fit$sigma[2]) + fit$sigma[3],
			source=source)

	scores <- scores %>%
		dplyr::mutate(
			squared_mmd =
				ref_cut_set_sum/(ref_set_size*ref_set_size) -
				2*cut_set_sum/(ref_set_size*query_set_size) +
				query_cut_set_sum/(query_set_size*query_set_size))

}

plot_scores_fit <- function(
	scores,
	output_fname_base
){
	point_size <- .5
	point_alpha <- .8
	plot_height <- 4
	plot_width <- 10

	if(!dir.exists(output_fname_base)){
		dir.create(output_fname_base)
	}

	scale_x_prod_set_size <- list(
		scale_x_continuous("Set Size", breaks=c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000)))

	save_plot <- function(plot, plot_tag){
		sources <- scores %>%
			dplyr::distinct(source) %>%
			magrittr::extract2("source") %>%
			paste0(collapse="_")

		dir.create(file.path(output_fname_base, sources), showWarnings = FALSE)
		fname_base <- paste0(
			output_fname_base, "/",
			sources, "/",
			plot_tag, "_",
			date_code())

		cat("saving plots to '", fname_base, ".[pdf,png]'\n", sep="")

#		ggplot2::ggsave(
#			filename=paste0(fname_base, ".pdf"),
#			plot=plot,
#			height=plot_height,
#			width=plot_width)

		ggplot2::ggsave(
			filename=paste0(fname_base, ".png"),
			plot=plot,
			height=plot_height,
			width=plot_width)
	}

	scores <- scores %>% filter(target1 < target2)
	scores <- scores %>% mutate(facet=ifelse(MaxTC==1, "MaxTC==1", "MaxTC < 1"))
	scores_model <- scores %>% sample_n(10000)

	# cut set sum
	p <- ggplot() +
		theme_bw() +
		ggtitle(
			label="Target Target Similarity: CutSetSum",
			subtitle="Blue: E-value < 1e-20, Red: E-value >= 1e-20") +
		geom_point(
			data=scores,
			mapping=aes(
				x=ref_set_size * query_set_size,
				y=cut_set_sum,
				color=(Evalue < 1e-20)),
			size=.1,
			alpha=.3) +
		geom_point(
			data=scores %>% filter(MaxTC==1),
			mapping=aes(
				x=ref_set_size * query_set_size,
				y=cut_set_sum,
				color=(Evalue < 1e-20)),
			size=.1,
			alpha=.3) +
		geom_line(
			data=scores_model,
			mapping=aes(
				x=ref_set_size * query_set_size,
				y=pmax(.28, cut_set_sum_fit_mu)),
			color="black") +
		geom_line(
			data=scores_model,
			mapping=aes(
				x=ref_set_size * query_set_size,
				y=pmax(.28, cut_set_sum_fit_mu - cut_set_sum_fit_stddev)),
			color="darkgreen") +
 		geom_line(
			data=scores_model,
			mapping=aes(
				x=ref_set_size * query_set_size,
				y=pmax(.28, cut_set_sum_fit_mu + cut_set_sum_fit_stddev)),
			color="darkgreen") +
		geom_line(
			data=scores_model,
			mapping=aes(
				x=ref_set_size * query_set_size,
				y=pmax(.28, cut_set_sum_fit_mu - 2*cut_set_sum_fit_stddev)),
			color="green") +
 		geom_line(
			data=scores_model,
			mapping=aes(
				x=ref_set_size * query_set_size,
				y=pmax(.28, cut_set_sum_fit_mu + 2*cut_set_sum_fit_stddev)),
			color="green") +
		facet_wrap(~facet) +
		coord_trans(x="log10", y="log10") +
		scale_y_continuous("CutSetSum", breaks=c(.1, 1,10,100,1000,10000, 100000)) +
		scale_x_continuous(
			"Product of Target Set Sizes",
			breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
		scale_color_discrete("Evalue", breaks=c("< 1e-20", ">= 1e-20"))
#	save_plot(p, "hits_CutSetSum")

	scores <- scores %>%
		dplyr::mutate(
			squared_mmd =
				ref_cut_set_sum/(ref_set_size*ref_set_size) -
				2*cut_set_sum/(ref_set_size*query_set_size) +
				query_cut_set_sum/(query_set_size*query_set_size))

	scores <- scores %>% arrange(desc(Evalue))

	p <- ggplot() +
		theme_bw() +
		ggtitle(
			label="Target Target Similarity: MMD By Product of Target Set Sizes",
			subtitle="Blue: E-value < 1e-20, Red: E-value >= 1e-20") +
		geom_point(
			data=scores,
			mapping=aes(
				x=log(ref_set_size * query_set_size),
				y=squared_mmd,
				color=(Evalue < 1e-20)),
			size=.1,
			alpha=.3) +
		geom_density_2d(
			data=scores,
			mapping=aes(
				x=log(ref_set_size * query_set_size),
				y=squared_mmd),
			color="black",
			size=.4) +
		facet_wrap(~facet) +
		scale_y_continuous("MMD^2") +
		scale_x_continuous(
			"Product of Target Set Sizes",
			expand=c(0,0),
			limits=log(c(10, 5000000)),
			breaks=log(c('10'=10, '100'=100, '1,000'=1000, '10,000'=10000, '100,000'=100000, '1,000,000'=1000000))) +
		scale_color_discrete("Evalue", breaks=c("< 1e-20", ">= 1e-20"))
	save_plot(p, "hits_squared_mmd")

	p <- ggplot() +
		theme_bw() +
		ggtitle(
			label="Target Target Similarity: MMD by CutSetSum",
			subtitle="Blue: E-value < 1e-20, Red: E-value >= 1e-20") +
		geom_point(
			data=scores,
			mapping=aes(
				x=log(cut_set_sum),
				y=squared_mmd,
				color=(Evalue < 1e-20)),
			size=.1,
			alpha=.3) +
		geom_density_2d(
			data=scores,
			mapping=aes(
				x=log(cut_set_sum),
				y=squared_mmd),
			color="black",
			size=.4)+
		facet_wrap(~facet) +
		scale_y_continuous("MMD^2") +
		scale_x_continuous(
			"CutSetSum",
			expand=c(0,0),
			limits=log(c(.1, 10000)),
			breaks=log(c('.1'=.1, '1'=1,'10'=10,'100'=100,'1,000'=1000,'10,000'=10000))) +
		scale_color_discrete("Evalue", breaks=c("< 1e-20", ">= 1e-20"))
	save_plot(p, "hits_squared_mmd_by_cut_set_sum")

}


plot_scores_overlap <- function(
	scores,
	output_fname_base
){
	plot_height <- 5
	plot_width <- 10

	if(!dir.exists(output_fname_base)){
		dir.create(output_fname_base)
	}


	save_plot <- function(plot, plot_tag){
		sources <- scores %>%
			dplyr::distinct(source) %>%
			magrittr::extract2("source") %>%
			paste0(collapse="_")

		dir.create(file.path(output_fname_base, sources), showWarnings = FALSE)
		fname_base <- paste0(
			output_fname_base, "/",
			sources, "/",
			plot_tag, "_",
			date_code())

		cat("saving plots to '", fname_base, ".[pdf,png]'\n", sep="")

#		ggplot2::ggsave(
#			filename=paste0(fname_base, ".pdf"),
#			plot=plot,
#			height=plot_height,
#			width=plot_width)

		ggplot2::ggsave(
			filename=paste0(fname_base, ".png"),
			plot=plot,
			height=plot_height,
			width=plot_width)
	}

	scores <- scores %>% filter(target1 < target2)
		dplyr::group_by(source) %>%
			dplyr::sample_n(50000) %>%
		dplyr::ungroup()

	scores <- scores %>% mutate(facet=ifelse(MaxTC==1, "MaxTC==1", "MaxTC < 1"))

	browser()
	p <- ggplot() +
		theme_bw() +
		ggtitle(
			label="Target Target Similarity: CutSetSum by product of set sizes",
			subtitle="Background: HitLead   Contour: ChEMBL23") +
		stat_density_2d(
			data=scores %>% filter(source=="HitLead2"),
			geom="raster",
			contour=FALSE,
			mapping=aes(
				x=log(ref_set_size * query_set_size),
				y=log(cut_set_sum),
				fill=..density..)) +
		geom_density_2d(
			data=scores %>% filter(source=="chembl23_rdkit_ecfp4"),
			mapping=aes(
				x=log(ref_set_size * query_set_size),
				y=log(cut_set_sum)),
			color="green",
			size=1)+
		facet_wrap(~facet)+
		scale_y_continuous(
			"CutSetSum",
			expand=c(0,0),
			limits=log(c(.1, 10000)),
			breaks=log(c('.1'=.1, '1'=1,'10'=10,'100'=100,'1,000'=1000,'10,000'=10000))) +
		scale_x_continuous(
			"Product of Target Set Sizes",
			expand=c(0,0),
			limits=log(c(10, 5000000)),
			breaks=log(c('10'=10, '100'=100, '1,000'=1000, '10,000'=10000, '100,000'=100000, '1,000,000'=1000000))) +
		scale_fill_gradient2(
			"Density",
			low="black",
			mid=muted("blue"),
			high="white")
	save_plot(p, "hit_CutSetSum_overlap")
}

##############################################################3
# Collect par-target statistics

stats_HitLead2 <- collect_per_target_stats(
		library_fname=paste0(staging_directory, "/data/hit_lead_background_2.sea"),
		tc_cutoff=.35,
		source="HitLead2")

stats_chembl23 <- collect_per_target_stats(
	library_fname=paste0(staging_directory, "/data/chembl23_rdkit_ecfp4.sea"),
	tc_cutoff=0.28,
	source="chembl23_rdkit_ecfp4")

# plot per-target statistics
plot_stats(
	stats=rbind(
		stats_HitLead2),
	output_fname_base=paste0(staging_directory, "/product/hit_lead_background_2"))

plot_stats(
	stats=rbind(
		stats_chembl23),
	output_fname_base=paste0(staging_directory, "/product/hit_lead_background_2"))

plot_stats(
	stats=rbind(
		stats_HitLead2,
		stats_chembl23),
	output_fname_base=paste0(staging_directory, "/product/hit_lead_background_model"))

# collect target vs target scores
scores_HitLead2 <- collect_target_target_scores(
	scores_fname=paste0(staging_directory, '/data/hit_lead_background_2.scores.all.csv'),
	library_fname=paste0(staging_directory, "/data/hit_lead_background_2.sea"),
	stats=stats_HitLead2,
	source="HitLead2")

scores_chembl23 <- collect_target_target_scores(
	scores_fname=paste0(staging_directory, '/data/chembl23_rdkit_ecfp4.scores.all.csv'),
	library_fname=paste0(staging_directory, "/data/chembl23_rdkit_ecfp4.sea"),
	stats=stats_chembl23,
	source="chembl23_rdkit_ecfp4")

# plot target vs. target scores
plot_scores_fit(
	scores=scores_HitLead2,
	output_fname_base=paste0(staging_directory, "/product/hit_lead_background_model"))

plot_scores_fit(
	scores=scores_chembl23,
	output_fname_base=paste0(staging_directory, "/product/hit_lead_background_model"))

plot_scores_overlap(
	scores=rbind(
		scores_HitLead2,
		scores_chembl23),
	output_fname_base=paste0(staging_directory, "/product/hit_lead_background_model"))

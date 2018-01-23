# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(data.table)
library(stringr)
library(sqldf)

source("~/work/sets/uniprot_idmapping/scripts/uniprot_id_map.R")
source("~/work/sea/scripts/sea.R")

#########################
guide_to_pharmacology_targets_table_fname <- "data/targets_and_families_140708.csv"
tt_sea_scores_fname <- "~/work/phenologs/sea_ChEMBL_filtered_131030_ChEMBL_filtered_131030/sea_predictions/scores_1e-5.csv"
tc_sea_scores_fname <- "~/work/phenologs/sea_HMDB_filtered_131030_ChEMBL_filtered_131030/sea_predictions/scores_1e-5.csv"


targets <- read.table(guide_to_pharmacology_targets_table_fname, header=T, sep=",", quote="\"")
gpcr <- targets %>%
	select(Type, Target.name, Human.Entrez.Gene, Human.SwissProt) %>%
	filter(Type == "gpcr") %>%
	mutate(Target.name = str_replace_all(Target.name, "<[^>]+>", "")) %>%
	mutate(Target.name = str_replace_all(Target.name, "&alpha;", "a")) %>%
	mutate(Target.name = str_replace_all(Target.name, "&beta;", "b")) %>%
	mutate(Target.name = str_replace_all(Target.name, "&gamma;", "g")) %>%
	mutate(Target.name = str_replace_all(Target.name, "&deta;", "d")) %>%
	mutate(uniprot_entry = uniprot_id_to_entry(Human.SwissProt))

#missing swissprot entries:
gpcr %>% filter(Human.SwissProt == "")
#   Type     Target.name Human.Entrez.Gene Human.SwissProt
#1  gpcr 5-ht5b receptor            645694                
#2  gpcr   AMY1 receptor                                  
#3  gpcr   AMY2 receptor                                  
#4  gpcr   AMY3 receptor                                  
#5  gpcr   CGRP receptor                                  
#6  gpcr    AM1 receptor                                  
#7  gpcr    AM2 receptor                                  
#8  gpcr  GABAB receptor                                  
#9  gpcr   TRH2 receptor                                  
#10 gpcr           GPR79             27200                
#11 gpcr         TAAR4P             503612                

tt_scores <- tt_sea_scores_fname %>% read.scores()
tc_scores <- tc_sea_scores_fname %>% read.scores()

gpcr_to_nongpcr_scores <- tt_scores %>%
	filter(
		target1 %in% gpcr$uniprot_entry,
		!(target2 %in% gpcr$uniprot_entry)) %>%
	arrange(EValue) %>%
	group_by(target1, target2) %>%
	filter(row_number(EValue) == 1)



tct <- sqldf("
SELECT
	tt.target1 AS target1,
	tc.target1 AS metabolite,
	tt.target2 AS target2,
	tt.EValue AS tt_EValue,
	tt.MaxTC AS tt_MaxTC,
	max(tc.EValue, ct.EValue) AS max_EValue,
	min(tc.MaxTC, ct.MaxTC) AS min_MaxTC
FROM
	gpcr_to_nongpcr_scores AS tt,
	tc_scores AS tc,
	tc_scores AS ct
WHERE
	tc.target2 = tt.target1 AND
	ct.target1 = tc.target1 AND
	ct.target2 = tt.target2;")

tct2 <- tct %>%
	arrange(min_MaxTC) %>%
	group_by(target1, target2) %>%
	filter(row_number(tt_EValue) == 1)

tct2 <- tct2 %>% mutate(target2_prefix = str_extract(target2, "^[A-Z0-9]+"))

tct3 <- tct2 %>% group_by(target1, target2_prefix) %>% filter(row_number(tt_EValue) == 1)



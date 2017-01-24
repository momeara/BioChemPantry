# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(curl)

curl::curl_download(
	"http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip",
	paste0("data/hmdb_metabolites.zip")

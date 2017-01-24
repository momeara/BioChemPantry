# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

require(data.table)
require(plyr)
library(dplyr)
require(stringr)
require(sqldf)


library(XML)

#####################################
# Analysis parameters

# input paths
hmdb_metabolites_xml_paths <- paste0(
	"data/hmdb_metabolites_xml/",
	list.files("data/hmdb_metabolites_xml", "HMDB[0-9]+.xml"))

hmdb_metabolite_activities_fname <- "data/hmdb_metabolite_activities_140908.csv"
####################################
# prepare primary data


xmlChild <- function(xml, name){
	child_index <- which(sapply(xmlChildren(xml), xmlName) == name)
	xmlChildren(xml)[[child_index]]
}



# this takes an hour or so
activities <- hmdb_metabolites_xml_paths %>% lapply(function(xml_path){
	t <- xmlInternalTreeParse(xml_path, getDTD=F) %>% xmlChild("metabolite")
	metabolite_accession <- xmlValue( t %>% xmlChild("accession") )

	targets <- t %>% xmlChild("protein_associations")
	if(is.null(xmlChildren(targets)$protein)){
		return( data.frame() )
	} else {
		ldply(.data=targets %>% xmlChildren, .fun=function(target) {
			data.frame(
				metabolite_accession = metabolite_accession,
				hmdb_protein_accession = target %>%
					xmlChild("protein_accession") %>% xmlValue,
				uniprot_id = target %>%
					xmlChild("uniprot_id") %>% xmlValue,
				protein_name = target %>%
					xmlChild("name") %>% xmlValue)
		}) %>% select(-.id)
	}
})
activities_tbl <- activities %>% rbind_all

activities_tbl %>% write.table(
		hmdb_metabolite_activities_fname, row.names=F, sep="\t")

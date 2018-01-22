# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(stringr)
library(plyr)
library(dplyr)
library(magrittr)
library(httr)

uniprot_host <- "http://www.uniprot.org"
user_agent_arg <- user_agent("httr http://github.com/momeara/BioChemPantry")


## helper functions
prepare_data <- function(data, col){
	if(!("data.frame" %in% class(data))){
		data <- data.frame(col=data)
		names(data) <- c(col)
	}
	if(nrow(data) == 0){
		cat("WARNING: No uniprot_entry values to convert\n")
	}
	data
}

uniprot_entry_web_lookup_debug <- function(
	request_bin_host,
	columns
){
	library(httr)
	library(plyr)
	library(dplyr)
	library(readr)
	column_arg <- paste0(columns, collapse=",")
	r <- GET(
		request_bin_host,
		user_agent_arg,
		query = list(
			query="DRD2_HUMAN",
			format = 'tab',
			columns = column_arg)) %>%
	content
}

uniprot_web_batch_debug <- function(
	request_bin_host = "http://www.uniprot.org/batch/",
	format='txt',
	columns=c("entry name")
){
	library(httr)
	library(plyr)
	library(dplyr)
	library(readr)
	uniprot_entries <- data.frame(uniprot_entry="DRD2_HUMAN")
	uniprot_entries_fname <- tempfile()
	write.table(uniprot_entries, uniprot_entries_fname, quote=F, col.names=F, row.names=F, sep="")
	uniprot_entries_f <- upload_file(uniprot_entries_fname, type="text/plain")

	column_arg <- paste0(columns, collapse=",")
	r <- POST(
		request_bin_host,
		user_agent_arg,
		query = list(
			query=uniprot_entries_f,
			format = format,
			columns = column_arg)) %>%
	content
}

#' Get lookup info for given uniprot accessions one at a time
#' @export
uniprot_web_lookup <- function(
	uniprot_accns,
	columns = c("entry name"),
	verbose=F
){
	library(httr)
	library(plyr)
	library(dplyr)
	library(readr)
	column_arg <- paste0(columns, collapse=",")
	if(verbose){
		cat("Getting data for ", length(uniprot_entries), " uniprot accn ... \n", sep="")
	}
	r2 <- adply(
		data_frame(uniprot_accn = uniprot_accns),
		1,
		function(df){
			if(verbose){
				cat("Looking up entry for '", df$uniprot_accn, "' ... ", sep="")
			}
			r <- GET(
				paste0(uniprot_host, "/uniprot/", df$uniprot_accn),
				user_agent_arg,
				query = list(
					format = 'tab',
					columns = column_arg)) %>%
			content
			if(r %>% is.null){
				if(verbose){
					cat(" MISSING\n")
				}
				return(data_frame())
			} else {
				if(verbose){
					cat(" GOT IT\n")
				}
				r %>%
					read_tsv %>%
					return
			}
		})
	if(nrow(r2) == 0){
		cat("WARNING: no entries were found\n")
	}
	r2
}

#' Get lookup info for given uniprot entries one at a time
#'
#'  #http://www.uniprot.org/help/uniprotkb_column_names
#'  given a vector of uniprot entries return information about them
#'
#' @export
uniprot_entry_web_lookup <- function(
	uniprot_entries,
	columns = c("entry name"),
	format='tab',
	verbose=FALSE){
	library(httr)
	library(plyr)
	library(dplyr)
	library(readr)
	column_arg <- paste0(columns, collapse=",")
	if(verbose){
		cat("Getting data for ", length(uniprot_entries), " uniprot entries ... \n", sep="")
	}
	r2 <- adply(
		data_frame(uniprot_entry = uniprot_entries),
		1,
		function(df){
			if(verbose){
				cat("Looking up entry for '", df$uniprot_entry, "' ... ", sep="")
			}
			r <- GET(
				paste0(uniprot_host, "/uniprot/"),
				user_agent_arg,
				query = list(
					query=df$uniprot_entry,
					format = format,
					columns = column_arg)) %>%
			content
			if(r %>% is.null){
				if(verbose){
					cat(" MISSING\n")
				}
				return(data_frame())
			} else {
				if(verbose){
					cat(" GOT IT\n")
				}
				r %>%
					read_tsv %>%
					filter(`Entry name` == df$uniprot_entry) %>%
					head(1) %>%
					return
			}
		})
	if(nrow(r2) == 0){
		cat("WARNING: no entries were found\n")
	}
	r2
}

#' Get all taxa below a given ancestor taxa
#' @export
uniprot_taxa <- function(
	ancestor=7711, # chordate
	reviewed=TRUE,
	format='tab',
	verbose=F
){
	library(httr)
	library(plyr)
	library(dplyr)
	library(readr)
	if(verbose){
		cat("Retriving taxa descending from '", ancestor, "' from uniprot...\n")
	}
	query_arg <- paste0(
		"rank:species ",
		ifelse(reviewed, "uniprot:(reviewed:yes) ", ""),
		"ancestor:", ancestor)

	r <- httr::GET(
		paste0(uniprot_host, "/taxonomy/"),
		user_agent_arg,
		query = list(
			query=query_arg,
			format = 'tab',
			compress='yes'))

	if(r %>% is.null){
		if(verbose){
			cat(" None retrieved\n")
		}
		return(data_frame())
	}

	r %>%
		httr::content() %>%
		rawConnection() %>%
		gzcon() %>%
		readr::read_tsv(
			col_types=readr::cols(
				Taxon = readr::col_double(),
				Mnemonic = readr::col_character(),
				`Scientific name` = readr::col_character(),
				`Common name` = readr::col_character(),
				Synonym = readr::col_character(),
				`Other Names` = readr::col_character(),
				Reviewed = readr::col_character(),
				Rank = readr::col_character(),
				Lineage = readr::col_character(),
				Parent = readr::col_double(),
				`Virus hosts` = readr::col_character())) %>%
		dplyr::select(
			taxon=Taxon,
			mnemonic=Mnemonic,
			scientific_name=`Scientific name`,
			common_name=`Common name`,
			synonym=Synonym,
			other_name=`Other Names`,
			reviewed=Reviewed,
			rank=Rank,
			lineage=Lineage,
			parent=Parent,
			virus_hosts=`Virus hosts`) %>%
		return
}

# don't try to get tab results, it doesn't work :/
uniprot_web_batch <- function(
	uniprot_entries,
	format='txt',
	columns=c("entry name"),
	url=paste0(uniprot_host, "/batch/"),
	verbose=F
){
	library(httr)
	library(plyr)
	library(dplyr)
	library(readr)
	uniprot_entries_fname <- tempfile()
	write.table(uniprot_entries, uniprot_entries_fname, quote=F, col.names=F, row.names=F, sep="")
	uniprot_entries_f <- upload_file(uniprot_entries_fname, type="text/plain")
	r <- httr::POST(
		url=url,
		user_agent_arg,
		body=list(
			file = uniprot_entries_f,
			format=format))
	file.remove(uniprot_entries_fname)
	if(r %>% is.null){
		if(verbose){
			cat(" None retrieved\n")
		}
		return(data_frame())
	}

	if(format=='tab'){
		z <- r %>% httr::content() %>%
		readr::read_tsv()
	} else if(format=="txt"){
		z <- r %>% httr::content() %>% str_split("//\n")
	}
	z
}

#parse_uniprot_entries <- function(
#	uniprot_entries_full,
#	parsers = list()
#) {
#	library(stringr)
#	str_split(uniprot_entries_full, "\n//\n")
#	



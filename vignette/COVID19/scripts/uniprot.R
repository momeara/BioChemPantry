
#' See http://www.uniprot.org/help/uniprotkb_column_names for available columns
uniprot_entry_web_lookup <- function(
	uniprot_entries,
	columns,
	user_agent_arg,
	verbose=F
){
	column_arg <- paste0(columns, collapse=",")
	if(verbose){
		cat("Getting data for ", length(uniprot_entries), " uniprot entries ... \n", sep="")
	}
	r2 <- plyr::adply(
		dplyr::data_frame(uniprot_entry = uniprot_entries),
		1,
		function(df){
			if(verbose){
				cat("Looking up entry for '", df$uniprot_entry, "' ... ", sep="")
			}
			r <- httr::GET(
				"http://www.uniprot.org/uniprot/",
				user_agent_arg,
				query = list(
					query=df$uniprot_entry,
					format = 'tab',
					columns = column_arg)) %>%
				httr::content()
			if(r %>% is.null){
				if(verbose){
					cat(" MISSING\n")
				}
				return(dplyr::data_frame())
			} else {
				if(verbose){
					cat(" GOT IT\n")
				}
				r %>%
					readr::read_tsv() %>%
					dplyr::filter(`Entry name` == df$uniprot_entry) %>%
					head(1) %>%
					return
			}
		})
	if(nrow(r2) == 0){
		cat("WARNING: no entries were found\n")
	}
	r2
}

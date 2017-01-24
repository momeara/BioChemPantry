
#' test if the columns names(ids) are symmetric to the columns ids
#' and for those that are symmetric, test equality over the value columns
#' @export
test_is_symmetric <- function(data, ids, values){
	library(plyr)
	library(dplyr)
	if(class(ids) != "character"){
		paste0("class(ids) == ", class(ids), " is not 'character'.") %>%
			stop
	}
	if(!all( names(ids) %in% names(data))){
		paste0("names(ids): [", paste0(names(ids), collapse=", "), "] is not a subset of names(data) = [", paste0(names(data), collapse=", "), "].")
	}
	if(!all( ids %in% names(data))){
		paste0("ids: [", paste0(ids, collapse=", "), "] is not a subset of names(data) = [", paste0(names(data), collapse=", "), "].")
	}
	anti_ids = array(names(ids), dimnames=list(ids))
	n_different <- rbind(
		data %>% anti_join(data, by=ids),
		data %>% anti_join(data, by=anti_ids)) %>%
		nrow
	return(n_different == 0)
}

#' make a dataset symmetric
#' @export
make_symmetric <- function(data, ids){
	data %>%
		rbind(data %>% rename_(.dots=ids)) %>%
		distinct_(.dots=c(ids, names(ids)))
}

#' test fraction null
#' @export
test_null_frac <- function(data, ids, frac){
	library(plyr)
	library(dplyr)
	(((data %>% select_(.dots=ids) %>% complete.cases) / nrow) < frac) %>% return
}

#' test number null
#' @export
test_null_count <- function(data, ids, count){
	library(plyr)
	library(dplyr)
	((data %>% select_(.dots=ids) %>% complete.cases) < count) %>% return
}

#' test distinct rows
#' e.g. for testing for network loops
#' @export
test_distinct <- function(data, ids){
	(nrow(data) == nrow(data %>% distinct_(.dots=ids))) %>% return
}

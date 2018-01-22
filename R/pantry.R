# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

#' Get and check login information for Postgres database
#'
#' @param pantry_config either the name of a json file with config information or
#'   the config info parsed into a list. The login key should correspond with the argument
#'   keys of the dplyr::src_postgres function.
#' @return list of login parameters for Postgres database
#' @export
get_pantry_config <- function(pantry_config="~/.pantry_config"){
	if(is.character(pantry_config)){
		pantry_config <- jsonlite::fromJSON(pantry_config)
	} else if(!is.list(config_config)){
		stop("Login must either be the name of a json file of login info, or a list of login info.\n")
	}
	pantry_config
}

#' Get the staging directory for the given dataset
#'
#' @param data_set A tag for the given dataset
#' @param pantry_config see get_pantry_config
#' @return '<pantry_config$staging_dir>/<data_set>', creating the directory if necessary
#' @export
get_staging_directory <- function(data_set, pantry_config="~/.pantry_config", verbose=TRUE){
	pantry_config <- BioChemPantry::get_pantry_config(pantry_config)
	staging_directory <- paste0(pantry_config$staging_directory, "/", data_set)
	if(!file.exists(staging_directory)){
		tryCatch({
			dir.create(staging_directory, recursive=TRUE)
		}, error=function(e){
			stop(paste0("Unable to get staging area for dataset '", data_set, "': '", staging_directory, "'."))
		})
	}
	if(verbose){
		cat("Staging directory: ", staging_directory, "\n")
	}
	staging_directory
}

#' Attach to the pantry database, optionally using a schema. If the
#' specified schema doesn't exist, create it.
#'
#' @param pantry_config see get_pantry_config
#' @return a connection to the the pantry
#' @export
get_pantry <- function(schema=NULL, pantry_config="~/.pantry_config"){
	pantry_config <- BioChemPantry::get_pantry_config(pantry_config)

	if(!is.null(schema)){
		schema <- tolower(schema)
	}

	# this two step process is to allow DBI::dbConnect to pick the right specialized method
	# while allowing arguments to be passed in as a list
	get_pantry_from_login <- function(...){DBI::dbConnect(RPostgres::Postgres(), ...)}
	pantry <- do.call(get_pantry_from_login, pantry_config$login)

	if(!is.null(schema)){
		schemas <- BioChemPantry::get_schemas(pantry)
		if(!(schema %in% schemas$schema_name)){
			stop(paste0("Schema '", schema, "' doesn't exist. To create it do:\n\n  pantry <- get_pantry()\n  pantry %>% create_schema('", schema, "')\n  pantry %>% set_schema('", schema, "')\n"))
		}
		pantry %>% BioChemPantry::set_schema(schema)
	}
	return(pantry)
}

#' Create a schema in the pantry
#' @export
create_schema <- function(pantry, schema_name){
	DBI::dbSendQuery(pantry, paste0("CREATE SCHEMA ", schema_name %>% dplyr::sql(), ";")) %T>%
		DBI::dbFetch() %T>%
		DBI::dbClearResult()

}

#' Drop a schema from the pantry
#' @export
drop_schema <- function(pantry, schema_name){
	DBI::dbSendQuery(
		pantry,
		paste0("DROP SCHEMA \"", schema_name, "\" CASCADE;")) %T>%
		DBI::dbFetch() %T>%
		DBI::dbClearResult()
}

#' Set a schema as the active schema
#' @export
set_schema <- function(pantry, schema_name, include_temp=TRUE){
	search_path <- c("\"$user\"")
	if(include_temp){
			# https://www.postgresql.org/message-id/il7g0n%24r4j%241%40dough.gmane.org
			temp_namespace <- pantry %>%
				DBI::dbGetQuery("SELECT nspname FROM pg_namespace WHERE oid = pg_my_temp_schema();")
			if (nrow(temp_namespace) > 0){
				search_path <- c(search_path, temp_namespace)
			}
	}
	if(!is.null(schema_name)){
		search_path <- c(search_path, schema_name)
	}
	search_path <- c(search_path, "public")

	search_path <- paste(search_path, collapse=",")

	pantry %>% DBI::dbSendQuery(paste0("SET search_path TO ", search_path, ";")) %T>%
		DBI::dbFetch() %T>%
		DBI::dbClearResult()
}

#' Get the search path for looking for namespaces
#' @export
get_search_path <- function(pantry){
	a <- DBI::dbGetQuery(pantry, "SHOW search_path")
	stringr::str_split(a$search_path, ", ")[[1]]
}

#' Alias for get_search_path
#' @export
get_schema <- function(pantry){
	get_search_path(pantry)
}

#' Get all available schemas
#' @export
get_schemas <- function(pantry, verbose=T){
	if(verbose){
		cat("Current schema: ", BioChemPantry::get_search_path(pantry), "\n", sep="")
	}
	pantry %>%
		dplyr::tbl(dbplyr::build_sql("SELECT schema_name FROM information_schema.schemata")) %>%
		dplyr::arrange(schema_name) %>%
		dplyr::collect() %>%
		data.frame
}

#' Get all tables in the active schema or search path
#' @export
get_tables <- function(pantry all_tables=FALSE){
	tables <- DBI::dbGetQuery(pantry, paste0(
		"SELECT schemaname, tablename FROM pg_tables ",
		"WHERE schemaname != 'information_schema' AND schemaname !='pg_catalog'"))
	if(is.null(schema)){
		schema <- BioChemPantry::get_search_path(pantry)
	}
	if(!all_tables){
		if(length(schema) > 0){
			tables <- tables[tables$schemaname %in% schema,]
		}
	}
	tables
}

#' Drop a table
#' @export
drop_table <- function(pantry, table){
	DBI::dbSendQuery(pantry, paste0("DROP TABLE \"", table, "\" CASCADE;")) %T>%
		DBI::dbFetch() %T>%
		DBI::dbClearResult()

}

#' Identify a table qualified by a schema. This works around issue 244 in dplyr
#' @export
schema_tbl <- function(pantry, schema_table){
	# see http://stackoverflow.com/questions/21592266/i-cannot-connect-postgresql-schema-table-with-dplyr-package
	# https://github.com/hadley/dplyr/issues/244
	pantry %>%
		dplyr::tbl(structure(paste0("SELECT * FROM ",  schema_table), class=c("sql", "character")))
}

# #' this adds the fast argument working around a known problem in dplyr and DBI
# #' https://github.com/hadley/dplyr/blob/master/R/tbl-sql.r#L294
# #'
# #' https://github.com/hadley/dplyr/issues/1471
# #' https://github.com/rstats-db/DBI/issues/62
# #' @export
# copy_to.src_postgres <- function(
# 	dest,
# 	df,
# 	name = deparse(substitute(df)),
# 	types = NULL,
# 	temporary = TRUE,
# 	unique_indexes = NULL,
# 	indexes = NULL,
# 	analyze = TRUE,
# 	fast=FALSE,
# 	...
# ) {
# 	assertthat::assert_that(
# 		is.data.frame(df),
# 		is.string(name),
# 		is.flag(temporary))
# 	class(df) <- "data.frame" # avoid S4 dispatch problem in dbSendPreparedQuery
# 
# 	if (isTRUE(dplyr::db_has_table(pantry, name))) {
# 		stop("Table ", name, " already exists.", call. = FALSE)
# 	}
# 
# 	types <- if(is.null(types)) dplyr::db_data_type(pantry, df) else types
# 	names(types) <- names(df)
# 
# 	dplyr::db_begin(pantry)
# 	tryCatch({
# 		dplyr::db_create_table(pantry, name, types, temporary = temporary)
# 
# 		if(fast){
# 			#http://james.hiebert.name/blog/work/2011/10/24/RPostgreSQL-and-COPY-IN/
# 			#warning this doesn't report any errors on failure!
# 			field_names <- dplyr::escape(dplyr::ident(names(types)), collapse = NULL, con = pantry)
# 			fields <- dplyr:::sql_vector(field_names, parens = TRUE, collapse = ", ", con = pantry)
# 			sql <- dbplyr::build_sql(
# 				"COPY ", dplyr::ident(name), " ", fields, " FROM STDIN;", con=pantry)
# 			DBI::dbSendQuery(pantry, sql)
# 			RPostgreSQL::postgresqlCopyInDataframe(con, df)
# 		} else {
# 			dplyr::db_insert_into(pantry, name, df)
# 		}
# 
# 		dplyr:::db_create_indexes(pantry, name, unique_indexes, unique = TRUE)
# 		dplyr:::db_create_indexes(pantry, name, indexes, unique = FALSE)
# 
# 		if (analyze) dplyr::db_analyze(pantry, name)
# 		dplyr::db_commit(pantry)
# 	}, error = function(err){
# 		dplyr::db_rollback(pantry)
# 	})
# 	dplyr::tbl(dest, name)
# }

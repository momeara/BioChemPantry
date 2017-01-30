# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


#library(plyr)
#library(dplyr)
#library(stringr)
#library(jsonlite)

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
get_staging_directory <- function(data_set, pantry_config="~/.pantry/config"){
	pantry_config <- get_pantry_config(pantry_config)
	staging_directory <- paste0(pantry_config$staging_directory, "/", data_set)
	if(!file.exists(staging_directory)){
		tryCatch({
			dir.create(staging_directory, recursive=TRUE)
		}, error=function(e){
			stop(paste0("Unable to get staging area for dataset '", data_set, "': '", staging_directory, "'."))
		})
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
	pantry_config <- get_pantry_config(pantry_config)

	if(!is.null(schema)){
		schema <- tolower(schema)
	}

	pantry <- do.call(dplyr::src_postgres, pantry_config$login)

	if(!is.null(schema)){
		schemas <- get_schemas(pantry)
		if(!(schema %in% schemas$schema_name)){
			stop(paste0("Schema '", schema, "' doesn't exist. To create it do:\n\n  pantry <- get_pantry()\n  pantry %>% create_schema('", schema, "')\n  pantry %>% set_schema('", schema, "')\n"))
		}
		set_schema(pantry, schema)
	}
	return(pantry)
}

#' Create a schema in the pantry
#' @export
create_schema <- function(pantry, schema_name){
	require(DBI)
	a <- pantry$con %>%
		dbGetQuery(paste0("CREATE SCHEMA ", schema_name %>% sql, ";"))
}

#' Drop a schema from the pantry
#' @export
drop_schema <- function(pantry, schema_name){
	require(DBI)
	a <- pantry$con %>%
		dbGetQuery(paste0("DROP SCHEMA ", schema_name %>% sql, " CASCADE;"))
}

#' Set a schema as the active schema
#' @export
set_schema <- function(pantry, schema_name){
	require(DBI)
	if(is.null(schema_name)){
		search_path <- paste("\"$user\"", "public", sep=",")
	} else {
		search_path <- paste("\"$user\"", schema_name, "public", sep=",")
	}
	a <- pantry$con %>%
		dbGetQuery(paste0("SET search_path TO ", search_path, ";"))
}

#' Get the search path for looking for namespaces
#' @export
get_search_path <- function(pantry){
	require(DBI)
	a <- pantry$con %>% dbGetQuery("SHOW search_path")
	str_split(a$search_path, ", ")[[1]]
}

#' Get all available schemas
#' @export
get_schemas <- function(pantry, verbose=T){
	if(verbose){
		cat("Current schema: ", get_search_path(pantry), "\n", sep="")
	}
	pantry %>%
		tbl(build_sql("SELECT schema_name FROM information_schema.schemata")) %>%
		collect %>%
		data.frame
}

#' Get all tables in the active schema or search path
#' @export
get_tables <- function(pantry, all_tables=FALSE){
	require(DBI)
	tables <- dbGetQuery(pantry$con, paste0(
		"SELECT schemaname, tablename FROM pg_tables ",
		"WHERE schemaname != 'information_schema' AND schemaname !='pg_catalog'"))
	schema <- get_search_path(pantry)
	if(!all_tables){
		if(length(schema) > 0){
			tables <- tables[tables$schemaname %in% schema,]
		}
	}
	return(tables);
}

#' Drop a table
#' @export
drop_table <- function(pantry, table){
	require(DBI)
	x <- dbGetQuery(pantry$con, paste0("DROP TABLE ", table %>% sql, " CASCADE;"))
}

#' Identify a table qualified by a schema. This works around issue 244 in dplyr
#' @export
schema_tbl <- function(pantry, schema_table){
	require(dplyr)
	# see http://stackoverflow.com/questions/21592266/i-cannot-connect-postgresql-schema-table-with-dplyr-package
	# https://github.com/hadley/dplyr/issues/244
	tbl(pantry, structure(paste0("SELECT * FROM ",  schema_table), class=c("sql", "character")))
}

#' this adds the fast argument working around a known problem in dplyr and DBI
#' https://github.com/hadley/dplyr/issues/1471
#' https://github.com/rstats-db/DBI/issues/62
copy_to.src_postgres <- function(
	dest, df, name = deparse(substitute(df)),
	types = NULL, temporary = TRUE, indexes = NULL,
	analyze = TRUE,
	fast=FALSE,
	...
) {
	require(assertthat)
	require(RPostgreSQL)
	assertthat::assert_that(is.data.frame(df), is.string(name), is.flag(temporary))
	class(df) <- "data.frame" # avoid S4 dispatch problem in dbSendPreparedQuery


	if (name %in% dbGetQuery(
		dest$con,
		"SELECT relname FROM pg_class WHERE relkind IN ('r', 'v') AND pg_table_is_visible(oid);")[,1]){
		stop("Table ", name, " already exists.", call. = FALSE)
	}

	types <- if(is.null(types)) db_data_type(dest$con, df) else types
	names(types) <- names(df)

	con <- dest$con
	db_begin(con)
	on.exit(db_rollback(con))

	db_create_table(con, name, types, temporary = temporary)
	if(fast){
		#http://james.hiebert.name/blog/work/2011/10/24/RPostgreSQL-and-COPY-IN/
		#warning this doesn't report any errors on failure!
		field_names <- escape(ident(names(types)), collapse = NULL, con = con)
		fields <- dplyr:::sql_vector(field_names, parens = TRUE, collapse = ", ", con = con)
		sql <- build_sql(
			"COPY ", ident(name), " ", fields, " FROM STDIN;",
			con=con)
	  dbSendQuery(con, sql)
	  postgresqlCopyInDataframe(con, df)
	} else {
		db_insert_into(con, name, df)
	}
	dplyr:::db_create_indexes(con, name, indexes)
	if (analyze) db_analyze(con, name)

	db_commit(con)
	on.exit(NULL)

	tbl(dest, name)
}

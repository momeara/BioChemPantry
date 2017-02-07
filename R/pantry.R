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
get_staging_directory <- function(data_set, pantry_config="~/.pantry_config"){
	pantry_config <- BioChemPantry::get_pantry_config(pantry_config)
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
	pantry_config <- BioChemPantry::get_pantry_config(pantry_config)

	if(!is.null(schema)){
		schema <- tolower(schema)
	}

	pantry <- do.call(dplyr::src_postgres, pantry_config$login)

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
	con <- dplyr::con_acquire(pantry)
	tryCatch({
		a <- DBI::dbGetQuery(con, paste0("CREATE SCHEMA ", schema_name %>% dplyr::sql(), ";"))
	}, finally = {
		dplyr::con_release(pantry, con)
	})
}

#' Drop a schema from the pantry
#' @export
drop_schema <- function(pantry, schema_name){
	con <- dplyr::con_acquire(pantry)
	tryCatch({
		a <- DBI::dbGetQuery(con, paste0("DROP SCHEMA ", schema_name %>% dplyr::sql(), " CASCADE;"))
	}, finally = {
		dplyr::con_release(pantry, con)
	})
}

#' Set a schema as the active schema
#' @export
set_schema <- function(pantry, schema_name){
	if(is.null(schema_name)){
		search_path <- paste("\"$user\"", "public", sep=",")
	} else {
		search_path <- paste("\"$user\"", schema_name, "public", sep=",")
	}

	con <- dplyr::con_acquire(pantry)
	tryCatch({
		DBI::dbGetQuery(con, paste0("SET search_path TO ", search_path, ";"))
	}, finally = {
		dplyr::con_release(pantry, con)
	})
}

#' Get the search path for looking for namespaces
#' @export
get_search_path <- function(pantry){
	con <- dplyr::con_acquire(pantry)
	tryCatch({
		a <- DBI::dbGetQuery(con, "SHOW search_path")
		search_path <- stringr::str_split(a$search_path, ", ")[[1]]
	}, finally = {
		dplyr::con_release(pantry, con)
	})
	search_path
}

#' Get all available schemas
#' @export
get_schemas <- function(pantry, verbose=T){
	if(verbose){
		cat("Current schema: ", BioChemPantry::get_search_path(pantry), "\n", sep="")
	}
	pantry %>%
		dplyr::tbl(dplyr::build_sql("SELECT schema_name FROM information_schema.schemata")) %>%
		dplyr::collect() %>%
		data.frame
}

#' Get all tables in the active schema or search path
#' @export
get_tables <- function(pantry, all_tables=FALSE){
	con <- dplyr::con_acquire(pantry)
	tryCatch({
		tables <- DBI::dbGetQuery(con, paste0(
			"SELECT schemaname, tablename FROM pg_tables ",
			"WHERE schemaname != 'information_schema' AND schemaname !='pg_catalog'"))
		schema <- BioChemPantry::get_search_path(pantry)
		if(!all_tables){
			if(length(schema) > 0){
				tables <- tables[tables$schemaname %in% schema,]
			}
		}
	}, finally = {
		dplyr::con_release(pantry, con)
	})
	tables
}

#' Drop a table
#' @export
drop_table <- function(pantry, table){
	con <- dplyr::con_acquire(pantry)
	tryCatch({
		x <- DBI::dbGetQuery(con, paste0("DROP TABLE ", table %>% sql, " CASCADE;"))
	}, finally = {
		dplyr::con_release(pantry, con)
	})
	x
}

#' Identify a table qualified by a schema. This works around issue 244 in dplyr
#' @export
schema_tbl <- function(pantry, schema_table){
	# see http://stackoverflow.com/questions/21592266/i-cannot-connect-postgresql-schema-table-with-dplyr-package
	# https://github.com/hadley/dplyr/issues/244
	pantry %>%
		dplyr::tbl(structure(paste0("SELECT * FROM ",  schema_table), class=c("sql", "character")))
}

#' this adds the fast argument working around a known problem in dplyr and DBI
#' https://github.com/hadley/dplyr/blob/master/R/tbl-sql.r#L294
#'
#' https://github.com/hadley/dplyr/issues/1471
#' https://github.com/rstats-db/DBI/issues/62
#' @export
copy_to.src_postgres <- function(
	dest,
	df,
	name = deparse(substitute(df)),
	types = NULL,
	temporary = TRUE,
	unique_indexes = NULL,
	indexes = NULL,
	analyze = TRUE,
	fast=FALSE,
	...
) {
	assertthat::assert_that(
		assertthat::is.data.frame(df),
		assertthat::is.string(name),
		assertthat::is.flag(temporary))
	class(df) <- "data.frame" # avoid S4 dispatch problem in dbSendPreparedQuery

	con <- dplyr::con_acquire(pantry)
	tryCatch({
		if (isTRUE(dplyr::db_has_table(con, name))) {
			stop("Table ", name, " already exists.", call. = FALSE)
		}

		types <- types %||% dplyr::db_data_type(con, df)
		names(types) <- names(df)

		dplyr::db_begin(con)
		tryCatch({
			dplyr::db_create_table(con, name, types, temporary = temporary)

			if(fast){
				#http://james.hiebert.name/blog/work/2011/10/24/RPostgreSQL-and-COPY-IN/
				#warning this doesn't report any errors on failure!
				field_names <- dplyr::escape(dplyr::ident(names(types)), collapse = NULL, con = con)
				fields <- dplyr:::sql_vector(field_names, parens = TRUE, collapse = ", ", con = con)
				sql <- dplyr::build_sql(
					"COPY ", dplyr::ident(name), " ", fields, " FROM STDIN;", con=con)
				DBI::dbSendQuery(con, sql)
				RPostgreSQL::postgresqlCopyInDataframe(con, df)
			} else {
				dplyr::db_insert_into(con, name, df)
			}

			dplyr:::db_create_indexes(con, name, unique_indexes, unique = TRUE)
			dplyr:::db_create_indexes(con, name, indexes, unique = FALSE)

			if (analyze) dplyr::db_analyze(con, name)
		}, error = function(err){
			dplyr::db_rollback(con)
		})
	}, finally = {
		dplyr::con_release(pantry, con)
	})

	dplyr::tbl(dest, name)
}

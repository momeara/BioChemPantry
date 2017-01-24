# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


#library(plyr)
#library(dplyr)
#library(stringr)
#library(jsonlite)


#' attach to the data repository, optionally using a schema. If the
#' specified schema doesn't exist, create it.
#'
#' @param login either the name of a json file with login information or a
#'   list of login info. The keys should correspond with the argument
#'   keys of the dplyr::src_postgres function.
#' @return a connection to the the data_repo
#' @export
get_data_repo <- function(schema=NULL, login="~/.data_repo_login"){
	if(is.character(login)){
		login <- fromJSON(login)
	} else if(!is.list(login)){
		stop("Login must either be the name of a json file of login info, or a list of login info.\n")
	}
	if(!is.null(schema)){
		schema <- tolower(schema)
	}
	data_repo <- do.call(dplyr::src_postgres, login)

	if(!is.null(schema)){
		schemas <- get_schemas(data_repo)
		if(!(schema %in% schemas$schema_name)){
			stop(paste0("Schema '", schema, "' doesn't exist. To create it do:\n\n  data_repo <- get_data_repo()\n  data_repo %>% create_schema('", schema, "')\n  data_repo %>% set_schema('", schema, "')\n"))
		}
		set_schema(data_repo, schema)
	}
	return(data_repo)
}

#' Create a schema in the data_repo
#' @export
create_schema <- function(data_repo, schema_name){
	require(DBI)
	a <- data_repo$con %>%
		dbGetQuery(paste0("CREATE SCHEMA ", schema_name %>% sql, ";"))
}

#' Drop a schema from the data_repo
#' @export
drop_schema <- function(data_repo, schema_name){
	require(DBI)
	a <- data_repo$con %>%
		dbGetQuery(paste0("DROP SCHEMA ", schema_name %>% sql, " CASCADE;"))
}

#' Set a schema as the active schema
#' @export
set_schema <- function(data_repo, schema_name){
	require(DBI)
	if(is.null(schema_name)){
		search_path <- paste("\"$user\"", "public", sep=",")
	} else {
		search_path <- paste("\"$user\"", schema_name, "public", sep=",")
	}
	a <- data_repo$con %>%
		dbGetQuery(paste0("SET search_path TO ", search_path, ";"))
}

#' Get the search path for looking for namespaces
#' @export
get_search_path <- function(data_repo){
	require(DBI)
	a <- data_repo$con %>% dbGetQuery("SHOW search_path")
	str_split(a$search_path, ", ")[[1]]
}

#' Get all available schemas
#' @export
get_schemas <- function(data_repo, verbose=T){
	if(verbose){
		cat("Current schema: ", get_search_path(data_repo), "\n", sep="")
	}
	data_repo %>%
		tbl(build_sql("SELECT schema_name FROM information_schema.schemata")) %>%
		collect %>%
		data.frame
}

#' Get all tables in the active schema or search path
#' @export
get_tables <- function(data_repo, all_tables=FALSE){
	require(DBI)
	tables <- dbGetQuery(data_repo$con, paste0(
		"SELECT schemaname, tablename FROM pg_tables ",
		"WHERE schemaname != 'information_schema' AND schemaname !='pg_catalog'"))
	schema <- get_search_path(data_repo)
	if(!all_tables){
		if(length(schema) > 0){
			tables <- tables[tables$schemaname %in% schema,]
		}
	}
	return(tables);
}

#' Drop a table
#' @export
drop_table <- function(data_repo, table){
	require(DBI)
	x <- dbGetQuery(data_repo$con, paste0("DROP TABLE ", table %>% sql, " CASCADE;"))
}

#' Identify a table qualified by a schema. This works around issue 244 in dplyr
#' @export
schema_tbl <- function(data_repo, schema_table){
	require(dplyr)
	# see http://stackoverflow.com/questions/21592266/i-cannot-connect-postgresql-schema-table-with-dplyr-package
	# https://github.com/hadley/dplyr/issues/244
	tbl(data_repo, structure(paste0("SELECT * FROM ",  schema_table), class=c("sql", "character")))
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

Once ingredients have been procured, they are shelved in the pantry
ready for use.

Bioinformatics datasets are provided by a wide range of groups, in a
wide range formats, and are rapidly changing. This is an opinionated
guide to collect and curate data sets. It sets up a local database of
datasets that can be used for an integrated analysis.


## to install:
1. Install a postgres database (http://www.postgresql.org/)

1a. Create the `json` file ~/.pantry_login that will be passed to `dplyr::src_postgres` to login ot the database. For example:

      "staging_directory" : "/mnt/nfs/work/momeara/pantry_sets",
      "login" : {
          "dbname" : "<database name",
          "host" : "<host>",
          "user" : "<user>",
          "password" : "<password>",
          "port" : <port>
      }     

2. Install the package in R

   install.packages("devtools")
   devtools::install_github("momeara/BioChemPantry")
  
# usage:

1. Install datasets

    library(BioChemPantry)
  
2. Load datasets by following the vignettes in `vignettes/sets/`
  
3. Use datasets

    library(plyr)
    library(dplyr)
    library(BioChemPantry)
    pantry <- get_pantry(schema=<dataset>)

    tbl <- pantry %>% schema_tbl("<tbl>")


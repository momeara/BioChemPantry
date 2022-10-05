Once ingredients have been procured, they are shelved in the pantry
ready for use.

Bioinformatics datasets are provided by a wide range of groups, in a
wide range formats, and are rapidly changing. This is an opinionated
guide to collect and curate data sets. It sets up a local database of
datasets that can be used for an integrated analysis.

ChEMBL is evolving each year. It might on ChEMBL25 does not mean it would with future releases

## to install:
1. Install a postgres database (http://www.postgresql.org/)

1a. Create the `json` file ~/.pantry_login that will be passed to `dplyr::src_postgres` to login ot the database. For example:

```json
{
  "staging_directory" : "/mnt/nfs/work/momeara/pantry_sets",
  "login" : {
    "dbname" : "<database name",
    "host" : "<host>",
    "user" : "<user>",
    "password" : "<password>",
    "port" : <port>
  }
}
```

2. Install the package in R

```R
install.packages("devtools")
devtools::install_github("momeara/BioChemPantry")
```

# usage:

1. Install datasets

```R
library(BioChemPantry)
```

2. Load datasets by following the vignettes in `vignettes/sets/`
  
3. Use datasets

```R
library(plyr)
library(dplyr)
library(BioChemPantry)
pantry <- get_pantry(schema=<dataset>)

tbl <- pantry %>% schema_tbl("<tbl>")
```

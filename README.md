Once ingredients have been procured, they are shelved in the pantry
ready for use.

Bioinformatics datasets are provided by a wide range of groups, in a
wide range formats, and are rapidly changing. This is an opinionated
guide to collect and curate data sets. It sets up a local database of
datasets that can be used for an integrated analysis.


## to install:
1. Install a postgres database (http://www.postgresql.org/)

For example, to install postgres on a mac using Homebrew

    brew install postgresql
    
then start the database with

    brew services start postgresql

Setup up the database

    psql postgres
    # now in the psql prompt
    
    CREATE ROLE pantry WITH LOGIN PASSWORD ‘password’;
    CREATE DATABASE pantry;
    GRANT CREATE ON DATABASE pantry TO pantry;

Enable password-free login (see (documentation](https://tableplus.com/blog/2019/09/how-to-use-pgpass-in-postgresql.html))

     echo "<host>:<port>:<dbname>:<user>:<password>" >> ${HOME}/.pgpass
     # where the values in <> are specific for your setup
     sudo chmod 600 ~/.pgpass
     echo "export PGPASSFILE=\'${HOME}.pgpass\'" >> ${HOME}/.zshrc


1a. Create the `json` file ~/.pantry_config that will be passed to
    `dplyr::src_postgres` to login to the database. For example:

```json
{
  "staging_directory" : "/Users/maom/opt/pantry_sets",
  "login" : {
    "dbname" : "<database name>",
    "host" : "<host>",
    "user" : "<user>",
    "port" : <port>
  }
}
```

2. Install the package in R

```R
install.packages("remotes")
remotes::install_github("momeara/BioChemPantry")
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

tbl <- pantry |> schema_tbl("<tbl>")
```

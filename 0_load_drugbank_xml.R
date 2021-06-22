library(curl)
library(BioChemPantry)
library(Zr)
library(plyr)
library(dplyr)

schema <- paste0("drugbank")
staging_directory <- BioChemPantry::get_staging_directory(schema)

dir.create(paste0(staging_directory, "/dump"), recursive=TRUE)

####### to download from website using curl
### if this doesn't work, download directly from drugbank website
dump_fname <- paste0(staging_directory, "/dump/drugbank_compounds.zip")
h <- curl::new_handle()
curl::handle_setopt(handle = h, httpauth = 1, userpwd = "user:pwd")
curl::curl_download(url = "https://go.drugbank.com/releases/5-1-8/downloads/all-full-database", 
                    destfile = dump_fname, quiet = FALSE, handle = h)
unzip(dump_fname, exdir=paste0(staging_directory, "/dump"))


drugbank_substances <- Zr::catalog_items(
  "dball",
  output_fields=c(
    "zinc_id",
    "supplier_code",
    "substance.preferred_name",
    "substance.smiles"),
  result_batch_size=10000,
  verbose=T) %>%
  dplyr::select(
    accession=supplier_code,
    zinc_id,
    zinc_name = substance.preferred_name,
    zinc_smiles = substance.smiles)

drugbank_substances %>% readr::write_tsv(
  paste0(staging_directory, "/dball_substances.tsv"))

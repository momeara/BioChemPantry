library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(assertr)
library(BioChemPantry)

source("scripts/parameters.R")

pantry <- BioChemPantry::get_pantry(schema)
staging_directory <- BioChemPantry::get_staging_directory(schema)

chembl_targets <- pantry |>
	dplyr::tbl("component_sequences") |>
	dplyr::rename(uniprot_accn = accession) |>
	dplyr::select(-tax_id, -organism) |> # get these from uniprot instead
	dplyr::collect(n = Inf) |>
	dplyr::filter(!is.na(uniprot_accn)) |> # ribosomal RNA
	dplyr::filter(!(uniprot_accn %in% c(
		"A0A090LW34", # deleted
		"A0A0I9CLG9", # deleted
		"C6EG73"))) # obsolete


# not sure why this REST query isn't working but just downloading it from the browser does
url <- "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cxref_veupathdb%2Cgene_primary%2Corganism_id%2Corganelle%2Cft_act_site%2Cft_binding%2Ccc_cofactor%2Cft_dna_bind%2Cec%2Ccc_function%2Ccc_pathway%2Cph_dependence%2Cft_site%2Ccc_caution%2Ckeyword%2Cprotein_existence%2Ccomment_count%2Cfeature_count%2Ccc_interaction%2Ccc_subunit%2Cgo_p%2Cgo_c%2Cgo_f%2Cgo_id%2Cft_transmem%2Clit_doi_id%2Cversion%2Cxref_generif%2Cxref_pdb%2Cxref_alphafolddb%2Cxref_complexportal%2Cxref_chembl%2Cxref_guidetopharmacology%2Cxref_drugbank%2Cxref_cazy%2Cxref_geneid%2Cxref_wbparasitetranscriptprotein%2Cxref_cgd%2Cxref_hgnc%2Cxref_opentargets%2Cxref_nextprot%2Cxref_orthodb&format=tsv&query=%28%28database%3Achembl%29%29"
uniprot_targets_full <- httr::GET(
  url = url,
  httr::user_agent("httr mattjomeara@gmail.com")) |>
  httr::content() |>
  rawConnection() |>
  gzcon() |>
  xml2::read_xml() |>
  xml2::xml_children() |>
  head(-1) # copyright notice

uniprot_targets <- readr::read_tsv(
  file = paste0(
    staging_directory,
    "/data/uniprotkb_database_chembl_2024_01_22.tsv"),
  show_col_types = FALSE)

uniprot_targets <- uniprot_targets |>
  dplyr::transmute(
	  uniprot_accn = Entry,
	  uniprot_entry = `Entry Name`,
    reviewed = Reviewed,
	  protein_name = `Protein names`,
	  gene_name = `Gene Names (primary)`,
	  gene_names = `Gene Names`,
    organism = Organism,
	  taxon_id = `Organism (ID)`,
	  length = Length,
	  veupathdb_id = VEuPathDB,
	  site = Site,
	  active_site = `Active site`,
    binding_site = `Binding site`,
	  cofactor = Cofactor,
	  dna_binding = `DNA binding`,
	  ec_number = `EC number`,
	  function_cc = `Function [CC]`,
	  pathway = Pathway,
	  pH_dependence = `pH dependence`,
    caution = Caution,
    keywords = Keywords,
    protein_existance = `Protein existence`,
    comment = Comments,
    interacts_with = `Interacts with`,
    subunit_structure = `Subunit structure`,
    go_bp = `Gene Ontology (biological process)`,
	  go_cc = `Gene Ontology (cellular component)`,
	  go_mf = `Gene Ontology (molecular function)`,
	  go_ids = `Gene Ontology IDs`,
	  transmembrane = Transmembrane,
	  doi_ids = `DOI ID`,
	  uniprot_entry_version = `Entry version`,
	  gene_rif = GeneRIF,
	  pdb_ids = PDB,
	  alphafolddb_ids = AlphaFoldDB,
	  complex_portal_ids = ComplexPortal,
	  chembl_ids = ChEMBL,
	  iuphar_ids = GuidetoPHARMACOLOGY,
	  drug_bank_ids = DrugBank,
	  cazy_ids = CAZy,
	  gene_ids = GeneID,
	  cgd_ids = CGD,
	  hgnc_ids = HGNC,
	  open_targets_ids = OpenTargets,
	  nextprot_ids = neXtProt,
	  orthodb_ids = OrthoDB)
	  
	  
	  
uniprot_targets |>
	readr::write_tsv(
	  file = paste0(staging_directory, "/data/uniprot_targets_raw.tsv"))

uniprot_targets <- readr::read_tsv(
	paste0(staging_directory, "/data/uniprot_targets_raw.tsv"))

problems <- dplyr::full_join(
  chembl_targets |>
    dplyr::transmute(
      uniprot_accn,
      chembl_uniprot_accn = uniprot_accn,
      description),
  uniprot_targets |>
    dplyr::transmute(
      uniprot_accn,
      uniprot_uniprot_accn = uniprot_accn),
  by=c("uniprot_accn")) |>
  summarize_map(
    x_cols = c("chembl_uniprot_accn"),
    y_cols = c("uniprot_uniprot_accn"))

target_classes <- readr::read_tsv(
  file = paste0(staging_directory, "/data/chembl_target_classes.tsv"),
  show_col_types = FALSE) |>
	dplyr::select(
	  uniprot_accn, class_1, class_2, class_3, class_4, class_5, class_6) |>
	dplyr::distinct(uniprot_accn, .keep_all = TRUE)

target_info <- chembl_targets |>
	dplyr::inner_join(uniprot_targets, by = "uniprot_accn") |>
	dplyr::left_join(target_classes, by = "uniprot_accn")


target_info |>
	readr::write_tsv(paste0(staging_directory, "/data/target_info.tsv"))

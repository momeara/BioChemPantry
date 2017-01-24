# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(assertr)


chembl_db <- get_data_repo("chembl21")

chembl_targets <- chembl_db %>% dplyr::tbl("component_sequences") %>%
	dplyr::rename(uniprot_accn = accession) %>%
	dplyr::select(-tax_id, -organism) %>% # get these from uniprot instead
	dplyr::collect() %>%
	dplyr::filter(!is.na(uniprot_accn)) %>% # ribosomal RNA
	dplyr::filter(!(uniprot_accn %in% c(
		"A0A090LW34", # deleted
		"A0A0I9CLG9", # deleted
		"C6EG73"))) # obsolete

uniprot_target_columns <- paste(
			"id",
			"entry name",
			"protein names",
			"genes(PREFERRED)",
			"keywords",
			"organism",
			"organism-id",
			"ec",
			"feature(BINDING SITE)",
			"comment(PATHWAY)",
			"comment(FUNCTION)",
			"feature(ACTIVE SITE)",
			"feature(SITE)",
			"features",
			"interactor",
			"comment(TISSUE SPECIFICITY)",
			"go(biological process)",
			"go(molecular function)",
			"go(cellular component)",
			"go-id",
			"comment(ALLERGEN)",
			"comment(BIOTECHNOLOGY)",
			"comment(DISRUPTION PHENOTYPE)",
			"comment(DISEASE)",
			"comment(PHARMACEUTICAL)",
			"comment(TOXIC DOSE)",
			"feature(MUTAGENESIS)",
			"comment(SUBCELLULAR LOCATION)",
			"3d",
			"comment(DOMAIN)",
			"families",
			"database(Ensembl)",
			"database(GeneID)",
			"database(PDB)",
			"database(BioGrid)",
			"database(ChEMBL)",
			"database(GuidetoPHARMACOLOGY)",
			"database(Reactome)",
			"database(BRENDA)",
			"database(BioCyc)",
			"database(Pfam)", sep=",")

uniprot_targets <- httr::GET(
	url="http://www.uniprot.org/uniprot/",
	httr::user_agent("httr mattjomeara@gmail.com"),
	query=list(
		query="database:(type:chembl)",
		format='tab',
		compress='yes',
		columns=uniprot_target_columns)) %>%
	httr::content() %>%
	rawConnection() %>%
	gzcon() %>%
	readr::read_tsv()

uniprot_targets <- uniprot_targets %>% rename(
	uniprot_accn = Entry,
	uniprot_entry = `Entry name`,
	gene_name = `Gene names  (primary )`,
	entrez_id_alt = `Cross-reference (GeneID)`,
	ensembl_id_alt = `Cross-reference (Ensembl)`,
	keywords = Keywords,
	organism = Organism,
	tax_id = `Organism ID`,
	ec_number =  `EC number`,
	binding_site = `Binding site`,
	uniprot_pathway = Pathway,
	uniprot_function = `Function [CC]`,
	active_site = `Active site`,
	site = Site,
	protein_features = Features,
	interactors = `Interacts with`,
	tissue_specificity = `Tissue specificity`,
	GO_bp =	`Gene ontology (biological process)`,
	GO_mf = `Gene ontology (molecular function)`,
	GO_cc = `Gene ontology (cellular component)`,
	GO_terms = `Gene ontology IDs`,
	allergenic_properties = `Allergenic properties`,
	biotech_use = `Biotechnological use`,
	disruption_phenotype = `Disruption phenotype`,
	disease_relevance = `Involvement in disease`,
	pharmaceutical_use = `Pharmaceutical use`,
	toxic_dose = `Toxic dose`,
	mutagenesis = `Mutagenesis`,
	uniprot_subcellular_location = `Subcellular location [CC]`,
	experimental_structures = `3D`,
	uniprot_domains = `Domain [CC]`,
	protein_families = `Protein families`,
	pdb_ids = `Cross-reference (PDB)`,
	biogrid_ids = `Cross-reference (BioGrid)`,
	chembl_ids = `Cross-reference (ChEMBL)`,
	iuphar_ids = `Cross-reference (GuidetoPHARMACOLOGY)`,
	reactom_ids = `Cross-reference (Reactome)`,
	brenda_ids = `Cross-reference (BRENDA)`,
	biocyc_ids = `Cross-reference (BioCyc)`,
	pfam_ids = `Cross-reference (Pfam)`) %>%
	dplyr::mutate(
		entrez_id = entrez_id_alt %>% str_extract("^[0-9]+")

uniprot_targets %>% write_tsv("data/uniprot_targets_raw_160627.tsv")

problems <- dplyr::full_join(
	chembl_targets %>%
		dplyr::transmute(uniprot_accn, chembl_uniprot_accn=uniprot_accn, description),
	uniprot_targets %>%
		dplyr::transmute(uniprot_accn, uniprot_uniprot_accn=uniprot_accn),
	by=c("uniprot_accn")) %>%
	summarize_map(
		xcols=c("chembl_uniprot_accn"),
		ycols=c("uniprot_uniprot_accn"))

# this is a half baked attemed at getting reasonable names for uniprot entry
# the next step is to parse the XML to get out the relevant information
uniprot_targets_full <- httr::GET(
	url="http://www.uniprot.org/uniprot/",
	user_agent("httr mattjomeara@gmail.com"),
	query=list(
		query="database:(type:chembl)",
		format='xml',
		limit=5,
		compress='yes')) %>%
	httr::content() %>%
	rawConnection() %>%
	gzcon() %>%
	xml2::read_xml() %>%
	xml2::xml_children() %>%
	head(-1) # copyright notice


target_classes <- read_tsv("data/chembl_target_classes.tsv") %>%
	dplyr::select(uniprot_accn, class_1, class_2, class_3, class_4, class_5, class_6) %>%
	dplyr::distinct(uniprot_accn, .keep_all=T)

target_info <- chembl_targets %>%
	dplyr::left_join(uniprot_targets, by="uniprot_accn") %>%
	dplyr::left_join(target_classes, by="uniprot_accn")


target_info %>%
	write_tsv("data/target_info.tsv")

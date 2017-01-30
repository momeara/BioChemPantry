# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(BioChemPantry)
source("~/work/sea/scripts/id_utils.R")


# the idea is to use HGNC's HCOP data to identify orthologs
# So far I've collected the data and started to look at it.

staging_directory <- get_staging_directory("hgnc/release_160324")

human_orthologs_fname <- paste0(staging_directory, "/data/human_all_hcop_sixteen_column.txt")
system(paste0(
   "wget 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_all_hcop_sixteen_column.txt.gz' ",
   "-O ", human_orthologs_fname))
system(past0("cd ", staging_directory, "/data && tar -xzvf ", human_orthologs_fname))


#prepend "row_id	" to the header row
system(paste0("echo 'row_id	' | cat -  ", human_orthologs_fname, " > /tmp/out && mv /tmp/out ", human_orthologs_fname))
human_orthologs <- read_tsv(human_orthologs_fname)

human_orthologs %>% summarize_map("human_symbol", "ortholog_species_ensembl_gene")
human_orthologs %>% summarize_map("human_symbol", "human_entrez_gene")

human_orthologs <- human_orthologs %>%
	mutate(ortholog_species_prefix = ortholog_species_ensembl_gene %>% str_extract("^[a-zA-Z]+"))

Chimp	ENSPTR	Pan troglodytes (Chimpanzee)
Macaque	ENSMMU	Macaca mulatta (Macaque)
Mouse	ENSMUS	Mus musculus (Mouse)
Rat	ENSRNO	Rattus norvegicus (Rat)
Dog	ENSCAF	Canis lupus familiaris (Dog)
Horse	ENSECA	Equus caballus (Horse)
Cow	ENSBTA	Bos taurus (Cow)
Pig	ENSSSC	Sus scrofa (Pig)
Opossum	ENSMOD	Monodelphis domestica (Opossum)
Platypus	ENSOAN	Ornithorhynchus anatinus (Platypus)
Chicken	ENSGAL	Gallus gallus (Chicken)
Anole lizard	ENSACA	Anolis carolinensis (Anole lizard)
Xenopus	ENSXET	Xenopus tropicalis (Xenopus)
Zebrafish	ENSDAR	Danio rerio (Zebrafish)
C.elegans	ENSCEL	Caenorhabditis elegans (Caenorhabditis elegans)
Fruitfly	FB	Drosophila melanogaster (Fruitfly)
S.cerevisiae	ENSSCE	Saccharomyces cerevisiae (Saccharomyces cerevisiae)



   ortho_species
1        ENSPTRG Chimp
2           <NA>
3        ENSMMUG Macaque
4        ENSMUSG Mouse
5        ENSRNOG Rat
6        ENSCAFG Dog
7        ENSECAG Horse
8        ENSBTAG Cow
9        ENSSSCG Pig
10       ENSMODG Opossum
11       ENSOANG Platypus
12       ENSGALG Chicken
13       ENSGALP Chicken (Protein)
14       ENSACAG Anole lizard
15       ENSXETG Xenopus
16       ENSXETT Xenopus (Transcript)
17       ENSDARG Zebrafish
18        WBGene WormBase Gene
19          FBgn 
20           YAL
21           YBR
22           YGL
23           YHR
24           YIL
25           YLR
26           YNL
27           YIR
28           YDR
29           YBL
30           YOR
31           YNR
32           YGR
33           YOL
34           YKL
35           YPL
36           YMR
37           YLL
38           YHL
39           YFR
40           YER
41           YCR
42           YDL
43           YJL
44           YJR
45           YPR
46           YFL
47           YEL
48           YML
49           YKR
50             Q
51           YCL
52           YAR
